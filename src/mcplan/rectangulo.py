"""Ejercicio de la pizarra: rectángulo 2D que se enfría por todo su borde.

    rho cp dT/dt - k lap(T) = Q   en un rectángulo Lx x Ly (4 m x 2 m)
    T(x, y, 0) = Ti (20 °C)
    T = Tb (10 °C) en todo el borde Γ, para todo t (Dirichlet, "C.B. esencial")

Sin hidratación (Q = 0) la solución exacta es el producto de dos losas 1D
(la misma serie de Fourier que usa `mcplan verificar`):

    (T - Tb) / (Ti - Tb) = X(x, t) * Y(y, t)
    X(x, t) = sum_{m impar} 4/(m pi) sin(m pi x / Lx) exp(-alpha (m pi / Lx)^2 t)

A tiempos largos manda el primer modo, que decae con exp(-alpha lambda1 t),
lambda1 = (pi/Lx)^2 + (pi/Ly)^2: eso da el tiempo para "estabilizar".

El calor Q puede ser:
- la curva de Ju y Lei (2019) con edad equivalente (mismo modelo del solver 1D):
  depende de T, no hay solución cerrada;
- el pulso coseno de la clase, Q(t) = Q_max/2 (1 - cos(2 pi t/T)) para t < T:
  sólo depende de t, el problema sigue lineal y la solución exacta es la de
  arriba más la respuesta de cada modo a la fuente.
Con calor, el tiempo para estabilizar se busca simulando.

Coordenadas como en la pizarra: origen en el centro del rectángulo
(x en [-Lx/2, Lx/2], y en [-Ly/2, Ly/2]). Con origen="esquina", x en [0, Lx].
"""

from dataclasses import dataclass
from math import log, pi, sqrt
from typing import Optional

import numpy as np

from .hidratacion import CalorPulsoCoseno, HidratacionJuLei2019
from .solver_nd import paso_critico, paso_critico_malla, simular_serie
from .verificacion import solucion_analitica

CARAS = ("x_min", "x_max", "y_min", "y_max")


def _a_esquina(punto, Lx, Ly, origen):
    """Convierte un punto del usuario a coordenadas con origen en la esquina."""
    if origen not in ("centro", "esquina"):
        raise ValueError(f"origen debe ser 'centro' o 'esquina', no {origen!r}")
    x, y = punto
    if origen == "centro":
        x, y = x + Lx / 2, y + Ly / 2
    tol = 1e-9 * max(Lx, Ly)
    if not (-tol <= x <= Lx + tol and -tol <= y <= Ly + tol):
        raise ValueError(f"El punto {tuple(punto)} está fuera del rectángulo")
    return min(max(x, 0.0), Lx), min(max(y, 0.0), Ly)


def en_el_borde(punto, Lx, Ly, origen="centro") -> bool:
    x, y = _a_esquina(punto, Lx, Ly, origen)
    tol = 1e-9 * max(Lx, Ly)
    return min(x, Lx - x, y, Ly - y) <= tol


def _pesos_bilineales(puntos_esquina, Lx, Ly, nx, ny) -> np.ndarray:
    """Matriz W tal que T(puntos) = W @ P.ravel(), con P el campo completado con Tb en las caras.

    Interpolación bilineal entre centros de celda; entre el último centro y
    el borde se interpola contra la temperatura de la cara (Dirichlet).
    """
    gx = np.concatenate(([0.0], (np.arange(nx) + 0.5) * Lx / nx, [Lx]))
    gy = np.concatenate(([0.0], (np.arange(ny) + 0.5) * Ly / ny, [Ly]))
    W = np.zeros((len(puntos_esquina), (nx + 2) * (ny + 2)))
    for f, (x, y) in enumerate(puntos_esquina):
        i = min(int(np.searchsorted(gx, x, side="right")) - 1, nx)
        j = min(int(np.searchsorted(gy, y, side="right")) - 1, ny)
        sx = (x - gx[i]) / (gx[i + 1] - gx[i])
        sy = (y - gy[j]) / (gy[j + 1] - gy[j])
        for di, wx in ((0, 1 - sx), (1, sx)):
            for dj, wy in ((0, 1 - sy), (1, sy)):
                W[f, (i + di) * (ny + 2) + (j + dj)] += wx * wy
    return W


def _completar(T, Tb) -> np.ndarray:
    """Agrega las caras (a Tb) a un campo (nx, ny) o a una serie (n, nx, ny) y aplana cada campo."""
    T = np.asarray(T, dtype=float)
    ancho = ((0, 0),) * (T.ndim - 2) + ((1, 1), (1, 1))
    P = np.pad(T, ancho, constant_values=Tb)
    return P.reshape(P.shape[:-2] + (-1,))


def temperatura_en(T, Lx, Ly, Tb, punto, origen="centro") -> float:
    """Temperatura en un punto, interpolando bilinealmente entre centros de celda."""
    W = _pesos_bilineales([_a_esquina(punto, Lx, Ly, origen)], Lx, Ly, *np.shape(T))
    return float((W @ _completar(T, Tb))[0])


def contorno_interior(Lx, Ly, nx, ny, profundidad) -> list[tuple[float, float]]:
    """Puntos a `profundidad` del borde (coordenadas desde la esquina), donde se mide Tmin.

    Sobre esa línea la interpolación es lineal por tramos entre centros de
    celda, así que basta muestrear las abscisas de los centros y las esquinas.
    """
    d = profundidad
    if not 0 < d < min(Lx, Ly) / 2:
        raise ValueError(f"profundidad = {d} m debe estar entre 0 y la mitad del lado menor")
    xs = np.concatenate(([d], [x for x in (np.arange(nx) + 0.5) * Lx / nx if d < x < Lx - d], [Lx - d]))
    ys = np.concatenate(([d], [y for y in (np.arange(ny) + 0.5) * Ly / ny if d < y < Ly - d], [Ly - d]))
    return ([(x, d) for x in xs] + [(x, Ly - d) for x in xs]
            + [(d, y) for y in ys] + [(Lx - d, y) for y in ys])


def _losa(s, t, L, alpha):
    """Serie 1D normalizada (1 adentro, 0 en las caras), con los términos necesarios."""
    if t <= 0:
        s = np.asarray(s, dtype=float)
        return np.where((s > 0) & (s < L), 1.0, 0.0)
    # se corta cuando exp(-alpha (n pi / L)^2 t) < 1e-15
    n = min(4001, int(L / pi * sqrt(35.0 / (alpha * t))) + 3)
    return solucion_analitica(s, t, L, alpha, 1.0, 0.0, n_terminos=n)


def _respuesta_pulso(xs, ys, t, Lx, Ly, alpha, rho, cp, pulso, n_terminos=401):
    """Aumento de T por el pulso coseno con el contorno a 0 °C, en la grilla xs × ys (desde la esquina).

    La fuente uniforme se descompone en los modos sin(m pi x/Lx) sin(n pi y/Ly)
    (m, n impares) con coeficientes c = 16/(m n pi^2). Cada modo cumple
        dA/dt = -a A + c Q(t)/(rho cp),   a = alpha pi^2 (m^2/Lx^2 + n^2/Ly^2)
    y con Q = Q_max/2 (1 - cos w t) la integral G = int_0^t e^{-a(t-s)} Q(s) ds es cerrada.
    """
    m = np.arange(1, n_terminos + 1, 2)
    a = alpha * pi**2 * ((m[:, None] / Lx) ** 2 + (m[None, :] / Ly) ** 2)
    c = 16.0 / (pi**2 * m[:, None] * m[None, :])
    G = pulso.integral_amortiguada(a, t)
    Sx = np.sin(pi * np.outer(xs, m) / Lx)
    Sy = np.sin(pi * np.outer(ys, m) / Ly)
    return Sx @ (c * G) @ Sy.T / (rho * cp)


def solucion_exacta(xs, ys, t, Lx, Ly, alpha, Ti, Tb, pulso=None, rho=None, cp=None):
    """Solución exacta en la grilla xs × ys (medidos desde la esquina): arreglo (len(xs), len(ys)).

    Sin calor es el producto de dos losas; con el pulso coseno (problema lineal)
    se le suma la respuesta a la fuente.
    """
    xs = np.atleast_1d(np.asarray(xs, dtype=float)).ravel()
    ys = np.atleast_1d(np.asarray(ys, dtype=float)).ravel()
    T = Tb + (Ti - Tb) * np.outer(_losa(xs, t, Lx, alpha), _losa(ys, t, Ly, alpha))
    if pulso is not None:
        T = T + _respuesta_pulso(xs, ys, t, Lx, Ly, alpha, rho, cp, pulso)
    return T


def tiempo_estabilizacion(Lx, Ly, alpha, Ti, Tb, tolerancia) -> float:
    """Tiempo para que |T - Tb| < tolerancia en todo el rectángulo, sin hidratación (primer modo).

    El primer modo vale (16/pi^2)(Ti - Tb) en el centro y decae con
    exp(-alpha lambda1 t): t = ln((16/pi^2)|Ti - Tb| / tol) / (alpha lambda1).
    """
    amplitud = 16 / pi**2 * abs(Ti - Tb)
    lam1 = (pi / Lx) ** 2 + (pi / Ly) ** 2
    return max(0.0, log(amplitud / tolerancia)) / (alpha * lam1)


def _estabilizacion_simulada(T_inicial, alpha, dx, rho, cp, bordes, hidratacion, Tb, tolerancia,
                             bloque=30 * 86400.0, tope=730 * 86400.0) -> float:
    """Busca, simulando por bloques de 30 días (cada 6 h, dt <= 600 s), cuándo |T - Tb| < tol en todo."""
    T, te, t0 = T_inicial, 0.0, 0.0
    while t0 < tope:
        s = simular_serie(T, alpha, (dx, dx), bloque, 6 * 3600.0, rho=rho, cp=cp, bordes=bordes,
                          dt_max=600.0, hidratacion=hidratacion, te_inicial=te)
        dentro = np.max(np.abs(s.T - Tb), axis=(1, 2)) < tolerancia
        if dentro.any():
            return t0 + float(s.t[np.argmax(dentro)])
        T, te, t0 = s.T[-1], s.te[-1], t0 + bloque
    raise ValueError(f"No se estabiliza a ±{tolerancia} °C en {tope / 86400:.0f} días: indique --t")


@dataclass
class ResultadoRectangulo:
    Lx: float
    Ly: float
    dx: float
    k: float
    rho: float
    cp: float
    Ti: float
    Tb: float
    origen: str
    puntos: list[tuple[float, float]]
    modelo: object  # calor de hidratación (HidratacionJuLei2019, CalorPulsoCoseno) o None
    t: np.ndarray  # [s]
    T: np.ndarray  # (n_cuadros, nx, ny)
    te: Optional[np.ndarray]  # edad equivalente [h] (None sin hidratación)
    T_puntos: np.ndarray  # (n_cuadros, n_puntos), numérica
    T_puntos_exacta: Optional[np.ndarray]  # (n_cuadros, n_puntos); None si no hay solución exacta (Ju y Lei)
    error_campo: Optional[np.ndarray]  # máx |T - exacta| por cuadro (0 en t = 0); None si no hay exacta
    Tmin_superficie: np.ndarray  # mínimo a `profundidad` del borde, por cuadro
    profundidad: float  # [m] (25 mm, criterio del curso)
    dt: float
    dt_critico: float  # real, con el borde Dirichlet
    dt_critico_interior: float  # dx^2 / (4 alpha)
    dt_max: Optional[float]  # tope de precisión (con hidratación)
    tolerancia: float
    n_cuadros: int  # cuadros de la animación (se toman de los instantes guardados)

    def indices_animacion(self) -> np.ndarray:
        return np.unique(np.linspace(0, len(self.t) - 1, self.n_cuadros + 1).round().astype(int))

    @property
    def hidratacion(self) -> bool:
        return self.te is not None

    @property
    def alpha(self) -> float:
        return self.k / (self.rho * self.cp)

    @property
    def tau(self) -> float:
        """Constante de tiempo del modo más lento de conducción [s]."""
        return 1.0 / (self.alpha * ((pi / self.Lx) ** 2 + (pi / self.Ly) ** 2))

    @property
    def t_estabilizacion_estimado(self) -> float:
        """Estimación analítica SIN hidratación (primer modo)."""
        return tiempo_estabilizacion(self.Lx, self.Ly, self.alpha, self.Ti, self.Tb, self.tolerancia)

    @property
    def t_estabilizacion(self):
        """Instante en que máx |T - Tb| baja de la tolerancia (None si no llega).

        Se interpola entre los dos cuadros guardados que lo encierran, en escala
        logarítmica (la desviación decae casi exponencialmente).
        """
        desvio = np.max(np.abs(self.T - self.Tb), axis=(1, 2))
        dentro = desvio < self.tolerancia
        if not dentro.any():
            return None
        i = int(np.argmax(dentro))
        if i == 0:
            return 0.0
        a, b = np.log(desvio[i - 1]), np.log(max(desvio[i], 1e-300))
        return float(self.t[i - 1] + (a - log(self.tolerancia)) / (a - b) * (self.t[i] - self.t[i - 1]))

    @property
    def Tmax_t(self) -> np.ndarray:
        return self.T.max(axis=(1, 2))

    @property
    def dT_t(self) -> np.ndarray:
        """ΔT = Tmax − Tmin a `profundidad` del borde, por cuadro (criterio del curso)."""
        return self.Tmax_t - self.Tmin_superficie

    def posicion_Tmax(self) -> tuple[int, float, float]:
        """(cuadro, x, y) de la temperatura máxima en todo el tiempo, en coordenadas del usuario.

        Si varias celdas empatan (simetría), se informa el promedio de sus centros.
        """
        n = int(np.argmax(self.Tmax_t))
        i, j = np.nonzero(self.T[n] >= self.T[n].max() - 1e-9)
        x, y = (i.mean() + 0.5) * self.dx, (j.mean() + 0.5) * self.dx
        if self.origen == "centro":
            x, y = x - self.Lx / 2, y - self.Ly / 2
        return int(n), float(x), float(y)

    @property
    def extension(self):
        """Límites (x0, x1, y0, y1) en las coordenadas del usuario."""
        if self.origen == "centro":
            return (-self.Lx / 2, self.Lx / 2, -self.Ly / 2, self.Ly / 2)
        return (0.0, self.Lx, 0.0, self.Ly)


def simular_rectangulo(Lx=4.0, Ly=2.0, dx=0.1, k=2.222, rho=2400.0, cp=970.0, Ti=20.0, Tb=10.0,
                       puntos=((0.0, 0.0), (1.0, 0.5)), origen="centro", t_final=None,
                       tolerancia=0.1, n_cuadros=150, factor_seguridad=0.9,
                       hidratacion=False, dt_max=None, profundidad=0.025) -> ResultadoRectangulo:
    """Simula el ejercicio (con o sin calor de hidratación).

    hidratacion: False (sin calor), True (Ju y Lei 2019) o un modelo ya armado,
    por ejemplo CalorPulsoCoseno(Q_max, T_h).

    Por defecto el material es el hormigón de Ju y Lei (2019). Si t_final es
    None, se simula 1,2 veces el tiempo de estabilización (analítico sin
    calor; buscado simulando con calor). Con calor dt se limita a
    dt_max = 150 s por precisión (estudio de Δt del modelo 1D).
    """
    nx, ny = round(Lx / dx), round(Ly / dx)
    if nx < 2 or ny < 2 or abs(nx * dx - Lx) > 1e-9 * Lx or abs(ny * dx - Ly) > 1e-9 * Ly:
        raise ValueError(f"dx = {dx} m debe dividir exacto a Lx = {Lx} y Ly = {Ly} (celdas cuadradas)")
    if tolerancia <= 0 or n_cuadros < 1:
        raise ValueError("tolerancia y n_cuadros deben ser positivos")
    puntos = [tuple(float(c) for c in p) for p in puntos]
    esquina = [_a_esquina(p, Lx, Ly, origen) for p in puntos]  # también valida que estén dentro

    alpha = k / (rho * cp)
    if hidratacion is True:
        modelo = HidratacionJuLei2019()
    elif hidratacion is False or hidratacion is None:
        modelo = None
    else:
        modelo = hidratacion  # un modelo ya armado (p. ej. CalorPulsoCoseno)
    if modelo is not None and dt_max is None:
        dt_max = 150.0
    T_inicial = np.full((nx, ny), float(Ti))
    bordes = {c: ("dirichlet", Tb) for c in CARAS}

    if t_final is None:
        if modelo is None:
            t_final = 1.2 * tiempo_estabilizacion(Lx, Ly, alpha, Ti, Tb, tolerancia)
        else:
            t_final = 1.2 * _estabilizacion_simulada(T_inicial, alpha, dx, rho, cp, bordes, modelo, Tb, tolerancia)
        if t_final <= 0:
            raise ValueError("Con esa tolerancia el rectángulo ya está estable en t = 0: indique --t")

    # Se guarda al menos cada 1 h (como el 1D) para no saltarse el peak de Tmax;
    # la animación toma n_cuadros de estos instantes.
    n_guardados = max(n_cuadros, int(np.ceil(t_final / 3600.0)))
    s = simular_serie(T_inicial, alpha, (dx, dx), t_final, t_final / n_guardados, rho=rho, cp=cp, bordes=bordes,
                      factor_seguridad=factor_seguridad, dt_max=dt_max, hidratacion=modelo)

    P = _completar(s.T, Tb)  # (n_cuadros, (nx+2)(ny+2))
    T_puntos = P @ _pesos_bilineales(esquina, Lx, Ly, nx, ny).T
    Tmin_sup = (P @ _pesos_bilineales(contorno_interior(Lx, Ly, nx, ny, profundidad), Lx, Ly, nx, ny).T).min(axis=1)

    # Solución exacta: sin calor, o con el pulso coseno (lineal). Con Ju y Lei no hay.
    T_exacta_pts = error_campo = None
    pulso = modelo if isinstance(modelo, CalorPulsoCoseno) else None
    if modelo is None or pulso is not None:
        xc = (np.arange(nx) + 0.5) * dx  # centros, desde la esquina
        yc = (np.arange(ny) + 0.5) * dx

        def exacta(xs, ys, tk):
            return solucion_exacta(xs, ys, tk, Lx, Ly, alpha, Ti, Tb, pulso, rho, cp)

        T_exacta_pts = np.array([[exacta(x, y, tk)[0, 0] for x, y in esquina] for tk in s.t])
        error_campo = np.array([0.0] + [float(np.max(np.abs(Tk - exacta(xc, yc, tk))))
                                        for tk, Tk in zip(s.t[1:], s.T[1:])])

    return ResultadoRectangulo(
        Lx=Lx, Ly=Ly, dx=dx, k=k, rho=rho, cp=cp, Ti=Ti, Tb=Tb, origen=origen, puntos=puntos, modelo=modelo,
        t=s.t, T=s.T, te=s.te, T_puntos=T_puntos, T_puntos_exacta=T_exacta_pts, error_campo=error_campo,
        Tmin_superficie=Tmin_sup, profundidad=profundidad,
        dt=s.dt, dt_critico=paso_critico_malla(T_inicial, alpha, (dx, dx), rho, cp, bordes),
        dt_critico_interior=paso_critico(alpha, (dx, dx)), dt_max=dt_max, tolerancia=tolerancia,
        n_cuadros=n_cuadros,
    )


def unidad_tiempo(t_max: float) -> tuple[float, str]:
    """(segundos por unidad, nombre): horas hasta 3 días, si no días."""
    return (3600.0, "h") if t_max <= 3 * 86400 else (86400.0, "días")
