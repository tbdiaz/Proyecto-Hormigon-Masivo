"""Ejercicio de la pizarra: bloque 3D con el pulso de calor y cuatro condiciones de borde.

    rho cp dT/dt - k lap(T) = Q(t)     en un bloque Lx x Ly x Lz (3 x 3 x 0,6 m)
    T(x, y, z, 0) = Ti (20 °C)
    Q(t) = Q_max/2 (1 - cos(2 pi t/T)) si t < T, 0 después (pulso de la clase)

Casos, con la misma condición en las seis caras:
    Dirichlet   T = 20 °C  y  T = 30 °C
    Neumann     dT/dn = 0 (adiabático)
    Robin       -k dT/dn = h (T - T_aire), es decir dT/dn = -beta (T - T_aire), beta = h/k [1/m]

Solución exacta (los cuatro casos). Como Q sólo depende de t, el problema es
lineal. Con el origen en el centro, cada eje (medio lado l) tiene modos
simétricos cos(mu s / l), con mu tan(mu) = Bi = h l / k:
    Dirichlet: Bi = inf -> mu = (2k - 1) pi/2
    Neumann:   Bi = 0   -> sólo mu = 0 (el campo queda uniforme)
    Robin:     raíces de mu tan(mu) = Bi
La temperatura inicial y la fuente son uniformes, y el 1 se descompone como
    1 = sum_k b_k cos(mu_k s / l),   b_k = 2 sin(mu_k) / (mu_k + sin(mu_k) cos(mu_k)).
En 3D cada modo (i, j, m) decae con a = alpha (mu_i^2/lx^2 + mu_j^2/ly^2 + mu_m^2/lz^2):
    T = T_ref + sum b_i b_j b_m cos(.) cos(.) cos(.) [(Ti - T_ref) e^{-a t} + G(a, t) / (rho cp)]
con T_ref la temperatura del borde o del aire y G(a, t) = int_0^t e^{-a (t-s)} Q(s) ds,
que es cerrada para el pulso (CalorPulsoCoseno.integral_amortiguada).

Coordenadas: origen en el centro del bloque (x, y en [-1,5, 1,5], z en [-0,3, 0,3]);
con origen="esquina", x, y, z parten de 0 en una esquina de la base.
"""

import time
from dataclasses import dataclass
from math import ceil, inf, pi
from typing import Optional

import numpy as np

from .hidratacion import CalorPulsoCoseno
from .solver_nd import _pasos_para, paso_critico, paso_critico_malla, simular_difusion

EJES = ("x", "y", "z")
H_HORMIGON_AIRE = 27.8  # W/(m2 K), Ju y Lei (2019) Tabla 2: 100 kJ/(m2 h °C)


@dataclass
class CasoBorde:
    """Condición de borde, la misma en las seis caras del bloque."""

    tipo: str  # "dirichlet" | "neumann" | "robin"
    T_borde: float = 20.0  # Dirichlet: T de la cara; Robin: T del aire; Neumann: no se usa
    h: float = 0.0  # Robin: coeficiente de convección [W/(m2 K)]

    def __post_init__(self):
        if self.tipo not in ("dirichlet", "neumann", "robin"):
            raise ValueError(f"Tipo de borde desconocido: {self.tipo}")
        if self.tipo == "robin" and self.h <= 0:
            raise ValueError("Robin requiere h > 0")

    @property
    def nombre(self) -> str:
        return f"dirichlet{self.T_borde:g}" if self.tipo == "dirichlet" else self.tipo

    def etiqueta(self, k: float) -> str:
        if self.tipo == "dirichlet":
            return f"Dirichlet T = {self.T_borde:g} °C"
        if self.tipo == "neumann":
            return "Neumann ∂T/∂n = 0 (adiabático)"
        return f"Robin ∂T/∂n = −β(T − {self.T_borde:g}), β = {self.h / k:.3g} 1/m"

    def bordes(self) -> dict:
        """Condiciones en el formato de solver_nd, iguales en las seis caras."""
        condicion = {"dirichlet": ("dirichlet", self.T_borde), "neumann": ("neumann", 0.0),
                     "robin": ("robin", (self.h, self.T_borde))}[self.tipo]
        return {f"{e}_{lado}": condicion for e in EJES for lado in ("min", "max")}

    def Bi(self, l: float, k: float) -> float:
        """Número de Biot con el medio lado l: cuánto frena la superficie frente a la conducción."""
        return inf if self.tipo == "dirichlet" else 0.0 if self.tipo == "neumann" else self.h * l / k

    def temperatura_cara(self, T_celda, k: float, d: float):
        """T de la cara a partir de la celda vecina (a media celda), con el mismo flujo que usa el solver."""
        if self.tipo == "dirichlet":
            return np.full_like(T_celda, self.T_borde)
        if self.tipo == "neumann":
            return T_celda
        g = 2.0 * k / d  # conductancia de media celda
        return (g * T_celda + self.h * self.T_borde) / (g + self.h)


def casos_pizarra(k=2.222, beta=None, T_dirichlet=(20.0, 30.0), T_aire=20.0) -> list[CasoBorde]:
    """Los casos de la pizarra. Sin beta se usa h = 27,8 W/(m2 K) de Ju y Lei: beta = h/k."""
    h = H_HORMIGON_AIRE if beta is None else beta * k
    return ([CasoBorde("dirichlet", T) for T in T_dirichlet]
            + [CasoBorde("neumann"), CasoBorde("robin", T_aire, h)])


# ---------------------------------------------------------------------------
# Solución exacta
# ---------------------------------------------------------------------------

def autovalores_simetricos(Bi: float, n: int) -> np.ndarray:
    """Primeros n autovalores de la losa simétrica: mu tan(mu) = Bi (modos cos(mu s/l)).

    Cada raíz está en [(k-1) pi, (k-1) pi + pi/2], donde f = mu sin(mu) - Bi cos(mu)
    cambia de signo: se encuentran todas a la vez por bisección.
    """
    k = np.arange(1, n + 1)
    if Bi == inf:
        return (2 * k - 1) * pi / 2
    if Bi == 0:
        return (k - 1) * pi
    a, b = (k - 1) * pi, (k - 1) * pi + pi / 2

    def f(mu):
        return mu * np.sin(mu) - Bi * np.cos(mu)

    fa = f(a)
    for _ in range(100):
        m = 0.5 * (a + b)
        fm = f(m)
        mismo_signo = np.sign(fm) == np.sign(fa)
        a, fa, b = np.where(mismo_signo, m, a), np.where(mismo_signo, fm, fa), np.where(mismo_signo, b, m)
    return 0.5 * (a + b)


def coeficientes_uno(mu: np.ndarray) -> np.ndarray:
    """b_k tales que 1 = sum b_k cos(mu_k s/l) en [-l, l] (con b = 1 para mu = 0)."""
    mu = np.asarray(mu, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        b = 2.0 * np.sin(mu) / (mu + np.sin(mu) * np.cos(mu))
    return np.where(mu == 0, 1.0, b)


class SolucionExacta:
    """Solución exacta del bloque con el pulso, para un caso de borde (coordenadas desde el centro)."""

    def __init__(self, L, caso: CasoBorde, k, rho, cp, Ti, pulso: CalorPulsoCoseno, n_modos=(100, 100, 60)):
        self.l = tuple(Li / 2 for Li in L)
        self.caso, self.rho, self.cp, self.Ti, self.pulso = caso, rho, cp, Ti, pulso
        self.mu, self.b = [], []
        for li, n in zip(self.l, n_modos):
            Bi = caso.Bi(li, k)
            mu = np.array([0.0]) if Bi == 0 else autovalores_simetricos(Bi, n)
            self.mu.append(mu)
            self.b.append(coeficientes_uno(mu))
        alpha = k / (rho * cp)
        mx, my, mz = (m / li for m, li in zip(self.mu, self.l))
        self.a = alpha * (mx[:, None, None] ** 2 + my[None, :, None] ** 2 + mz[None, None, :] ** 2)
        self.B = self.b[0][:, None, None] * self.b[1][None, :, None] * self.b[2][None, None, :]
        # En Neumann el único modo es el uniforme y la referencia es la temperatura inicial.
        self.T_ref = Ti if caso.tipo == "neumann" else caso.T_borde

    def amplitudes(self, t) -> np.ndarray:
        """Amplitud de cada modo (i, j, m) en el instante t [s]: (Ti - T_ref) e^{-a t} + G(a, t)/(rho cp)."""
        return self.B * ((self.Ti - self.T_ref) * np.exp(-self.a * t)
                         + self.pulso.integral_amortiguada(self.a, t) / (self.rho * self.cp))

    def _cosenos(self, s, eje):
        return np.cos(np.outer(np.asarray(s, dtype=float), self.mu[eje]) / self.l[eje])

    def campo(self, t, xs, ys, zs, A=None) -> np.ndarray:
        """T en la grilla xs × ys × zs: arreglo (len(xs), len(ys), len(zs)). A: amplitudes ya calculadas."""
        A = self.amplitudes(t) if A is None else A
        R = np.tensordot(A, self._cosenos(zs, 2), axes=([2], [1]))  # (Kx, Ky, nz)
        R = np.tensordot(R, self._cosenos(ys, 1), axes=([1], [1]))  # (Kx, nz, ny)
        R = np.tensordot(R, self._cosenos(xs, 0), axes=([0], [1]))  # (nz, ny, nx)
        return self.T_ref + R.transpose(2, 1, 0)

    def puntos(self, t, P, A=None) -> np.ndarray:
        """T en los puntos P (n, 3). A: amplitudes ya calculadas."""
        A = self.amplitudes(t) if A is None else A
        P = np.atleast_2d(np.asarray(P, dtype=float))
        Fx, Fy, Fz = (self._cosenos(P[:, e], e) for e in range(3))
        return self.T_ref + np.einsum("ijk,pi,pj,pk->p", A, Fx, Fy, Fz, optimize=True)

    @property
    def T_equilibrio(self) -> float:
        """Estado final: la T del borde o del aire; adiabático, Ti + toda la energía del pulso."""
        if self.caso.tipo == "neumann":
            return self.Ti + float(self.pulso.energia(self.pulso.T_h)) / (self.rho * self.cp)
        return self.caso.T_borde

    def tiempo_estabilizacion(self, tolerancia: float, paso: float = 900.0) -> float:
        """Primer instante desde el cual |T - T_eq| < tolerancia en el centro (ahí manda el modo más lento).

        En el centro todos los cosenos valen 1: T = T_ref + suma de amplitudes.
        Durante el pulso se evalúa instante a instante; después, cada modo sólo
        decae desde su valor en t = T y se evalúan todos los instantes juntos.
        Se recorre cada 15 min hasta 25 constantes de tiempo del modo más lento
        y se interpola el cruce. Con 30 x 30 x 20 modos basta para el centro.
        """
        T_s = self.pulso.T_h * 3600.0
        a, B = self.a.ravel(), self.B.ravel()
        a_min = a[a > 0].min() if (a > 0).any() else None
        horizonte = T_s if a_min is None else T_s + 25.0 / a_min
        durante = np.arange(0.0, T_s, paso)
        despues = np.arange(T_s, horizonte + paso, paso)
        T_centro = [self.T_ref + self.amplitudes(tk).sum() for tk in durante]
        G_T = self.pulso.integral_amortiguada(a, T_s) / (self.rho * self.cp)
        for bloque in np.array_split(despues, max(1, len(despues) // 200)):
            decae = np.exp(-np.outer(a, bloque - T_s))
            T_centro.extend(self.T_ref + (B[:, None] * ((self.Ti - self.T_ref) * np.exp(-a * T_s)[:, None] * decae
                                                        + G_T[:, None] * decae)).sum(axis=0))
        t = np.concatenate((durante, despues))
        desvio = np.abs(np.array(T_centro) - self.T_equilibrio)
        fuera = np.nonzero(desvio >= tolerancia)[0]
        if len(fuera) == 0:
            return 0.0
        i = fuera[-1]
        if i + 1 >= len(t):
            return float(t[-1])
        return float(t[i] + (desvio[i] - tolerancia) / (desvio[i] - desvio[i + 1]) * (t[i + 1] - t[i]))


# ---------------------------------------------------------------------------
# Interpolación en el campo numérico
# ---------------------------------------------------------------------------

def completar(T, caso: CasoBorde, k: float, d: float) -> np.ndarray:
    """Agrega una capa con la temperatura de cada cara, según su condición de borde."""
    P = np.pad(np.asarray(T, dtype=float), 1, mode="edge")
    for eje in range(3):
        for lado, vecina in ((0, 1), (-1, -2)):
            cara = [slice(None)] * 3
            junto = [slice(None)] * 3
            cara[eje], junto[eje] = lado, vecina
            P[tuple(cara)] = caso.temperatura_cara(P[tuple(junto)], k, d)
    return P


class Interpolador:
    """Interpolación trilineal en el campo completado (coordenadas desde la esquina).

    Entre el último centro de celda y la cara se interpola contra la
    temperatura de la cara, así un punto en la superficie da la T superficial.
    """

    def __init__(self, puntos, L, n):
        P = np.atleast_2d(np.asarray(puntos, dtype=float))
        indices, pesos = [], []
        for e, (Li, ni) in enumerate(zip(L, n)):
            g = np.concatenate(([0.0], (np.arange(ni) + 0.5) * Li / ni, [Li]))
            i = np.clip(np.searchsorted(g, P[:, e], side="right") - 1, 0, ni)
            s = np.clip((P[:, e] - g[i]) / (g[i + 1] - g[i]), 0.0, 1.0)
            indices.append(i)
            pesos.append(s)
        forma = tuple(ni + 2 for ni in n)
        idx, w = [], []
        for di in (0, 1):
            for dj in (0, 1):
                for dk in (0, 1):
                    idx.append(np.ravel_multi_index((indices[0] + di, indices[1] + dj, indices[2] + dk), forma))
                    w.append((pesos[0] if di else 1 - pesos[0]) * (pesos[1] if dj else 1 - pesos[1])
                             * (pesos[2] if dk else 1 - pesos[2]))
        self.idx, self.w = np.stack(idx, axis=1), np.stack(w, axis=1)

    def __call__(self, P_completo) -> np.ndarray:
        return (P_completo.ravel()[self.idx] * self.w).sum(axis=1)


def puntos_cascara(L, n, profundidad) -> np.ndarray:
    """Puntos a `profundidad` de las seis caras (desde la esquina), donde se mide Tmin.

    Sobre esas superficies la interpolación es lineal por tramos entre centros
    de celda: basta con las coordenadas de los centros y los bordes.
    """
    if not 0 < profundidad < min(L) / 2:
        raise ValueError(f"profundidad = {profundidad} m debe estar entre 0 y la mitad del lado menor")
    ejes = []
    for Li, ni in zip(L, n):
        c = (np.arange(ni) + 0.5) * Li / ni
        ejes.append(np.concatenate(([profundidad], c[(c > profundidad) & (c < Li - profundidad)],
                                    [Li - profundidad])))
    bloques = []
    for e in range(3):
        o1, o2 = (o for o in range(3) if o != e)
        A, B = np.meshgrid(ejes[o1], ejes[o2], indexing="ij")
        for valor in (profundidad, L[e] - profundidad):
            p = np.empty((A.size, 3))
            p[:, e], p[:, o1], p[:, o2] = valor, A.ravel(), B.ravel()
            bloques.append(p)
    return np.vstack(bloques)


def _corte_central(T, eje):
    """Promedio de las dos capas centrales a lo largo de `eje` (o la central si son impares)."""
    n = T.shape[eje]
    capa = np.take(T, n // 2, axis=eje)
    return capa if n % 2 else 0.5 * (capa + np.take(T, n // 2 - 1, axis=eje))


# ---------------------------------------------------------------------------
# Simulación
# ---------------------------------------------------------------------------

@dataclass
class ResultadoBloque:
    caso: CasoBorde
    L: tuple
    dx: float
    k: float
    rho: float
    cp: float
    Ti: float
    pulso: CalorPulsoCoseno
    origen: str
    puntos: list  # coordenadas del usuario
    t: np.ndarray  # [s]
    T_centro: np.ndarray
    T_centro_exacta: Optional[np.ndarray]
    T_puntos: np.ndarray  # (n_t, n_puntos)
    T_puntos_exacta: Optional[np.ndarray]
    Tmax_t: np.ndarray  # máximo del campo en cada instante
    Tmin_cascara_t: np.ndarray  # mínimo a `profundidad` de las caras
    desvio_t: np.ndarray  # máx |T - T_eq| en el campo
    error_campo: Optional[np.ndarray]  # máx |T - exacta| en el campo
    corte_xy: np.ndarray  # (n_t, nx, ny), plano z = centro
    corte_xz: np.ndarray  # (n_t, nx, nz), plano y = centro
    T_eq: float
    Tmax: float
    t_Tmax: float
    posicion_Tmax: tuple
    dt: float
    dt_critico: float  # real, con la condición de borde
    dt_critico_interior: float  # dx^2 / (6 alpha)
    dt_max: Optional[float]
    tolerancia: float
    profundidad: float
    t_estabilizacion_exacta: Optional[float]
    computo_s: float

    @property
    def alpha(self) -> float:
        return self.k / (self.rho * self.cp)

    @property
    def dT_t(self) -> np.ndarray:
        return self.Tmax_t - self.Tmin_cascara_t

    @property
    def t_estabilizacion(self):
        """Instante desde el cual máx |T - T_eq| < tolerancia en todo el bloque, para siempre (interpolado).

        Se busca la última vez que estuvo fuera: un caso que parte en su
        equilibrio (Dirichlet 20 °C, Robin con aire a 20 °C) está "dentro" en
        t = 0 y después sale por el calor del pulso. None si no se estabiliza.
        """
        fuera = np.nonzero(self.desvio_t >= self.tolerancia)[0]
        if len(fuera) == 0:
            return 0.0
        i = int(fuera[-1])
        if i + 1 >= len(self.t):
            return None
        a, b = self.desvio_t[i], self.desvio_t[i + 1]
        return float(self.t[i] + (a - self.tolerancia) / (a - b) * (self.t[i + 1] - self.t[i]))

    def extension(self, ejes) -> tuple:
        """Límites para imshow de un corte en los ejes dados (coordenadas del usuario)."""
        lim = []
        for e in ejes:
            lim += [-self.L[e] / 2, self.L[e] / 2] if self.origen == "centro" else [0.0, self.L[e]]
        return tuple(lim)


def _a_esquina(punto, L, origen):
    if origen not in ("centro", "esquina"):
        raise ValueError(f"origen debe ser 'centro' o 'esquina', no {origen!r}")
    p = np.asarray(punto, dtype=float) + (np.asarray(L) / 2 if origen == "centro" else 0.0)
    if np.any(p < -1e-9) or np.any(p > np.asarray(L) + 1e-9):
        raise ValueError(f"El punto {tuple(punto)} está fuera del bloque")
    return np.clip(p, 0.0, L)


def simular_bloque(caso: CasoBorde, L=(3.0, 3.0, 0.6), dx=0.05, k=2.222, rho=2400.0, cp=970.0, Ti=20.0,
                   pulso: Optional[CalorPulsoCoseno] = None, puntos=((0.0, 0.0, 0.0),), origen="centro",
                   t_final=None, tolerancia=0.1, profundidad=0.025, dt_max=150.0, factor_seguridad=0.9,
                   n_modos=(100, 100, 60), con_exacta=True, observador=None) -> ResultadoBloque:
    """Simula un caso con el solver 3D y lo compara, instante a instante, con la solución exacta.

    Se guarda cada 1 h (o más seguido si t_final es corto): T en el centro y en
    los puntos, Tmax, Tmin a 25 mm de las caras, la desviación respecto del
    equilibrio, el error contra la exacta y los cortes centrales del campo.
    observador(t, P), si se indica, recibe en cada instante guardado el campo
    completado con las caras (para leer otras superficies, p. ej. la animación 3D).
    """
    inicio = time.perf_counter()
    pulso = CalorPulsoCoseno() if pulso is None else pulso
    L = tuple(float(Li) for Li in L)
    n = tuple(round(Li / dx) for Li in L)
    if min(n) < 2 or any(abs(ni * dx - Li) > 1e-9 * Li for ni, Li in zip(n, L)):
        raise ValueError(f"dx = {dx} m debe dividir exacto a los lados {L} (celdas cúbicas)")
    alpha = k / (rho * cp)
    esquina = [_a_esquina(p, L, origen) for p in puntos]
    puntos = [tuple(float(c) for c in p) for p in puntos]

    exacta = SolucionExacta(L, caso, k, rho, cp, Ti, pulso, n_modos) if con_exacta else None
    exacta_rapida = SolucionExacta(L, caso, k, rho, cp, Ti, pulso, (30, 30, 20))
    t_estab_exacta = exacta_rapida.tiempo_estabilizacion(tolerancia)
    if t_final is None:
        t_final = 1.2 * max(t_estab_exacta, pulso.T_h * 3600.0)

    T = np.full(n, float(Ti))
    bordes = caso.bordes()
    n_guardados = max(24, ceil(t_final / 3600.0 - 1e-9))
    guardar_cada = t_final / n_guardados
    pasos = _pasos_para(guardar_cada, T, alpha, (dx,) * 3, rho, cp, bordes, None, factor_seguridad, dt_max)
    dt = guardar_cada / pasos
    edad = np.zeros(n)

    centros = [(np.arange(ni) + 0.5) * dx - Li / 2 for ni, Li in zip(n, L)]  # desde el centro
    medio = np.asarray(L) / 2
    interp = Interpolador([medio] + esquina, L, n)
    interp_cascara = Interpolador(puntos_cascara(L, n, profundidad), L, n)
    P_exacta = np.array([np.zeros(3)] + [p - medio for p in esquina])
    T_eq = exacta_rapida.T_equilibrio

    reg = {c: [] for c in ("t", "Tpts", "Tpts_ex", "Tmax", "Tmin", "desvio", "error", "xy", "xz")}
    Tmax, t_Tmax, pos_Tmax = -inf, 0.0, (0.0, 0.0, 0.0)
    for i in range(n_guardados + 1):
        if i > 0:
            T = simular_difusion(T, alpha, (dx,) * 3, dt, pasos, rho=rho, cp=cp, bordes=bordes,
                                 hidratacion=pulso, edad=edad)
        t = i * guardar_cada
        reg["t"].append(t)
        P = completar(T, caso, k, dx)
        if observador is not None:
            observador(t, P)
        reg["Tpts"].append(interp(P))
        reg["Tmin"].append(interp_cascara(P).min())
        reg["Tmax"].append(T.max())
        reg["desvio"].append(np.abs(T - T_eq).max())
        reg["xy"].append(_corte_central(T, 2))
        reg["xz"].append(_corte_central(T, 1))
        if T.max() > Tmax + 1e-9:  # el 1e-9 evita que un redondeo mueva el instante de un máximo que ya no sube
            # Si varias celdas empatan (simetría), se informa el promedio de sus centros.
            idx = np.nonzero(T >= T.max() - 1e-9)
            Tmax, t_Tmax = float(T.max()), t
            pos_Tmax = tuple(round(float(np.mean(c[j])), 9) + 0.0 for c, j in zip(centros, idx))
        if exacta is not None:
            A = exacta.amplitudes(t)
            reg["Tpts_ex"].append(exacta.puntos(t, P_exacta, A))
            # En t = 0 la exacta es la condición inicial; con un salto en el borde (Dirichlet 30)
            # la malla no puede representarlo: el error se mide desde el primer instante guardado.
            reg["error"].append(0.0 if i == 0 else float(np.abs(T - exacta.campo(t, *centros, A=A)).max()))

    if origen == "esquina":
        pos_Tmax = tuple(round(p + m, 9) for p, m in zip(pos_Tmax, medio))
    Tpts = np.array(reg["Tpts"])
    Tpts_ex = np.array(reg["Tpts_ex"]) if exacta is not None else None
    return ResultadoBloque(
        caso=caso, L=L, dx=dx, k=k, rho=rho, cp=cp, Ti=Ti, pulso=pulso, origen=origen, puntos=puntos,
        t=np.array(reg["t"]), T_centro=Tpts[:, 0], T_centro_exacta=None if Tpts_ex is None else Tpts_ex[:, 0],
        T_puntos=Tpts[:, 1:], T_puntos_exacta=None if Tpts_ex is None else Tpts_ex[:, 1:],
        Tmax_t=np.array(reg["Tmax"]), Tmin_cascara_t=np.array(reg["Tmin"]), desvio_t=np.array(reg["desvio"]),
        error_campo=np.array(reg["error"]) if exacta is not None else None,
        corte_xy=np.array(reg["xy"]), corte_xz=np.array(reg["xz"]), T_eq=T_eq,
        Tmax=Tmax, t_Tmax=t_Tmax, posicion_Tmax=pos_Tmax,
        dt=dt, dt_critico=paso_critico_malla(np.full(n, Ti), alpha, (dx,) * 3, rho, cp, bordes),
        dt_critico_interior=paso_critico(alpha, (dx,) * 3), dt_max=dt_max, tolerancia=tolerancia,
        profundidad=profundidad, t_estabilizacion_exacta=t_estab_exacta,
        computo_s=time.perf_counter() - inicio,
    )


def t_final_comun(casos, L=(3.0, 3.0, 0.6), k=2.222, rho=2400.0, cp=970.0, Ti=20.0, pulso=None,
                  tolerancia=0.1) -> float:
    """1,2 veces el mayor tiempo de estabilización exacto, para comparar los casos en la misma ventana."""
    pulso = CalorPulsoCoseno() if pulso is None else pulso
    t = max(SolucionExacta(L, c, k, rho, cp, Ti, pulso, (30, 30, 20)).tiempo_estabilizacion(tolerancia)
            for c in casos)
    return 1.2 * max(t, pulso.T_h * 3600.0)
