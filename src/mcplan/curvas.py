"""Curvas de hidratación del curso: fuente volumétrica Q(t) medida en probetas (E13/F02).

El archivo data/hidratacion/source_terms.txt (entregado por el profesor) trae 8 curvas
Q(t) [W/m³] de cilindros de hormigón (r = 3,75 cm, h = 15 cm): probetas M1-M4 en dos
campañas, SS (13-20 abr 2025, 240,5 h) y 22A (22 abr 2025, 163,1 h). Fueron recuperadas
con una PINN a la conductancia β de menor RMSE de cada probeta. Ya son fuente de calor:
no hay que derivar una curva adiabática ni hacer el balance del ensayo.

Cómo se usan como fuente por edad, q_h(a) [W/m³], con a = t - t_activación:
- **Inicio:** el tiempo 0 del archivo es "inicio del registro + 2 h", así que
  a = t_archivo + 2 h. En 0 <= a < 2 h no hay datos: rampa lineal de 0 a Q(2 h)
  (energía < 1 MJ/m³, declarada).
- **Suavizado:** promedio por bloques de `paso_h` (10 min). Conserva la energía del
  tramo y quita el ruido de muestreo de 60 s (σ ≈ 0,1-0,3 W/m³).
- **Cola:** después del registro, decaimiento exponencial desde el último valor con τ
  ajustado a las últimas 48 h (acotado a 6-240 h). Alternativas para sensibilidad:
  "cero" o "constante".
- **Escala:** `escala` multiplica Q (1 = sin cambio). Sólo sirve para una sensibilidad
  declarada, por ejemplo si el profesor indica otra energía total.
- **Temperatura:** no se aplica madurez ni Arrhenius. La fuente queda prescrita por edad
  (limitación declarada: la cinética real se acelera en el núcleo caliente).
"""

import hashlib
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

ARCHIVO = Path("data/hidratacion/source_terms.txt")


def leer_source_terms(ruta=ARCHIVO) -> dict:
    """Curvas crudas: nombre -> {t_h, Q, beta, probeta, campaña}. Verifica el formato del encabezado."""
    ruta = Path(ruta)
    texto = ruta.read_text(encoding="utf-8")
    if "W/m^3" not in texto or "time_h" not in texto:
        raise ValueError(f"{ruta}: no tiene el formato esperado (columnas time_h y Q en W/m^3)")
    betas = dict(re.findall(r"#\s+\d+\s+(\w+)_W_m3 \[W/m\^3\]\s+beta = ([\d.]+)", texto))
    encabezado = next(l for l in texto.splitlines() if l.startswith("# time_h"))
    nombres = encabezado.lstrip("#").split()[1:]
    datos = np.loadtxt(ruta, comments="#")
    curvas = {}
    for i, n in enumerate(nombres, start=1):
        ok = np.isfinite(datos[:, i])
        probeta, campana = n.split("_")
        curvas[n] = {"t_h": datos[ok, 0], "Q": datos[ok, i], "beta": float(betas.get(n, "nan")),
                     "probeta": probeta, "campana": campana}
    return curvas


def huella(ruta=ARCHIVO) -> str:
    return hashlib.sha256(Path(ruta).read_bytes()).hexdigest()


def _energia(t_h, Q):
    """Energía acumulada [J/m³] por trapecios."""
    return np.concatenate([[0.0], np.cumsum((Q[1:] + Q[:-1]) / 2 * np.diff(t_h) * 3600.0)])


@dataclass
class FuenteTabulada:
    """Fuente de hidratación prescrita por edad desde una curva medida. Interfaz igual a HidratacionJuLei2019."""

    nombre: str
    edad_h: np.ndarray  # grilla uniforme de edades [h], desde 0
    Q_tabla: np.ndarray  # W/m³ en esa grilla
    info: dict = field(default_factory=dict)

    @property
    def paso_h(self) -> float:
        return float(self.edad_h[1] - self.edad_h[0])

    def Q(self, a_h):
        """Q(a) [W/m³] por interpolación lineal en la grilla uniforme (rápida: índice directo)."""
        a = np.maximum(np.asarray(a_h, dtype=float), 0.0) / self.paso_h
        i = np.minimum(a.astype(np.int64), len(self.Q_tabla) - 2)
        f = np.minimum(a - i, 1.0)
        return self.Q_tabla[i] * (1 - f) + self.Q_tabla[i + 1] * f

    def energia(self, a_h):
        """Calor liberado hasta la edad a [J/m³]."""
        E = _energia(self.edad_h, self.Q_tabla)
        return np.interp(a_h, self.edad_h, E)

    def factor_arrhenius(self, T_C):
        """Sin madurez: la edad avanza con el tiempo real (fuente prescrita por edad)."""
        return np.ones_like(np.asarray(T_C, dtype=float))

    def potencia(self, te_h, T_C, rho, cp):
        """Q [W/m³]; rho, cp y T no se usan (la curva ya está en W/m³ y no depende de T)."""
        return self.Q(te_h)


def fuente_desde_curva(nombre, curva, edad_max_h=600.0, paso_h=1 / 6, desfase_h=2.0, cola="exponencial",
                       ventana_cola_h=48.0, escala=1.0) -> FuenteTabulada:
    """Convierte una curva cruda en una fuente por edad (ver el docstring del módulo)."""
    t, Q = curva["t_h"] + desfase_h, curva["Q"]
    edad = np.arange(0.0, edad_max_h + paso_h / 2, paso_h)
    # Promedio por bloques de paso_h centrado en cada edad (integral del tramo / paso): conserva la energía.
    E = _energia(t, Q)
    bordes = np.clip(edad + paso_h / 2, t[0], t[-1])
    E_b = np.interp(bordes, t, E)
    tabla = np.empty_like(edad)
    tabla[1:] = np.diff(E_b) / (paso_h * 3600.0)
    tabla[0] = 0.0
    dentro = (edad - paso_h / 2 >= t[0]) & (edad + paso_h / 2 <= t[-1])
    tabla = np.where(dentro, tabla, np.nan)
    # Antes del registro: rampa 0 -> Q(desfase).
    antes = edad < t[0] + paso_h / 2
    tabla[antes] = Q[0] * edad[antes] / t[0]
    # Después del registro: cola.
    despues = ~antes & ~dentro
    tau = None
    if cola == "exponencial":
        m = t >= t[-1] - ventana_cola_h
        pendiente = np.polyfit(t[m], np.log(np.maximum(Q[m], 1e-6)), 1)[0]
        tau = float(np.clip(-1.0 / pendiente if pendiente < 0 else np.inf, 6.0, 240.0))
        tabla[despues] = Q[-1] * np.exp(-(edad[despues] - t[-1]) / tau)
    elif cola == "cero":
        tabla[despues] = 0.0
    elif cola == "constante":
        tabla[despues] = Q[-1]
    else:
        raise ValueError(f"Cola desconocida: {cola}")
    tabla *= escala  # 1 = la curva tal como viene; otro valor sólo como sensibilidad declarada
    info = {"desfase_h": desfase_h, "cola": cola, "escala": escala, "tau_cola_h": tau, "fin_registro_h": float(t[-1]),
            "Q_max": float(Q.max()), "a_Qmax_h": float(t[np.argmax(Q)]), "beta": curva["beta"],
            "probeta": curva["probeta"], "campana": curva["campana"],
            "energia_registro_MJ_m3": float(E[-1] / 1e6)}
    return FuenteTabulada(nombre, edad, tabla, info)


def cargar_fuentes(ruta=ARCHIVO, nombres=None, **opciones) -> dict:
    curvas = leer_source_terms(ruta)
    nombres = nombres or list(curvas)
    return {n: fuente_desde_curva(n, curvas[n], **opciones) for n in nombres}


def tabla_curvas(fuentes: dict, rho_cp: float, edades=(24, 72, 168, 336)) -> list[dict]:
    """Resumen por curva: pico, energía a varias edades y aumento adiabático equivalente E/(ρ cp)."""
    filas = []
    for n, f in fuentes.items():
        fila = {"curva": n, **{k: f.info[k] for k in ("probeta", "campana", "beta", "Q_max", "a_Qmax_h",
                                                       "fin_registro_h", "tau_cola_h")}}
        for a in edades:
            E = float(f.energia(a))
            fila[f"E{a}_MJ_m3"] = E / 1e6
            fila[f"dTad{a}_C"] = E / rho_cp
        filas.append(fila)
    return filas
