"""Métricas térmicas según los criterios del curso (Introducción General, punto 3).

    Tmax(t)  = temperatura máxima dentro del hormigón en el instante t
    Tmin(t)  = mínima entre los puntos ubicados a 25 mm bajo las superficies
               exteriores del hormigón (en 1D: cara superior e inferior)
    dT(t)    = Tmax(t) - Tmin(t)

    Tmax = max_t Tmax(t)  <= 70 °C
    dTmax = max_t dT(t)   <= 20 °C
"""

from dataclasses import dataclass

import numpy as np

from .solver1d import Resultado


def temperatura_en(res: Resultado, z: float) -> np.ndarray:
    """Historia temporal de la temperatura en la cota z [m] (interpolación lineal).

    Se interpola sobre los centros de celda más las dos caras del dominio,
    cuyas temperaturas reconstruye el solver con la condición de borde. Así
    un punto a 25 mm de la superficie se puede evaluar aunque quede entre la
    cara y el primer centro de celda.
    """
    m = res.malla
    z_ext = np.concatenate(([m.z_caras[0]], m.z, [m.z_caras[-1]]))
    T_ext = np.column_stack([res.T_cara_inf, res.T, res.T_cara_sup])
    return np.array([np.interp(z, z_ext, fila) for fila in T_ext])


@dataclass
class MetricasTermicas:
    Tmax: float  # [°C]
    t_Tmax_h: float  # instante en que ocurre [h]
    z_Tmax: float  # altura donde ocurre, medida desde la base del hormigón [m]
    dTmax: float  # [°C]
    t_dTmax_h: float  # [h]
    limite_Tmax: float
    limite_dT: float
    # Series temporales (una por instante guardado)
    Tmax_t: np.ndarray
    Tmin_t: np.ndarray
    dT_t: np.ndarray
    T_25mm_sup_t: np.ndarray
    T_25mm_inf_t: np.ndarray
    z_25mm_sup: float
    z_25mm_inf: float

    @property
    def cumple_Tmax(self) -> bool:
        return self.Tmax <= self.limite_Tmax

    @property
    def cumple_dT(self) -> bool:
        return self.dTmax <= self.limite_dT

    def resumen(self) -> str:
        ok = lambda c: "CUMPLE" if c else "NO CUMPLE"  # noqa: E731
        return (
            f"Tmax  = {self.Tmax:6.2f} °C  (t = {self.t_Tmax_h:6.1f} h, a {self.z_Tmax:.2f} m sobre la base)"
            f"   límite {self.limite_Tmax:.0f} °C -> {ok(self.cumple_Tmax)}\n"
            f"dTmax = {self.dTmax:6.2f} °C  (t = {self.t_dTmax_h:6.1f} h)"
            f"                               límite {self.limite_dT:.0f} °C -> {ok(self.cumple_dT)}"
        )


def calcular_metricas(
    res: Resultado,
    profundidad: float = 0.025,
    limite_Tmax: float = 70.0,
    limite_dT: float = 20.0,
) -> MetricasTermicas:
    malla = res.malla
    es_h = malla.es_hormigon
    if not es_h.any():
        raise ValueError("El dominio no tiene celdas de hormigón")
    idx = np.flatnonzero(es_h)
    if np.any(np.diff(idx) != 1):
        raise ValueError("Se asume un único bloque continuo de hormigón en la columna")

    z_base = malla.z_caras[idx[0]]  # cara inferior del hormigón
    z_tope = malla.z_caras[idx[-1] + 1]  # cara superior del hormigón

    T_h = res.T[:, es_h]
    Tmax_t = T_h.max(axis=1)
    z_max_t = malla.z[es_h][T_h.argmax(axis=1)]

    z_sup, z_inf = z_tope - profundidad, z_base + profundidad
    T_sup = temperatura_en(res, z_sup)
    T_inf = temperatura_en(res, z_inf)
    Tmin_t = np.minimum(T_sup, T_inf)
    dT_t = Tmax_t - Tmin_t

    i_T = int(np.argmax(Tmax_t))
    i_dT = int(np.argmax(dT_t))
    return MetricasTermicas(
        Tmax=float(Tmax_t[i_T]),
        t_Tmax_h=float(res.t_h[i_T]),
        z_Tmax=float(z_max_t[i_T] - z_base),
        dTmax=float(dT_t[i_dT]),
        t_dTmax_h=float(res.t_h[i_dT]),
        limite_Tmax=limite_Tmax,
        limite_dT=limite_dT,
        Tmax_t=Tmax_t,
        Tmin_t=Tmin_t,
        dT_t=dT_t,
        T_25mm_sup_t=T_sup,
        T_25mm_inf_t=T_inf,
        z_25mm_sup=z_sup,
        z_25mm_inf=z_inf,
    )
