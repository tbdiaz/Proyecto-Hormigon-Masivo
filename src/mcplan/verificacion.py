"""Verificación del solver contra una solución analítica conocida.

Problema: losa de espesor L, temperatura inicial uniforme T0 y ambas caras
llevadas a Tb en t = 0 (Dirichlet). Sin fuente de calor. Solución en serie de
Fourier (separación de variables):

    T(z,t) = Tb + (T0 - Tb) * sum_{n impar} 4/(n pi) sin(n pi z / L) exp(-alpha (n pi / L)^2 t)

Si el código está bien, el error debe bajar al refinar la malla. Con dt
proporcional a dz^2 (paso automático) esperamos orden ~2: dividir dz por 2
divide el error por ~4.
"""

from dataclasses import dataclass
from math import log
from typing import Optional

import numpy as np

from .solver1d import Borde, Capa, Material, crear_malla, simular

HORMIGON_PRUEBA = Material("hormigon", k=2.222, rho=2400.0, cp=970.0, es_hormigon=True)


def solucion_analitica(z, t, L, alpha, T0, Tb, n_terminos=2001):
    z = np.asarray(z, dtype=float)
    suma = np.zeros_like(z)
    for n in range(1, n_terminos + 1, 2):
        suma += 4.0 / (n * np.pi) * np.sin(n * np.pi * z / L) * np.exp(-alpha * (n * np.pi / L) ** 2 * t)
    return Tb + (T0 - Tb) * suma


@dataclass
class FilaConvergencia:
    dz: float
    dt: float
    dt_critico: float
    n_pasos: int
    error_max: float
    orden: Optional[float]  # log(e_anterior/e) / log(dz_anterior/dz)


def error_contra_analitica(dz, dt=None, forzar=False, L=1.0, T0=50.0, Tb=20.0, Fo=0.05,
                           material=HORMIGON_PRUEBA):
    """Simula el problema de verificación y devuelve (resultado, error máximo [°C])."""
    t_final = Fo * L**2 / material.alpha
    malla = crear_malla([Capa(material, L, T0)], dz)
    borde = Borde.dirichlet(Tb)
    res = simular(malla, borde, borde, t_final=t_final, dt=dt, guardar_cada=t_final, forzar=forzar)
    exacta = solucion_analitica(malla.z, res.t[-1], L, material.alpha, T0, Tb)
    return res, float(np.max(np.abs(res.T[-1] - exacta)))


def estudio_convergencia(dzs, dt=None, forzar=False) -> list[FilaConvergencia]:
    filas = []
    for dz in sorted(dzs, reverse=True):
        res, err = error_contra_analitica(dz, dt=dt, forzar=forzar)
        orden = None
        if filas and err > 0 and filas[-1].error_max > 0:
            orden = log(filas[-1].error_max / err) / log(filas[-1].dz / dz)
        filas.append(FilaConvergencia(dz, res.dt, res.dt_critico, res.n_pasos, err, orden))
    return filas
