"""Verificación del solver 2D/3D contra soluciones exactas (batería E05).

Casos de la guía *Difusión de calor 1D/2D/3D* (secciones 16, 19, 21 y 22):

- Sinusoidal con todo el contorno a T0 (Dirichlet):

      T = T0 + A * prod(sin(pi x_i / L_i)) * exp(-alpha * lambda * t)
      lambda = sum((pi / L_i)^2)

- Fuente uniforme en dominio adiabático: dT/dt = q0 / (rho cp).
- Dominio adiabático sin fuente: E = sum(rho cp V T) constante.

Todos los parámetros son argumentos, así se pueden cambiar desde los tests o
con `mcplan sinusoidal | fuente | energia`. El dt se elige solo (bajo el paso
crítico real de la malla), salvo que se indique uno.
"""

from dataclasses import dataclass
from math import log, prod
from typing import Optional

import numpy as np

from .solver_nd import paso_critico_malla, simular_difusion, simular_hasta

EJES = ("x", "y", "z")

# Propiedades de Ju y Lei (2019), Tabla 2, en SI
HORMIGON = {"k": 2.222, "rho": 2400.0, "cp": 970.0}
SUELO = {"k": 1.389, "rho": 1800.0, "cp": 1273.0}


def centros(L, n) -> list[np.ndarray]:
    """Coordenadas de los centros de celda de una malla uniforme (indexing ij)."""
    return np.meshgrid(*[(np.arange(ni) + 0.5) * Li / ni for Li, ni in zip(L, n)], indexing="ij")


def solucion_sinusoidal(coords, L, t, alpha, A, T0):
    lam = sum((np.pi / Li) ** 2 for Li in L)
    forma = prod(np.sin(np.pi * c / Li) for c, Li in zip(coords, L))
    return T0 + A * forma * np.exp(-alpha * lam * t)


@dataclass
class FilaSinusoidal:
    n: int  # celdas por eje
    d: tuple[float, ...]
    dt: float
    dt_critico: float
    pasos: int
    error_max: float  # inf si divergió
    orden: Optional[float] = None  # log(e_anterior/e) / log(d_anterior/d)


def caso_sinusoidal(L, n, alpha=1e-4, A=50.0, T0=20.0, t_final=1200.0, dt=None, forzar=False,
                    factor_seguridad=0.9) -> tuple[FilaSinusoidal, np.ndarray]:
    """Caso sinusoidal 2D (L = (Lx, Ly)) o 3D (L = (Lx, Ly, Lz)) con n celdas por eje.

    Con dt = None el paso es automático; con un dt fijo se hacen
    round(t_final/dt) pasos y la solución exacta se evalúa en pasos*dt.
    """
    L = tuple(float(Li) for Li in L)
    d = tuple(Li / n for Li in L)
    coords = centros(L, (n,) * len(L))
    T_inicial = solucion_sinusoidal(coords, L, 0.0, alpha, A, T0)
    bordes = {f"{e}_{lado}": ("dirichlet", T0) for e in EJES[:len(L)] for lado in ("min", "max")}

    dt_critico = paso_critico_malla(T_inicial, alpha, d, bordes=bordes)
    if dt is None:
        T, dt = simular_hasta(T_inicial, alpha, d, t_final, bordes=bordes, factor_seguridad=factor_seguridad)
        pasos = round(t_final / dt)
    else:
        pasos = max(1, round(t_final / dt))
        T = simular_difusion(T_inicial, alpha, d, dt, pasos, bordes=bordes, forzar=forzar)

    exacta = solucion_sinusoidal(coords, L, pasos * dt, alpha, A, T0)
    error = np.max(np.abs(T - exacta))
    error = float(error) if np.isfinite(error) and np.max(np.abs(T)) < 1e4 else float("inf")
    return FilaSinusoidal(n, d, dt, dt_critico, pasos, error), T


def convergencia_sinusoidal(L, ns, **kwargs) -> list[FilaSinusoidal]:
    """Repite caso_sinusoidal refinando la malla. Con dt ∝ d², se espera orden ≈ 2."""
    filas = []
    for n in sorted(ns):
        fila, _ = caso_sinusoidal(L, n, **kwargs)
        if filas and 0 < fila.error_max < float("inf") and filas[-1].error_max > 0:
            fila.orden = log(filas[-1].error_max / fila.error_max) / log(filas[-1].d[0] / fila.d[0])
        filas.append(fila)
    return filas


@dataclass
class ResultadoFuente:
    T_media: float
    T_exacta: float
    error_max: float  # máx |T - T_exacta| en todas las celdas
    dispersion: float  # Tmax - Tmin: debe ser 0 (el campo sigue uniforme)
    energia_ganada: float  # E_final - E_inicial [J]
    energia_aportada: float  # q0 * V_total * t [J]
    dt: float
    pasos: int

    @property
    def error_balance(self) -> float:
        return abs(self.energia_ganada - self.energia_aportada) / self.energia_aportada


def caso_fuente_uniforme(dim=3, q0=10000.0, rho=2400.0, cp=900.0, T0=20.0, t_final=8 * 3600.0,
                         n=5, L=1.0, k=2.0) -> ResultadoFuente:
    """Bloque adiabático con fuente uniforme q0 [W/m³]: T(t) = T0 + q0 t / (rho cp)."""
    if q0 <= 0:
        raise ValueError("q0 debe ser positivo (sin fuente, ver caso_energia_adiabatica)")
    d = (L / n,) * dim
    T_inicial = np.full((n,) * dim, float(T0))
    T, dt = simular_hasta(T_inicial, k / (rho * cp), d, t_final, fuente=q0, rho=rho, cp=cp)

    V = prod(d)
    exacta = T0 + q0 * t_final / (rho * cp)
    return ResultadoFuente(
        T_media=float(T.mean()), T_exacta=exacta,
        error_max=float(np.max(np.abs(T - exacta))), dispersion=float(T.max() - T.min()),
        energia_ganada=float(rho * cp * V * np.sum(T - T_inicial)),
        energia_aportada=q0 * V * T.size * t_final,
        dt=dt, pasos=round(t_final / dt),
    )


@dataclass
class ResultadoEnergia:
    E_inicial: float  # J
    E_final: float
    Tmax_inicial: float
    Tmax_final: float
    dt: float
    pasos: int

    @property
    def error_relativo(self) -> float:
        return abs(self.E_final - self.E_inicial) / abs(self.E_inicial)


def caso_energia_adiabatica(dim=3, n=10, L=1.0, t_final=3 * 86400.0, dos_materiales=False) -> ResultadoEnergia:
    """Campo inicial con un punto caliente, todas las caras adiabáticas y sin fuente.

    Con dos_materiales=True la mitad x < L/2 es hormigón y el resto suelo:
    la interfaz no debe crear ni destruir energía.
    """
    forma = (n,) * dim
    d = (L / n,) * dim
    coords = centros((L,) * dim, forma)
    r2 = (coords[0] - 0.35 * L) ** 2 + sum((c - 0.5 * L) ** 2 for c in coords[1:])
    T_inicial = 20.0 + 30.0 * np.exp(-r2 / (0.15 * L) ** 2)

    k, rho, cp = (np.full(forma, HORMIGON[p]) for p in ("k", "rho", "cp"))
    if dos_materiales:
        suelo = coords[0] > 0.5 * L
        for arreglo, p in ((k, "k"), (rho, "rho"), (cp, "cp")):
            arreglo[suelo] = SUELO[p]

    alpha = HORMIGON["k"] / (HORMIGON["rho"] * HORMIGON["cp"])  # sólo se valida: manda conductividad
    T, dt = simular_hasta(T_inicial, alpha, d, t_final, rho=rho, cp=cp, conductividad=k)

    V = prod(d)
    return ResultadoEnergia(
        E_inicial=float(np.sum(rho * cp * V * T_inicial)), E_final=float(np.sum(rho * cp * V * T)),
        Tmax_inicial=float(T_inicial.max()), Tmax_final=float(T.max()),
        dt=dt, pasos=round(t_final / dt),
    )
