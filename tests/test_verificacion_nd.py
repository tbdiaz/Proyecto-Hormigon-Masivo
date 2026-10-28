"""Batería E05 en 2D/3D contra soluciones exactas (guía de difusión, secciones 16 a 22).

Cada caso toma sus parámetros de una lista de parametrize: para probar otros
valores basta agregar o editar una fila (el dt se calcula solo bajo el paso
crítico real de la malla).
"""

import numpy as np
import pytest

from mcplan.solver_nd import paso_critico, paso_critico_malla, simular_hasta
from mcplan.verificacion_nd import (
    caso_energia_adiabatica,
    caso_fuente_uniforme,
    caso_sinusoidal,
    convergencia_sinusoidal,
)

# (L por eje [m], celdas por eje, alpha [m²/s], A [°C], T0 [°C], t [s])
CASOS_2D = [
    ((1.0, 1.0), 30, 1e-4, 50.0, 20.0, 1200.0),  # guía del curso
    ((2.0, 1.0), 30, 5e-5, 30.0, 15.0, 3000.0),  # rectángulo
    ((1.0, 1.0), 20, 1e-6, 30.0, 20.0, 1e5),  # difusividad del hormigón
]
CASOS_3D = [
    ((1.0, 1.0, 1.0), 16, 1e-4, 50.0, 20.0, 600.0),  # guía del curso
    ((1.0, 2.0, 1.5), 16, 2e-4, 30.0, 15.0, 900.0),  # caja
]


@pytest.mark.parametrize("L, n, alpha, A, T0, t", CASOS_2D + CASOS_3D)
def test_sinusoidal_exacta(L, n, alpha, A, T0, t):
    fila, T = caso_sinusoidal(L, n, alpha=alpha, A=A, T0=T0, t_final=t)
    assert fila.dt <= fila.dt_critico
    assert fila.error_max < 0.01 * A  # menos de 1 % de la amplitud
    assert T.min() > T0 - 1e-9  # principio del máximo: nada baja de T0 ni sube de T0 + A
    assert T.max() < T0 + A


@pytest.mark.parametrize("L, ns", [((1.0, 1.0), (10, 20, 40)), ((1.0, 1.0, 1.0), (6, 12, 24))])
def test_sinusoidal_converge_con_orden_2(L, ns):
    filas = convergencia_sinusoidal(L, ns, alpha=1e-4, A=50.0, T0=20.0, t_final=600.0)
    assert filas[0].error_max > filas[1].error_max > filas[2].error_max
    assert filas[-1].orden == pytest.approx(2.0, abs=0.15)


@pytest.mark.parametrize("dim", [2, 3])
@pytest.mark.parametrize("q0, rho, cp, T0, t", [
    (10000.0, 2400.0, 900.0, 20.0, 8 * 3600.0),  # guía del curso
    (2500.0, 2300.0, 1000.0, 10.0, 72 * 3600.0),
])
def test_fuente_uniforme_adiabatica_exacta(dim, q0, rho, cp, T0, t):
    r = caso_fuente_uniforme(dim, q0, rho, cp, T0, t)
    assert r.error_max < 1e-9
    assert r.dispersion < 1e-9  # el campo sigue uniforme: no hay flujos internos
    assert r.error_balance < 1e-12  # energía ganada = q0 * V * t


@pytest.mark.parametrize("dim", [2, 3])
@pytest.mark.parametrize("dos_materiales", [False, True])
def test_energia_constante_en_dominio_adiabatico(dim, dos_materiales):
    r = caso_energia_adiabatica(dim, dos_materiales=dos_materiales)
    assert r.Tmax_final < r.Tmax_inicial - 10  # el calor sí se redistribuyó
    assert r.error_relativo < 1e-12


def test_neumann_cero_es_adiabatico():
    T_inicial = np.full((8, 8), 20.0)
    T_inicial[2:4, 3:6] = 50.0
    caras = {f"{e}_{lado}": ("neumann", 0.0) for e in "xy" for lado in ("min", "max")}
    explicito, _ = simular_hasta(T_inicial, 1e-6, (0.1, 0.1), 1e5, bordes=caras)
    por_defecto, _ = simular_hasta(T_inicial, 1e-6, (0.1, 0.1), 1e5)
    assert np.array_equal(explicito, por_defecto)
    assert np.sum(explicito) == pytest.approx(np.sum(T_inicial), rel=1e-13)


@pytest.mark.parametrize("forma", [(10, 3), (10, 3, 3)])
def test_robin_estacionario_exacto(forma):
    """Dirichlet 50 °C en x = 0, convección (h = 10, 20 °C) en x = L, resto adiabático.

    Estacionario: un solo flujo q = (50 - 20) / (L/k + 1/h) cruza dos resistencias
    en serie, y T(x) = 50 - q x / k.
    """
    k, rho, cp, L, h = 2.0, 1000.0, 1000.0, 0.5, 10.0
    d = L / forma[0]
    alpha = k / (rho * cp)
    T, _ = simular_hasta(
        np.full(forma, 20.0), alpha, (d,) * len(forma), 15 * L**2 / alpha, rho=rho, cp=cp,
        bordes={"x_min": ("dirichlet", 50.0), "x_max": ("robin", (h, 20.0))},
    )
    q = (50.0 - 20.0) / (L / k + 1 / h)
    x = (np.arange(forma[0]) + 0.5) * d
    exacta = (50.0 - q * x / k).reshape((-1,) + (1,) * (len(forma) - 1))
    assert np.allclose(T, exacta, atol=1e-8)


@pytest.mark.parametrize("forma", [(10, 3), (10, 3, 3)])
def test_dos_materiales_en_serie_estacionario_exacto(forma):
    """Hormigón (k = 2) y suelo (k = 0.5) en serie entre 50 °C y 20 °C.

    Estacionario: el mismo flujo cruza ambos, q = (50 - 20) / (L1/k1 + L2/k2), y el
    perfil es lineal por tramos con más pendiente en el material que conduce menos.
    Sólo la media armónica de k en la interfaz lo reproduce exacto.
    """
    L1 = L2 = 0.25
    k1, k2 = 2.0, 0.5
    n = forma[0]
    d = (L1 + L2) / n
    x = (np.arange(n) + 0.5) * d
    k = np.where(x < L1, k1, k2).reshape((-1,) + (1,) * (len(forma) - 1)) * np.ones(forma)
    alpha_min = k2 / 1e6
    T, _ = simular_hasta(
        np.full(forma, 20.0), alpha_min, (d,) * len(forma), 3 * (L1 + L2) ** 2 / alpha_min,
        rho=1000.0, cp=1000.0, conductividad=k,
        bordes={"x_min": ("dirichlet", 50.0), "x_max": ("dirichlet", 20.0)},
    )
    q = (50.0 - 20.0) / (L1 / k1 + L2 / k2)
    exacta = np.where(x < L1, 50.0 - q * x / k1, 50.0 - q * L1 / k1 - q * (x - L1) / k2)
    assert np.allclose(T, exacta.reshape((-1,) + (1,) * (len(forma) - 1)), atol=1e-8)


def test_robin_el_calor_va_de_caliente_a_frio():
    caliente, _ = simular_hasta(np.full((5, 5, 5), 60.0), 1e-6, (0.1,) * 3, 1e4, bordes={"z_max": ("robin", (10.0, 20.0))})
    frio, _ = simular_hasta(np.full((5, 5, 5), 10.0), 1e-6, (0.1,) * 3, 1e4, bordes={"z_max": ("robin", (10.0, 20.0))})
    assert caliente.max() <= 60.0 and caliente.mean() < 60.0
    assert frio.min() >= 10.0 and frio.mean() > 10.0


@pytest.mark.parametrize("dim", [2, 3])
def test_paso_critico_con_y_sin_dirichlet(dim):
    """Sin bordes rige la celda interior; con Dirichlet, la celda de la esquina.

    Por eje, la celda junto a la cara Dirichlet suma 1 vecina + 2 (media celda
    hasta la cara) = 3 en vez de 2: el límite baja a 2/3 del clásico.
    """
    alpha, d = 1e-4, 0.05
    T = np.full((10,) * dim, 20.0)
    clasico = paso_critico(alpha, (d,) * dim)
    assert clasico == pytest.approx(d**2 / (2 * dim * alpha))
    assert paso_critico_malla(T, alpha, (d,) * dim) == pytest.approx(clasico)

    caras = {f"{e}_{lado}": ("dirichlet", 0.0) for e in "xyz"[:dim] for lado in ("min", "max")}
    assert paso_critico_malla(T, alpha, (d,) * dim, bordes=caras) == pytest.approx(d**2 / (3 * dim * alpha))


def test_dt_inestable_se_rechaza_y_forzado_diverge():
    fila, _ = caso_sinusoidal((1.0, 1.0), 20)
    with pytest.raises(ValueError, match="paso critico"):
        caso_sinusoidal((1.0, 1.0), 20, dt=2.0 * fila.dt_critico)
    divergido, T = caso_sinusoidal((1.0, 1.0), 20, dt=2.0 * fila.dt_critico, forzar=True)
    assert divergido.error_max == float("inf")
    assert np.max(np.abs(T)) > 1e4
