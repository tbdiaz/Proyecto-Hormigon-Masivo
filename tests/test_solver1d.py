"""Tests del solver 1D: cada uno compara contra un resultado que se conoce exacto."""

import numpy as np
import pytest

from mcplan.metricas import calcular_metricas
from mcplan.solver1d import Borde, Capa, Material, crear_malla, paso_critico, simular
from mcplan.verificacion import estudio_convergencia

HORMIGON = Material("hormigon", k=2.222, rho=2400.0, cp=970.0, es_hormigon=True)
SUELO = Material("suelo", k=1.389, rho=1800.0, cp=1273.0)
DIA = 86400.0


def test_malla_respeta_espesores_e_interfaces():
    malla = crear_malla([Capa(SUELO, 5.0, 15.0), Capa(HORMIGON, 3.3, 20.0)], 0.1)
    assert malla.dz.sum() == pytest.approx(8.3)
    assert np.any(np.isclose(malla.z_caras, 5.0))  # la interfaz cae sobre una cara
    assert malla.es_hormigon.sum() == 33


def test_paso_critico_equivale_a_fourier_un_medio():
    """Lejos de los bordes, dt_crit = dz^2 / (2 alpha), es decir Fo = 1/2."""
    dz = 0.1
    malla = crear_malla([Capa(HORMIGON, 1.0, 20.0)], dz)
    adiabatico = Borde.neumann(0.0)
    assert paso_critico(malla, adiabatico, adiabatico) == pytest.approx(dz**2 / (2 * HORMIGON.alpha))


def test_borde_dirichlet_es_mas_restrictivo():
    """La celda junto a un Dirichlet conduce por media celda (2k/dz): dt_crit = dz^2/(3 alpha)."""
    dz = 0.1
    malla = crear_malla([Capa(HORMIGON, 1.0, 20.0)], dz)
    b = Borde.dirichlet(20.0)
    assert paso_critico(malla, b, b) == pytest.approx(dz**2 / (3 * HORMIGON.alpha))


def test_estacionario_dos_materiales_es_lineal_por_tramos():
    """Con T fija en ambos extremos el flujo es igual en ambos materiales:
    q = (Tb - Ta) / (L1/k1 + L2/k2), y el perfil es lineal en cada material."""
    L1, L2, Ta, Tb = 0.5, 0.5, 10.0, 40.0
    malla = crear_malla([Capa(SUELO, L1, 25.0), Capa(HORMIGON, L2, 25.0)], 0.05)
    res = simular(malla, Borde.dirichlet(Ta), Borde.dirichlet(Tb), t_final=40 * DIA, guardar_cada=40 * DIA)
    q = (Tb - Ta) / (L1 / SUELO.k + L2 / HORMIGON.k)
    T_interfaz = Ta + q * L1 / SUELO.k
    z = malla.z
    exacta = np.where(z < L1, Ta + q * z / SUELO.k, T_interfaz + q * (z - L1) / HORMIGON.k)
    np.testing.assert_allclose(res.T[-1], exacta, atol=1e-6)


def test_estacionario_robin():
    """Base a 30 °C, superficie con convección a 0 °C: resistencias en serie L/k + 1/h."""
    h, T_inf, T_base, L = 10.0, 0.0, 30.0, 0.5
    malla = crear_malla([Capa(HORMIGON, L, 15.0)], 0.025)
    res = simular(malla, Borde.dirichlet(T_base), Borde.robin(h, T_inf), t_final=30 * DIA, guardar_cada=30 * DIA)
    q = (T_base - T_inf) / (L / HORMIGON.k + 1.0 / h)
    np.testing.assert_allclose(res.T[-1], T_base - q * malla.z / HORMIGON.k, atol=1e-6)
    assert res.T_cara_sup[-1] == pytest.approx(T_inf + q / h, abs=1e-6)


def test_conservacion_de_energia_con_bordes_adiabaticos():
    """Sin flujo por los bordes ni fuente, la energía total sum(C_i T_i) no cambia."""
    malla = crear_malla([Capa(SUELO, 1.0, 10.0), Capa(HORMIGON, 1.0, 60.0)], 0.05)
    adiabatico = Borde.neumann(0.0)
    res = simular(malla, adiabatico, adiabatico, t_final=5 * DIA)
    energia = (malla.capacidad * res.T).sum(axis=1)
    np.testing.assert_allclose(energia, energia[0], rtol=1e-12)


def test_flujo_neumann_entra_completo():
    """Todo el flujo impuesto queda almacenado: ganancia = q * t."""
    q = 50.0
    malla = crear_malla([Capa(HORMIGON, 1.0, 20.0)], 0.05)
    res = simular(malla, Borde.neumann(q), Borde.neumann(0.0), t_final=2 * DIA)
    ganancia = (malla.capacidad * (res.T[-1] - res.T[0])).sum()
    assert ganancia == pytest.approx(q * res.t[-1], rel=1e-10)


class FuenteUniforme:
    """Fuente constante q [W/m3] que no depende de la temperatura ni de la edad."""

    def __init__(self, q):
        self.q = q

    def potencia(self, te_h, T_C, rho, cp):
        return np.full_like(np.asarray(T_C, dtype=float), self.q)

    def factor_arrhenius(self, T_C):
        return np.ones_like(np.asarray(T_C, dtype=float))


def test_fuente_uniforme_adiabatica_sube_linealmente():
    """Guía G05, test 5: sin pérdidas y con fuente uniforme, dT/dt = q / (rho cp) en todo el dominio."""
    q = 500.0
    mat = Material("hormigon", k=2.222, rho=2400.0, cp=970.0, es_hormigon=True, hidratacion=FuenteUniforme(q))
    malla = crear_malla([Capa(mat, 1.0, 20.0)], 0.05)
    adiabatico = Borde.neumann(0.0)
    res = simular(malla, adiabatico, adiabatico, t_final=2 * DIA)
    esperado = 20.0 + q / (mat.rho * mat.cp) * res.t
    np.testing.assert_allclose(res.T, np.repeat(esperado[:, None], malla.n, axis=1), rtol=1e-12)


def test_converge_a_la_solucion_analitica_con_orden_dos():
    filas = estudio_convergencia([0.1, 0.05, 0.025])
    assert filas[-1].error_max < 5e-3
    for f in filas[1:]:
        assert 1.7 < f.orden < 2.5


def test_dt_estable_respeta_principio_del_maximo():
    """Sin fuente, la temperatura nunca sale del rango de los datos (no aparecen
    temperaturas mayores que la inicial ni menores que la del borde)."""
    malla = crear_malla([Capa(HORMIGON, 1.0, 50.0)], 0.05)
    b = Borde.dirichlet(20.0)
    res = simular(malla, b, b, t_final=2 * DIA)
    assert res.T.max() <= 50.0 + 1e-9
    assert res.T.min() >= 20.0 - 1e-9


def test_paso_automatico_respeta_dt_max_y_el_intervalo_de_guardado():
    malla = crear_malla([Capa(HORMIGON, 1.0, 20.0)], 0.2)  # malla gruesa: dt crítico > 1 h
    b = Borde.neumann(0.0)
    res = simular(malla, b, b, t_final=DIA, guardar_cada=3600.0, dt_max=600.0)
    assert res.dt <= 600.0
    assert 3600.0 / res.dt == pytest.approx(round(3600.0 / res.dt))


def test_dt_inestable_lanza_error():
    malla = crear_malla([Capa(HORMIGON, 1.0, 50.0)], 0.05)
    b = Borde.dirichlet(20.0)
    dt_crit = paso_critico(malla, b, b)
    with pytest.raises(ValueError, match="paso crítico"):
        simular(malla, b, b, t_final=DIA, dt=1.5 * dt_crit)


def test_dt_inestable_forzado_diverge():
    """Con 2x el paso crítico, las celdas interiores quedan en Fo = 0.67 > 0.5 y el error
    se amplifica en cada paso. (Con un exceso pequeño la solución puede sólo oscilar:
    el criterio C/sum(G) es suficiente, no el límite exacto de divergencia.)"""
    malla = crear_malla([Capa(HORMIGON, 1.0, 50.0)], 0.05)
    b = Borde.dirichlet(20.0)
    dt_crit = paso_critico(malla, b, b)
    res = simular(malla, b, b, t_final=10 * DIA, dt=2.0 * dt_crit, forzar=True)
    assert res.diverge


def test_temperatura_a_25mm_usa_la_superficie_reconstruida():
    """Con dz = 0.1 el punto a 25 mm queda entre la cara y el primer centro (a 50 mm).
    En estado estacionario el perfil es lineal, así que la interpolación debe ser exacta."""
    h, T_inf, T_base, L = 10.0, 0.0, 30.0, 0.5
    malla = crear_malla([Capa(HORMIGON, L, 15.0)], 0.1)
    res = simular(malla, Borde.dirichlet(T_base), Borde.robin(h, T_inf), t_final=30 * DIA, guardar_cada=30 * DIA)
    met = calcular_metricas(res)
    q = (T_base - T_inf) / (L / HORMIGON.k + 1.0 / h)
    assert met.T_25mm_sup_t[-1] == pytest.approx(T_base - q * (L - 0.025) / HORMIGON.k, abs=1e-6)
