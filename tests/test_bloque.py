"""Ejercicio de la pizarra: bloque 3 x 3 x 0,6 m con el pulso de calor y cuatro condiciones de borde."""

from math import inf, pi

import numpy as np
import pytest

from mcplan.bloque import (
    CasoBorde,
    Interpolador,
    SolucionExacta,
    autovalores_simetricos,
    casos_pizarra,
    coeficientes_uno,
    completar,
    puntos_cascara,
    simular_bloque,
)
from mcplan.cli import main
from mcplan.hidratacion import CalorPulsoCoseno

K, RHO, CP = 2.222, 2400.0, 970.0
ALPHA = K / (RHO * CP)
L = (3.0, 3.0, 0.6)
PULSO = CalorPulsoCoseno()
D20, D30, NEUMANN, ROBIN = casos_pizarra(K)


# --- Solución exacta ---------------------------------------------------------

def test_autovalores_de_los_tres_bordes():
    assert autovalores_simetricos(inf, 3) == pytest.approx([pi / 2, 3 * pi / 2, 5 * pi / 2])
    assert autovalores_simetricos(0.0, 3) == pytest.approx([0.0, pi, 2 * pi])
    Bi = ROBIN.h * 0.3 / K
    mu = autovalores_simetricos(Bi, 50)
    assert np.allclose(mu * np.tan(mu), Bi, atol=1e-9)
    assert np.all((mu > np.arange(50) * pi) & (mu < np.arange(50) * pi + pi / 2))  # una raíz por intervalo


@pytest.mark.parametrize("Bi", [inf, 3.75, 0.0])
def test_los_modos_reconstruyen_la_temperatura_uniforme(Bi):
    mu = autovalores_simetricos(Bi, 4000) if Bi else np.array([0.0])
    b = coeficientes_uno(mu)
    for s in (0.0, 0.5, -0.8):
        assert np.sum(b * np.cos(mu * s)) == pytest.approx(1.0, abs=1e-3)


def test_exacta_robin_cumple_la_condicion_de_borde():
    """∂T/∂n = −β (T − T_aire) en las caras: el signo correcto (sale calor si T > T_aire)."""
    ex = SolucionExacta(L, ROBIN, K, RHO, CP, 20.0, PULSO)
    beta = ROBIN.h / K
    e = 1e-6
    for t in (10 * 3600.0, 30 * 3600.0, 80 * 3600.0):
        for cara, adentro in (((0, 0, 0.3), (0, 0, 0.3 - e)), ((1.5, 0.4, 0.1), (1.5 - e, 0.4, 0.1))):
            T_cara, T_adentro = ex.puntos(t, [cara, adentro])
            assert T_cara > 20.0
            assert (T_cara - T_adentro) / e == pytest.approx(-beta * (T_cara - 20.0), rel=2e-3)


def test_exacta_dirichlet_vale_la_temperatura_de_la_cara():
    ex = SolucionExacta(L, D30, K, RHO, CP, 20.0, PULSO)
    assert ex.puntos(15 * 3600.0, [(0, 0, 0.3), (-1.5, 0.2, 0.0)]) == pytest.approx([30.0, 30.0], abs=1e-9)


def test_robin_con_h_enorme_es_dirichlet():
    casi_dirichlet = CasoBorde("robin", 20.0, 1e9)
    for t in (12 * 3600.0, 30 * 3600.0, 60 * 3600.0):
        a = SolucionExacta(L, casi_dirichlet, K, RHO, CP, 20.0, PULSO).puntos(t, [(0, 0, 0), (1, 0.5, 0.2)])
        b = SolucionExacta(L, D20, K, RHO, CP, 20.0, PULSO).puntos(t, [(0, 0, 0), (1, 0.5, 0.2)])
        assert a == pytest.approx(b, abs=1e-3)


def test_exacta_neumann_es_la_elevacion_adiabatica():
    ex = SolucionExacta(L, NEUMANN, K, RHO, CP, 20.0, PULSO)
    for t in (5 * 3600.0, 21 * 3600.0, 100 * 3600.0):
        assert ex.puntos(t, [(0.3, -1.0, 0.1)])[0] == pytest.approx(20.0 + PULSO.energia(t / 3600) / (RHO * CP))
    assert ex.T_equilibrio == pytest.approx(20.0 + 1800.0 * 42 * 3600 / (2 * RHO * CP))


# --- Herramientas numéricas ----------------------------------------------------

def test_temperatura_de_cara_robin_conserva_el_flujo():
    d, T_celda = 0.05, np.array([45.0])
    T_cara = ROBIN.temperatura_cara(T_celda, K, d)
    assert 2 * K / d * (T_celda - T_cara) == pytest.approx(ROBIN.h * (T_cara - 20.0))  # lo que llega = lo que sale


def test_interpolador_es_exacto_en_un_campo_lineal_y_lee_la_cara():
    n = (6, 6, 4)
    paso = [Li / ni for Li, ni in zip(L, n)]
    X, Y, Z = np.meshgrid(*[(np.arange(ni) + 0.5) * p for ni, p in zip(n, paso)], indexing="ij")
    T = 10 + 2 * X - 3 * Y + 5 * Z
    P = completar(T, D30, K, paso[0])
    interior = [(0.9, 1.1, 0.3), (2.0, 0.7, 0.25)]
    assert Interpolador(interior, L, n)(P) == pytest.approx([10 + 2 * x - 3 * y + 5 * z for x, y, z in interior])
    assert Interpolador([(1.0, 1.0, 0.0), (3.0, 2.0, 0.4)], L, n)(P) == pytest.approx([30.0, 30.0])  # en la cara


def test_puntos_cascara_a_25_mm_de_alguna_cara():
    pts = puntos_cascara(L, (60, 60, 12), 0.025)
    distancia = np.minimum(pts, np.asarray(L) - pts).min(axis=1)
    assert np.allclose(distancia, 0.025)
    assert np.all(pts >= 0.025 - 1e-12) and np.all(pts <= np.asarray(L) - 0.025 + 1e-12)


# --- Simulación contra la exacta ----------------------------------------------

def test_neumann_numerico_es_uniforme_y_sube_lo_del_ensayo():
    r = simular_bloque(NEUMANN, dx=0.2, t_final=60 * 3600.0, puntos=[(1.2, 1.2, 0.2)])
    assert np.allclose(r.corte_xy[-1], r.corte_xy[-1, 0, 0], atol=1e-9)  # adiabático: el campo sigue uniforme
    assert r.dT_t.max() < 1e-9
    assert r.T_centro[-1] == pytest.approx(r.T_eq, abs=1e-6)  # terminó el pulso: toda la energía quedó adentro
    # Durante el pulso, Euler explícito atrasa medio paso: error ~ (dt/2) Q/(rho cp)
    assert np.abs(r.T_centro - r.T_centro_exacta).max() < 0.6 * r.dt * PULSO.Q_max / (RHO * CP)


@pytest.mark.lento
@pytest.mark.parametrize("caso", [D20, D30, ROBIN], ids=lambda c: c.nombre)
def test_numerica_coincide_con_la_exacta_y_converge_con_orden_2(caso):
    """Con dx = 0,1 y 0,05 m (6 y 12 celdas en el espesor): el error baja ~4 veces, orden ~2.

    Con dx = 0,2 m (3 celdas en 0,6 m) la malla todavía no está en el régimen asintótico.
    """
    errores, centro = [], []
    for dx in (0.1, 0.05):
        r = simular_bloque(caso, dx=dx, t_final=36 * 3600.0, puntos=[(0, 0, 0.275)], n_modos=(60, 60, 40))
        desde_6h = r.t >= 6 * 3600.0  # deja pasar el salto inicial de Dirichlet 30 °C en el borde
        errores.append(r.error_campo[desde_6h].max())
        centro.append(np.abs(r.T_centro - r.T_centro_exacta)[desde_6h].max())
    assert centro[1] < 0.15  # °C, con peaks de 47 a 56 °C
    assert np.log2(errores[0] / errores[1]) > 1.7
    assert np.log2(centro[0] / centro[1]) > 1.7


def test_paso_critico_segun_el_borde():
    """Dirichlet: manda la esquina (dx²/9α). Neumann y Robin con Biot de celda chico: el interior (dx²/6α)."""
    dx = 0.1
    t = 2 * 3600.0
    assert simular_bloque(D20, dx=dx, t_final=t, con_exacta=False).dt_critico == pytest.approx(dx**2 / (9 * ALPHA))
    assert simular_bloque(NEUMANN, dx=dx, t_final=t, con_exacta=False).dt_critico == pytest.approx(dx**2 / (6 * ALPHA))
    assert simular_bloque(ROBIN, dx=dx, t_final=t, con_exacta=False).dt_critico == pytest.approx(dx**2 / (6 * ALPHA))


def test_simetria_y_tmax_en_el_centro():
    r = simular_bloque(ROBIN, dx=0.2, t_final=40 * 3600.0, con_exacta=False)
    i = int(np.argmax(r.Tmax_t))
    xy, xz = r.corte_xy[i], r.corte_xz[i]
    assert np.allclose(xy, xy[::-1, :]) and np.allclose(xy, xy[:, ::-1]) and np.allclose(xy, xy.T)
    assert np.allclose(xz, xz[:, ::-1])
    assert r.posicion_Tmax == pytest.approx((0.0, 0.0, 0.0))


def test_estabilizacion_simulada_igual_a_la_exacta():
    r = simular_bloque(D20, dx=0.1, con_exacta=False)  # t_final = 1,2 × la estabilización exacta
    assert r.t_estabilizacion == pytest.approx(r.t_estabilizacion_exacta, rel=0.03)
    assert r.t_estabilizacion > PULSO.T_h * 3600  # parte en su equilibrio, sale con el pulso y vuelve


@pytest.mark.lento
def test_comando_bloque(tmp_path, capsys):
    codigo = main(["bloque", "--dx", "0.2", "--t", "43200", "--cuadros", "4", "--salida", str(tmp_path)])
    assert codigo == 0
    salida = capsys.readouterr().out
    for texto in ("Dirichlet T = 20", "Dirichlet T = 30", "Neumann", "Robin", "manda la esquina", "manda el interior"):
        assert texto in salida
    for nombre in ("comparacion.png", "evolucion.gif", "resumen.json", "puntos_robin.png", "cortes_dirichlet30.png",
                   "serie_neumann.csv"):
        assert (tmp_path / nombre).stat().st_size > 0
