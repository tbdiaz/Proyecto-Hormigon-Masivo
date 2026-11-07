"""Condiciones de borde (E10): convección con viento, onda larga, sol, moldaje en serie y desmolde."""

import numpy as np
import pytest

from mcplan.bordes import (
    LAT,
    SIGMA,
    SUP,
    Ambiente,
    Capa,
    Clima,
    Cubierta,
    PlanCubiertas,
    T_cielo,
    coeficiente_U,
    h_conveccion,
)
from mcplan.hidratacion import HidratacionJuLei2019
from mcplan.simulacion3d import Material, construir_dominio, simular_3d
from mcplan.solver_nd import simular_difusion
from mcplan.voxelizacion import PlanVertidas

from test_simulacion3d import HORM, SUEL, _robin, _vox_cajas

TERCIADO = Capa("terciado", 0.018, 0.12)


def _clima(T=15.0, v=3.0, G=0.0, c=0.0, horas=48):
    t = np.arange(horas + 1, dtype=float)
    return Clima(t, np.full(t.shape, T), np.full(t.shape, v), np.full(t.shape, G), np.full(t.shape, c))


def test_conveccion_con_viento():
    assert h_conveccion(0.0) == pytest.approx(5.6)
    assert h_conveccion(5.0) == pytest.approx(25.6)
    assert h_conveccion(10.0) == pytest.approx(7.2 * 10**0.78)
    assert np.all(np.diff(h_conveccion(np.linspace(0, 15, 50))) > 0)  # más viento, más convección


def test_cielo_mas_frio_que_el_aire_salvo_nublado():
    assert T_cielo(15.0, 0.0) < 15.0 - 10  # cielo despejado: varios grados bajo el aire
    assert T_cielo(15.0, 1.0) == pytest.approx(15.0)  # nube gruesa: emite como el aire


def test_modelo_combinado_no_agrega_radiacion():
    amb = Ambiente(12.0, h_superior=27.8, h_lateral=8.0)
    h, T_eq = amb.exterior(0.0, np.array([SUP, LAT]), np.array([0.65, 0.65]), np.array([0.9, 0.9]))
    assert list(h) == [27.8, 8.0] and np.all(T_eq == 12.0)  # el h fijo ya incluye la radiación
    assert coeficiente_U(27.8, 0.0, 0.1, 2.222) == pytest.approx(1 / (1 / 27.8 + 0.05 / 2.222))
    assert coeficiente_U(0.0, 0.15, 0.1, 2.222) == 0.0  # cara adiabática


def test_temperatura_sol_aire_suma_cada_mecanismo_una_vez():
    """q = h_ext (T_eq - T_s) reproduce convección + onda larga + sol, para cualquier T_s."""
    amb = Ambiente(modelo="fisico", clima=_clima(T=14.0, v=4.0, G=800.0, c=0.3))
    orient, alfa, eps = np.array([SUP, LAT]), np.array([0.65, 0.65]), np.array([0.9, 0.9])
    h, T_eq = amb.exterior(0.0, orient, alfa, eps)
    T_rad = np.array([1.0, 0.5]) * T_cielo(14.0, 0.3) + np.array([0.0, 0.5]) * 14.0
    h_c = h_conveccion(4.0)
    h_r = 4 * 0.9 * SIGMA * ((14.0 + T_rad) / 2 + 273.15) ** 3
    for T_s in (5.0, 30.0, 60.0):
        q = h * (T_eq - T_s)
        assert q == pytest.approx(h_c * (14.0 - T_s) + h_r * (T_rad - T_s) + alfa * 800.0 * np.array([1, 0]))
    # Sin sol y con cielo a la temperatura del aire, T_eq es exactamente la del aire: no hay doble conteo.
    amb_n = Ambiente(modelo="fisico", clima=_clima(T=14.0, v=4.0, G=0.0, c=1.0))
    _, T_eq_n = amb_n.exterior(0.0, orient, alfa, eps)
    assert np.allclose(T_eq_n, 14.0)


def test_plan_de_moldaje_y_desmolde():
    plan = PlanCubiertas([Cubierta("moldaje", TERCIADO, (LAT,), (1, 2), {1: 72.0, 2: 72.0})])
    inicio = np.array([0.0, 6.0, np.inf])  # V3 oculta
    tabla, combos = plan.tabla(10.0, inicio, 3)
    assert combos[tabla[1, LAT]] == (TERCIADO,) and combos[tabla[2, LAT]] == (TERCIADO,)
    assert combos[tabla[1, SUP]] == () and tabla[3, LAT] == 0
    tabla, combos = plan.tabla(75.0, inicio, 3)  # V1 desmoldada a las 72 h; V2 recién a las 78 h
    assert combos[tabla[1, LAT]] == () and combos[tabla[2, LAT]] == (TERCIADO,)
    assert [t for t, _ in plan.eventos(inicio)] == [72.0, 78.0]


def test_moldaje_es_una_resistencia_en_serie_y_el_desmolde_la_quita():
    """Caja con moldaje lateral que se retira a las 5 h, contra solver_nd en dos tramos con h' = 1/(1/h + R)."""
    vox = _vox_cajas([4])
    dom = construir_dominio(vox, profundidad_suelo=0.0)
    amb = Ambiente(12.0, h_superior=27.8, h_lateral=15.0)
    cubiertas = PlanCubiertas([Cubierta("moldaje", TERCIADO, (LAT,), (1,), {1: 5.0})])
    hid = HidratacionJuLei2019()
    res = simular_3d(dom, PlanVertidas(vox), HORM, SUEL, hid, amb, 9.0, dt_max_s=600.0, instantes_3d_h=[5, 9],
                     cubiertas=cubiertas)
    assert [t for t, _ in res.retiros] == [5.0]
    dt = res.dts[0][1]
    alpha = HORM.k / HORM.capacidad
    h_mol = 1 / (1 / 15.0 + TERCIADO.R)
    bordes = _robin(27.8, h_mol, 12.0)
    edad = np.zeros((6, 5, 4))
    T1 = simular_difusion(np.full((6, 5, 4), 18.0), alpha, (0.25,) * 3, dt, round(5 * 3600 / dt), rho=HORM.rho,
                          cp=HORM.cp, bordes=bordes, hidratacion=hid, edad=edad)
    assert np.allclose(res.campos[5.0][1:-1, 1:-1, 1:-1], T1, atol=1e-4)
    T2 = simular_difusion(T1, alpha, (0.25,) * 3, dt, round(4 * 3600 / dt), rho=HORM.rho, cp=HORM.cp,
                          bordes=_robin(27.8, 15.0, 12.0), hidratacion=hid, edad=edad)
    assert np.allclose(res.campos[9.0][1:-1, 1:-1, 1:-1], T2, atol=1e-4)


def test_tmin_a_25_mm_entre_la_superficie_y_el_centro_de_la_celda():
    vox = _vox_cajas([4], dx=0.2)
    dom = construir_dominio(vox, profundidad_suelo=0.0)
    amb = Ambiente(5.0, h_superior=27.8, h_lateral=27.8)
    hid = HidratacionJuLei2019()
    kw = dict(dt_max_s=600.0)
    centro = simular_3d(dom, PlanVertidas(vox), HORM, SUEL, hid, amb, 24.0, **kw)
    mitad = simular_3d(dom, PlanVertidas(vox), HORM, SUEL, hid, amb, 24.0, profundidad_tmin=0.1, **kw)
    p25 = simular_3d(dom, PlanVertidas(vox), HORM, SUEL, hid, amb, 24.0, profundidad_tmin=0.025, **kw)
    sup = simular_3d(dom, PlanVertidas(vox), HORM, SUEL, hid, amb, 24.0, profundidad_tmin=0.0, **kw)
    assert np.allclose(mitad.Tmin_sup[1:], centro.Tmin_sup[1:])  # dx/2 = el centro de la celda
    assert np.all(sup.Tmin_sup[1:] < p25.Tmin_sup[1:]) and np.all(p25.Tmin_sup[1:] < centro.Tmin_sup[1:])
    assert np.all(p25.dT[1:] > centro.dT[1:])  # el criterio a 25 mm es más exigente que el centro de la celda
    with pytest.raises(ValueError):
        simular_3d(dom, PlanVertidas(vox), HORM, SUEL, hid, amb, 1.0, profundidad_tmin=0.2)


def test_tmin_en_caras_planas_excluye_aristas():
    """Las aristas y vértices se enfrían por dos o tres caras: el criterio 'plana' los deja fuera y los informa aparte."""
    vox = _vox_cajas([4], dx=0.2)
    dom = construir_dominio(vox, profundidad_suelo=0.0)
    amb = Ambiente(5.0, h_superior=27.8, h_lateral=27.8)
    hid = HidratacionJuLei2019()
    plana = simular_3d(dom, PlanVertidas(vox), HORM, SUEL, hid, amb, 24.0, dt_max_s=600.0, profundidad_tmin=0.025,
                       superficie_tmin="plana")
    todas = simular_3d(dom, PlanVertidas(vox), HORM, SUEL, hid, amb, 24.0, dt_max_s=600.0, profundidad_tmin=0.025)
    assert np.all(plana.Tmin_sup[1:] > todas.Tmin_sup[1:])
    assert np.allclose(plana.Tmin_todas, todas.Tmin_sup)  # la arista más fría queda registrada igual
    assert "caras planas" in plana.criterio_Tmin
