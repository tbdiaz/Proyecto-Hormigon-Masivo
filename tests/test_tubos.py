"""Tubos de enfriamiento y hormigón distinto por vertida."""

import numpy as np
import pytest

from mcplan.bordes import Ambiente
from mcplan.hidratacion import HidratacionJuLei2019
from mcplan.simulacion3d import Material, construir_dominio, simular_3d
from mcplan.solver_nd import simular_difusion
from mcplan.tubos import Circuito, agua_a_lo_largo, circuito_desde_dict, recorrer, serpentin
from mcplan.voxelizacion import PlanVertidas

from test_simulacion3d import HORM, SUEL, _robin, _vox_cajas


def test_resistencias_y_material_del_tubo():
    pvc, cobre = Circuito("a", np.zeros((2, 3)), "pvc"), Circuito("b", np.zeros((2, 3)), "cobre")
    _, R_pvc, Re = pvc.resistencias()
    assert R_pvc == pytest.approx(np.log(32 / 28) / (2 * np.pi * 0.19))  # pared de PVC: ~0,11 m K/W
    assert Re > 2300  # 20 L/min en Ø28 mm interior: turbulento
    assert cobre.U_L(2.222, 0.275) > pvc.U_L(2.222, 0.275)  # el cobre casi no opone resistencia
    assert pvc.U_L(4.0, 0.275) > pvc.U_L(2.222, 0.275)  # un hormigón más conductor ayuda
    assert pvc.m_c == pytest.approx(998 * 20 / 60000 * 4186)


def test_agua_a_lo_largo_del_tubo_exacta():
    a = np.full(500, np.exp(-0.01))
    Tw = agua_a_lo_largo(a, np.full(500, 50.0), 15.0)
    n = np.arange(1, 501)
    assert np.allclose(Tw, 50 + (15 - 50) * np.exp(-0.01 * n))  # se acerca a T del hormigón exponencialmente
    a[100:200] = 1.0  # tramo sin intercambio: el agua pasa sin cambiar
    Tw = agua_a_lo_largo(a, np.full(500, 50.0), 15.0)
    assert np.all(Tw[100:200] == Tw[99])


def _caja():
    vox = _vox_cajas([4], dx=0.25, n=(8, 8))
    return vox, construir_dominio(vox, profundidad_suelo=0.0)


def test_recorrido_en_celdas():
    vox, dom = _caja()
    celdas, largos = recorrer(dom, [[0.1, 0.6, 0.5], [1.9, 0.6, 0.5]])
    assert largos.sum() == pytest.approx(1.8) and len(celdas) == 8  # cruza 8 celdas de 0,25 m
    P = serpentin(dom, 1, 0.5, separacion=0.5, margen=0.2)
    c = circuito_desde_dict(dom, {"vertida": 1, "puntos": P.tolist()})
    assert np.all(dom.material.reshape(-1)[c.celdas] == 2)  # todo el tubo dentro del hormigón
    with pytest.raises(ValueError, match="fuera del hormigón"):
        circuito_desde_dict(dom, {"puntos": [[0.1, 0.6, 0.5], [0.1, 0.6, 1.1]]})  # sube al aire


def test_tubos_extraen_calor_con_balance_exacto_y_se_apagan():
    vox, dom = _caja()
    tubo = circuito_desde_dict(dom, {"nombre": "t", "vertida": 1, "material": "cobre", "durante_h": 12,
                                     "serpentin": {"z": 0.5, "separacion": 0.5, "margen": 0.2}})
    hid = HidratacionJuLei2019()
    amb = Ambiente(12.0)
    sin = simular_3d(dom, PlanVertidas(vox), HORM, SUEL, hid, amb, 24.0)
    con = simular_3d(dom, PlanVertidas(vox), HORM, SUEL, hid, amb, 24.0, tubos=[tubo])
    assert np.nanmax(con.Tmax) < np.nanmax(sin.Tmax)
    b = con.balance
    assert b["residuo_relativo"] < 1e-9 and b["tubos"] < 0  # lo que pierde el hormigón lo gana el agua
    s = con.tubos["t"]
    assert s["energia_J"][-1] == pytest.approx(-b["tubos"], rel=1e-12)  # energía acumulada por el agua = balance
    assert np.all(s["potencia_W"][s["t_h"] > 12] == 0)  # apagado después de durante_h
    assert np.nanmax(s["T_salida_C"]) < np.nanmax(con.Tmax)


def test_hormigon_distinto_por_vertida_igual_a_solver_nd():
    """Dos vertidas con k distinto (p. ej. con grafito): media armónica en la junta, como solver_nd."""
    vox = _vox_cajas([2, 2])
    dom = construir_dominio(vox, profundidad_suelo=0.0)
    grafito = Material(4.5, 2400, 970)
    amb = Ambiente(12.0, h_superior=27.8, h_lateral=8.0)
    res = simular_3d(dom, PlanVertidas(vox), HORM, SUEL, HidratacionJuLei2019(), amb, 6.0, dt_max_s=600.0,
                     instantes_3d_h=[6], materiales_vertida={2: grafito})
    k = np.full((6, 5, 4), HORM.k)
    k[:, :, 2:] = grafito.k
    dt = res.dts[0][1]
    edad = np.zeros((6, 5, 4))
    ref = simular_difusion(np.full((6, 5, 4), 18.0), HORM.k / HORM.capacidad, (0.25,) * 3, dt, round(6 * 3600 / dt),
                           rho=HORM.rho, cp=HORM.cp, bordes=_robin(27.8, 8.0, 12.0), conductividad=k,
                           hidratacion=HidratacionJuLei2019(), edad=edad)
    assert np.allclose(res.campos[6.0][1:-1, 1:-1, 1:-1], ref, atol=1e-4)
