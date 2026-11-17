"""Programa de agua por tramos en los tubos y horizonte que incluye el retiro de las protecciones."""

import numpy as np
import pytest

from mcplan.bordes import Ambiente
from mcplan.estudio import _max_retiro
from mcplan.hidratacion import HidratacionJuLei2019
from mcplan.simulacion3d import construir_dominio, simular_3d
from mcplan.tubos import circuito_desde_dict, leer_programa
from mcplan.voxelizacion import PlanVertidas

from test_simulacion3d import HORM, SUEL, _vox_cajas

SERPENTIN = {"z": 0.5, "separacion": 0.5, "margen": 0.2}


def _caja():
    vox = _vox_cajas([4], dx=0.25, n=(8, 8))
    return vox, construir_dominio(vox, profundidad_suelo=0.0)


def test_programa_de_un_tramo_igual_a_temperatura_constante():
    vox, dom = _caja()
    fijo = circuito_desde_dict(dom, {"nombre": "t", "vertida": 1, "T_agua_C": 15, "durante_h": 12,
                                     "serpentin": SERPENTIN})
    prog = circuito_desde_dict(dom, {"nombre": "t", "vertida": 1, "programa": [{"hasta_h": 12, "T_agua_C": 15}],
                                     "serpentin": SERPENTIN})
    assert prog.tramos() == fijo.tramos() == [(0.0, 12.0, 15.0)]
    hid, amb = HidratacionJuLei2019(), Ambiente(12.0)
    a = simular_3d(dom, PlanVertidas(vox), HORM, SUEL, hid, amb, 24.0, tubos=[fijo])
    b = simular_3d(dom, PlanVertidas(vox), HORM, SUEL, hid, amb, 24.0, tubos=[prog])
    assert np.array_equal(a.Tmax, b.Tmax) and a.balance["tubos"] == b.balance["tubos"]


def test_el_agua_cambia_de_temperatura_en_cada_tramo():
    """Bloque a 18 °C sin hidratación: con agua a 18 °C no hay intercambio; desde las 6 h, a 5 °C, sí."""
    vox, dom = _caja()
    tubo = circuito_desde_dict(dom, {"nombre": "t", "vertida": 1, "material": "cobre", "serpentin": SERPENTIN,
                                     "programa": [{"hasta_h": 6, "T_agua_C": 18}, {"hasta_h": 12, "T_agua_C": 5}]})
    assert tubo.durante_h == 12 and tubo.T_agua_C == 18

    class SinCalor(HidratacionJuLei2019):
        def dtheta_dte(self, te):
            return np.zeros_like(np.asarray(te, dtype=float))

    res = simular_3d(dom, PlanVertidas(vox), HORM, SUEL, SinCalor(), Ambiente(18.0, h_superior=0.0, h_lateral=0.0), 18.0,
                     tubos=[tubo])
    s = res.tubos["t"]
    t, P = s["t_h"], s["potencia_W"]
    assert np.all(np.abs(P[(t > 0) & (t <= 6)]) < 1e-6)  # agua a la misma temperatura que el hormigón
    assert np.all(P[(t > 6.5) & (t < 12)] > 1.0)  # agua fría: el hormigón le entrega calor
    assert np.all(P[t > 12] == 0)  # el programa terminó
    assert res.balance["residuo_relativo"] < 1e-9 and res.balance["tubos"] < 0


def test_programa_mal_escrito_se_rechaza():
    assert leer_programa({}) == []
    assert leer_programa({"programa": [{"hasta_h": 168, "T_agua_C": 25}, {"hasta_h": 336, "T_agua_C": 20}]}) == \
        [(168.0, 25.0), (336.0, 20.0)]
    for malo in ([{"hasta_h": 100, "T_agua_C": 25}, {"hasta_h": 50, "T_agua_C": 20}],  # horas que no crecen
                 [{"hasta_h": 100}],  # falta la temperatura
                 [{"hasta_h": 100, "T_agua_C": 90}]):  # fuera de rango
        with pytest.raises(ValueError):
            leer_programa({"programa": malo})
    with pytest.raises(ValueError, match="salen del programa"):
        leer_programa({"programa": [{"hasta_h": 100, "T_agua_C": 25}], "durante_h": 50})


def test_mayor_retiro_de_las_protecciones():
    assert _max_retiro({}) == 0.0
    assert _max_retiro({"moldaje": {"retiro_h": 600}, "aislacion": {"capa": "manta_10cm", "retiro_h": 480}}) == 600
    assert _max_retiro({"moldaje": {"retiro_h": "pset"}, "aislacion": {"retiro_h": {1: 300, 2: 350}}}) == 350
