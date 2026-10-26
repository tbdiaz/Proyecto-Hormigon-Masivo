"""Los casos YAML del repositorio cargan y producen resultados físicamente razonables."""

from pathlib import Path

import numpy as np
import pytest
import yaml

from mcplan.caso import _borde, cargar_caso, temperatura_ambiente
from mcplan.metricas import calcular_metricas
from mcplan.solver1d import crear_malla, simular

RAIZ = Path(__file__).resolve().parents[1]
CONFIG = RAIZ / "config"


def _correr(caso, horas=48, dz=0.1):
    malla = crear_malla(caso.capas, dz)
    return simular(malla, caso.borde_inf, caso.borde_sup, t_final=horas * 3600.0)


def test_caso_base_corre_y_se_calienta():
    caso = cargar_caso(CONFIG / "caso_base_1d.yaml")
    res = _correr(caso)
    assert np.all(np.isfinite(res.T))
    met = calcular_metricas(res)
    assert met.Tmax > 40.0  # el calor de hidratación eleva el núcleo
    assert met.dTmax > 0.0


def test_caso_ifc_lee_las_vertidas_del_modelo():
    pytest.importorskip("ifcopenshell")
    caso = cargar_caso(CONFIG / "caso_ifc_1d.yaml")
    hormigon = [c for c in caso.capas if c.material.es_hormigon]
    assert len(hormigon) >= 1
    assert all(c.espesor > 0 for c in hormigon)


def test_puente_ifc_reproduce_el_caso_base(tmp_path):
    """El IFC oficial leído en el eje debe dar exactamente lo mismo que el caso base
    escrito a mano (misma columna: 4.30 m a 18 °C). Si no, la lectura del IFC falla."""
    pytest.importorskip("ifcopenshell")
    from mcplan.ifc_io import crear_modelo_base

    ruta_ifc = crear_modelo_base(tmp_path / "foundation.ifc", RAIZ / "geometry" / "course_geometry.json")
    cfg = yaml.safe_load((CONFIG / "caso_ifc_1d.yaml").read_text(encoding="utf-8"))
    cfg["hormigon_desde_ifc"]["archivo"] = str(ruta_ifc)
    (tmp_path / "caso.yaml").write_text(yaml.safe_dump(cfg), encoding="utf-8")

    desde_ifc = calcular_metricas(_correr(cargar_caso(tmp_path / "caso.yaml"), horas=72, dz=0.05))
    a_mano = calcular_metricas(_correr(cargar_caso(CONFIG / "caso_base_1d.yaml"), horas=72, dz=0.05))
    assert desde_ifc.Tmax == pytest.approx(a_mano.Tmax, abs=1e-6)
    assert desde_ifc.dTmax == pytest.approx(a_mano.dTmax, abs=1e-6)


def test_ambiente_sinusoidal_es_maximo_a_la_hora_indicada():
    T = temperatura_ambiente({"tipo": "sinusoidal", "media": 16, "amplitud": 7, "hora_max": 15}, hora_inicio=8)
    assert T(7 * 3600.0) == pytest.approx(23.0)  # 8 h + 7 h = 15 h
    assert T(19 * 3600.0) == pytest.approx(9.0)  # 3 h de la madrugada


def test_aislacion_reduce_el_coeficiente_equivalente():
    sin = _borde({"tipo": "robin", "h": 27.8, "T_amb": 15.0}, 0)
    con = _borde({"tipo": "robin", "h": 27.8, "T_amb": 15.0, "aislacion": {"espesor": 0.05, "k": 0.04}}, 0)
    assert con.h == pytest.approx(1.0 / (1.0 / 27.8 + 0.05 / 0.04))
    assert con.h < sin.h
