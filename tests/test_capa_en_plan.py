"""Un plan puede declarar su propia capa de aislación (espesor y k) sin tocar el catálogo del caso."""

from types import SimpleNamespace

import pytest

from mcplan.bordes import Capa
from mcplan.simulacion3d import armar_cubiertas


def _caso(aislacion):
    capas = {"manta_10cm": Capa("manta_10cm", 0.10, 0.04, 0.4, 0.9)}
    return SimpleNamespace(capas=capas, cubiertas_cfg={"aislacion": aislacion})


VOX = SimpleNamespace(n_vertidas=2, propiedades=[{}, {}])


def test_capa_del_catalogo_y_capa_declarada_en_el_plan():
    del_catalogo = armar_cubiertas(_caso({"capa": "manta_10cm", "caras": "todas", "retiro_h": 100}), VOX)
    propia = armar_cubiertas(_caso({"capa": {"espesor": 0.10, "k": 0.05, "nombre": "manta_R2",
                                             "absortividad": 0.4},
                                    "caras": "todas", "retiro_h": 100}), VOX)
    c0, c1 = del_catalogo.cubiertas[0], propia.cubiertas[0]
    assert c0.capa.R == pytest.approx(2.5)
    assert c1.capa.R == pytest.approx(2.0)
    assert c1.capa.nombre == "manta_R2" and c1.capa.absortividad == 0.4 and c1.capa.emisividad == 0.9
    assert (c1.orientaciones, c1.vertidas, c1.retiro_h) == (c0.orientaciones, c0.vertidas, c0.retiro_h)


def test_capa_declarada_incompleta_avisa():
    with pytest.raises(ValueError, match="le falta k"):
        armar_cubiertas(_caso({"capa": {"espesor": 0.10}}), VOX)
