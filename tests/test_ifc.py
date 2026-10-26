"""IFC del curso: lectura, propiedades, detección de cambios y columna IFC -> solver.

Los tests generan su propio IFC (misma geometría oficial) en una carpeta temporal,
así no dependen del archivo geometry/foundation.ifc que el grupo edita en Bonsai.
"""

from pathlib import Path

import numpy as np
import pytest

ifcopenshell = pytest.importorskip("ifcopenshell")
import ifcopenshell.api.geometry  # noqa: E402
import ifcopenshell.api.pset  # noqa: E402
import ifcopenshell.util.element  # noqa: E402

from mcplan.ifc_io import (  # noqa: E402
    PROPIEDADES_VERTIDA,
    PSET_VERTIDA,
    agregar_propiedades_curso,
    columna_desde_ifc,
    comparar_resumenes,
    crear_modelo_base,
    leer_vertidas,
    resumir_ifc,
)

GEOMETRIA = Path(__file__).resolve().parents[1] / "geometry" / "course_geometry.json"
# Tabla de control de volúmenes de la guía (sección 19). V5 es un círculo en la
# guía y un polígono de 48 lados en el modelo: difiere 0.3 %.
VOLUMENES_GUIA = {"Vertida_01": 147.63, "Vertida_02": 147.63, "Vertida_03": 214.22,
                  "Vertida_04": 88.92, "Vertida_05": 38.48}


@pytest.fixture
def ruta_ifc(tmp_path):
    return crear_modelo_base(tmp_path / "foundation.ifc", GEOMETRIA)


def test_cinco_vertidas_con_las_propiedades_del_curso(ruta_ifc):
    vertidas = leer_vertidas(ruta_ifc)
    assert [v.nombre for v in vertidas] == [f"Vertida_0{i}" for i in range(1, 6)]
    for v in vertidas:
        assert set(PROPIEDADES_VERTIDA) <= set(v.propiedades)


def test_volumenes_y_dimensiones_de_la_guia(ruta_ifc):
    elementos = {e["nombre"]: e for e in resumir_ifc(ruta_ifc)["elementos"].values()}
    for nombre, volumen in VOLUMENES_GUIA.items():
        assert elementos[nombre]["volumen_m3"] == pytest.approx(volumen, rel=5e-3)
    lo = np.min([e["bbox_min"] for e in elementos.values()], axis=0)
    hi = np.max([e["bbox_max"] for e in elementos.values()], axis=0)
    np.testing.assert_allclose(hi - lo, [18.0, 18.0, 4.30], atol=1e-3)  # 18 m entre caras
    assert lo[2] == pytest.approx(0.0, abs=1e-9)


def test_agregar_propiedades_actualiza_sin_duplicar(ruta_ifc):
    agregar_propiedades_curso(ruta_ifc, valores={"Vertida_03": {"TemperaturaColocacion_C": 12.0}})
    f = ifcopenshell.open(str(ruta_ifc))
    v3 = next(p for p in f.by_type("IfcFooting") if p.Name == "Vertida_03")
    psets = ifcopenshell.util.element.get_psets(v3)
    assert psets[PSET_VERTIDA]["TemperaturaColocacion_C"] == 12.0
    assert psets[PSET_VERTIDA]["VertidaID"] == 3  # las demás propiedades se mantienen
    nombres = [r.RelatingPropertyDefinition.Name for r in v3.IsDefinedBy
               if r.is_a("IfcRelDefinesByProperties")]
    assert nombres.count(PSET_VERTIDA) == 1


def test_detecta_cambio_de_propiedad_y_de_geometria(ruta_ifc):
    """Simula lo que se hace en Bonsai: editar una propiedad y mover el pedestal."""
    antes = resumir_ifc(ruta_ifc)

    f = ifcopenshell.open(str(ruta_ifc))
    pedestal = next(p for p in f.by_type("IfcFooting") if p.Name == "Vertida_05")
    pset_id = ifcopenshell.util.element.get_psets(pedestal)[PSET_VERTIDA]["id"]
    ifcopenshell.api.pset.edit_pset(f, pset=f.by_id(pset_id), properties={"TemperaturaColocacion_C": 15.0})
    matriz = np.eye(4)
    matriz[0, 3] = 5.0
    ifcopenshell.api.geometry.edit_object_placement(f, product=pedestal, matrix=matriz)
    f.write(str(ruta_ifc))

    cambios = comparar_resumenes(antes, resumir_ifc(ruta_ifc))
    assert any("TemperaturaColocacion_C 18.0 -> 15.0" in c for c in cambios)
    assert any("Vertida_05: bbox_max" in c for c in cambios)
    assert comparar_resumenes(antes, antes) == []


def test_columna_en_el_eje_atraviesa_las_cinco_vertidas(ruta_ifc):
    columna = columna_desde_ifc(ruta_ifc)
    assert [c["nombre"] for c in columna] == [f"Vertida_0{i}" for i in range(1, 6)]
    assert columna[0]["z0"] == pytest.approx(0.0, abs=1e-6)
    assert columna[-1]["z1"] == pytest.approx(4.30, abs=1e-6)
    assert [c["inicio_h"] for c in columna] == [0.0, 6.0, 12.0, 18.0, 24.0]


@pytest.mark.parametrize("x, vertidas, espesor", [
    (5.0, 4, 2.70),  # el tronco V4 tiene apotema 5 m a z = 2.70 m
    (8.5, 3, 1.30),  # el tronco V3 tiene apotema 8.5 m a z = 1.30 m
])
def test_columna_sobre_los_troncos_usa_la_geometria_real(ruta_ifc, x, vertidas, espesor):
    """Con el bounding box saldría el espesor completo de cada tronco; el rayo vertical no."""
    columna = columna_desde_ifc(ruta_ifc, x=x, y=0.0)
    assert len(columna) == vertidas
    assert columna[-1]["z1"] == pytest.approx(espesor, abs=1e-6)


def test_columna_fuera_de_la_fundacion(ruta_ifc):
    with pytest.raises(ValueError, match="Ninguna vertida"):
        columna_desde_ifc(ruta_ifc, x=20.0, y=0.0)


def test_sin_propiedades_avisa_como_agregarlas(tmp_path):
    ruta = crear_modelo_base(tmp_path / "sin_props.ifc", GEOMETRIA, valores={})
    with pytest.raises(ValueError, match="mcplan ifc-propiedades"):
        columna_desde_ifc(ruta)
