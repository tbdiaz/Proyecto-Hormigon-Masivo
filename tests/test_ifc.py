"""Round-trip IFC: crear, leer, modificar y detectar el cambio desde Python."""

import pytest

ifcopenshell = pytest.importorskip("ifcopenshell")
import ifcopenshell.api.pset  # noqa: E402
import ifcopenshell.util.element  # noqa: E402

from mcplan.ifc_io import (  # noqa: E402
    PROPIEDADES_VERTIDA,
    PSET_VERTIDA,
    columna_desde_ifc,
    comparar_resumenes,
    crear_modelo_base,
    resumir_ifc,
)


@pytest.fixture
def ruta_ifc(tmp_path):
    return crear_modelo_base(tmp_path / "fundacion.ifc")


def test_tres_vertidas_identificables_con_sus_propiedades(ruta_ifc):
    elementos = list(resumir_ifc(ruta_ifc)["elementos"].values())
    assert len(elementos) == 3
    assert sorted(e["psets"][PSET_VERTIDA]["VertidaID"] for e in elementos) == ["V1", "V2", "V3"]
    for e in elementos:
        assert set(PROPIEDADES_VERTIDA) <= set(e["psets"][PSET_VERTIDA])


def test_dimensiones_globales_del_caso(ruta_ifc):
    elementos = resumir_ifc(ruta_ifc)["elementos"].values()
    xmin = min(e["bbox_min"][0] for e in elementos)
    xmax = max(e["bbox_max"][0] for e in elementos)
    zmax = max(e["bbox_max"][2] for e in elementos)
    pedestal = next(e for e in elementos if e["psets"][PSET_VERTIDA]["VertidaID"] == "V3")
    assert xmax - xmin == pytest.approx(18.0, abs=1e-3)
    assert zmax == pytest.approx(3.30, abs=1e-3)
    assert pedestal["bbox_max"][0] - pedestal["bbox_min"][0] == pytest.approx(7.0, abs=1e-2)


def test_detecta_cambio_de_propiedad_y_de_geometria(ruta_ifc):
    """Simula lo que se hace en Bonsai: editar una propiedad y la altura del pedestal."""
    antes = resumir_ifc(ruta_ifc)

    f = ifcopenshell.open(str(ruta_ifc))
    pedestal = next(e for e in f.by_type("IfcFooting") if e.Name.startswith("V3"))
    pset_id = ifcopenshell.util.element.get_psets(pedestal)[PSET_VERTIDA]["id"]
    ifcopenshell.api.pset.edit_pset(f, pset=f.by_id(pset_id), properties={"TemperaturaColocacion_C": 15.0})
    pedestal.Representation.Representations[0].Items[0].Depth = 1.5
    f.write(str(ruta_ifc))

    cambios = comparar_resumenes(antes, resumir_ifc(ruta_ifc))
    assert any("TemperaturaColocacion_C 20.0 -> 15.0" in c for c in cambios)
    assert any("bbox_max" in c for c in cambios)
    assert comparar_resumenes(antes, antes) == []


def test_columna_en_el_eje(ruta_ifc):
    columna = columna_desde_ifc(ruta_ifc)
    assert [c["VertidaID"] for c in columna] == ["V1", "V2", "V3"]
    assert columna[0]["z0"] == pytest.approx(0.0, abs=1e-6)
    assert columna[-1]["z1"] == pytest.approx(3.30, abs=1e-6)


def test_columna_fuera_del_pedestal_no_lo_incluye(ruta_ifc):
    columna = columna_desde_ifc(ruta_ifc, x=6.0, y=6.0)  # dentro de la losa, fuera del cilindro
    assert [c["VertidaID"] for c in columna] == ["V1", "V2"]
