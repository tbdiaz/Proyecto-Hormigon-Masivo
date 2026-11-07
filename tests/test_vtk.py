"""Exportación VTK (E11): el .vti tiene los mismos valores que la simulación, en el orden de VTK."""

import numpy as np
import pytest

from mcplan.vtk_io import escribir_vti, exportar_vtk, instantes_vtk, nucleo_y_superficie, superficie_expuesta

pv = pytest.importorskip("pyvista")


def test_vti_respeta_el_orden_de_vtk(tmp_path):
    a = np.arange(4 * 3 * 2, dtype=np.float32).reshape(4, 3, 2)  # (nx, ny, nz)
    escribir_vti(tmp_path / "a.vti", (1.0, 2.0, -3.0), 0.5, {"T": a, "material": np.full(a.shape, 2, np.int8)})
    m = pv.read(str(tmp_path / "a.vti"))
    assert m.dimensions == (5, 4, 3) and m.n_cells == 24
    assert np.allclose(m.bounds, (1.0, 3.0, 2.0, 3.5, -3.0, -2.0))
    centros = m.cell_centers().points
    i, j, k = ((centros - (1.0, 2.0, -3.0)) / 0.5 - 0.5).round().astype(int).T
    assert np.array_equal(m.cell_data["T"], a[i, j, k])  # cada celda trae el valor de su (i, j, k)


def test_superficie_expuesta():
    material = np.zeros((5, 5, 5), dtype=np.int8)
    material[1:4, 1:4, 1:4] = 2
    sup = superficie_expuesta(material, material > 0)
    assert sup[2, 2, 2] == 0 and sup[1, 2, 2] == 1 and sup.sum() == 26  # cubo 3x3x3: todo menos el centro


def test_exportar_resultado_y_ubicar_nucleo(tmp_path):
    from mcplan.hidratacion import HidratacionJuLei2019
    from mcplan.simulacion3d import Ambiente, construir_dominio, simular_3d
    from mcplan.voxelizacion import PlanVertidas

    from test_simulacion3d import HORM, SUEL, _vox_cajas

    vox = _vox_cajas([3, 2], dx=0.25)
    vox.propiedades[1]["Inicio_h"] = 4.0
    dom = construir_dominio(vox, profundidad_suelo=1.0, margen_suelo=0.5)
    res = simular_3d(dom, PlanVertidas(vox), HORM, SUEL, HidratacionJuLei2019(), Ambiente(10.0), 12.0,
                     instantes_3d_h=[2, 12])
    pvd = exportar_vtk(res, tmp_path / "vtk")
    instantes = instantes_vtk(pvd.parent)
    assert [t for t, _ in instantes] == sorted(res.campos)
    t, ruta = instantes[-1]
    m = pv.read(str(ruta))
    T = np.asarray(m.cell_data["T"])
    assert np.nanmax(T) == pytest.approx(np.nanmax(res.campos[t]), abs=1e-4)
    info = nucleo_y_superficie(m.threshold(0.5, scalars="material"))
    assert info["nucleo"][1] == pytest.approx(res.Tmax[list(res.t_h).index(t)], abs=1e-3)
    assert info["superficie"][1] < info["nucleo"][1]
    v2 = np.asarray(m.cell_data["vertida_id"]) == 2
    temprano = pv.read(str(instantes[0][1]))
    assert np.isnan(np.asarray(temprano.cell_data["T"])[v2]).all()  # a las 2 h la V2 todavía no existe
