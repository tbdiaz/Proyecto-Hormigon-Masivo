"""Voxelización del IFC (E07): pertenencia por paridad, vertida_id, volumen, superficies y activación."""

import numpy as np
import pytest

from mcplan.ifc_io import VertidaIFC, crear_modelo_base, leer_vertidas, malla_prisma
from mcplan.voxelizacion import (
    VACIO,
    PlanVertidas,
    area_planta_octagono,
    caras,
    dentro_del_solido,
    volumen_malla,
    voxelizar,
)


def _caja(x0, x1, y0, y1, z0, z1):
    v = np.array([[x, y, z] for x in (x0, x1) for y in (y0, y1) for z in (z0, z1)], dtype=float)
    q = [(0, 2, 3, 1), (4, 5, 7, 6), (0, 1, 5, 4), (2, 6, 7, 3), (0, 4, 6, 2), (1, 3, 7, 5)]
    t = [tri for a, b, c, d in q for tri in ((a, b, c), (a, c, d))]
    return v, np.array(t)


def _prisma_L():
    """Prisma en forma de L (no convexo): planta [0,2]x[0,2] sin el cuadrado [1,2]x[1,2], alto 1."""
    planta = [(0, 0), (2, 0), (2, 1), (1, 1), (1, 2), (0, 2)]
    n = len(planta)
    v = np.array([(x, y, z) for z in (0.0, 1.0) for x, y in planta], dtype=float)
    tapa = [(0, 1, 2), (0, 2, 3), (0, 3, 4), (0, 4, 5)]  # triangulación en abanico desde (0, 0): válida aquí
    t = [tuple(reversed(f)) for f in tapa] + [(a + n, b + n, c + n) for a, b, c in tapa]
    for i in range(n):
        j = (i + 1) % n
        t += [(i, j, j + n), (i, j + n, i + n)]
    return v, np.array(t)


@pytest.fixture(scope="module")
def ifc(tmp_path_factory):
    ruta = tmp_path_factory.mktemp("ifc") / "fundacion.ifc"
    crear_modelo_base(ruta, "geometry/course_geometry.json")
    return leer_vertidas(ruta)


def test_volumen_exacto_de_mallas():
    assert volumen_malla(*_caja(0, 2, 0, 3, 0, 4)) == pytest.approx(24.0)
    assert volumen_malla(*_prisma_L()) == pytest.approx(3.0)
    v, poligonos = malla_prisma(8, 9.74153, 9.74153, 0.0, 0.55, 22.5)
    t = [(p[0], p[i], p[i + 1]) for p in poligonos for i in range(1, len(p) - 1)]  # abanico (caras convexas)
    assert volumen_malla(v, t) == pytest.approx(area_planta_octagono(18.0) * 0.55, rel=1e-5)


def test_paridad_en_un_solido_no_convexo():
    """Puntos genéricos: voxelizar() corre el rayo para no pasar justo por una arista compartida."""
    v, t = _prisma_L()
    xs, ys = np.array([0.43, 1.57]), np.array([0.61, 1.39])
    zs = np.array([0.25, 0.75, 1.5])
    dentro = dentro_del_solido(v, t, xs, ys, zs)
    assert dentro[0, 0, :2].all() and dentro[1, 0, :2].all() and dentro[0, 1, :2].all()
    assert not dentro[1, 1].any()  # el hueco de la L: un método "entre z mínimo y máximo" lo marcaría adentro
    assert not dentro[:, :, 2].any()  # sobre el sólido


def test_caja_alineada_con_la_grilla_se_voxeliza_exacta():
    v, t = _caja(0, 1, 0, 2, 0, 0.5)
    vox = voxelizar([VertidaIFC("g", "Vertida_01", v, t)], 0.1)
    assert vox.volumen_voxel(1) == pytest.approx(1.0)
    assert vox.superposiciones == 0


def test_fundacion_cinco_vertidas_apiladas(ifc):
    vox = voxelizar(ifc, 0.1)
    assert vox.n_vertidas == 5 and vox.superposiciones == 0
    j = int(np.argmin(np.abs(vox.centros(1))))
    i0 = int(np.argmin(np.abs(vox.centros(0))))
    columna = vox.vertida_id[i0, j][vox.vertida_id[i0, j] != VACIO]
    assert list(np.unique(columna)) == [1, 2, 3, 4, 5]
    assert np.all(np.diff(columna) >= 0)  # de abajo hacia arriba: V1, V2, ..., V5
    assert len(columna) * vox.dx == pytest.approx(4.30)
    i_borde = int(np.argmin(np.abs(vox.centros(0) - 8.5)))
    assert set(np.unique(vox.vertida_id[i_borde, j])) - {VACIO} == {1, 2, 3}  # cerca del borde: 1,30 m


def test_error_de_volumen_baja_con_dx_y_es_minimo_si_dx_divide_las_juntas(ifc):
    gruesa, fina = voxelizar(ifc, 0.25), voxelizar(ifc, 0.05)
    assert abs(fina.error_volumen_total) < abs(gruesa.error_volumen_total)
    assert abs(fina.error_volumen_total) < 5e-4
    assert all(abs(fina.error_volumen(v)) < 1e-3 for v in range(1, 6))  # juntas en múltiplos de 0,05 m
    media = voxelizar(ifc, 0.1)  # 0,55 m no es múltiplo de 0,1: una capa entera cambia de vertida
    assert abs(media.error_volumen(1)) > 0.05 and abs(media.error_volumen_total) < 0.005


def test_superficies_exteriores_e_interfaces(ifc):
    vox = voxelizar(ifc, 0.05)
    c = caras(vox.vertida_id)
    area = vox.dx**2
    assert c["exterior"]["-z"] * area == pytest.approx(area_planta_octagono(18.0), rel=5e-3)  # base sobre el suelo
    assert c["interfaz"][(1, 2)] * area == pytest.approx(area_planta_octagono(18.0), rel=5e-3)
    assert c["interfaz"][(4, 5)] * area == pytest.approx(np.pi * 3.5**2, rel=1e-2)  # pedestal sobre V4
    assert set(c["interfaz"]) == {(1, 2), (2, 3), (3, 4), (4, 5)}  # sólo vertidas consecutivas se tocan


def test_activacion_temporal(ifc):
    vox = voxelizar(ifc, 0.1)
    plan = PlanVertidas(vox)
    assert [e[0] for e in plan.eventos()] == [0.0, 6.0, 12.0, 18.0, 24.0]
    assert plan.vertidas_activas(0.0) == [1] and plan.vertidas_activas(13.0) == [1, 2, 3]
    T = np.full(vox.forma, np.nan)
    edad = np.full(vox.forma, -1.0)
    assert plan.activar(T, edad, -1.0, 0.0) == [1]
    assert np.all(T[vox.vertida_id == 1] == 18.0) and np.all(edad[vox.vertida_id == 1] == 0.0)
    assert np.isnan(T[vox.vertida_id == 2]).all()  # todavía no existe
    T[vox.vertida_id == 1] = 40.0  # el hormigón ya colocado se calentó
    assert plan.activar(T, edad, 0.0, 6.0) == [2]
    assert np.all(T[vox.vertida_id == 1] == 40.0) and np.all(T[vox.vertida_id == 2] == 18.0)
    # con sólo V1 activa, su cara superior completa está expuesta; al poner V2 pasa a ser interfaz
    solo_v1 = caras(vox.vertida_id, plan.celdas_activas(0.0))
    assert solo_v1["interfaz"] == {} and solo_v1["exterior"]["+z"] == caras(vox.vertida_id)["interfaz"][(1, 2)]
