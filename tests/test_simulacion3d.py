"""Solver 3D por vertidas (E08/E09): equivalencia con solver_nd, activación, hidratación e inputs."""

import numpy as np
import pytest

from mcplan.caso import temperatura_ambiente
from mcplan.hidratacion import HidratacionJuLei2019
from mcplan.ifc_io import crear_modelo_base, leer_vertidas
from mcplan.simulacion3d import (
    HORMIGON,
    SUELO,
    Ambiente,
    Material,
    aplicar_cambios,
    cargar_caso_3d,
    construir_dominio,
    preparar,
    simular_3d,
)
from mcplan.solver_nd import simular_difusion
from mcplan.voxelizacion import PlanVertidas, Voxelizacion

HORM = Material(2.222, 2400, 970)
SUEL = Material(1.389, 1800, 1273)
HID = HidratacionJuLei2019()


def _vox_cajas(capas, dx=0.25, n=(6, 5)):
    """Voxelización a mano: cajas apiladas de n[0] x n[1] celdas, `capas` = celdas de alto de cada vertida."""
    nz = sum(capas)
    vid = np.zeros((n[0] + 2, n[1] + 2, nz + 2), dtype=np.int8)
    z = 1
    for v, c in enumerate(capas, start=1):
        vid[1:-1, 1:-1, z:z + c] = v
        z += c
    props = [{"Inicio_h": 0.0, "TemperaturaColocacion_C": 18.0} for _ in capas]
    return Voxelizacion(dx=dx, origen=np.array([-dx, -dx, -dx]), vertida_id=vid,
                        nombres=[f"Vertida_0{v}" for v in range(1, len(capas) + 1)], propiedades=props,
                        volumen_ifc=[n[0] * n[1] * c * dx**3 for c in capas])


def _robin(h_sup, h_lat, T_amb):
    return {f"{e}_{l}": ("robin", ((h_sup if e == "z" else h_lat), T_amb)) for e in "xyz" for l in ("min", "max")}


def test_caja_igual_que_solver_nd():
    """Una vertida sin suelo, Robin en todas las caras: mismo balance que solver_nd (con hidratación)."""
    vox = _vox_cajas([4])
    dom = construir_dominio(vox, profundidad_suelo=0.0)
    amb = Ambiente(12.0, h_superior=27.8, h_lateral=8.0)
    res = simular_3d(dom, PlanVertidas(vox), HORM, SUEL, HID, amb, 10.0, dt_max_s=600.0, instantes_3d_h=[10])
    dt = res.dts[0][1]
    edad = np.zeros((6, 5, 4))
    ref = simular_difusion(np.full((6, 5, 4), 18.0), HORM.k / HORM.capacidad, (0.25,) * 3, dt, round(36000 / dt),
                           rho=HORM.rho, cp=HORM.cp, bordes=_robin(27.8, 8.0, 12.0), hidratacion=HID, edad=edad)
    assert np.allclose(res.campos[10.0][1:-1, 1:-1, 1:-1], ref, atol=1e-4)  # campos en float32
    assert np.allclose(res.edades[10.0][1:-1, 1:-1, 1:-1], edad, atol=1e-4)
    assert np.isnan(res.campos[10.0][0]).all()  # el aire no tiene temperatura


def test_activacion_por_etapas_igual_que_solver_nd_en_dos_tramos():
    """V2 se coloca a las 5 h sobre V1: hasta ahí la cara superior de V1 es Robin; después es junta."""
    vox = _vox_cajas([3, 2])
    vox.propiedades[1].update(Inicio_h=5.0, TemperaturaColocacion_C=25.0)
    dom = construir_dominio(vox, profundidad_suelo=0.0)
    amb = Ambiente(12.0, h_superior=27.8, h_lateral=8.0)
    res = simular_3d(dom, PlanVertidas(vox), HORM, SUEL, HID, amb, 9.0, dt_max_s=600.0, instantes_3d_h=[4, 9])
    dt = res.dts[0][1]
    alpha = HORM.k / HORM.capacidad
    bordes = _robin(27.8, 8.0, 12.0)
    e1 = np.zeros((6, 5, 3))
    T1 = simular_difusion(np.full((6, 5, 3), 18.0), alpha, (0.25,) * 3, dt, round(5 * 3600 / dt), rho=HORM.rho,
                          cp=HORM.cp, bordes=bordes, hidratacion=HID, edad=e1)
    T0 = np.concatenate([T1, np.full((6, 5, 2), 25.0)], axis=2)
    e = np.concatenate([e1, np.zeros((6, 5, 2))], axis=2)
    T2 = simular_difusion(T0, alpha, (0.25,) * 3, dt, round(4 * 3600 / dt), rho=HORM.rho, cp=HORM.cp,
                          bordes=bordes, hidratacion=HID, edad=e)
    assert np.isnan(res.campos[4.0][1:-1, 1:-1, 4:6]).all()  # V2 todavía no existe
    assert np.allclose(res.campos[9.0][1:-1, 1:-1, 1:-1], T2, atol=1e-4)
    assert res.puntos == {} and [ev[:2] for ev in res.eventos] == [(0.0, 1), (5.0, 2)]
    assert np.isnan(res.Tmax_vertida[4, 1]) and res.Tmax_vertida[5, 1] == pytest.approx(25.0)  # recién colocada


def test_bloque_adiabatico_sigue_la_curva_del_ensayo():
    vox = _vox_cajas([3])
    dom = construir_dominio(vox, profundidad_suelo=0.0)
    res = simular_3d(dom, PlanVertidas(vox), HORM, SUEL, HID, Ambiente(0.0, 0.0, 0.0), 72.0, dt_max_s=60.0,
                     instantes_3d_h=[72])
    T, te = res.campos[72.0][1:-1, 1:-1, 1:-1], res.edades[72.0][1:-1, 1:-1, 1:-1]
    assert np.ptp(T) < 1e-3  # sin pérdidas el campo queda uniforme
    assert T.mean() == pytest.approx(18.0 + float(HID.theta(te.mean())), abs=0.1)
    assert te.mean() > 72.0  # más caliente que 20 °C: envejece más rápido (Arrhenius)


@pytest.fixture(scope="module")
def ifc(tmp_path_factory):
    ruta = tmp_path_factory.mktemp("ifc") / "fundacion.ifc"
    crear_modelo_base(ruta, "geometry/course_geometry.json")
    return leer_vertidas(ruta)


def test_dominio_con_suelo_y_aire(ifc):
    from mcplan.voxelizacion import voxelizar

    vox = voxelizar(ifc, 0.5)
    dom = construir_dominio(vox, profundidad_suelo=2.0, margen_suelo=1.0)
    z = dom.centros(2)
    assert (dom.material[:, :, z < 0] == SUELO).all() and z[0] == pytest.approx(-1.75)
    assert np.count_nonzero(dom.material == HORMIGON) == np.count_nonzero(vox.vertida_id)
    assert dom.forma[0] == vox.forma[0] + 4


def test_la_secuencia_esta_en_los_inputs(ifc):
    """Cambiar inicio, temperatura de colocación u ocultar una vertida es un input: el solver no cambia."""
    caso = cargar_caso_3d("config/caso_3d.yaml")
    caso.dx, caso.duracion_h, caso.instantes_3d_h = 0.5, 14.0, [5, 7, 14]
    caso.profundidad_suelo, caso.margen_suelo = 1.0, 0.5
    caso.cambios = {"Vertida_02": {"Inicio_h": 7, "TemperaturaColocacion_C": 25}, "Vertida_05": {"activa": False}}
    vox, dom, plan, registro = preparar(caso, ifc)
    assert len(registro) == 3 and [e[0] for e in plan.eventos()][:3] == [0.0, 7.0, 12.0]
    assert plan.vertidas_activas(1e6) == [1, 2, 3, 4]
    res = simular_3d(dom, plan, caso.hormigon, caso.suelo, caso.hidratacion, caso.ambiente, caso.duracion_h,
                     guardar_cada_h=1.0, puntos={"V2": (0.0, 0.0, 0.8)}, instantes_3d_h=caso.instantes_3d_h)
    v2 = dom.vertida_id == 2
    assert np.isnan(res.campos[5.0][v2]).all()
    assert np.allclose(res.campos[7.0][v2], 25.0)  # recién colocada a la T pedida
    assert np.isnan(res.puntos["V2"][:7]).all() and res.puntos["V2"][7] == pytest.approx(25.0)
    assert np.isnan(res.Tmax_vertida[:, 4]).all()  # V5 oculta: nunca existe
    assert res.campos[14.0][v2].min() > 18.0  # V1 caliente le pasa calor y además hidrata
    assert np.nanmax(res.Tmax) > 18.0 and len(res.eventos) == 3
    with pytest.raises(ValueError):
        aplicar_cambios(vox, {"Vertida_01": {"Moldaje": "x"}})


def test_panel_interactivo(ifc, tmp_path):
    """El panel lleva beta = h/k de cada borde, JSON válido y un campo igual al de la simulación."""
    import base64
    import gzip
    import json

    from mcplan.panel import Recolector, datos_panel, escribir_panel

    caso = cargar_caso_3d("config/caso_3d.yaml")
    caso.dx, caso.duracion_h, caso.instantes_3d_h = 0.5, 12.0, [12]
    caso.cambios = {"Vertida_05": {"activa": False}}
    vox, dom, plan, registro = preparar(caso, ifc)
    rec = Recolector(dom, caso.duracion_h, [e[0] for e in plan.eventos()], dx_objetivo=0.5)
    res = simular_3d(dom, plan, caso.hormigon, caso.suelo, caso.hidratacion, caso.ambiente, caso.duracion_h,
                     puntos=caso.puntos, instantes_3d_h=caso.instantes_3d_h, observador=rec)
    ruta = tmp_path / "panel.html"
    escribir_panel(ruta, datos_panel(res, caso, vox, plan, registro, rec, "mcplan simular3d"))
    html = ruta.read_text(encoding="utf-8")
    assert "__DATOS_PANEL__" not in html
    texto = html.split('id="datos">')[1].split("</script>")[0]
    d = json.loads(texto, parse_constant=lambda c: pytest.fail(f"JSON inválido: {c}"))

    bordes = {b["superficie"]: b for b in d["bordes"]}
    sup = bordes["Hormigón: caras horizontales al aire (±z)"]
    assert sup["beta"] == pytest.approx(27.8 / 2.222) and sup["h_ef"] < sup["h"]
    assert bordes["Suelo: caras laterales"]["beta"] == 0.0 and bordes["Suelo: fondo"]["tipo"] == "Dirichlet"
    assert d["malla"]["dx"] == d["malla"]["dy"] == d["malla"]["dz"] == 0.5
    assert d["vertidas"][4]["inicio_h"] is None and d["vertidas"][4]["Tmax"] is None  # oculta
    assert [e["activas"] for e in d["etapas"]][-1] == [1, 2, 3]

    v = d["vista"]
    assert v["tiempos"][0] == 0.0 and v["tiempos"][-1] == 12.0
    campos = np.frombuffer(gzip.decompress(base64.b64decode(v["T"])), "<i2").reshape(len(v["tiempos"]), *v["forma"])
    ultimo = np.where(campos[-1] == v["sin_dato"], np.nan, campos[-1] / v["escala"])
    assert np.allclose(ultimo, res.campos[12.0], atol=0.03, equal_nan=True)  # s = 1: la misma malla, en pasos de 0,05 °C
