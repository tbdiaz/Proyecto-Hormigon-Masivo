"""Relleno de tierra: entierra la fundación en parte o entera (simulacion3d.Relleno)."""

import numpy as np
import pytest

from mcplan.ifc_io import crear_modelo_base, leer_vertidas
from mcplan.simulacion3d import (
    HORMIGON,
    SUELO,
    VACIO_M,
    cargar_caso_3d,
    construir_dominio,
    preparar,
    simular_3d,
)
from mcplan.voxelizacion import voxelizar


@pytest.fixture(scope="module")
def ifc(tmp_path_factory):
    ruta = tmp_path_factory.mktemp("ifc") / "fundacion.ifc"
    crear_modelo_base(ruta, "geometry/course_geometry.json")
    return leer_vertidas(ruta)


def _tope(ifc):
    return max(float(v.vertices[:, 2].max()) for v in ifc)


def _correr(ifc, relleno, duracion=48.0, instantes=(48.0,)):
    caso = cargar_caso_3d("config/caso_3d.yaml")
    caso.dx, caso.duracion_h = 0.5, duracion
    caso.relleno_cfg = relleno
    vox, dom, plan, _ = preparar(caso, ifc)
    res = simular_3d(dom, plan, caso.hormigon, caso.suelo, caso.hidratacion, caso.ambiente, caso.duracion_h,
                     puntos=caso.puntos, instantes_3d_h=list(instantes), profundidad_tmin=0.025,
                     superficie_tmin="plana", relleno=caso.relleno())
    return caso, dom, plan, res


def test_relleno_parcial_ocupa_el_aire_bajo_z(ifc):
    vox = voxelizar(ifc, 0.5)
    sin = construir_dominio(vox)
    dom = construir_dominio(vox, relleno_z=1.0)
    assert dom.forma == sin.forma  # bajo la cara superior: la grilla no crece
    r, zc = dom.relleno, dom.centros(2)
    assert (dom.material[r] == SUELO).all() and (dom.vertida_id[r] == 0).all()
    assert np.allclose(zc[r.any(axis=(0, 1))], [0.25, 0.75])  # dos capas de 0,5 m
    bajo = (zc > 0) & (zc < 1.0)
    assert ((sin.material[:, :, bajo] == VACIO_M) == r[:, :, bajo]).all()  # todo el aire bajo 1 m, nada más
    assert np.count_nonzero(dom.material == HORMIGON) == np.count_nonzero(sin.material == HORMIGON)


def test_relleno_entero_tapa_la_fundacion_y_deja_aire_encima(ifc):
    vox = voxelizar(ifc, 0.5)
    tope = _tope(ifc)
    dom = construir_dominio(vox, relleno_z=tope + 0.6)
    zc = dom.centros(2)
    assert dom.relleno[:, :, zc > tope].any()  # hay tierra sobre la cara superior
    assert (dom.material[:, :, -1] == VACIO_M).all()  # y una capa de aire encima: su cara superior ve el clima
    h = dom.material == HORMIGON
    k_alto = np.where(h.any(axis=2), dom.forma[2] - 1 - np.argmax(h[:, :, ::-1], axis=2), -1)
    i, j = np.nonzero(k_alto >= 0)
    assert dom.relleno[i, j, k_alto[i, j] + 1].all()  # cada columna de hormigón queda tapada
    with pytest.raises(ValueError):
        construir_dominio(vox, profundidad_suelo=0, relleno_z=1.0)


def test_relleno_aparece_a_su_hora_y_conserva_energia(ifc):
    caso, dom, plan, res = _correr(ifc, {"z": 1.0, "desde_h": 12, "T": 10.0}, instantes=(6.0, 12.0, 48.0))
    r = dom.relleno
    assert np.isnan(res.campos[6.0][r]).all()  # todavía no existe
    assert np.allclose(res.campos[12.0][r], 10.0)  # recién colocado a su temperatura
    assert np.isfinite(res.campos[48.0][r]).all()
    info = res.relleno
    assert info["colocado"] and info["t_h"] == 12 and info["T_C"] == 10 and not info["entera"]
    assert any(t == 12 and "relleno" in d for t, d in res.retiros)
    assert res.balance["residuo_relativo"] < 1e-9
    assert "enterradas" in res.criterio_Tmin


def test_hora_y_temperatura_automaticas(ifc):
    caso, dom, plan, res = _correr(ifc, {"z": 1.0}, duracion=24.0, instantes=(24.0,))
    info = res.relleno
    tocadas = [int(n[-2:]) for n in info["vertidas"]]
    assert tocadas and info["t_h"] == max(plan.inicio_h[v - 1] for v in tocadas)  # sin capas: al colocar la última
    horas = np.arange(max(info["t_h"] - 24, 0), info["t_h"] + 1e-9, 1.0)
    assert info["T_C"] == pytest.approx(np.mean([caso.ambiente.temperatura(h * 3600) for h in horas]))


def test_enterrar_protege_el_borde_y_tmin_sigue_definido(ifc):
    """Con la fundación enterrada entera ya no hay caras al aire, pero el criterio de ΔT sigue midiendo la
    superficie del hormigón (bajo las caras enterradas). Y el borde de la losa queda más tibio."""
    _, _, _, sin = _correr(ifc, {})
    caso, dom, plan, con = _correr(ifc, {"z": _tope(ifc) + 0.6, "desde_h": 24, "T": 15.0})
    assert con.relleno["entera"]
    assert np.isfinite(con.Tmin_sup[-1]) and np.isfinite(con.dT[-1])
    assert con.puntos["borde_mitad"][-1] > sin.puntos["borde_mitad"][-1]
    assert con.Tmin_sup[-1] > sin.Tmin_sup[-1]
    assert con.balance["residuo_relativo"] < 1e-9
