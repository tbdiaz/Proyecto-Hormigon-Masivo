"""Excavación: la fundación en un hoyo recibe menos viento y ve menos cielo bajo el terreno natural."""

import math

import numpy as np
import pytest

from mcplan.bordes import LAT, SUP
from mcplan.ifc_io import crear_modelo_base, leer_vertidas
from mcplan.simulacion3d import (
    HORMIGON,
    SUELO,
    VACIO_M,
    Excavacion,
    _caras_expuestas,
    _distancia_huella,
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


def test_perfil_de_viento_en_1991():
    e = Excavacion(2.0)  # categoría II: z0 = 0,05 m, z_min = 2 m
    assert e.factor_altura(10.0) == pytest.approx(1.0)
    assert e.factor_altura(0.5) == pytest.approx(math.log(2 / 0.05) / math.log(10 / 0.05))  # bajo z_min: el de z_min
    assert e.factor_altura(0.5) == pytest.approx(0.696, abs=1e-3)
    assert Excavacion(2.0, perfil_viento=False).factor_altura(0.5) == 1.0
    with pytest.raises(ValueError):
        Excavacion(2.0, categoria="V")


def test_distancia_al_borde_de_la_huella():
    huella = np.zeros((11, 11), dtype=bool)
    huella[3:8, 3:8] = True
    s = _distancia_huella(huella, 0.5)
    assert s[2, 5] == pytest.approx(0.25) and s[0, 5] == pytest.approx(1.25)  # afuera: al borde
    assert s[3, 5] == pytest.approx(-0.25) and s[5, 5] == pytest.approx(-1.25)  # adentro


def test_hoyo_con_paredes_de_suelo_natural(ifc):
    vox = voxelizar(ifc, 0.5)
    e = Excavacion(2.0, sobreancho_m=1.0)
    dom = construir_dominio(vox, excavacion=e)
    assert dom.cota_terreno == pytest.approx(2.0)
    zc = dom.centros(2)
    s = dom.distancia_plano
    nat = dom.natural
    assert (dom.material[nat] == SUELO).all()
    assert (s[nat.any(axis=2)] >= 1.0 - 1e-9).all()  # las paredes están a 1 m de la base o más
    hoyo = (s > 0) & (s < 1.0)
    capas = (zc > 0) & (zc < 2.0)
    assert (dom.material[hoyo][:, capas] == VACIO_M).all()  # dentro del hoyo: aire hasta el terreno
    assert (dom.material[:, :, -1] == VACIO_M).all()
    assert s.max() >= 1.0 + 1.0 - 1e-9  # queda al menos 1 m de suelo natural alrededor


def test_factores_de_las_caras(ifc):
    vox = voxelizar(ifc, 0.5)
    e = Excavacion(2.0, sobreancho_m=1.0, factor_abrigo=0.5)
    dom = construir_dominio(vox, excavacion=e)
    activas = dom.material != VACIO_M
    celda, direccion, orientacion = _caras_expuestas(activas)
    h_ = dom.material.reshape(-1)[celda] == HORMIGON
    fv, fc, fs = e.factores(dom, celda, direccion, orientacion)
    z = dom.centros(2)[np.unravel_index(celda, dom.forma)[2]]
    bajo = h_ & (orientacion == LAT) & (z < 2.0)
    sobre = h_ & (orientacion == LAT) & (z > 2.0)
    assert bajo.any() and sobre.any()
    assert np.allclose(fv[bajo], 0.696 * 0.5, atol=1e-3)  # en el hoyo: perfil bajo z_min × abrigo
    assert np.allclose(fv[sobre], e.factor_altura(z[sobre] - 2.0))  # sobre el terreno: perfil a su altura
    assert fv[sobre].min() == pytest.approx(0.696, abs=1e-3) and fv[sobre].max() > 0.7  # el pedestal pasa z_min
    assert (fc[bajo] < 0.5).all() and np.allclose(fc[sobre], 0.5)  # en el hoyo ven menos cielo
    # Más hondo, menos cielo (misma distancia a la pared): la cara del borde de la losa
    borde = bajo & (dom.distancia_plano[np.unravel_index(celda, dom.forma)[:2]] > -0.3)
    zs, cs = z[borde], fc[borde]
    assert cs[zs == zs.min()].max() < cs[zs == zs.max()].min()
    # Fórmula 2D (cuerdas cruzadas): a 0,75 m bajo el terreno y 1 m de la pared
    d, G = 2.0 - 1.25, 1.0
    assert 0.5 * (1 - d / math.hypot(d, G)) == pytest.approx(0.2, abs=1e-12)
    sup_bajo = h_ & (orientacion == SUP) & (z < 2.0)
    assert (fs[sup_bajo] == fc[sup_bajo]).all() and (fc[sup_bajo] <= 1).all()


def _correr(ifc, excavacion):
    caso = cargar_caso_3d("config/caso_14d.yaml")
    caso.dx, caso.duracion_h = 0.5, 36.0
    caso.excavacion_cfg = excavacion
    vox, dom, plan, _ = preparar(caso, ifc)
    res = simular_3d(dom, plan, caso.hormigon, caso.suelo, caso.hidratacion, caso.ambiente, caso.duracion_h,
                     T_suelo_inicial=caso.T_suelo_inicial, puntos=caso.puntos, cubiertas=caso.cubiertas,
                     profundidad_tmin=0.025, superficie_tmin="plana")
    return res


def test_el_hoyo_abriga_la_superficie(ifc):
    sin = _correr(ifc, {})
    con = _correr(ifc, {"profundidad": 2.2, "sobreancho": 1.0})
    assert con.balance["residuo_relativo"] < 1e-9
    assert con.puntos["borde_mitad"][-1] > sin.puntos["borde_mitad"][-1]  # junto a la cara lateral de V1/V2
    assert con.Tmin_sup[-1] > sin.Tmin_sup[-1]
