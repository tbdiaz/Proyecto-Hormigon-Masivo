"""Decisión probabilística (E13/F02): curvas, planes, activación por capas, balance de energía y agregación."""

import numpy as np
import pytest

from mcplan.curvas import ARCHIVO, cargar_fuentes, leer_source_terms, tabla_curvas
from mcplan.estudio import agregar, decision
from mcplan.planes import Etapa, PlanObra, activacion_por_capas, cronograma, expandir_familias

from test_simulacion3d import HORM, SUEL, _vox_cajas


def test_ocho_curvas_y_energia_conservada():
    crudas = leer_source_terms(ARCHIVO)
    assert len(crudas) == 8 and {c["campana"] for c in crudas.values()} == {"SS", "22A"}
    fuentes = cargar_fuentes(ARCHIVO)
    for n, c in crudas.items():
        f = fuentes[n]
        E_crudo = np.sum((c["Q"][1:] + c["Q"][:-1]) / 2 * np.diff(c["t_h"])) * 3600
        E_tabla = float(f.energia(c["t_h"][-1] + 2.0) - f.energia(2.0))
        assert E_tabla == pytest.approx(E_crudo, rel=1e-4)  # el suavizado por bloques conserva la energía
        assert f.Q(0.0) == 0.0 and f.Q(2.0) == pytest.approx(c["Q"][0], rel=0.05)  # rampa en las 2 h sin datos
        assert np.all(f.factor_arrhenius(np.array([10.0, 80.0])) == 1.0)  # sin madurez
    assert all(r["E336_MJ_m3"] > r["E72_MJ_m3"] > 0 for r in tabla_curvas(fuentes, 2400 * 970))
    cero = cargar_fuentes(ARCHIVO, ["M1_SS"], cola="cero")["M1_SS"]
    assert cero.Q(300.0) == 0.0 and fuentes["M1_SS"].Q(300.0) > 0.0


def _dom_cajas():
    from mcplan.simulacion3d import construir_dominio

    vox = _vox_cajas([2, 2, 2], dx=0.25)  # tres vertidas de 0,5 m, 6 x 5 celdas en planta
    return vox, construir_dominio(vox, profundidad_suelo=0.0)


def test_cronograma_D_y_particion():
    vox, dom = _dom_cajas()
    V = sum(vox.volumen_ifc)
    unica = PlanObra("u", [Etapa((1, 2, 3))])
    tres = PlanObra("t", [Etapa((1,)), Etapa((2,), espera_h=10), Etapa((3,), espera_h=5)])
    assert cronograma(unica, dom, vox, 2.0)["D_h"] == pytest.approx(V / 2.0)
    c = cronograma(tres, dom, vox, 2.0)
    assert c["D_h"] == pytest.approx(V / 2.0 + 15)  # mismas horas de bomba + esperas
    assert [e["inicio_h"] for e in c["etapas"]] == pytest.approx([0, vox.volumen_ifc[0] / 2 + 10,
                                                                  (vox.volumen_ifc[0] + vox.volumen_ifc[1]) / 2 + 15])
    with pytest.raises(ValueError, match="no cubre"):
        cronograma(PlanObra("x", [Etapa((1, 2))]), dom, vox, 2.0)
    with pytest.raises(ValueError, match="no cubre"):
        cronograma(PlanObra("x", [Etapa((1, 2, 3)), Etapa((3,))]), dom, vox, 2.0)
    planes = expandir_familias([{"familia": "F", "etapas": [[1], [2], [3]], "esperas_h": [0, 24]}], 18.0)
    assert [p.nombre for p in planes] == ["F-e0", "F-e24"] and planes[1].etapas[2].espera_h == 24


def test_activacion_por_capas_sube_a_la_tasa():
    vox, dom = _dom_cajas()
    act, crono = activacion_por_capas(PlanObra("u", [Etapa((1, 2, 3), T_colocacion=12.0)]), dom, vox, 1.0)
    capas = sorted({float(t) for t in np.unique(act.t_act[np.isfinite(act.t_act)])})
    assert len(capas) == 6 and capas[0] == 0.0  # 6 capas de celdas, la primera en t0
    assert capas[-1] == pytest.approx(crono["D_h"] * 5 / 6)  # la última empieza cuando falta una capa
    assert np.all(act.T_col[np.isfinite(act.t_act)] == 12.0)
    assert act.inicio_h[1] == pytest.approx(crono["D_h"] / 3) and np.all(act.fin_h == pytest.approx(crono["D_h"]))
    dividida = PlanObra("d", [Etapa((1, 2)), Etapa((3,), z=(0.0, 1.25)), Etapa((3,), z=(1.25, 9.0), espera_h=4)])
    act2, c2 = activacion_por_capas(dividida, dom, vox, 1.0)
    assert c2["etapas"][1]["volumen_m3"] == pytest.approx(vox.volumen_ifc[2] / 2)


def test_balance_de_energia_y_activacion_en_el_solver():
    from mcplan.bordes import Ambiente
    from mcplan.simulacion3d import simular_3d
    from mcplan.voxelizacion import PlanVertidas

    vox, dom = _dom_cajas()
    fuente = cargar_fuentes(ARCHIVO, ["M4_22A"])["M4_22A"]
    act, crono = activacion_por_capas(PlanObra("t", [Etapa((1, 2)), Etapa((3,), espera_h=6)]), dom, vox, 1.0)
    res = simular_3d(dom, PlanVertidas(vox), HORM, SUEL, fuente, Ambiente(10.0), 48.0, activacion=act,
                     profundidad_tmin=0.025, superficie_tmin="plana", instantes_3d_h=[1, 48])
    b = res.balance
    assert b["residuo_relativo"] < 1e-9  # Euler explícito conserva la energía a precisión de máquina
    assert b["hidratacion"] > 0 > b["bordes"] and b["colocacion"] > 0
    assert np.isnan(res.campos[1.0][dom.vertida_id == 3]).all()  # la tercera vertida todavía no existe
    assert [e[1] for e in res.eventos] == [1, 2, 3]


def _esc(planes, curvas=("A", "B"), climas=(("k1", 0.5), ("k2", 0.5))):
    return {"planes": planes, "p_curva": {c: 1 / len(curvas) for c in curvas},
            "historias": [{"clave": k, "peso": q} for k, q in climas],
            "cronos": {p.nombre: {"D_h": float(i)} for i, p in enumerate(planes)}}


def test_probabilidad_conjunta_sin_doble_conteo():
    """Ejemplo de la pauta: 24 escenarios, uno excede sólo T y otro sólo Δ -> P_T = P_Δ = 4,17 %, P_f = 8,33 %."""
    plan = PlanObra("p", [Etapa((1,))])
    curvas = [f"c{i}" for i in range(8)]
    esc = _esc([plan], curvas, (("k1", 1 / 3), ("k2", 1 / 3), ("k3", 1 / 3)))
    filas = [{"plan": "p", "curva": c, "clima": k, "I_T": 0, "I_D": 0, "I_F": 0, "M_T": 50, "M_D": 10,
              "pico_al_final": 0, "residuo_energia": 0} for c in curvas for k in ("k1", "k2", "k3")]
    filas[0].update(I_T=1, I_F=1)
    filas[1].update(I_D=1, I_F=1)
    r = agregar(filas, esc)[0]
    assert r["P_T"] == pytest.approx(1 / 24) and r["P_D"] == pytest.approx(1 / 24)
    assert r["P_f"] == pytest.approx(2 / 24) and not r["cumple"]
    filas[1].update(I_D=0, I_F=0)
    filas[0].update(I_D=1)  # una sola corrida excede ambos: cuenta una vez
    r = agregar(filas, esc)[0]
    assert r["P_f"] == pytest.approx(1 / 24) and r["cumple"]


def test_decision_elige_el_mas_rapido_que_cumple():
    planes = [PlanObra(n, [Etapa((1,))]) for n in ("rapido", "medio", "lento")]
    esc = _esc(planes)
    filas = []
    for p, falla in (("rapido", 1), ("medio", 0), ("lento", 0)):
        filas += [{"plan": p, "curva": c, "clima": k, "I_T": falla, "I_D": 0, "I_F": falla, "M_T": 50, "M_D": 10,
                   "pico_al_final": 0, "residuo_energia": 0} for c in ("A", "B") for k in ("k1", "k2")]
    d = decision(agregar(filas, esc))
    assert d == {"propuesto": "medio", "mas_rapido": "rapido", "mas_lento": "lento"}
    incompleto = agregar([f for f in filas if f["plan"] != "medio" or f["curva"] == "A"], esc)
    assert not next(r for r in incompleto if r["plan"] == "medio")["cumple"]  # sin todas las corridas no hay decisión


def test_historias_reducidas_y_pesos():
    from mcplan.bordes import Clima
    from mcplan.clima import reducir_historias

    t = np.arange(24 * 8 + 1, dtype=float)
    historias = [(f"m{i}", Clima(t, 10 + i + 5 * np.sin(2 * np.pi * t / 24), t * 0 + 3, t * 0, t * 0))
                 for i in range(9)]
    elegidas = reducir_historias(historias, 3)
    assert [h["peso"] for h in elegidas] == pytest.approx([1 / 3] * 3)
    assert [h["representante"] for h in elegidas] == ["m1", "m4", "m7"]  # la mediana de cada tercil
    assert elegidas[0]["etiqueta"] == "fría" and elegidas[0]["miembros"] == ["m0", "m1", "m2"]


def test_sin_planes_factibles_informa_el_de_menor_temperatura():
    planes = [PlanObra(n, [Etapa((1,))]) for n in ("rapido", "lento")]
    esc = _esc(planes)
    filas = [{"plan": p, "curva": c, "clima": k, "I_T": 1, "I_D": 1, "I_F": 1, "M_T": mt, "M_D": 80,
              "pico_al_final": 0, "residuo_energia": 0}
             for p, mt in (("rapido", 120), ("lento", 90)) for c in ("A", "B") for k in ("k1", "k2")]
    d = decision(agregar(filas, esc))
    assert d["propuesto"] is None and d["menor_riesgo"] == "lento"  # empate en P_f = 100 %: decide M_T
