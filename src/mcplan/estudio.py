"""Estudio de decisión (E13/F02): planes × curvas de hidratación × historias climáticas.

Para cada plan u y escenario s = (curva m, historia k), con peso w = p_m q_k:
    M_T(u,s) = max_t max_{x en Ω_c(t)} T            (sólo hormigón colocado)
    M_Δ(u,s) = max_t [max_{Ω_c(t)} T - min_{S25(t)} T]  (mismo instante; S25: 25 mm bajo superficies existentes)
    F = {M_T > 70} ∪ {M_Δ > 20};   P_f(u) = Σ w 1[F],  y también P_T, P_Δ marginales.
El plan recomendado es el de menor D(u) con P_f <= 0,05 entre los planes evaluados.

Las corridas se guardan fila a fila en escenarios.csv (y su serie en series/), así que el
estudio se puede interrumpir y retomar: sólo corre lo que falta para la misma configuración.
"""

import csv
import hashlib
import json
import math
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path

import numpy as np
import yaml

LIMITE_T, LIMITE_D = 70.0, 20.0
EXTENSION_H = 168.0  # si un máximo cae al final del horizonte, se extiende una vez esta cantidad de horas
COLUMNAS = ["plan", "D_h", "curva", "clima", "w", "p_curva", "q_clima", "M_T", "t_MT_h", "x_MT", "y_MT", "z_MT",
            "M_D", "t_MD_h", "Tmax_en_MD", "Tmin_en_MD", "x_Tmin", "y_Tmin", "z_Tmin", "M_T_14d", "M_D_14d",
            "M_D_aristas", "I_T", "I_D", "I_F", "horizonte_h", "pico_al_final", "residuo_energia", "dx", "dt_s",
            "segundos", "firma"]


# ---------------------------------------------------------------------------
# Configuración
# ---------------------------------------------------------------------------

def cargar_estudio(ruta) -> dict:
    ruta = Path(ruta)
    cfg = yaml.safe_load(ruta.read_text(encoding="utf-8"))
    base = ruta.parent
    cfg["_ruta"] = str(ruta)
    cfg["_caso"] = str((base / cfg["caso"]).resolve())
    cfg["_curvas"] = str((base / cfg["curvas"]["archivo"]).resolve())
    cfg["_salida"] = str((base / cfg.get("salida", f"../resultados/{ruta.stem}")).resolve())
    cfg["_texto"] = ruta.read_text(encoding="utf-8")
    return cfg


def _pesos_curvas(cfg, nombres) -> dict:
    p = cfg["curvas"].get("pesos", "iguales")
    if p == "iguales":
        return {n: 1.0 / len(nombres) for n in nombres}
    pesos = {n: float(p[n]) for n in nombres}
    total = sum(pesos.values())
    return {n: w / total for n, w in pesos.items()}


def preparar_escenarios(cfg, informar=print) -> dict:
    """Curvas, historias climáticas (con pesos), planes y su cronograma. Se calcula una vez, en el proceso principal."""
    from .clima import historias_ensemble, historias_era5, leer_snapshot, perfil_suelo, reducir_historias, ultimo_snapshot
    from .curvas import cargar_fuentes, huella
    from .planes import cronograma, expandir_familias
    from .simulacion3d import cargar_caso_3d, construir_dominio
    from .ifc_io import leer_vertidas
    from .voxelizacion import voxelizar

    caso = cargar_caso_3d(cfg["_caso"])
    c = cfg["curvas"]
    fuentes = cargar_fuentes(cfg["_curvas"], c.get("usar") if c.get("usar") not in (None, "todas") else None,
                             cola=c.get("cola", "exponencial"), desfase_h=float(c.get("desfase_h", 2.0)),
                             escala=float(c.get("escala", 1.0)))
    p_curva = _pesos_curvas(cfg, list(fuentes))

    T_col = float(cfg.get("T_colocacion", 18.0))
    planes = expandir_familias(cfg["planes"], T_col)
    tasa = float(cfg["recursos"]["tasa_colocacion_m3_h"])
    dx_ref = float(cfg.get("dx", caso.dx))
    vox = voxelizar(leer_vertidas(caso.ifc), dx_ref)
    dom = construir_dominio(vox, caso.profundidad_suelo, caso.margen_suelo)
    cronos = {p.nombre: cronograma(p, dom, vox, tasa) for p in planes}
    H14 = float(cfg.get("horizonte_h", 336))
    post = float(cfg.get("despues_ultima_colocacion_h", 240))
    # El clima cubre además una extensión de una semana por si algún máximo cae al final del horizonte.
    horas = max([H14] + [cr["D_h"] + post for cr in cronos.values()]) + EXTENSION_H

    cl = cfg["clima"]
    carpeta = Path(cfg["_caso"]).parent / "../data/clima/snapshots"
    era5 = leer_snapshot(ultimo_snapshot("reanalisis", carpeta) if cl.get("climatologia", "ultimo") == "ultimo"
                         else Path(cfg["_ruta"]).parent / cl["climatologia"])
    t0 = cfg["t0"]
    if cl.get("fuente", "ensemble") == "ensemble":
        ruta_ens = ultimo_snapshot("ensemble", carpeta) if cl.get("ensemble", "ultimo") == "ultimo" else (
            Path(cfg["_ruta"]).parent / cl["ensemble"])
        ens = leer_snapshot(ruta_ens)
        todas = historias_ensemble(ens, era5, t0, horas, float(cl.get("transicion_h", 24)), int(cl.get("ventana_dias", 3)))
        origen_clima = {"fuente": "ensemble", "archivo": Path(ens["archivo"]).name,
                        "descargado_utc": ens["metadata"]["descargado_utc"], "modelo": ens["metadata"]["modelo"],
                        "sha256": ens["metadata"]["sha256_respuesta"], "miembros": len(todas),
                        "cubre_hasta": ens["metadata"]["horizonte"]["ultima_hora_con_dato"]["temperature_2m"]}
    else:
        todas = historias_era5(era5, t0, horas)
        origen_clima = {"fuente": "era5_bloques", "anios": era5["metadata"]["anios"]}
    origen_clima["climatologia"] = Path(era5["archivo"]).name
    elegidas = reducir_historias(todas, int(cl.get("historias", 3)), float(cl.get("dias_indice", 7)))
    for h in elegidas:
        h["clave"] = h["etiqueta"]
    T_z, T_fondo, info_suelo = perfil_suelo(era5, t0, int(cl.get("ventana_dias", 3)))
    suelo = {"prof": info_suelo["profundidades_m"], "T": info_suelo["T_capas_C"], "fondo": T_fondo}
    informar(f"Curvas: {len(fuentes)} ({', '.join(fuentes)}); historias: {len(elegidas)} de {len(todas)} "
             f"({origen_clima['fuente']}); planes: {len(planes)}; horizonte de simulación {horas:g} h")
    # Firma de todo lo que afecta a una corrida, salvo el plan y la malla (que van en la firma de cada tarea):
    # agregar o quitar planes del YAML no obliga a repetir las corridas ya hechas.
    afecta = {k: cfg.get(k) for k in ("t0", "horizonte_h", "despues_ultima_colocacion_h", "guardar_cada_h", "dt_max_s",
                                      "profundidad_superficie_m", "T_colocacion")}
    afecta.update(recursos=cfg["recursos"], curvas={k: v for k, v in cfg["curvas"].items() if k != "archivo"},
                  clima=cfg["clima"])
    firma_base = hashlib.sha256("|".join([json.dumps(afecta, sort_keys=True, default=str),
                                          Path(cfg["_caso"]).read_text(encoding="utf-8"), huella(cfg["_curvas"]),
                                          json.dumps(origen_clima, sort_keys=True)]).encode()).hexdigest()[:16]
    return {"caso": caso, "fuentes": fuentes, "p_curva": p_curva, "planes": planes, "cronos": cronos, "tasa": tasa,
            "historias": elegidas, "todas": todas, "origen_clima": origen_clima, "horas": horas, "H14": H14,
            "suelo": suelo, "info_suelo": info_suelo, "firma": firma_base, "t0": t0}


# ---------------------------------------------------------------------------
# Una corrida (se ejecuta en un proceso de trabajo)
# ---------------------------------------------------------------------------

_BASE = {}


def _base(caso_ruta, dx):
    """Caso, voxelización y dominio, una vez por proceso y por dx."""
    clave = (caso_ruta, dx)
    if clave not in _BASE:
        from .ifc_io import leer_vertidas
        from .simulacion3d import cargar_caso_3d, construir_dominio
        from .voxelizacion import voxelizar

        caso = cargar_caso_3d(caso_ruta)
        vox = voxelizar(leer_vertidas(caso.ifc), dx)
        dom = construir_dominio(vox, caso.profundidad_suelo, caso.margen_suelo)
        _BASE[clave] = (caso, vox, dom)
    return _BASE[clave]


_FUENTES = {}


def simular_tarea(tarea: dict, instantes_3d_h=()):
    """La simulación de una tarea: devuelve (res, cronograma). res.plan queda con la activación del plan."""
    import copy

    from .bordes import Ambiente, Clima
    from .curvas import cargar_fuentes
    from .planes import PlanObra, activacion_por_capas, Etapa
    from .simulacion3d import armar_cubiertas, simular_3d
    from .voxelizacion import PlanVertidas

    caso, vox, dom = _base(tarea["caso"], tarea["dx"])
    caso = copy.deepcopy(caso)
    plan = PlanObra(tarea["plan"]["nombre"], [Etapa(**e) for e in tarea["plan"]["etapas"]],
                    cubiertas=tarea["plan"]["cubiertas"])
    for tipo, cambios in plan.cubiertas.items():
        caso.cubiertas_cfg[tipo] = {**caso.cubiertas_cfg.get(tipo, {}), **cambios}
    cubiertas = armar_cubiertas(caso, vox)
    act, crono = activacion_por_capas(plan, dom, vox, tarea["tasa"])
    clave = (tarea["curvas"], tarea["curva"], json.dumps(tarea["opciones_curvas"], sort_keys=True))
    if clave not in _FUENTES:
        _FUENTES[clave] = cargar_fuentes(tarea["curvas"], [tarea["curva"]], **tarea["opciones_curvas"])[tarea["curva"]]
    fuente = _FUENTES[clave]
    c = tarea["clima"]
    clima = Clima(np.array(c["t_h"]), np.array(c["T_aire"]), np.array(c["viento"]), np.array(c["G"]),
                  np.array(c["nubosidad"]))
    sup = caso.superficie_cfg
    amb = Ambiente(20.0, float(sup.get("h_superior", 27.8)), float(sup.get("h_lateral", 27.8)), tarea["suelo"]["fondo"],
                   modelo=sup.get("modelo", "fisico"), clima=clima, absortividad=float(sup.get("absortividad", 0.65)),
                   emisividad=float(sup.get("emisividad", 0.9)))
    prof, Tz = np.array(tarea["suelo"]["prof"]), np.array(tarea["suelo"]["T"])
    H = max(tarea["H14"], crono["D_h"] + tarea["post_h"]) + tarea.get("extra_h", 0.0)
    res = simular_3d(dom, PlanVertidas(vox), caso.hormigon, caso.suelo, fuente, amb, H,
                     T_suelo_inicial=lambda z: np.interp(-np.asarray(z), prof, Tz), guardar_cada_h=tarea["guardar_h"],
                     dt_max_s=tarea["dt_s"], cubiertas=cubiertas, profundidad_tmin=tarea["profundidad"],
                     superficie_tmin="plana", activacion=act, puntos=caso.puntos, instantes_3d_h=instantes_3d_h)
    res.plan = act  # para las figuras: qué vertidas existen en cada instante según este plan
    return res, crono


def _pico_al_final(res) -> bool:
    n = len(res.t_h)
    return bool(int(np.nanargmax(res.Tmax)) >= n - 2 or int(np.nanargmax(res.dT)) >= n - 2)


def _firma_plan(plan) -> str:
    texto = json.dumps({"etapas": [asdict(e) for e in plan.etapas], "cubiertas": plan.cubiertas}, sort_keys=True,
                       default=str)
    return hashlib.sha256(texto.encode()).hexdigest()[:8]


def correr_escenario(tarea: dict) -> dict:
    """Simula un plan con una curva y una historia climática; devuelve la fila de métricas y la serie."""
    t_ini = time.perf_counter()
    res, crono = simular_tarea(tarea)
    if _pico_al_final(res):  # no esconder un pico fuera de la ventana: se extiende y se informan ambos horizontes
        res, crono = simular_tarea({**tarea, "extra_h": EXTENSION_H})
    H = res.t_h[-1]
    t = res.t_h
    en14 = t <= tarea["H14"] + 1e-9
    i = int(np.nanargmax(res.Tmax))
    j = int(np.nanargmax(res.dT))
    pico_final = _pico_al_final(res)
    MT, MD = float(res.Tmax[i]), float(res.dT[j])
    fila = {
        "plan": tarea["plan"]["nombre"], "D_h": crono["D_h"], "curva": tarea["curva"], "clima": tarea["clima_clave"],
        "w": tarea["w"], "p_curva": tarea["p"], "q_clima": tarea["q"],
        "M_T": MT, "t_MT_h": float(t[i]), "x_MT": res.Tmax_pos[i][0], "y_MT": res.Tmax_pos[i][1],
        "z_MT": res.Tmax_pos[i][2], "M_D": MD, "t_MD_h": float(t[j]), "Tmax_en_MD": float(res.Tmax[j]),
        "Tmin_en_MD": float(res.Tmin_sup[j]), "x_Tmin": res.Tmin_pos[j][0], "y_Tmin": res.Tmin_pos[j][1],
        "z_Tmin": res.Tmin_pos[j][2], "M_T_14d": float(np.nanmax(res.Tmax[en14])),
        "M_D_14d": float(np.nanmax(res.dT[en14])), "M_D_aristas": float(np.nanmax(res.Tmax - res.Tmin_todas)),
        "I_T": int(MT > LIMITE_T), "I_D": int(MD > LIMITE_D), "I_F": int(MT > LIMITE_T or MD > LIMITE_D),
        "horizonte_h": float(H), "pico_al_final": int(pico_final), "residuo_energia": res.balance["residuo_relativo"],
        "dx": tarea["dx"], "dt_s": tarea["dt_s"], "segundos": time.perf_counter() - t_ini, "firma": tarea["firma"],
    }
    serie = {"t_h": t, "Tmax": res.Tmax, "Tmin": res.Tmin_sup, "dT": res.dT, "T_amb": res.T_amb,
             "nucleo": res.puntos.get("nucleo_eje", np.full(t.shape, np.nan))}
    return {"fila": fila, "serie": serie}


def _clima_dict(clima) -> dict:
    return {"t_h": clima.t_h.tolist(), "T_aire": clima.T_aire.tolist(), "viento": clima.viento.tolist(),
            "G": clima.G.tolist(), "nubosidad": clima.nubosidad.tolist()}


def tareas(cfg, esc, dx, planes=None) -> list[dict]:
    cl = cfg["curvas"]
    opciones = {"cola": cl.get("cola", "exponencial"), "desfase_h": float(cl.get("desfase_h", 2.0))}
    if float(cl.get("escala", 1.0)) != 1.0:
        opciones["escala"] = float(cl["escala"])
    lista = []
    for plan in esc["planes"]:
        if planes and plan.nombre not in planes:
            continue
        for curva in esc["fuentes"]:
            for h in esc["historias"]:
                lista.append({
                    "caso": cfg["_caso"], "dx": dx, "plan": {"nombre": plan.nombre,
                                                             "etapas": [asdict(e) for e in plan.etapas],
                                                             "cubiertas": plan.cubiertas},
                    "tasa": esc["tasa"], "curvas": cfg["_curvas"], "curva": curva, "opciones_curvas": opciones,
                    "clima": _clima_dict(h["clima"]), "clima_clave": h["clave"], "p": esc["p_curva"][curva],
                    "q": h["peso"], "w": esc["p_curva"][curva] * h["peso"], "suelo": esc["suelo"],
                    "H14": esc["H14"], "post_h": float(cfg.get("despues_ultima_colocacion_h", 240)),
                    "guardar_h": float(cfg.get("guardar_cada_h", 1.0)), "dt_s": float(cfg.get("dt_max_s", 600)),
                    "profundidad": float(cfg.get("profundidad_superficie_m", 0.025)),
                    "firma": f"{esc['firma']}-{_firma_plan(plan)}-{dx:g}",
                })
    return lista


def _clave(f) -> tuple:
    return (f["plan"], f["curva"], f["clima"], f["firma"])


def leer_filas(ruta) -> list[dict]:
    if not Path(ruta).exists():
        return []
    with open(ruta, encoding="utf-8") as f:
        filas = list(csv.DictReader(f))
    for r in filas:
        for k in COLUMNAS:
            if k not in ("plan", "curva", "clima", "firma"):
                r[k] = float(r[k])
    return filas


def ejecutar(cfg, esc, dx, procesos=None, planes=None, informar=print) -> list[dict]:
    """Corre las tareas que faltan (en paralelo) y agrega cada fila a escenarios.csv apenas termina."""
    salida = Path(cfg["_salida"])
    (salida / "series").mkdir(parents=True, exist_ok=True)
    ruta = salida / "escenarios.csv"
    hechas = {_clave(f) for f in leer_filas(ruta)}
    pendientes = [t for t in tareas(cfg, esc, dx, planes)
                  if (t["plan"]["nombre"], t["curva"], t["clima_clave"], t["firma"]) not in hechas]
    total = len(tareas(cfg, esc, dx, planes))
    informar(f"{total - len(pendientes)} de {total} corridas ya hechas; faltan {len(pendientes)} (dx = {dx:g} m)")
    if not pendientes:
        return leer_filas(ruta)
    procesos = procesos or max(1, (os.cpu_count() or 2) - 1)
    nuevo = not ruta.exists()
    t0 = time.perf_counter()
    with open(ruta, "a", newline="", encoding="utf-8") as f, ProcessPoolExecutor(procesos) as pool:
        w = csv.DictWriter(f, COLUMNAS)
        if nuevo:
            w.writeheader()
        futuros = {pool.submit(correr_escenario, t): t for t in pendientes}
        for n, fut in enumerate(as_completed(futuros), start=1):
            r = fut.result()
            fila = r["fila"]
            w.writerow({k: fila[k] for k in COLUMNAS})
            f.flush()
            nombre = f"{fila['plan']}__{fila['curva']}__{fila['clima']}__{fila['firma']}".replace(" ", "_")
            np.savez_compressed(salida / "series" / f"{nombre}.npz", **r["serie"])
            resta = (time.perf_counter() - t0) / n * (len(pendientes) - n)
            informar(f"  [{n}/{len(pendientes)}] {fila['plan']:<12} {fila['curva']:<7} {fila['clima']:<7} "
                     f"M_T {fila['M_T']:6.1f}  M_Δ {fila['M_D']:6.1f}  falla {fila['I_F']}  "
                     f"({fila['segundos']:.0f} s; quedan ~{resta / 60:.0f} min)")
    return leer_filas(ruta)


# ---------------------------------------------------------------------------
# Agregación y decisión
# ---------------------------------------------------------------------------

def agregar(filas, esc, firma=None, pesos_curva=None, pesos_clima=None) -> list[dict]:
    """P_T, P_Δ y P_f por plan (suma de pesos de las corridas que fallan), con la opción de cambiar pesos."""
    if firma:
        filas = [f for f in filas if f["firma"] == firma]
    pc = pesos_curva or esc["p_curva"]
    qc = pesos_clima or {h["clave"]: h["peso"] for h in esc["historias"]}
    salida = []
    for plan in esc["planes"]:
        fs = [f for f in filas if f["plan"] == plan.nombre]
        esperadas = len(pc) * len(qc)
        completo = len({(f["curva"], f["clima"]) for f in fs}) == esperadas
        w = {(f["curva"], f["clima"]): pc.get(f["curva"], 0.0) * qc.get(f["clima"], 0.0) for f in fs}
        suma = sum(w.values())
        P = {k: sum(w[(f["curva"], f["clima"])] * f[c] for f in fs) for k, c in (("P_T", "I_T"), ("P_D", "I_D"),
                                                                                  ("P_f", "I_F"))}
        salida.append({"plan": plan.nombre, "descripcion": plan.descripcion or plan.resumen(),
                       "D_h": esc["cronos"][plan.nombre]["D_h"], **P, "suma_pesos": suma, "corridas": len(fs),
                       "completo": completo and abs(suma - 1) < 1e-9, "cumple": completo and P["P_f"] <= 0.05,
                       "M_T_max": max((f["M_T"] for f in fs), default=math.nan),
                       "M_D_max": max((f["M_D"] for f in fs), default=math.nan),
                       "M_T_min": min((f["M_T"] for f in fs), default=math.nan),
                       "M_D_min": min((f["M_D"] for f in fs), default=math.nan),
                       "pico_al_final": any(f["pico_al_final"] for f in fs),
                       "residuo_max": max((f["residuo_energia"] for f in fs), default=math.nan)})
    return salida


def decision(resumen) -> dict:
    """Propuesto = menor D que cumple; más rápido = el de D inmediatamente menor; más lento = inmediatamente mayor."""
    por_D = sorted(resumen, key=lambda r: r["D_h"])
    factibles = [r for r in por_D if r["cumple"]]
    if not factibles:
        # Sin planes factibles: el de menor P_f; si empatan (p. ej. todos al 100 %), el de menor M_T y M_Δ.
        return {"propuesto": None, "mas_rapido": por_D[0]["plan"],
                "menor_riesgo": min(por_D, key=lambda r: (r["P_f"], r["M_T_max"], r["M_D_max"], r["D_h"]))["plan"]}
    p = factibles[0]
    i = por_D.index(p)
    return {"propuesto": p["plan"], "mas_rapido": por_D[i - 1]["plan"] if i > 0 else None,
            "mas_lento": next((r["plan"] for r in por_D[i + 1:]), None)}
