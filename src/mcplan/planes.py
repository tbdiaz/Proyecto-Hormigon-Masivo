"""Planes de obra (E13/F02): partición en etapas, cronograma, duración D y activación por capas.

Un plan u define:
- **Etapas:** partición del hormigón en vertidas de obra. Cada etapa agrupa vertidas del IFC
  (p. ej. todo junto = vertida única) o toma una franja de altura de una de ellas.
- **Esperas:** el tiempo entre el fin de una etapa y el inicio de la siguiente. La primera
  espera se cuenta desde t0: postergar el inicio también suma a D.
- **Colocación:** temperatura de colocación, moldaje y aislación.

**Recursos (iguales para todos los planes).** Una sola cuadrilla y bomba: las etapas no se
traslapan. La duración física de una etapa es d_j = V_j / r, con V_j del IFC y r la tasa de
colocación.

**Objetivo.** D(u) = max_j t_j,fin - t0 = fin de la última etapa, sin el enfriamiento posterior.

**Activación por capas.** El modelo no hace aparecer la etapa entera de golpe: sube de abajo
hacia arriba a la tasa r. Cada capa de celdas aparece cuando el volumen colocado llega a su
altura, a la temperatura de colocación y con edad 0. Así la vertida única (todo el volumen en
una operación continua) es una columna que crece durante d = V/r, y no un bloque instantáneo.
"""

from dataclasses import dataclass, field

import numpy as np

from .simulacion3d import HORMIGON, Activacion


@dataclass
class Etapa:
    vertidas: tuple  # vertidas del IFC que forman la etapa
    espera_h: float = 0.0  # desde el fin de la etapa anterior (la primera: desde t0)
    T_colocacion: float = 18.0
    z: tuple | None = None  # (z_min, z_max): sólo esa franja de altura de las vertidas (dividir una vertida)
    nombre: str = ""


@dataclass
class PlanObra:
    nombre: str
    etapas: list
    descripcion: str = ""
    cubiertas: dict = field(default_factory=dict)  # cambios a `moldaje`/`aislacion` del caso para este plan

    def resumen(self) -> str:
        partes = []
        for e in self.etapas:
            v = "+".join(f"V{x}" for x in e.vertidas)
            if e.z:
                v += f"[{e.z[0]:g}-{e.z[1]:g} m]"
            partes.append(f"{v}" + (f" (espera {e.espera_h:g} h)" if e.espera_h else ""))
        return " → ".join(partes)


def _celdas_etapa(dom, etapa: Etapa) -> np.ndarray:
    m = np.isin(dom.vertida_id, etapa.vertidas) & (dom.material == HORMIGON)
    if etapa.z is not None:
        z = dom.centros(2)[None, None, :]
        m &= (z >= etapa.z[0]) & (z < etapa.z[1])
    return m


def volumen_etapa(dom, vox, etapa: Etapa) -> float:
    """Volumen del IFC de la etapa [m³]. Con una franja z se usa la fracción de celdas de cada vertida."""
    V = 0.0
    for v in etapa.vertidas:
        if etapa.z is None:
            V += vox.volumen_ifc[v - 1]
        else:
            de_v = (dom.vertida_id == v) & (dom.material == HORMIGON)
            z = dom.centros(2)[None, None, :]
            en_z = de_v & (z >= etapa.z[0]) & (z < etapa.z[1])
            V += vox.volumen_ifc[v - 1] * np.count_nonzero(en_z) / max(np.count_nonzero(de_v), 1)
    return V


def cronograma(plan: PlanObra, dom, vox, tasa_m3_h: float) -> dict:
    """Inicio, duración y fin de cada etapa [h desde t0], y D. Verifica que la partición cubra todo una vez."""
    if tasa_m3_h <= 0:
        raise ValueError("La tasa de colocación debe ser positiva")
    cubierto = np.zeros(dom.forma, dtype=np.int16)
    filas, t = [], 0.0
    for n, e in enumerate(plan.etapas, start=1):
        celdas = _celdas_etapa(dom, e)
        if not celdas.any():
            raise ValueError(f"{plan.nombre}: la etapa {n} no tiene celdas")
        cubierto += celdas
        V = volumen_etapa(dom, vox, e)
        inicio = t + e.espera_h
        d = V / tasa_m3_h
        filas.append({"etapa": n, "nombre": e.nombre or f"E{n}", "vertidas": list(e.vertidas), "z": e.z,
                      "volumen_m3": V, "espera_h": e.espera_h, "inicio_h": inicio, "duracion_h": d,
                      "fin_h": inicio + d, "T_colocacion": e.T_colocacion})
        t = inicio + d
    hormigon = dom.material == HORMIGON
    if np.any(cubierto[hormigon] != 1):
        faltan, dobles = int(np.sum(cubierto[hormigon] == 0)), int(np.sum(cubierto[hormigon] > 1))
        raise ValueError(f"{plan.nombre}: la partición no cubre el hormigón una sola vez ({faltan} celdas sin etapa, "
                         f"{dobles} en dos etapas)")
    return {"etapas": filas, "D_h": t, "volumen_m3": sum(f["volumen_m3"] for f in filas), "tasa_m3_h": tasa_m3_h}


def activacion_por_capas(plan: PlanObra, dom, vox, tasa_m3_h: float) -> tuple[Activacion, dict]:
    """Hora de colocación de cada celda: cada etapa sube por capas de celdas a la tasa r."""
    crono = cronograma(plan, dom, vox, tasa_m3_h)
    t_act = np.full(dom.forma, np.inf)
    T_col = np.full(dom.forma, np.nan)
    for e, f in zip(plan.etapas, crono["etapas"]):
        celdas = _celdas_etapa(dom, e)
        por_capa = celdas.sum(axis=(0, 1))  # celdas de la etapa en cada capa k
        acumulado = np.concatenate([[0], np.cumsum(por_capa)[:-1]]) / por_capa.sum()  # fracción bajo la capa
        hora_capa = f["inicio_h"] + f["duracion_h"] * acumulado
        t_act = np.where(celdas, hora_capa[None, None, :], t_act)
        T_col = np.where(celdas, e.T_colocacion, T_col)
    nv = vox.n_vertidas
    inicio = np.full(nv, np.inf)
    fin = np.full(nv, np.inf)
    T_v = np.full(nv, np.nan)
    for v in range(1, nv + 1):
        m = dom.vertida_id == v
        if np.isfinite(t_act[m]).any():
            inicio[v - 1] = float(t_act[m].min())
            T_v[v - 1] = float(T_col[m][np.argmin(t_act[m])])
            # La vertida queda colocada cuando termina la etapa que contiene su última celda.
            fin[v - 1] = max(f["fin_h"] for e, f in zip(plan.etapas, crono["etapas"]) if v in e.vertidas)
    return Activacion(t_act=t_act, T_col=T_col, inicio_h=inicio, fin_h=fin, T_col_v=T_v), crono


def plan_desde_dict(d: dict, T_colocacion: float) -> PlanObra:
    etapas = []
    for e in d["etapas"]:
        if isinstance(e, (list, tuple)):
            e = {"vertidas": e}
        etapas.append(Etapa(vertidas=tuple(int(v) for v in e["vertidas"]), espera_h=float(e.get("espera_h", 0.0)),
                            T_colocacion=float(e.get("T_colocacion", T_colocacion)),
                            z=None if e.get("z") is None else tuple(float(x) for x in e["z"]),
                            nombre=str(e.get("nombre", ""))))
    return PlanObra(nombre=str(d["nombre"]), etapas=etapas, descripcion=str(d.get("descripcion", "")),
                    cubiertas=dict(d.get("cubiertas") or {}))


def expandir_familias(planes_cfg: list, T_colocacion: float) -> list[PlanObra]:
    """Planes explícitos y familias: {familia: nombre, etapas: [...], esperas_h: [0, 6, ...]} da un plan por espera."""
    planes = []
    for d in planes_cfg:
        if "familia" in d:
            for w in d["esperas_h"]:
                etapas = []
                for n, e in enumerate(d["etapas"]):
                    e = {"vertidas": e} if isinstance(e, (list, tuple)) else dict(e)
                    if n > 0:
                        e["espera_h"] = w
                    etapas.append(e)
                planes.append(plan_desde_dict({"nombre": f"{d['familia']}-e{w:g}", "etapas": etapas,
                                               "descripcion": f"{d.get('descripcion', d['familia'])}, espera {w:g} h",
                                               "cubiertas": d.get("cubiertas")}, d.get("T_colocacion", T_colocacion)))
        else:
            planes.append(plan_desde_dict(d, d.get("T_colocacion", T_colocacion)))
    nombres = [p.nombre for p in planes]
    if len(set(nombres)) != len(nombres):
        raise ValueError(f"Hay planes con el mismo nombre: {nombres}")
    return planes
