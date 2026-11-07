"""Panel interactivo de una corrida 3D: un solo archivo HTML, sin internet ni dependencias.

Muestra la discretización (dx, dy, dz, dt, dt crítico), las condiciones de borde con
su beta = h/k, las etapas de colocación, el error de volumen, una vista 3D que se
gira, cortes con lectura de T bajo el puntero y las curvas T(t).

Guardar los campos 3D completos de cada instante no cabe en memoria con mallas
finas, así que `Recolector` se engancha a simular_3d (observador) y se queda sólo
con una copia submuestreada. Las métricas (Tmax, dT, curvas) salen del campo completo.
"""

import base64
import gzip
import json
from pathlib import Path

import numpy as np

from .simulacion3d import HORMIGON, SUELO
from .bordes import filas_bordes

_ESCALA_T = 20  # el campo se guarda como int16 en pasos de 0,05 °C
_SIN_DATO = -32768  # aire, o vertida que todavía no se coloca
_MAX_VALORES = 6_000_000  # tope de valores (celdas x instantes) para que el archivo no pese demasiado


def instantes_panel(duracion_h, colocaciones, celdas):
    """Instantes del deslizador: más densos al principio, que es cuando pasa casi todo."""
    tramos = ((24.0, 1.0), (48.0, 2.0), (120.0, 6.0), (float("inf"), 24.0))
    t, instantes = 0.0, [0.0]
    for hasta, paso in tramos:
        while t + paso <= min(hasta, duracion_h) + 1e-9:
            t += paso
            instantes.append(t)
    salto = max(1, int(np.ceil(len(instantes) * celdas / _MAX_VALORES)))
    elegidos = set(instantes[::salto])
    elegidos |= {float(c) for c in colocaciones if 0 <= c <= duracion_h}
    elegidos.add(float(duracion_h))
    return sorted(elegidos)


class Recolector:
    """Guarda el campo T submuestreado (una celda de cada `s` por eje) en los instantes del panel."""

    def __init__(self, dom, duracion_h, colocaciones, dx_objetivo=0.3):
        self.s = max(1, int(round(dx_objetivo / dom.dx)))
        self.corte = (slice(self.s // 2, None, self.s),) * 3
        self.vertida_id = dom.vertida_id[self.corte].copy()
        self.forma = self.vertida_id.shape
        self.instantes = instantes_panel(duracion_h, colocaciones, int(np.prod(self.forma)))
        self.origen = [float(dom.origen[e] + (self.s // 2) * dom.dx) for e in range(3)]  # esquina de la 1.ª muestra
        self.d = dom.dx * self.s
        self.tiempos, self.campos = [], []

    def __call__(self, t_h, T, activas):
        if not any(abs(t_h - ti) < 1e-6 for ti in self.instantes):
            return
        muestra = np.where(activas[self.corte], T[self.corte], np.nan)
        self.tiempos.append(float(t_h))
        self.campos.append(np.where(np.isfinite(muestra), np.round(np.nan_to_num(muestra) * _ESCALA_T),
                                    _SIN_DATO).astype("<i2"))


def _gz64(arreglo) -> str:
    return base64.b64encode(gzip.compress(arreglo.tobytes(), 6)).decode("ascii")


def _lista(a):
    """Arreglo a lista JSON: NaN e infinito pasan a null (JSON no los admite)."""
    return [None if not np.isfinite(x) else round(float(x), 3) for x in np.asarray(a, dtype=float)]


def _areas(dom, vox, plan, caso):
    """Área [m²] de cada tipo de superficie con todas las vertidas del plan colocadas (sin contar capas)."""
    from .simulacion3d import _caras_expuestas
    from .bordes import LAT, SUP

    final = (dom.material == SUELO) | np.isin(dom.vertida_id, plan.vertidas_activas(caso.duracion_h))
    celda, _, orientacion = _caras_expuestas(final)
    mat = dom.material.reshape(-1)[celda]
    vid = dom.vertida_id.reshape(-1)[celda]
    a = dom.dx**2
    con_moldaje = set()
    for c in caso.cubiertas.cubiertas:
        if c.tipo == "moldaje":
            con_moldaje |= set(c.vertidas)
    lat_h = (mat == HORMIGON) & (orientacion == LAT)
    moldeadas = lat_h & np.isin(vid, list(con_moldaje))
    n_suelo = int(np.count_nonzero((dom.material == SUELO).any(axis=(0, 1))))
    return {
        "sup": float(np.count_nonzero((mat == HORMIGON) & (orientacion != LAT)) * a),
        "lat_mol": float(np.count_nonzero(moldeadas) * a), "lat": float(np.count_nonzero(lat_h & ~moldeadas) * a),
        "suelo_sup": float(np.count_nonzero((mat == SUELO) & (orientacion == SUP)) * a),
        "fondo": float(dom.forma[0] * dom.forma[1] * a) if n_suelo else None,
        "suelo_lat": float(2 * (dom.forma[0] + dom.forma[1]) * n_suelo * a) if n_suelo else None,
    }


def _limpio(x):
    """Listas y números para JSON: NaN e infinito pasan a null."""
    if isinstance(x, dict):
        return {k: _limpio(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_limpio(v) for v in x]
    if isinstance(x, (float, np.floating)):
        return float(x) if np.isfinite(x) else None
    if isinstance(x, np.integer):
        return int(x)
    return x


def datos_panel(res, caso, vox, plan, registro, rec, comando) -> dict:
    dom = res.dominio
    i, j = int(np.nanargmax(res.Tmax)), int(np.nanargmax(res.dT))
    k, h_c = caso.hormigon, caso.hidratacion
    vertidas = []
    for n, nombre in enumerate(vox.nombres):
        col = res.Tmax_vertida[:, n]
        hay = bool(np.isfinite(col).any())
        ini = plan.inicio_h[n]
        vertidas.append({
            "id": n + 1, "nombre": nombre, "inicio_h": float(ini) if np.isfinite(ini) else None,
            "T_colocacion": float(plan.T_colocacion[n]), "celdas": vox.celdas(n + 1), "V_vox": vox.volumen_voxel(n + 1),
            "V_ifc": vox.volumen_ifc[n], "error": vox.error_volumen(n + 1),
            "Tmax": float(np.nanmax(col)) if hay else None,
            "t_Tmax": float(res.t_h[int(np.nanargmax(col))]) if hay else None,
        })
    dt_de = lambda t: next((d[1] for d in reversed(res.dts) if d[0] <= t + 1e-9), None)  # noqa: E731
    return {
        "caso": caso.nombre,
        "comando": comando,
        "cambios": registro,
        "malla": {"dx": dom.dx, "dy": dom.dx, "dz": dom.dx, "nx": dom.forma[0], "ny": dom.forma[1], "nz": dom.forma[2],
                  "celdas": int(np.prod(dom.forma)), "hormigon": int(np.count_nonzero(dom.material == HORMIGON)),
                  "suelo": int(np.count_nonzero(dom.material == SUELO)),
                  "error_total": vox.error_volumen_total, "superposiciones": vox.superposiciones},
        "tiempo": {"duracion_h": caso.duracion_h, "hora_inicio": caso.hora_inicio, "dt": res.dts[-1][1],
                   "dt_max": caso.dt_max_s, "seguridad": caso.factor_seguridad, "guardar_cada_h": caso.guardar_cada_h},
        "bordes": _limpio(filas_bordes(caso.ambiente, caso.cubiertas, dom.dx, caso.hormigon.k, caso.suelo.k,
                                       _areas(dom, vox, plan, caso))),
        "ambiente": {"descripcion": caso.ambiente.describir(), "T_amb": caso.T_amb_texto,
                     "cubiertas": caso.cubiertas.describir(plan.inicio_h), "clima": _limpio(caso.clima_info or None),
                     "fecha_inicio": caso.fecha_inicio, "criterio_Tmin": res.criterio_Tmin},
        "retiros": [{"t": t, "evento": d} for t, d in res.retiros],
        "materiales": [
            {"nombre": nombre, "k": m.k, "rho": m.rho, "cp": m.cp, "alpha": m.k / m.capacidad}
            for nombre, m in (("Hormigón", k), ("Suelo", caso.suelo))
        ],
        "hidratacion": {"Eh_R": h_c.Eh_R, "T_ref": h_c.T_ref_C, "te_quiebre": h_c.te_quiebre_h},
        "vertidas": vertidas,
        "etapas": [{"t": t, "activas": a, "dt_critico": c, "dt": dt_de(t),
                    "V": float(np.count_nonzero(np.isin(vox.vertida_id, a))) * dom.dx**3}
                   for t, a, c in res.etapas],
        "limites": {"Tmax": caso.limite_Tmax, "dT": caso.limite_dT},
        "resultado": {"Tmax": float(res.Tmax[i]), "t_Tmax": float(res.t_h[i]), "pos_Tmax": res.Tmax_pos[i].tolist(),
                      "dT": float(res.dT[j]), "t_dT": float(res.t_h[j]), "pos_Tmin": _limpio(res.Tmin_pos[j].tolist()),
                      "Tmin": _limpio(float(res.Tmin_sup[j]))},
        "series": {"t": _lista(res.t_h), "Tmax": _lista(res.Tmax), "dT": _lista(res.dT), "T_amb": _lista(res.T_amb),
                   "Tmin": _lista(res.Tmin_sup),
                   "vertidas": [_lista(res.Tmax_vertida[:, n]) for n in range(vox.n_vertidas)],
                   "puntos": [{"nombre": n, "pos": list(res.posiciones[n]), "T": _lista(T)}
                              for n, T in res.puntos.items()]},
        "colocaciones": [{"t": t, "v": v} for t, v, *_ in res.eventos],
        "vista": {"forma": list(rec.forma), "origen": rec.origen, "d": rec.d, "s": rec.s, "tiempos": rec.tiempos,
                  "escala": _ESCALA_T, "sin_dato": _SIN_DATO, "vid": _gz64(rec.vertida_id.astype(np.int8)),
                  "T": _gz64(np.stack(rec.campos))},
    }


def escribir_panel(ruta, datos: dict) -> None:
    plantilla = (Path(__file__).parent / "panel.html").read_text(encoding="utf-8")
    texto = json.dumps(datos, ensure_ascii=False, allow_nan=False).replace("</", "<\\/")
    Path(ruta).write_text(plantilla.replace("__DATOS_PANEL__", texto), encoding="utf-8")
