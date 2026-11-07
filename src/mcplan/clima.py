"""Clima real para las condiciones de borde (E10/E11): Open-Meteo -> snapshots -> serie horaria.

Tres tipos de dato, que no son lo mismo:

- **Pronóstico (forecast)**: lo que un modelo (aquí ECMWF IFS 0,25°) predice desde la hora
  en que se emitió. Cambia con cada emisión, por eso se guarda cada descarga.
- **Reanálisis (ERA5)**: el pasado reconstruido por un modelo que asimila observaciones;
  grilla regular y sin huecos. Se usa como historia y para la climatología.
- **Observación**: lo que mide una estación o un sensor en un punto. No la usamos todavía
  (no hay estación en el sitio); serviría para calibrar o validar.

Snapshots: cada descarga se guarda tal como llegó, en un archivo nuevo de
data/clima/snapshots/ con su metadata (fuente, modelo, fecha y hora de descarga, coordenadas
pedidas y del punto de grilla, unidades, variables, horizonte y SHA-256 de la respuesta).
Nunca se sobrescribe un snapshot: si el nombre ya existe, la descarga falla. El índice
data/clima/indice.csv crece con una línea por descarga.

Horizonte: un pronóstico llega a 15-16 días desde que se emite. Para la ventana de
la construcción (14 días desde la primera vertida) se usa el pronóstico donde exista; la
anomalía de su última hora decae linealmente a cero en `transicion_h`, y el resto es
climatología horaria ERA5 de los mismos días de años anteriores.
"""

import csv
import hashlib
import json
import time
import urllib.parse
import urllib.request
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np

from .bordes import Clima

API_PRONOSTICO = "https://api.open-meteo.com/v1/forecast"
API_ERA5 = "https://archive-api.open-meteo.com/v1/archive"
ZONA = "America/Santiago"
VARIABLES_PRONOSTICO = ("temperature_2m", "wind_speed_10m", "shortwave_radiation", "cloud_cover",
                        "relative_humidity_2m", "soil_temperature_54cm")
CAPAS_SUELO_ERA5 = {  # variable -> profundidad media de la capa [m]
    "soil_temperature_0_to_7cm": 0.035, "soil_temperature_7_to_28cm": 0.175,
    "soil_temperature_28_to_100cm": 0.64, "soil_temperature_100_to_255cm": 1.775,
}
VARIABLES_ERA5 = ("temperature_2m", "wind_speed_10m", "shortwave_radiation", "cloud_cover", *CAPAS_SUELO_ERA5)
CARPETA = Path("data/clima/snapshots")
INDICE = Path("data/clima/indice.csv")
_COLUMNAS_INDICE = ("archivo", "tipo", "fuente", "modelo", "descargado_utc", "horizonte_desde", "horizonte_hasta",
                    "latitud", "longitud", "sha256")


# ---------------------------------------------------------------------------
# Descarga y snapshots
# ---------------------------------------------------------------------------

def _huella(obj) -> str:
    """SHA-256 del JSON canónico: detecta si alguien editó la respuesta guardada."""
    texto = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(texto.encode("utf-8")).hexdigest()


def _pedir(url, params, intentos=3):
    completa = f"{url}?{urllib.parse.urlencode(params)}"
    pedido = urllib.request.Request(completa, headers={"User-Agent": "mcplan (curso de hormigon masivo)"})
    for n in range(intentos):
        try:
            with urllib.request.urlopen(pedido, timeout=60) as r:
                datos = json.loads(r.read().decode("utf-8"))
            if datos.get("error"):
                raise RuntimeError(f"Open-Meteo respondió con error: {datos.get('reason')}")
            return datos, completa
        except OSError:
            if n == intentos - 1:
                raise
            time.sleep(2 * (n + 1))


def _horizonte(hourly) -> dict:
    t = hourly["time"]
    con_dato = {v: next((t[i] for i in range(len(t) - 1, -1, -1) if hourly[v][i] is not None), None)
                for v in hourly if v != "time"}
    return {"desde": t[0], "hasta": t[-1], "horas": len(t), "ultima_hora_con_dato": con_dato}


def _guardar(carpeta, nombre, metadata, respuesta, indice) -> Path:
    carpeta = Path(carpeta)
    carpeta.mkdir(parents=True, exist_ok=True)
    ruta = carpeta / nombre
    documento = {"metadata": {**metadata, "sha256_respuesta": _huella(respuesta)}, "respuesta": respuesta}
    with open(ruta, "x", encoding="utf-8") as f:  # "x": falla si ya existe; un snapshot nunca se reemplaza
        json.dump(documento, f, ensure_ascii=False, indent=1)
    indice = Path(indice)
    nuevo = not indice.exists()
    indice.parent.mkdir(parents=True, exist_ok=True)
    with open(indice, "a", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        if nuevo:
            w.writerow(_COLUMNAS_INDICE)
        m = documento["metadata"]
        w.writerow([ruta.name, m["tipo"], m["fuente"], m["modelo"], m["descargado_utc"], m["horizonte"]["desde"],
                    m["horizonte"]["hasta"], m["sitio"]["latitud"], m["sitio"]["longitud"], m["sha256_respuesta"]])
    return ruta


def _sitio(ruta_sitio) -> dict:
    import yaml

    s = yaml.safe_load(Path(ruta_sitio).read_text(encoding="utf-8"))
    return {"grupo": s.get("grupo"), "nombre": s.get("sitio"), "latitud": float(s["latitud"]),
            "longitud": float(s["longitud"])}


def descargar_pronostico(ruta_sitio="config/sitio.yaml", modelo="ecmwf_ifs025", dias=16, carpeta=CARPETA,
                         indice=INDICE) -> Path:
    """Descarga el pronóstico horario vigente y lo guarda como snapshot nuevo."""
    sitio = _sitio(ruta_sitio)
    params = {"latitude": sitio["latitud"], "longitude": sitio["longitud"], "hourly": ",".join(VARIABLES_PRONOSTICO),
              "models": modelo, "forecast_days": dias, "timezone": ZONA, "wind_speed_unit": "ms"}
    ahora = datetime.now(timezone.utc)
    datos, url = _pedir(API_PRONOSTICO, params)
    metadata = {
        "tipo": "pronostico", "fuente": "Open-Meteo Forecast API", "endpoint": API_PRONOSTICO, "url": url,
        "parametros": params, "modelo": modelo, "descargado_utc": ahora.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "sitio": sitio,
        "punto_modelo": {"latitud": datos["latitude"], "longitud": datos["longitude"], "elevacion_m": datos["elevation"]},
        "zona_horaria": datos["timezone"], "utc_offset_s": datos["utc_offset_seconds"],
        "unidades": datos["hourly_units"], "variables": list(VARIABLES_PRONOSTICO),
        "horizonte": _horizonte(datos["hourly"]),
    }
    nombre = f"pronostico_{modelo}_{ahora.strftime('%Y%m%dT%H%MZ')}.json"
    return _guardar(carpeta, nombre, metadata, datos, indice)


def descargar_era5(fecha_inicio, ruta_sitio="config/sitio.yaml", dias=14, margen_dias=7, anios=range(2016, 2026),
                   carpeta=CARPETA, indice=INDICE) -> Path:
    """Reanálisis ERA5 de la misma ventana del calendario (± margen) en varios años: base de la climatología."""
    sitio = _sitio(ruta_sitio)
    inicio = datetime.fromisoformat(fecha_inicio)
    desde, hasta = inicio - timedelta(days=margen_dias), inicio + timedelta(days=dias + margen_dias)
    ahora = datetime.now(timezone.utc)
    respuestas, urls = {}, []
    for anio in anios:
        a = desde.replace(year=anio + (desde.year - inicio.year))
        b = hasta.replace(year=anio + (hasta.year - inicio.year))
        params = {"latitude": sitio["latitud"], "longitude": sitio["longitud"], "start_date": a.date().isoformat(),
                  "end_date": b.date().isoformat(), "hourly": ",".join(VARIABLES_ERA5), "models": "era5",
                  "timezone": ZONA, "wind_speed_unit": "ms"}
        datos, url = _pedir(API_ERA5, params)
        respuestas[str(anio)] = datos
        urls.append(url)
    primero = respuestas[str(anios[0])]
    metadata = {
        "tipo": "reanalisis", "fuente": "Open-Meteo Historical Weather API (ERA5)", "endpoint": API_ERA5, "urls": urls,
        "modelo": "era5", "descargado_utc": ahora.strftime("%Y-%m-%dT%H:%M:%SZ"), "sitio": sitio,
        "punto_modelo": {"latitud": primero["latitude"], "longitud": primero["longitude"],
                         "elevacion_m": primero["elevation"]},
        "zona_horaria": primero["timezone"], "utc_offset_s": primero["utc_offset_seconds"],
        "unidades": primero["hourly_units"], "variables": list(VARIABLES_ERA5), "anios": list(map(int, anios)),
        "ventana": {"fecha_inicio_obra": fecha_inicio, "dias": dias, "margen_dias": margen_dias},
        "horizonte": {"desde": primero["hourly"]["time"][0], "hasta": respuestas[str(anios[-1])]["hourly"]["time"][-1],
                      "horas": sum(len(r["hourly"]["time"]) for r in respuestas.values())},
    }
    nombre = f"era5_{inicio.strftime('%m%d')}_{anios[0]}-{anios[-1]}_{ahora.strftime('%Y%m%dT%H%MZ')}.json"
    return _guardar(carpeta, nombre, metadata, {"por_anio": respuestas}, indice)


def leer_snapshot(ruta) -> dict:
    documento = json.loads(Path(ruta).read_text(encoding="utf-8"))
    if _huella(documento["respuesta"]) != documento["metadata"]["sha256_respuesta"]:
        raise ValueError(f"{ruta}: la respuesta no coincide con su SHA-256 (el snapshot fue modificado)")
    documento["archivo"] = str(ruta)
    return documento


def ultimo_snapshot(tipo, carpeta=CARPETA) -> Path | None:
    """El snapshot más reciente (por hora de descarga) de un tipo: "pronostico" o "reanalisis"."""
    candidatos = []
    for ruta in Path(carpeta).glob("*.json"):
        m = json.loads(ruta.read_text(encoding="utf-8"))["metadata"]
        if m["tipo"] == tipo:
            candidatos.append((m["descargado_utc"], ruta))
    return max(candidatos)[1] if candidatos else None


# ---------------------------------------------------------------------------
# Series horarias: pronóstico + transición + climatología
# ---------------------------------------------------------------------------

def _serie(hourly, variables):
    t = np.array(hourly["time"], dtype="datetime64[m]")
    return t, {v: np.array([np.nan if x is None else x for x in hourly[v]], dtype=float) for v in variables}


def climatologia(era5: dict, tiempos, variables, ventana_dias=3):
    """Media de cada hora del calendario sobre todos los años y ± ventana_dias (misma hora del día).

    Conserva el ciclo diario típico de esos días; suaviza los eventos (frentes, días nublados).
    """
    tabla = {}
    for datos in era5["respuesta"]["por_anio"].values():
        t, vals = _serie(datos["hourly"], variables)
        for i, ti in enumerate(t):
            tabla[ti] = {v: vals[v][i] for v in variables}
    anios = sorted(int(a) for a in era5["respuesta"]["por_anio"])
    salida = {v: np.full(len(tiempos), np.nan) for v in variables}
    muestras = np.zeros(len(tiempos), dtype=int)
    for n, ti in enumerate(np.asarray(tiempos, dtype="datetime64[m]")):
        base = ti.astype(datetime)
        acumulado = {v: [] for v in variables}
        for anio in anios:
            for d in range(-ventana_dias, ventana_dias + 1):
                try:
                    clave = np.datetime64(base.replace(year=anio) + timedelta(days=d), "m")
                except ValueError:  # 29 de febrero
                    continue
                fila = tabla.get(clave)
                if fila is not None:
                    for v in variables:
                        acumulado[v].append(fila[v])
        for v in variables:
            x = np.array(acumulado[v], dtype=float)
            if np.isfinite(x).any():
                salida[v][n] = np.nanmean(x)
        muestras[n] = len(acumulado[variables[0]])
    return salida, muestras


def ensamblar(fecha_inicio, duracion_h, pronostico=None, era5=None, transicion_h=24.0, ventana_dias=3) -> Clima:
    """Serie horaria de T_aire, viento, G y nubosidad desde fecha_inicio (hora local), con su origen hora a hora."""
    if era5 is None and pronostico is None:
        raise ValueError("Se necesita al menos un snapshot (pronóstico o ERA5)")
    t0 = np.datetime64(fecha_inicio, "m")
    horas = np.arange(0, int(np.ceil(duracion_h)) + 1)
    tiempos = t0 + horas.astype("timedelta64[h]")
    nombres = {"T_aire": "temperature_2m", "viento": "wind_speed_10m", "G": "shortwave_radiation",
               "nubosidad": "cloud_cover"}
    if era5 is not None:
        clim, muestras = climatologia(era5, tiempos, tuple(nombres.values()), ventana_dias)
    else:
        clim, muestras = {v: np.full(len(tiempos), np.nan) for v in nombres.values()}, np.zeros(len(tiempos), int)
    valores = {k: clim[v].copy() for k, v in nombres.items()}
    origen = np.array(["climatología ERA5"] * len(tiempos), dtype=object)
    ultima = None
    if pronostico is not None:
        t_p, vals_p = _serie(pronostico["respuesta"]["hourly"], tuple(nombres.values()))
        pos = {ti: i for i, ti in enumerate(t_p)}
        cubiertas = np.array([ti in pos and np.isfinite(vals_p["temperature_2m"][pos[ti]]) for ti in tiempos])
        if cubiertas.any():
            ultima = int(np.flatnonzero(cubiertas)[-1])
            for k, v in nombres.items():
                p = np.array([vals_p[v][pos[ti]] if c else np.nan for ti, c in zip(tiempos, cubiertas)])
                ambos = np.flatnonzero(np.isfinite(p) & np.isfinite(clim[v]))
                anomalia = float(p[ambos[-1]] - clim[v][ambos[-1]]) if len(ambos) else 0.0
                valores[k] = np.where(cubiertas & np.isfinite(p), p, valores[k])
                # Después del pronóstico, su anomalía respecto de la climatología decae linealmente a 0.
                despues = horas > horas[ultima]
                peso = np.clip(1 - (horas - horas[ultima]) / transicion_h, 0, 1) if transicion_h > 0 else 0 * horas
                valores[k] = np.where(despues, valores[k] + peso * anomalia, valores[k])
            origen[cubiertas] = "pronóstico"
            origen[(horas > horas[ultima]) & (horas <= horas[ultima] + transicion_h)] = "transición"
    for k in valores:
        if not np.isfinite(valores[k]).all():
            raise ValueError(f"La serie de {k} tiene horas sin dato: falta un snapshot que cubra {fecha_inicio}")
    metadata = {
        "fecha_inicio": str(t0), "zona_horaria": ZONA, "horas": len(tiempos),
        "pronostico": None if pronostico is None else {
            "archivo": Path(pronostico["archivo"]).name, "modelo": pronostico["metadata"]["modelo"],
            "descargado_utc": pronostico["metadata"]["descargado_utc"],
            "horizonte": [pronostico["metadata"]["horizonte"]["desde"], pronostico["metadata"]["horizonte"]["hasta"]],
            "sha256": pronostico["metadata"]["sha256_respuesta"]},
        "climatologia": None if era5 is None else {
            "archivo": Path(era5["archivo"]).name, "anios": era5["metadata"]["anios"], "ventana_dias": ventana_dias,
            "muestras_por_hora": int(np.median(muestras)), "descargado_utc": era5["metadata"]["descargado_utc"],
            "sha256": era5["metadata"]["sha256_respuesta"]},
        "horas_pronostico": int(np.sum(origen == "pronóstico")), "horas_transicion": int(np.sum(origen == "transición")),
        "horas_climatologia": int(np.sum(origen == "climatología ERA5")), "transicion_h": transicion_h,
    }
    return Clima(t_h=horas.astype(float), T_aire=valores["T_aire"], viento=valores["viento"], G=valores["G"],
                 nubosidad=valores["nubosidad"] / 100.0, origen=list(origen), metadata=metadata)


def perfil_suelo(era5: dict, fecha_inicio, ventana_dias=3):
    """Perfil inicial del suelo T(z) y temperatura profunda, desde las capas de suelo de ERA5.

    Cada capa se ubica en su profundidad media (0,035; 0,175; 0,64; 1,775 m). El perfil
    inicial es la climatología de esas capas a la hora de inicio; bajo 1,775 m se
    mantiene el valor de la capa más profunda, que es también el Dirichlet del fondo
    (promedio de toda la ventana: a 2-3 m el ciclo diario ya no se siente).
    """
    t0 = np.datetime64(fecha_inicio, "m")
    capas = tuple(CAPAS_SUELO_ERA5)
    al_inicio, _ = climatologia(era5, np.array([t0]), capas, ventana_dias)
    horas = t0 + np.arange(0, 14 * 24 + 1, 6).astype("timedelta64[h]")
    en_ventana, _ = climatologia(era5, horas, capas[-1:], ventana_dias)
    prof = np.array(list(CAPAS_SUELO_ERA5.values()))
    T_capas = np.array([al_inicio[c][0] for c in capas])
    T_profunda = float(np.nanmean(en_ventana[capas[-1]])) if np.isfinite(en_ventana[capas[-1]]).any() else np.nan
    nota = "climatología a la hora de inicio"
    if not np.isfinite(T_capas).all() or not np.isfinite(T_profunda):
        # La fecha de inicio cae fuera de la ventana del snapshot (p. ej. una demo con fechas del pronóstico):
        # se usa el promedio de cada capa en todo el snapshot y se deja constancia.
        todas = {c: [] for c in capas}
        for datos in era5["respuesta"]["por_anio"].values():
            for c in capas:
                todas[c] += [x for x in datos["hourly"][c] if x is not None]
        T_capas = np.array([np.mean(todas[c]) for c in capas])
        T_profunda = float(T_capas[-1])
        nota = "fecha fuera de la ventana ERA5: promedio de cada capa en todo el snapshot"

    def T_inicial(z):
        return np.interp(-np.asarray(z, dtype=float), prof, T_capas)  # z < 0 bajo la base

    info = {"profundidades_m": prof.tolist(), "T_capas_C": [round(float(x), 2) for x in T_capas],
            "T_fondo_C": round(T_profunda, 2), "nota": nota}
    return T_inicial, T_profunda, info
