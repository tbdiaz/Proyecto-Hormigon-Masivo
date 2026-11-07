"""Pipeline meteorológico (E10/E11): snapshots trazables e inmutables, y serie pronóstico + climatología.

Funciona sin internet: usa snapshots sintéticos y los snapshots reales ya guardados en data/clima.
"""

import json
from pathlib import Path

import numpy as np
import pytest

from mcplan.clima import _guardar, climatologia, ensamblar, leer_snapshot, perfil_suelo, ultimo_snapshot

SNAPSHOTS = Path("data/clima/snapshots")
CAMPOS_METADATA = ("tipo", "fuente", "modelo", "descargado_utc", "sitio", "punto_modelo", "zona_horaria",
                   "unidades", "variables", "horizonte", "sha256_respuesta")


def _horas(desde, n):
    return [str(np.datetime64(desde, "m") + np.timedelta64(h, "h")) for h in range(n)]


def _hourly(desde, n, T, extra=None):
    h = {"time": _horas(desde, n), "temperature_2m": T, "wind_speed_10m": [3.0] * n,
         "shortwave_radiation": [0.0] * n, "cloud_cover": [50.0] * n}
    h.update(extra or {})
    return h


def _era5(tmp_path):
    capas = {c: [12.0 + i] * 24 * 20 for i, c in enumerate(("soil_temperature_0_to_7cm", "soil_temperature_7_to_28cm",
                                                            "soil_temperature_28_to_100cm",
                                                            "soil_temperature_100_to_255cm"))}
    por_anio = {str(a): {"hourly": _hourly(f"{a}-11-20T00:00", 24 * 20, [10.0 + (a - 2020)] * 24 * 20, capas)}
                for a in (2020, 2021)}
    meta = {"tipo": "reanalisis", "fuente": "prueba", "modelo": "era5", "descargado_utc": "2026-01-01T00:00:00Z",
            "sitio": {"latitud": 0, "longitud": 0}, "horizonte": {"desde": "", "hasta": ""}, "anios": [2020, 2021]}
    ruta = _guardar(tmp_path, "era5.json", meta, {"por_anio": por_anio}, tmp_path / "indice.csv")
    return leer_snapshot(ruta)


def _pronostico(tmp_path, desde, n, T):
    meta = {"tipo": "pronostico", "fuente": "prueba", "modelo": "m", "descargado_utc": "2026-11-07T00:00:00Z",
            "sitio": {"latitud": 0, "longitud": 0}, "horizonte": {"desde": desde, "hasta": _horas(desde, n)[-1]}}
    ruta = _guardar(tmp_path, "pron.json", meta, {"hourly": _hourly(desde, n, T)}, tmp_path / "indice.csv")
    return leer_snapshot(ruta)


def test_snapshot_inmutable_y_con_huella(tmp_path):
    snap = _pronostico(tmp_path, "2026-11-27T00:00", 10, [15.0] * 10)
    with pytest.raises(FileExistsError):  # nunca se reemplaza un snapshot
        _guardar(tmp_path, "pron.json", snap["metadata"], snap["respuesta"], tmp_path / "indice.csv")
    doc = json.loads((tmp_path / "pron.json").read_text(encoding="utf-8"))
    doc["respuesta"]["hourly"]["temperature_2m"][0] = 99.0  # alguien "corrige" el pronóstico a mano
    (tmp_path / "pron.json").write_text(json.dumps(doc), encoding="utf-8")
    with pytest.raises(ValueError, match="SHA-256"):
        leer_snapshot(tmp_path / "pron.json")
    assert len((tmp_path / "indice.csv").read_text(encoding="utf-8").splitlines()) == 2  # encabezado + 1 descarga


def test_climatologia_promedia_anios_y_dias():
    era5 = {"respuesta": {"por_anio": {str(a): {"hourly": _hourly(f"{a}-11-20T00:00", 24 * 20, [10.0 + (a - 2020)] * 480)}
                                       for a in (2020, 2021)}}}
    valores, muestras = climatologia(era5, np.array(["2026-11-27T08:00"], dtype="datetime64[m]"),
                                     ("temperature_2m",), ventana_dias=3)
    assert valores["temperature_2m"][0] == pytest.approx(10.5) and muestras[0] == 14  # 2 años x 7 días


def test_pronostico_transicion_y_climatologia(tmp_path):
    era5 = _era5(tmp_path)  # climatología: 10,5 °C
    pron = _pronostico(tmp_path, "2026-11-27T08:00", 30, [20.5] * 30)  # cubre las primeras 30 h
    c = ensamblar("2026-11-27T08:00", 72, pron, era5, transicion_h=24)
    assert c.metadata["horas_pronostico"] == 30 and c.metadata["horas_transicion"] == 24
    assert c.metadata["horas_climatologia"] == 73 - 54
    assert np.allclose(c.T_aire[:30], 20.5)
    assert c.T_aire[29 + 12] == pytest.approx(10.5 + 10.0 * 0.5)  # a mitad de la transición, media anomalía
    assert np.allclose(c.T_aire[29 + 24:], 10.5)
    assert c.origen[0] == "pronóstico" and c.origen[35] == "transición" and c.origen[-1] == "climatología ERA5"
    assert np.allclose(c.nubosidad, 0.5)  # % -> fracción
    sin_pron = ensamblar("2026-11-27T08:00", 24, None, era5)
    assert sin_pron.metadata["horas_climatologia"] == 25
    with pytest.raises(ValueError, match="sin dato"):
        ensamblar("2026-06-01T00:00", 24, None, era5)  # nadie cubre junio


def test_perfil_de_suelo_desde_era5(tmp_path):
    T_z, T_fondo, info = perfil_suelo(_era5(tmp_path), "2026-11-27T08:00")
    assert T_z(-0.035) == pytest.approx(12.0) and T_z(-3.0) == pytest.approx(15.0) and T_fondo == pytest.approx(15.0)
    assert info["profundidades_m"] == [0.035, 0.175, 0.64, 1.775]


@pytest.mark.skipif(not SNAPSHOTS.exists(), reason="no hay snapshots reales")
def test_snapshots_reales_tienen_la_metadata_del_curso():
    """Fuente, modelo, fecha de descarga, coordenadas, unidades, horizonte y variables (enunciado, sección 8)."""
    for tipo in ("pronostico", "reanalisis"):
        ruta = ultimo_snapshot(tipo, SNAPSHOTS)
        assert ruta is not None
        m = leer_snapshot(ruta)["metadata"]  # además verifica el SHA-256
        assert all(campo in m for campo in CAMPOS_METADATA)
        assert m["sitio"]["latitud"] == pytest.approx(-34.105457) and m["sitio"]["longitud"] == pytest.approx(-71.649539)
        assert m["unidades"]["temperature_2m"] == "°C" and m["unidades"]["wind_speed_10m"] == "m/s"
