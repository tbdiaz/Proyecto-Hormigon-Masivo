"""Lectura de un caso de simulación desde YAML.

La idea (punto 16 del enunciado) es que una decisión constructiva se cambie
editando el YAML, sin tocar el código. Las rutas dentro del YAML son relativas
a la carpeta donde está el propio YAML.
"""

import os
from dataclasses import dataclass
from math import cos, pi
from pathlib import Path
from typing import Optional

import yaml

from .hidratacion import HidratacionJuLei2019
from .solver1d import Borde, Capa, Material, Valor


@dataclass
class Caso:
    nombre: str
    capas: list[Capa]
    borde_inf: Borde
    borde_sup: Borde
    dz: float
    duracion_h: float
    dt_s: Optional[float]  # None = automático
    dt_max_s: Optional[float]  # tope del paso automático (precisión)
    guardar_cada_h: float
    factor_seguridad: float
    limite_Tmax: float
    limite_dT: float
    profundidad_superficie: float
    salida: Path


def _material(nombre: str, d: dict) -> Material:
    hid = None
    if "hidratacion" in d:
        h = dict(d["hidratacion"])
        modelo = h.pop("modelo")
        if modelo != "ju_lei_2019":
            raise ValueError(f"Modelo de hidratación no implementado: {modelo}")
        hid = HidratacionJuLei2019(**h)
    return Material(nombre, k=d["k"], rho=d["rho"], cp=d["cp"],
                    es_hormigon=d.get("es_hormigon", False), hidratacion=hid)


def temperatura_ambiente(d, hora_inicio: float) -> Valor:
    """Constante (número) o sinusoidal diaria {media, amplitud, hora_max}.

    T(t) = media + amplitud * cos(2*pi*(hora_del_dia - hora_max)/24)
    con hora_del_dia = hora_inicio + t/3600. Así la hora de inicio de la
    vertida (una decisión del plan) cambia el clima que "ve" el hormigón.
    """
    if isinstance(d, (int, float)):
        return float(d)
    if d.get("tipo") != "sinusoidal":
        raise ValueError(f"Temperatura ambiente no reconocida: {d}")
    media, amp, hmax = d["media"], d["amplitud"], d["hora_max"]
    return lambda t: media + amp * cos(2 * pi * (hora_inicio + t / 3600.0 - hmax) / 24.0)


def _borde(d: dict, hora_inicio: float) -> Borde:
    tipo = d["tipo"]
    if tipo == "dirichlet":
        return Borde.dirichlet(temperatura_ambiente(d["T"], hora_inicio))
    if tipo == "neumann":
        return Borde.neumann(d.get("q", 0.0))
    if tipo == "robin":
        h = d["h"]
        if "aislacion" in d:  # manta sobre la superficie: resistencia en serie con la convección
            a = d["aislacion"]
            h = 1.0 / (1.0 / h + a["espesor"] / a["k"])
        return Borde.robin(h, temperatura_ambiente(d["T_amb"], hora_inicio))
    raise ValueError(f"Tipo de borde desconocido: {tipo}")


def cargar_caso(ruta) -> Caso:
    ruta = Path(ruta)
    cfg = yaml.safe_load(ruta.read_text(encoding="utf-8"))
    base = ruta.parent

    materiales = {n: _material(n, d) for n, d in cfg["materiales"].items()}
    capas = [Capa(materiales[c["material"]], c["espesor"], c["T_inicial"], c.get("nombre", ""))
             for c in cfg["capas"]]

    # Opcional: el hormigón se arma desde el IFC (una capa por vertida en el eje).
    if "hormigon_desde_ifc" in cfg:
        from .ifc_io import columna_desde_ifc

        d = cfg["hormigon_desde_ifc"]
        mat = materiales[d["material"]]
        for v in columna_desde_ifc(base / d["archivo"], d.get("x", 0.0), d.get("y", 0.0)):
            capas.append(Capa(mat, v["z1"] - v["z0"], v["T_colocacion"], v["VertidaID"] or v["nombre"]))

    tiempo = cfg["tiempo"]
    hora_inicio = tiempo.get("hora_inicio", 0.0)
    dt = tiempo.get("dt_s", "auto")
    crit = cfg.get("criterios", {})
    return Caso(
        nombre=cfg.get("nombre", ruta.stem),
        capas=capas,
        borde_inf=_borde(cfg["bordes"]["inferior"], hora_inicio),
        borde_sup=_borde(cfg["bordes"]["superior"], hora_inicio),
        dz=cfg["malla"]["dz"],
        duracion_h=tiempo["duracion_h"],
        dt_s=None if dt == "auto" else float(dt),
        dt_max_s=tiempo.get("dt_max_s"),
        guardar_cada_h=tiempo.get("guardar_cada_h", 1.0),
        factor_seguridad=tiempo.get("factor_seguridad", 0.9),
        limite_Tmax=crit.get("Tmax_C", 70.0),
        limite_dT=crit.get("dT_C", 20.0),
        profundidad_superficie=crit.get("profundidad_superficie_m", 0.025),
        salida=Path(os.path.normpath(base / cfg.get("salida", f"../resultados/{ruta.stem}"))),
    )
