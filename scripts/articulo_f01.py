"""Artículo académico F01 (Word y PDF), armado desde los resultados de los estudios.

Uso, desde la raíz del repositorio y después de correr los estudios (ver docs/articulo_f01/LEEME.md):

    python scripts/articulo_f01.py          # docs/articulo_f01/Articulo_F01_Grupo4.docx + figuras + IFC del plan
    python scripts/articulo_f01.py --pdf    # además lo exporta a PDF con Word (Windows) y cuenta las páginas

Todo número del artículo sale de los resultados (resultados/estudio_era5, confiabilidad_r3, confiabilidad_era5,
estres_era5, retiro_tubos_era5, retiro_tubos2_era5, retiro_tubos3_era5, sensibilidad_era5, verificacion_ensemble): si
se vuelven a correr los estudios, basta con volver a correr este script. La propuesta es el primer plan de
config/confiabilidad_r3.yaml; la anterior, el de config/confiabilidad_era5.yaml.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import re
import subprocess
import sys
from collections import Counter
from datetime import datetime, timedelta
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))

from mcplan.estudio import cargar_estudio, leer_filas  # noqa: E402

RES = RAIZ / "resultados"
SALIDA = RAIZ / "docs" / "articulo_f01"
IMG = SALIDA / "img"
DOCX = SALIDA / "Articulo_F01_Grupo4.docx"
# La propuesta es el primer plan de CONF_PROP; CONF_ANT guarda la propuesta anterior (descartada con malla fina).
CONF_PROP, CONF_ANT = "confiabilidad_r3", "confiabilidad_era5"
PROPUESTA = cargar_estudio(RAIZ / "config" / f"{CONF_PROP}.yaml")["planes"][0]["nombre"]
LIM_T, LIM_D = 70.0, 20.0
RHO_CP = 2400 * 970  # J/m3K, el del modelo (Ju y Lei 2019)

# Paleta (skill de visualización, validada sobre blanco): categórica en orden fijo + tintas neutras.
AZUL, NARANJA, AQUA, AMARILLO, MAGENTA = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"
TINTA, TINTA2, GRIS, GRILLA, EJE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
GRIS_CLARO = "#cfcdc6"
COLOR_V = {1: AZUL, 2: NARANJA, 3: AQUA, 4: AMARILLO, 5: MAGENTA}
TEXTO_EN_V = {1: "white", 2: "white", 3: "white", 4: TINTA, 5: TINTA}

plt.rcParams.update({
    "font.family": ["Calibri", "Segoe UI", "DejaVu Sans"], "font.size": 8.5, "axes.titlesize": 9,
    "axes.titleweight": "bold", "axes.titlelocation": "left", "axes.labelsize": 8.5, "axes.edgecolor": EJE,
    "axes.labelcolor": TINTA2, "xtick.color": TINTA2, "ytick.color": TINTA2, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": GRILLA, "grid.linewidth": 0.6,
    "legend.frameon": False, "legend.fontsize": 8, "savefig.dpi": 220, "mathtext.fontset": "cm",
    "figure.facecolor": "white", "axes.facecolor": "white",
})
CM = 1 / 2.54  # pulgadas por cm


# ---------------------------------------------------------------------------
# Datos
# ---------------------------------------------------------------------------

def filas_de(carpeta: Path) -> list[dict]:
    """Filas de escenarios.csv sin repetidos (si una corrida se repitió, vale la última)."""
    unicas = {}
    for f in leer_filas(carpeta / "escenarios.csv"):
        unicas[(f["plan"], f["curva"], f["clima"])] = f
    return list(unicas.values())


def resumir(filas: list[dict]) -> dict:
    planes = {}
    for f in filas:
        planes.setdefault(f["plan"], []).append(f)
    out = {}
    for p, fs in planes.items():
        w = np.array([f["w"] for f in fs])
        out[p] = {"n": len(fs), "D": fs[0]["D_h"], "P_T": float(w @ [f["I_T"] for f in fs]),
                  "P_D": float(w @ [f["I_D"] for f in fs]), "P_F": float(w @ [f["I_F"] for f in fs]),
                  "MD_min": min(f["M_D"] for f in fs), "MD_max": max(f["M_D"] for f in fs),
                  "MT_max": max(f["M_T"] for f in fs), "peor": max(fs, key=lambda f: f["M_D"]),
                  "peor_T": max(fs, key=lambda f: f["M_T"]),
                  "MD14_max": max(f["M_D_14d"] for f in fs)}
    return out


def sensibilidades(carpetas: list[Path], filas_por_carpeta: dict) -> list[dict]:
    """Corridas con malla fina (mcplan estudio --sensibilidad): M_Δ con Δx fino frente al de la malla del estudio."""
    out = []
    for carpeta in carpetas:
        for ruta in sorted(carpeta.glob("sensibilidad_*.csv")):
            partes = ruta.stem[len("sensibilidad_"):].split("_")
            plan, curva, clima = partes[0], "_".join(partes[1:3]), "_".join(partes[3:])
            filas = list(csv.DictReader(open(ruta, encoding="utf-8")))
            base = [r for r in filas if float(r["dt_s"]) == 600 and abs(float(r["profundidad_m"]) - 0.025) < 1e-6]
            gruesa = next((r for r in base if abs(float(r["dx"]) - 0.275) < 1e-6), None)
            finas_r = sorted((r for r in base if float(r["dx"]) < 0.275 - 1e-6), key=lambda r: float(r["dx"]))
            fina = finas_r[0] if finas_r else None  # la malla más fina corrida en ese escenario
            intermedias = finas_r[1:]
            esc = next((f for f in filas_por_carpeta.get(carpeta, []) if (f["plan"], f["curva"], f["clima"]) == (plan, curva, clima)), None)
            MD_g = float(gruesa["M_D"]) if gruesa else (esc["M_D"] if esc else None)
            MT_g = float(gruesa["M_T"]) if gruesa else (esc["M_T"] if esc else None)
            if fina is None or MD_g is None:
                continue
            out.append({"plan": plan, "curva": curva, "clima": clima, "dx_fina": float(fina["dx"]),
                        "MD_gruesa": MD_g, "MD_fina": float(fina["M_D"]), "MT_gruesa": MT_g, "MT_fina": float(fina["M_T"]),
                        "delta": float(fina["M_D"]) - MD_g,
                        "otras": [r for r in filas if r is not fina and r is not gruesa and r not in intermedias],
                        "intermedias": [{"dx": float(r["dx"]), "M_D": float(r["M_D"]), "M_T": float(r["M_T"]),
                                         "delta": float(r["M_D"]) - MD_g} for r in intermedias],
                        "segundos_fina": float(fina["segundos"]), "t_MD_fina": float(fina["t_MD_h"])})
    return out


def _brecha_clima(X: dict) -> float:
    finas = [x["MD_fina"] for x in X.get("sens", []) if x["plan"] == PROPUESTA]
    return (max(finas) - min(finas)) if len(finas) > 1 else float("nan")


def corregir(filas: list[dict], delta: float, finos: dict | None = None) -> list[dict]:
    """M_Δ con la corrección de malla y los indicadores de falla recalculados con ella.

    Los escenarios corridos con malla fina (`finos`: (plan, curva, clima) -> M_Δ fino) usan ese valor directamente;
    el resto suma `delta`, la mayor diferencia fina - gruesa medida.
    """
    out = []
    finos = finos or {}
    for f in filas:
        g = dict(f)
        g["M_D_gruesa"] = f["M_D"]
        clave = (f["plan"], f["curva"], f["clima"])
        g["malla_fina"] = clave in finos
        g["M_D"] = finos[clave] if clave in finos else f["M_D"] + delta
        g["M_D_14d"] = f["M_D_14d"] + delta
        g["I_D"] = float(g["M_D"] > LIM_D)
        g["I_F"] = float(g["I_D"] or f["I_T"])
        out.append(g)
    return out


def tabla_delta(filas: list[dict], plan: str, delta: float, finos: dict) -> dict:
    """¿Puede el error de malla cambiar la clasificación? P_f de la propuesta en función de la corrección δ.

    Dos variantes: δ uniforme sumado a los escenarios de la malla del estudio, y la que usa el artículo (malla fina
    directa donde se corrió y δ en el resto). δ* es la menor corrección (paso 0,01 °C) con la que P_f supera 5 %.
    """
    fs = [f for f in filas if f["plan"] == plan]
    if not fs:
        return {}

    def pf(dl, usar_finos):
        k, p = 0, 0.0
        for f in fs:
            clave = (f["plan"], f["curva"], f["clima"])
            md = finos[clave] if usar_finos and clave in finos else f["M_D"] + dl
            if md > LIM_D or f["M_T"] > LIM_T:
                k, p = k + 1, p + f["w"]
        return p, k

    rejilla = [round(0.01 * i, 2) for i in range(0, 801)]
    crit = {m: next((x for x in rejilla if pf(x, m)[0] > 0.05), None) for m in (False, True)}
    valores = sorted({round(delta, 2), round(delta + 0.25, 2), 2.0, *(v for v in crit.values() if v is not None)})
    out = []
    for v in valores:
        (pu, ku), (pm, km) = pf(v, False), pf(v, True)
        out.append({"delta": v, "pf_unif": pu, "k_unif": ku, "pf_mixto": pm, "k_mixto": km})
    n_finos = sum(1 for f in fs if (f["plan"], f["curva"], f["clima"]) in finos)
    return {"filas": out, "crit_unif": crit[False], "crit_mixto": crit[True], "n": len(fs), "n_finos": n_finos}


def control_curvas(ruta: Path, desfase_h: float) -> dict:
    """Integral de cada curva frente a la tabla del encabezado (β del archivo) y energía de la rampa inicial."""
    texto = ruta.read_text(encoding="utf-8").splitlines()
    datos = np.genfromtxt(ruta, comments="#")
    cab = next((ln for ln in texto if ln.startswith("# time_h")), None)
    if cab is None:
        return {}
    nombres = cab[1:].split()[1:]
    tabla = {}
    for ln in texto:
        m = re.match(r"#\s+(\S+)\s+(.+)$", ln)
        if m and m.group(1) in nombres:
            grupos = [g.strip() for g in re.split(r"\s{2,}", m.group(2).strip()) if "/" in g]
            if len(grupos) >= 2:  # columnas β bajo / β del archivo / β alto: la del medio es la de este archivo
                tabla[m.group(1)] = float(grupos[1].split("/")[0])
    t = datos[:, 0]
    errores, rampa, total = [], [], []
    for i, n in enumerate(nombres):
        q = datos[:, i + 1]
        ok = np.isfinite(q)
        E = float(np.trapezoid(q[ok], t[ok])) * 3600 / 1e6
        if n in tabla:
            errores.append(abs(E / tabla[n] - 1))
        rampa.append(0.5 * float(q[ok][0]) * desfase_h * 3600 / 1e6)  # 0 -> Q(desfase) en las primeras horas
        total.append(E)
    return {"n": len(nombres), "n_tabla": len(tabla), "err_max": max(errores, default=float("nan")),
            "rampa_max": max(rampa), "rampa_frac": max(r / e for r, e in zip(rampa, total)), "desfase_h": desfase_h}


def hardware() -> str:
    """Procesador, núcleos y memoria del computador que corrió los estudios (Windows; si no, lo que se pueda)."""
    import os
    import platform

    nombre, nucleos, memoria = platform.processor() or platform.machine(), None, None
    if sys.platform == "win32":
        try:
            import ctypes
            import winreg

            with winreg.OpenKey(winreg.HKEY_LOCAL_MACHINE, r"HARDWARE\DESCRIPTION\System\CentralProcessor\0") as k:
                nombre = winreg.QueryValueEx(k, "ProcessorNameString")[0].strip()

            class _Mem(ctypes.Structure):
                _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                            ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                            ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                            ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                            ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]
            m = _Mem()
            m.dwLength = ctypes.sizeof(_Mem)
            ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))
            memoria = round(m.ullTotalPhys / 2 ** 30)
            r = subprocess.run(["powershell", "-NoProfile", "-Command",
                                "(Get-CimInstance Win32_Processor | Measure-Object NumberOfCores -Sum).Sum"],
                               capture_output=True, text=True, timeout=60)
            nucleos = int(r.stdout.strip()) if r.stdout.strip().isdigit() else None
        except (OSError, ValueError, subprocess.SubprocessError):
            pass
    hilos = os.cpu_count()
    partes = [re.sub(r"\(R\)|\(TM\)", "", nombre).replace("  ", " ")]
    if nucleos:
        partes.append(f"{nucleos} núcleos y {hilos} hilos")
    elif hilos:
        partes.append(f"{hilos} hilos")
    if memoria:
        partes.append(f"{memoria} GB de RAM")
    return ", ".join(partes)


def serie(carpeta: Path, plan, curva, clima):
    archivos = sorted(glob.glob(str(carpeta / "series" / f"{plan}__{curva}__{clima}__*.npz")))
    return dict(np.load(archivos[-1])) if archivos else None


def confiabilidad(filas: list[dict], plan: str = PROPUESTA) -> dict:
    """Enumeración ponderada de los escenarios de un plan: probabilidades, percentiles ponderados y extremos."""
    fs = [f for f in filas if f["plan"] == plan]
    if not fs:
        return {}
    curvas = sorted({f["curva"] for f in fs}, key=lambda c: (c.split("_")[1] != "SS", c))
    anios = sorted({f["clima"] for f in fs})
    MD = np.full((len(curvas), len(anios)), np.nan)
    MT = np.full_like(MD, np.nan)
    for f in fs:
        i, j = curvas.index(f["curva"]), anios.index(f["clima"])
        MD[i, j], MT[i, j] = f["M_D"], f["M_T"]
    w = np.array([f["w"] for f in fs])
    w = w / w.sum()

    def cuantil(valores, q):  # cuantil de la distribución ponderada de escenarios
        orden = np.argsort(valores)
        acum = np.cumsum(w[orden])
        return float(np.asarray(valores)[orden][np.searchsorted(acum, q)])

    md, mt = [f["M_D"] for f in fs], [f["M_T"] for f in fs]
    return {
        "n": len(fs), "curvas": curvas, "anios": anios, "MD": MD, "MT": MT,
        "k_T": sum(f["I_T"] for f in fs), "k_D": sum(f["I_D"] for f in fs), "k_F": sum(f["I_F"] for f in fs),
        "P_T": float(w @ [f["I_T"] for f in fs]), "P_D": float(w @ [f["I_D"] for f in fs]),
        "P_F": float(w @ [f["I_F"] for f in fs]),
        "MD_p50": cuantil(md, 0.5), "MD_p95": cuantil(md, 0.95), "MD_max": max(md),
        "MT_p50": cuantil(mt, 0.5), "MT_p95": cuantil(mt, 0.95), "MT_max": max(mt),
        "peor": max(fs, key=lambda f: f["M_D"]), "peor_T": max(fs, key=lambda f: f["M_T"]), "filas": fs,
    }





def representantes(carpeta: Path) -> dict:
    """Tercil climático -> temporada que lo representa ("fría" -> "ERA5_2013"), desde definicion.json del estudio."""
    ruta = carpeta / "definicion.json"
    if not ruta.exists():
        return {}
    clima = json.loads(ruta.read_text(encoding="utf-8")).get("clima", {})
    return {h["clave"]: h["representante"].replace(" ", "_") for h in clima.get("historias", [])
            if isinstance(h, dict) and h.get("representante")}


def margen_material(filas_base: list[dict], filas_estres: list[dict], fuentes: dict, plan=PROPUESTA,
                    escala=1.2, mapa: dict | None = None) -> dict:
    """Sensibilidad determinista: ¿con cuánta más energía llegaría M_Δ a 20 °C?

    Se ajusta M_Δ (peor clima por curva) contra la energía a 14 días de cada curva SS, con y sin el +20 % de la prueba
    de estrés, y se despeja la energía crítica E_crít. No es una probabilidad.
    """
    E = {n: fuentes[n].energia(336.0) / 1e6 for n in fuentes}  # MJ/m3
    base, estres = {}, {}
    for f in filas_estres:
        if f["plan"] == plan:
            estres[f["curva"]] = max(estres.get(f["curva"], -np.inf), f["M_D"])
    # Mismos climas en los dos lados: la base puede tener 30 temporadas y el estrés sólo los 3 terciles, cada uno
    # representado por una temporada (mapa: tercil -> nombre de esa temporada en la base).
    mapa = mapa or {}
    climas_estres = {mapa.get(f["clima"], f["clima"]) for f in filas_estres if f["plan"] == plan}
    for f in filas_base:
        if f["plan"] == plan and f["clima"] in climas_estres:
            base[f["curva"]] = max(base.get(f["curva"], -np.inf), f["M_D"])
    ss = [c for c in base if c.endswith("_SS") and c in estres]
    if len(ss) < 2:
        return {}
    x = np.array([E[c] for c in ss] + [escala * E[c] for c in ss])
    y = np.array([base[c] for c in ss] + [estres[c] for c in ss])
    b, a = np.polyfit(x, y, 1)
    r2 = 1 - np.sum((y - (a + b * x)) ** 2) / np.sum((y - y.mean()) ** 2)
    E_crit = (LIM_D - a) / b
    e_max = max(E.values())
    return {"a": float(a), "b": float(b), "r2": float(r2), "E_crit": float(E_crit), "E_max": float(e_max),
            "curva_max": max(E, key=E.get), "exceso": float(E_crit / e_max - 1)}


def geometria(dx=0.11):
    from mcplan.ifc_io import leer_vertidas
    from mcplan.voxelizacion import voxelizar

    vertidas = leer_vertidas(RAIZ / "geometry" / "foundation.ifc")
    vox = voxelizar(vertidas, dx)
    info = []
    for i, v in enumerate(vertidas):
        V = v.vertices
        info.append({"id": i + 1, "nombre": v.nombre, "z0": float(V[:, 2].min()), "z1": float(V[:, 2].max()),
                     "ancho": float(V[:, 0].max() - V[:, 0].min()), "V": float(vox.volumen_ifc[i])})
    return vox, info


def cronograma_vertidas(cfg_ruta: Path, plan=PROPUESTA):
    """Hora de inicio de cada vertida en la colocación continua por capas (la misma lógica del estudio)."""
    from mcplan.estudio import _base
    from mcplan.planes import activacion_por_capas, expandir_familias

    cfg = cargar_estudio(cfg_ruta)
    planes = {p.nombre: p for p in expandir_familias(cfg["planes"], float(cfg.get("T_colocacion", 18)))}
    _, vox, dom = _base(cfg["_caso"], float(cfg["dx"]))
    act, crono = activacion_por_capas(planes[plan], dom, vox, float(cfg["recursos"]["tasa_colocacion_m3_h"]))
    vid = dom.vertida_id
    inicios = [float(np.min(act.t_act[vid == v])) for v in range(1, vox.n_vertidas + 1)]
    fines = inicios[1:] + [float(crono["D_h"])]
    return inicios, fines, float(crono["D_h"]), cfg


def historias_reales(cfg) -> tuple[list, dict]:
    """Las temporadas ERA5 alineadas a t0 (una historia por año)."""
    from mcplan.clima import historias_era5, leer_snapshot

    era5 = leer_snapshot(Path(cfg["_ruta"]).parent / cfg["clima"]["climatologia"])
    return historias_era5(era5, cfg["t0"], 24 * 30), era5["metadata"]


def historias_ensemble_de(cfg) -> tuple[list, dict]:
    """Los 51 miembros del ensemble más reciente para t0, completados con la climatología de las temporadas ERA5."""
    from mcplan.clima import historias_ensemble, leer_snapshot, ultimo_snapshot

    carpeta = RAIZ / "data" / "clima" / "snapshots"
    ruta = ultimo_snapshot("ensemble", carpeta)
    if ruta is None:
        return [], {}
    ens = leer_snapshot(ruta)
    era5 = leer_snapshot(Path(cfg["_ruta"]).parent / cfg["clima"]["climatologia"])
    return historias_ensemble(ens, era5, cfg["t0"], 24 * 30, 24.0, 3), ens["metadata"]


def _guardar(fig, ruta):
    fig.savefig(ruta, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    return ruta


def _seccion_y_planta(ax_sec, ax_pla, vox, info, alfa=1.0, rotular=True):
    """Sección en y = 0 y planta (vertida más alta en cada x, y), con el color de cada vertida."""
    from matplotlib.colors import ListedColormap

    vid = vox.vertida_id
    x, y, z = vox.centros(0), vox.centros(1), vox.centros(2)
    j0 = int(np.argmin(np.abs(y)))
    cmap = ListedColormap(["white"] + [COLOR_V[v] for v in range(1, 6)])
    sec = np.ma.masked_equal(vid[:, j0, :].T, 0)
    ext = [x[0] - vox.dx / 2, x[-1] + vox.dx / 2, z[0] - vox.dx / 2, z[-1] + vox.dx / 2]
    ax_sec.imshow(sec, origin="lower", extent=ext, cmap=cmap, vmin=0, vmax=5, interpolation="nearest", alpha=alfa)
    ax_sec.fill_between([-11, 11], -1.0, 0, color="#ece6dc", zorder=0, lw=0)  # suelo
    ax_sec.set_xlim(-11, 11)
    ax_sec.set_ylim(-1.0, 4.8)
    ax_sec.set_aspect("equal")
    ax_sec.set_xlabel("x [m]")
    ax_sec.set_ylabel("z [m]")
    ax_sec.grid(False)
    alto = np.argmax(np.where(vid > 0, np.arange(vid.shape[2])[None, None, :], -1), axis=2)
    top = np.take_along_axis(vid, alto[:, :, None], axis=2)[:, :, 0]
    top = np.ma.masked_equal(np.where(vid.max(axis=2) > 0, top, 0).T, 0)
    ext_p = [x[0] - vox.dx / 2, x[-1] + vox.dx / 2, y[0] - vox.dx / 2, y[-1] + vox.dx / 2]
    ax_pla.imshow(top, origin="lower", extent=ext_p, cmap=cmap, vmin=0, vmax=5, interpolation="nearest", alpha=alfa)
    ax_pla.set_xlim(-10, 10)
    ax_pla.set_ylim(-10, 10)
    ax_pla.set_aspect("equal")
    ax_pla.set_xlabel("x [m]")
    ax_pla.set_ylabel("y [m]  (norte ↑)")
    ax_pla.grid(False)
    ax_sec.set_xticks(range(-10, 11, 5))
    ax_sec.set_yticks(range(0, 5))
    ax_pla.set_xticks(range(-10, 11, 5))
    ax_pla.set_yticks(range(-10, 11, 5))
    if rotular:
        for g in info:
            zc = (g["z0"] + g["z1"]) / 2
            ax_sec.text(0, zc, f"V{g['id']} · {g['V']:.0f} m³".replace(".", ","), ha="center", va="center",
                        fontsize=7.5, color=TEXTO_EN_V[g["id"]], fontweight="bold")
        ax_sec.text(-10.8, -0.6, "suelo (3 m modelados)", fontsize=7, color=TINTA2)


def fig_geometria(vox, info, ruta):
    fig = plt.figure(figsize=(16.5 * CM, 6.2 * CM))
    a1 = fig.add_axes([0.0, 0.1, 0.60, 0.85])
    a2 = fig.add_axes([0.71, 0.1, 0.29, 0.85])
    _seccion_y_planta(a1, a2, vox, info)
    a1.set_title("Sección por el eje (y = 0)")
    a2.set_title("Planta (vertida superior)")
    return _guardar(fig, ruta)


def fig_curvas(fuentes, ruta):
    edades = np.linspace(0, 336, 1345)
    fig, ejes = plt.subplots(2, 2, figsize=(16.5 * CM, 9.2 * CM), sharex="col", sharey="row")
    colores = {"M1": AZUL, "M2": NARANJA, "M3": AQUA, "M4": AMARILLO}
    finales = {}
    for c, campania in enumerate(("SS", "22A")):
        for nombre, f in fuentes.items():
            m, camp = nombre.split("_")
            if camp != campania:
                continue
            Q = np.array([f.Q(a) for a in edades])
            E = np.array([f.energia(a) for a in edades]) / RHO_CP
            ejes[0, c].plot(edades, Q, color=colores[m], lw=1.6, label=m)
            ejes[1, c].plot(edades, E, color=colores[m], lw=1.6)
            finales.setdefault(campania, []).append(f"{m} {E[-1]:.0f}")
        ejes[1, c].text(330, 3, "a 14 días: " + " · ".join(finales[campania]) + " °C", ha="right", fontsize=7.5,
                        color=TINTA)
        ejes[0, c].set_title(f"Campaña {campania}" + ("  (menor confianza)" if campania == "22A" else ""))
        ejes[1, c].set_xlabel("edad del hormigón [h]")
        ejes[1, c].set_xlim(0, 336)
        ejes[1, c].set_xticks([0, 48, 96, 168, 240, 336])
    ejes[0, 0].set_ylabel("fuente Q [W/m³]")
    ejes[1, 0].set_ylabel("aumento adiabático [°C]")
    ejes[0, 1].legend(ncol=2, loc="upper right")
    fig.subplots_adjust(wspace=0.08, hspace=0.12)
    return _guardar(fig, ruta)


def fig_clima(historias, ensemble, ruta, dias=14):
    T = np.array([c.T_aire for _, c in historias])[:, : 24 * dias + 1]
    v = np.array([c.viento for _, c in historias])[:, : 24 * dias + 1]
    t = np.arange(T.shape[1]) / 24
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(16.5 * CM, 8.4 * CM), sharex=True,
                                 gridspec_kw={"height_ratios": [1.6, 1]})
    for ax, X, clave, unidad in ((a1, T, "T_aire", "T del aire [°C]"), (a2, v, "viento", "viento a 10 m [m/s]")):
        lo, mid, hi = np.percentile(X, [10, 50, 90], axis=0)
        ax.fill_between(t, lo, hi, color=GRIS_CLARO, alpha=0.6, lw=0, label="ERA5 1996–2025: 10–90 %")
        ax.plot(t, mid, color=GRIS, lw=1.0, label="ERA5 1996–2025: mediana")
        if ensemble:
            E = np.array([getattr(c, clave) for _, c in ensemble])[:, : 24 * dias + 1]
            pron = np.array([o == "pronóstico" for o in ensemble[0][1].origen[: 24 * dias + 1]])
            em = np.where(pron, np.median(E, axis=0), np.nan)
            ax.plot(t, em, color=AZUL, lw=1.4, label="ensemble ECMWF 2026: mediana (pronóstico)")
        ax.set_ylabel(unidad)
    a1.legend(ncol=3, loc="upper left", bbox_to_anchor=(0, 1.22), fontsize=7.5)
    a2.set_xlabel("días desde t₀ (jue 19 nov 2026, 08:00)")
    a2.set_xlim(0, dias)
    a2.set_xticks(range(0, dias + 1, 2))
    fig.subplots_adjust(hspace=0.12)
    return _guardar(fig, ruta)


def fig_cronograma(inicios, fines, D, ruta, retiro_d=26.0, tubos=(), monitoreo_d=29.0):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(16.5 * CM, 5.4 * CM), gridspec_kw={"width_ratios": [1, 1.6]})
    for v in range(1, 6):
        a1.barh(6 - v, fines[v - 1] - inicios[v - 1], left=inicios[v - 1], height=0.62, color=COLOR_V[v], lw=0)
        a1.text(fines[v - 1] + 0.4, 6 - v, f"{inicios[v - 1]:.1f}–{fines[v - 1]:.1f} h".replace(".", ","),
                va="center", fontsize=7.5, color=TINTA)
    a1.set_yticks(range(1, 6), [f"V{v}" for v in range(5, 0, -1)])
    a1.set_xlim(0, 30)
    a1.set_xlabel("horas desde t₀")
    a1.set_title("Colocación (27 m³/h)")
    a1.grid(axis="y", visible=False)
    filas = [("Monitoreo", 0, monitoreo_d, GRIS, ""), ("Moldaje + manta 10 cm", 0, retiro_d, TINTA2, "")]
    for t in sorted(tubos, key=lambda t: -t["vertida"]):
        v = t["vertida"]
        t_on = inicios[v - 1] / 24 + t.get("desde_h", 0) / 24
        tramos = [(0.0, t["durante_h"], t.get("T_agua_C", 25))] if not t.get("programa") else \
            [(a["hasta_h"], b["hasta_h"], b["T_agua_C"]) for a, b in zip([{"hasta_h": 0.0}] + t["programa"], t["programa"])]
        if not t.get("programa"):
            tramos = [(0.0, t["durante_h"], t.get("T_agua_C", 25))]
        filas.append((f"Tubos V{v}", t_on, t_on + tramos[-1][1] / 24, COLOR_V[v],
                      " / ".join(f"{T:g}" for _, _, T in tramos) + " °C"))
    filas.append(("Horizonte de evaluación", 0, retiro_d + 7, EJE, ""))
    for k, (nombre, a, b, c, nota) in enumerate(filas):
        a2.barh(k, b - a, left=a, height=0.6, color=c, lw=0)
        a2.text(b + 0.3, k, nombre + (f" ({nota})" if nota else ""), va="center", fontsize=7.2, color=TINTA)
    a2.axvline(retiro_d, color=TINTA, lw=0.8)
    a2.text(retiro_d, len(filas) - 0.35, f" retiro día {retiro_d:.0f}", fontsize=7.5, color=TINTA, va="bottom")
    a2.set_yticks([])
    a2.set_xlim(0, retiro_d + 22)
    a2.set_xticks(range(0, int(retiro_d + 8), 5))
    a2.set_xlabel("días desde t₀")
    a2.set_title("Protecciones, agua de los tubos y monitoreo")
    a2.grid(axis="y", visible=False)
    a2.spines["left"].set_visible(False)
    fig.subplots_adjust(wspace=0.22)
    return _guardar(fig, ruta)


def fig_temperaturas(s_peor, rotulo_peor, series_anios, anio_peor, ruta, t_ret=None, delta=0.0):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(16.5 * CM, 6.4 * CM))
    t = s_peor["t_h"] / 24
    fin = float(t[-1])
    a1.plot(t, s_peor["Tmax"], color=NARANJA, lw=1.6, label="T máx. (núcleo)")
    a1.plot(t, s_peor["Tmin"], color=AZUL, lw=1.6, label="T mín. a 25 mm")
    a1.plot(t, s_peor["T_amb"], color=GRIS, lw=0.6, label="aire")
    a1.set_xlim(0, fin)
    a1.set_xlabel("días desde t₀")
    a1.set_ylabel("°C")
    a1.set_title(f"Peor escenario: {rotulo_peor}")
    a1.legend(loc="upper right")
    for nombre, s in series_anios.items():
        if nombre != anio_peor:
            a2.plot(s["t_h"] / 24, s["dT"], color=GRIS_CLARO, lw=0.6)
    if anio_peor in series_anios:
        s = series_anios[anio_peor]
        a2.plot(s["t_h"] / 24, s["dT"], color=NARANJA, lw=1.5, label=f"peor temporada ({anio_peor.split('_')[-1]})")
    a2.axhline(LIM_D - delta, color=TINTA, lw=1.0)
    a2.text(0.3, LIM_D - delta + 0.4, f"20 °C − corrección de malla ({delta:.1f})".replace(".", ","), fontsize=7.2,
            color=TINTA)
    for ax in (a1, a2):
        if t_ret is not None:
            ax.axvline(t_ret / 24, color=TINTA2, lw=0.8, ls=":")
    a2.set_xlim(0, fin)
    a2.set_ylim(0, 22)
    a2.set_xlabel("días desde t₀")
    a2.set_ylabel("ΔT [°C] (malla del estudio)")
    a2.set_title("ΔT(t) con la misma curva: 30 temporadas")
    a2.legend(loc="lower right")
    fig.subplots_adjust(wspace=0.25)
    return _guardar(fig, ruta)


def fig_confiabilidad(C, ruta):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(16.5 * CM, 6.4 * CM), sharey=True)
    curvas, anios = C["curvas"], C["anios"]
    rng = np.random.default_rng(3)
    for i, c in enumerate(curvas):
        y = len(curvas) - 1 - i
        jit = rng.uniform(-0.22, 0.22, len(anios))
        for ax, M in ((a1, C["MD"]), (a2, C["MT"])):
            ax.scatter(M[i], y + jit, s=9, color=AZUL, alpha=0.6, lw=0)
    a1.axvline(LIM_D, color=TINTA, lw=1.0)
    a1.text(LIM_D - 0.15, -0.45, "límite 20 °C", ha="right", fontsize=7.5, color=TINTA)
    a2.axvline(LIM_T, color=TINTA, lw=1.0)
    a2.text(LIM_T - 0.5, -0.45, "límite 70 °C", ha="right", fontsize=7.5, color=TINTA)
    a1.set_yticks(range(len(curvas)), list(reversed(curvas)))
    a1.set_ylim(-0.7, len(curvas) - 0.4)
    a1.set_xlabel(r"M$_\Delta$ en todo el horizonte (con corrección de malla) [°C]")
    a2.set_xlabel(r"M$_T$: temperatura máxima [°C]")
    a1.set_xlim(min(10, np.nanmin(C["MD"]) - 1), 21)
    a2.set_xlim(min(30, np.nanmin(C["MT"]) - 2), 72)
    a1.set_title(f"Diferencial ({len(curvas)} curvas × {len(anios)} temporadas)")
    a2.set_title("Temperatura máxima")
    a1.grid(axis="y", visible=False)
    a2.grid(axis="y", visible=False)
    fig.subplots_adjust(wspace=0.08)
    return _guardar(fig, ruta)


def fig_D_pf(fam: dict, prop: str, ruta, etiquetas: dict):
    """Duración de la colocación frente a P_f de cada plan de la familia, con la línea del 5 %."""
    fig, ax = plt.subplots(figsize=(16.5 * CM, 6.4 * CM))
    grupos = {}  # sólo los planes con la misma D y el mismo P_f se separan (horizontalmente, para que se vean)
    for p, r in fam.items():
        grupos.setdefault((round(r["D"], 1), round(100 * r["P_F"], 1)), []).append(p)
    Ds = sorted({k[0] for k in grupos})
    paso = 0.012 * (max(Ds) - min(Ds) + 10)
    for (D_, pf_), ps in grupos.items():
        for i, p in enumerate(sorted(ps, key=lambda p: p != prop)):
            x = D_ + (i - (len(ps) - 1) / 2) * paso
            es = p == prop
            ax.scatter(x, pf_, s=42 if es else 20, color=NARANJA if es else AZUL, edgecolor="white", lw=0.8,
                       zorder=3 if es else 2)
            if es:
                ax.annotate("propuesta final", (x, pf_), xytext=(14, 16), textcoords="offset points", fontsize=7.5,
                            color=TINTA, arrowprops={"arrowstyle": "-", "color": TINTA2, "lw": 0.6})
    for D_ in Ds:
        ps = [p for (dd, _), v in grupos.items() if dd == D_ for p in v]
        if len(ps) == 1 and ps[0] != prop:
            ax.annotate(etiquetas.get(ps[0], ps[0]), (D_, 100 * fam[ps[0]]["P_F"]), xytext=(0, -13),
                        textcoords="offset points", fontsize=7, color=TINTA2, ha="center")
        elif len(ps) > 1:
            ax.annotate(f"vertida única: {len(ps)} planes", (D_, 100 * max(fam[p]["P_F"] for p in ps)),
                        xytext=(8, -4), textcoords="offset points", fontsize=7, color=TINTA2, ha="left", va="top")
    ax.set_xlim(min(Ds) - 6, max(Ds) + 12)
    ax.axhline(5, color=TINTA, lw=1.0)
    ax.text(max(Ds) + 11, 7, "P$_f$ = 5 %", fontsize=7.5, color=TINTA, ha="right")
    ax.set_ylim(-5, 106)
    ax.set_xlabel("D: duración de la colocación [h]")
    ax.set_ylabel("P$_f$ [%]")
    ax.set_title(f"Familia de planes ({len(fam)}): duración frente a probabilidad de falla")
    return _guardar(fig, ruta)


def fig_sensores(vox, info, calientes, frios, sensores, ruta):
    fig = plt.figure(figsize=(16.5 * CM, 6.6 * CM))
    a1 = fig.add_axes([0.0, 0.1, 0.60, 0.8])
    a2 = fig.add_axes([0.71, 0.1, 0.29, 0.8])
    _seccion_y_planta(a1, a2, vox, info, alfa=0.28, rotular=False)
    a1.scatter(calientes[:, 0], calientes[:, 2], s=7, color=NARANJA, alpha=0.5, lw=0, label="punto más caliente (240 esc.)")
    a1.scatter(frios[:, 0], frios[:, 2], s=7, color=AZUL, alpha=0.5, lw=0, label="punto más frío a 25 mm (240 esc.)")
    a2.scatter(calientes[:, 0], calientes[:, 1], s=7, color=NARANJA, alpha=0.5, lw=0)
    a2.scatter(frios[:, 0], frios[:, 1], s=7, color=AZUL, alpha=0.5, lw=0)
    visibles = [s for s in sensores if s["tipo"] != "aire"]
    for s in visibles:
        mk = {"núcleo": "o", "25 mm": "D", "suelo": "s"}[s["tipo"]]
        a1.scatter(s["x"], s["z"], marker=mk, s=34, facecolor="white", edgecolor=TINTA, lw=1.1, zorder=4)
        a2.scatter(s["x"], s["y"], marker=mk, s=30, facecolor="white", edgecolor=TINTA, lw=1.1, zorder=4)
    for ax, (i, j), tam in ((a1, ("x", "z"), 7.5), (a2, ("x", "y"), 7)):
        grupos = []
        for s in visibles:  # un rótulo por grupo de sensores a menos de 0,6 m entre sí
            for g in grupos:
                if math.hypot(s[i] - g[0][i], s[j] - g[0][j]) < 0.6:
                    g.append(s)
                    break
            else:
                grupos.append([s])
        for g in grupos:
            ax.annotate(", ".join(s["id"] for s in g), (g[0][i], g[0][j]), xytext=(5, 4), textcoords="offset points",
                        fontsize=tam, color=TINTA)
    a1.scatter([], [], marker="o", s=30, facecolor="white", edgecolor=TINTA, label="sensor de núcleo")
    a1.scatter([], [], marker="D", s=30, facecolor="white", edgecolor=TINTA, label="sensor a 25 mm")
    a1.scatter([], [], marker="s", s=30, facecolor="white", edgecolor=TINTA, label="sensor en el suelo")
    a1.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, fontsize=7.5)
    a1.set_title("")
    a2.set_title("Planta", pad=4)
    return _guardar(fig, ruta)


def fig_retiro(rondas: list, ruta, delta=0.0):
    """Un panel por ronda: rondas = [(título, {nombre: (serie, día del retiro)}), ...]. A lo más 7 variantes por
    panel, para que cada una conserve un color distinguible (los colores no se repiten dentro de un panel)."""
    colores = [AZUL, NARANJA, AQUA, AMARILLO, MAGENTA, "#008300", "#6250d6"]
    rondas = [r for r in rondas if r[1]]
    fig, ejes = plt.subplots(len(rondas), 1, figsize=(16.5 * CM, 5.4 * CM * len(rondas)), sharex=True, squeeze=False)
    t_fin = max(float(s["t_h"][-1]) for _, series in rondas for s, _ in series.values()) / 24
    for ax, (titulo, series) in zip(ejes[:, 0], rondas):
        if len(series) > len(colores):
            raise ValueError(f"{titulo}: {len(series)} variantes en un panel; el máximo es {len(colores)}")
        for (nombre, (s, d)), c in zip(series.items(), colores):
            ax.plot(s["t_h"] / 24, s["dT"], color=c, lw=1.1, label=nombre)
            ax.axvline(d, color=c, lw=0.7, alpha=0.7)
        ax.axhline(LIM_D - delta, color=TINTA, lw=1.0)
        ax.text(0.3, LIM_D - delta + 0.4, f"20 °C − corrección de malla ({delta:.1f} °C)".replace(".", ","),
                fontsize=7, color=TINTA)
        ax.set_xlim(0, t_fin)
        ax.set_ylim(0, 26)
        ax.set_title(titulo, fontsize=8.5, loc="left")
        ax.set_ylabel("ΔT [°C]")
        ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=6.5, frameon=False)
    ejes[-1, 0].set_xlabel("días desde t₀ (cada línea vertical marca un retiro; malla del estudio)")
    fig.tight_layout()
    return _guardar(fig, ruta)


def ecuacion(tex: str, ruta: Path, tam=11.5):
    fig = plt.figure(figsize=(0.1, 0.1))
    fig.text(0, 0, tex, fontsize=tam, color=TINTA)
    fig.savefig(ruta, dpi=300, bbox_inches="tight", pad_inches=0.02, transparent=False)
    plt.close(fig)
    return ruta


# ---------------------------------------------------------------------------
# Sensores: ubicados donde el modelo dice que están el punto caliente y los puntos fríos
# ---------------------------------------------------------------------------

def _zona(z):
    if z < 1.10:
        return "losa (V1–V2)"
    if z < 2.20:
        return "tronco V3"
    if z < 3.30:
        return "tronco V4"
    return "pedestal V5"


def plan_sensores(C) -> tuple[list[dict], np.ndarray, np.ndarray, Counter]:
    fs = C["filas"]
    cal = np.array([[f["x_MT"], f["y_MT"], f["z_MT"]] for f in fs])
    fri = np.array([[f["x_Tmin"], f["y_Tmin"], f["z_Tmin"]] for f in fs])
    zonas = Counter(_zona(p[2]) for p in fri)
    sensores = []
    c = np.median(cal, axis=0)
    sensores.append({"id": "N1", "tipo": "núcleo", "x": float(c[0]), "y": float(c[1]), "z": float(c[2]),
                     "donde": "punto más caliente previsto (mediana de los 240 escenarios)", "n": 2})
    sensores.append({"id": "N2", "tipo": "núcleo", "x": 0.0, "y": 0.0, "z": 2.15,
                     "donde": "eje, mitad del espesor (control y comparación con Ju y Lei)", "n": 2})
    k = 1
    for zona, cuenta in zonas.most_common():
        pts = fri[[ _zona(p[2]) == zona for p in fri]]
        m = np.median(pts, axis=0)
        sensores.append({"id": f"S{k}", "tipo": "25 mm", "x": float(m[0]), "y": float(m[1]), "z": float(m[2]),
                         "donde": f"punto frío en {zona} ({cuenta} de {len(fs)} escenarios)", "n": 2})
        k += 1
        # La cara opuesta (sombra o sol distinto): misma altura, posición reflejada.
        sensores.append({"id": f"S{k}", "tipo": "25 mm", "x": float(-m[0]), "y": float(-m[1]), "z": float(m[2]),
                         "donde": f"cara opuesta de {zona} (asimetría de sol y viento)", "n": 1})
        k += 1
    for nombre, (x, y, z) in (("superficie superior de la losa", (0.0, -7.5, 1.075)),
                              ("superficie superior del pedestal", (0.0, 0.0, 4.275))):
        if not any(abs(s["z"] - z) < 0.2 and abs(s["y"] - y) < 1.5 for s in sensores):
            sensores.append({"id": f"S{k}", "tipo": "25 mm", "x": x, "y": y, "z": z, "donde": nombre, "n": 1})
            k += 1
    sensores.append({"id": "G1", "tipo": "suelo", "x": 0.0, "y": 0.0, "z": -0.5,
                     "donde": "suelo bajo el eje (interacción con el terreno)", "n": 1})
    sensores.append({"id": "A1", "tipo": "aire", "x": 0.0, "y": -12.0, "z": 1.5,
                     "donde": "aire a la sombra, 1,5 m sobre el terreno, con anemómetro", "n": 2})
    return sensores, cal, fri, zonas


# ---------------------------------------------------------------------------
# Documento Word
# ---------------------------------------------------------------------------

class Informe:
    def __init__(self):
        from docx import Document
        from docx.enum.text import WD_ALIGN_PARAGRAPH
        from docx.shared import Cm, Pt, RGBColor

        self.Cm, self.Pt, self.RGB, self.ALIGN = Cm, Pt, RGBColor, WD_ALIGN_PARAGRAPH
        self.doc = Document()
        self.nfig = 0
        self.ntab = 0
        self.refs = {}
        s = self.doc.sections[0]
        s.page_width, s.page_height = Cm(21.59), Cm(27.94)  # carta
        s.top_margin, s.bottom_margin, s.left_margin, s.right_margin = Cm(2.0), Cm(1.8), Cm(2.2), Cm(2.0)
        st = self.doc.styles
        st["Normal"].font.name = "Calibri"
        st["Normal"].font.size = Pt(10.5)
        st["Normal"]._element.rPr.rFonts.set(self._qn("w:eastAsia"), "Calibri")
        pf = st["Normal"].paragraph_format
        pf.space_after, pf.space_before, pf.line_spacing = Pt(4), Pt(0), 1.07
        for nombre, tam, antes, despues in (("Heading 1", 13.5, 12, 4), ("Heading 2", 11, 8, 3)):
            h = st[nombre]
            h.font.name, h.font.size, h.font.bold = "Calibri", Pt(tam), True
            h.font.color.rgb = RGBColor(0x18, 0x4F, 0x95)
            rf = h.element.rPr.rFonts
            for att in ("w:asciiTheme", "w:hAnsiTheme", "w:eastAsiaTheme", "w:cstheme"):
                rf.attrib.pop(self._qn(att), None)
            h.paragraph_format.space_before, h.paragraph_format.space_after = Pt(antes), Pt(despues)
            h.paragraph_format.keep_with_next = True
        for nombre in ("List Bullet", "List Number"):
            st[nombre].font.name, st[nombre].font.size = "Calibri", Pt(10.5)
            st[nombre].paragraph_format.space_after = Pt(2)
        cap = st["Caption"]
        cap.font.name, cap.font.size, cap.font.italic, cap.font.bold = "Calibri", Pt(9), False, False
        cap.font.color.rgb = RGBColor(0x52, 0x51, 0x4E)
        cap.paragraph_format.space_after = Pt(8)
        self._pie()

    @staticmethod
    def _qn(t):
        from docx.oxml.ns import qn

        return qn(t)

    def _campo(self, par, instr):
        from docx.oxml import OxmlElement

        fld = OxmlElement("w:fldSimple")
        fld.set(self._qn("w:instr"), instr)
        r = OxmlElement("w:r")
        t = OxmlElement("w:t")
        t.text = "1"
        r.append(t)
        fld.append(r)
        par._p.append(fld)

    def _pie(self):
        s = self.doc.sections[0]
        s.different_first_page_header_footer = True
        p = s.header.paragraphs[0]
        r = p.add_run("Grupo 4 · Proyecto 2 Hormigón Masivo · Artículo F01")
        r.font.size, r.font.color.rgb = self.Pt(8), self.RGB(0x89, 0x87, 0x81)
        p.alignment = self.ALIGN.RIGHT
        p = s.footer.paragraphs[0]
        p.alignment = self.ALIGN.CENTER
        r = p.add_run("Página ")
        r.font.size = self.Pt(8)
        self._campo(p, "PAGE")
        r = p.add_run(" de ")
        r.font.size = self.Pt(8)
        self._campo(p, "NUMPAGES")

    # Texto con marcas mínimas: **negrita**, *cursiva*, ~subíndice~, ^superíndice^ (anidables dentro de ** y *).
    def _runs(self, par, texto, tam=None, color=None, negrita=False, cursiva=False):
        for trozo in re.split(r"(\*\*.+?\*\*|\*[^*\s][^*]*?\*|~[^~]+?~|\^[^^]+?\^)", texto):
            if not trozo:
                continue
            if len(trozo) > 4 and trozo.startswith("**") and trozo.endswith("**"):
                self._runs(par, trozo[2:-2], tam, color, True, cursiva)
                continue
            if len(trozo) > 2 and trozo[0] == "*" and trozo[-1] == "*":
                self._runs(par, trozo[1:-1], tam, color, negrita, True)
                continue
            sub = len(trozo) > 2 and trozo[0] == "~" and trozo[-1] == "~"
            sup = len(trozo) > 2 and trozo[0] == "^" and trozo[-1] == "^"
            r = par.add_run(trozo[1:-1] if (sub or sup) else trozo)
            r.bold = True if negrita else None
            r.italic = True if cursiva else None
            if sub:
                r.font.subscript = True
            if sup:
                r.font.superscript = True
            if tam:
                r.font.size = self.Pt(tam)
            if color:
                r.font.color.rgb = self.RGB.from_string(color.lstrip("#"))
        return par

    def h1(self, texto):
        self.doc.add_heading(texto, level=1)

    def h2(self, texto):
        self.doc.add_heading(texto, level=2)

    def p(self, texto, tam=None, alinear=None, despues=None, color=None):
        par = self._runs(self.doc.add_paragraph(), texto, tam, color)
        if alinear == "centro":
            par.alignment = self.ALIGN.CENTER
        elif alinear == "justificado":
            par.alignment = self.ALIGN.JUSTIFY
        if despues is not None:
            par.paragraph_format.space_after = self.Pt(despues)
        return par

    def vinetas(self, items, estilo="List Bullet"):
        for it in items:
            self._runs(self.doc.add_paragraph(style=estilo), it)

    def ecuacion(self, ruta):
        par = self.doc.add_paragraph()
        par.alignment = self.ALIGN.CENTER
        par.add_run().add_picture(str(ruta))
        par.paragraph_format.space_after = self.Pt(4)

    def imagen(self, ruta, ancho_cm):
        par = self.doc.add_paragraph()
        par.alignment = self.ALIGN.CENTER
        par.add_run().add_picture(str(ruta), width=self.Cm(ancho_cm))

    def figura(self, ruta, leyenda, ancho_cm=16.5, ref=None):
        par = self.doc.add_paragraph()
        par.alignment = self.ALIGN.CENTER
        par.paragraph_format.keep_with_next = True
        par.paragraph_format.space_after = self.Pt(2)
        par.add_run().add_picture(str(ruta), width=self.Cm(ancho_cm))
        self.nfig += 1
        if ref:
            self.refs["f:" + ref] = self.nfig
        cap = self.doc.add_paragraph(style="Caption")
        self._runs(cap, f"**Figura {self.nfig}.** {leyenda}")
        return self.nfig

    def tabla(self, titulo, encabezado, filas, anchos, alinear=None, destacar=(), tam=8.5, ref=None):
        from docx.enum.table import WD_TABLE_ALIGNMENT
        from docx.oxml import OxmlElement

        self.ntab += 1
        if ref:
            self.refs["t:" + ref] = self.ntab
        cap = self.doc.add_paragraph(style="Caption")
        cap.paragraph_format.keep_with_next = True
        cap.paragraph_format.space_after = self.Pt(2)
        self._runs(cap, f"**Tabla {self.ntab}.** {titulo}")
        t = self.doc.add_table(rows=1 + len(filas), cols=len(encabezado))
        t.alignment = WD_TABLE_ALIGNMENT.CENTER
        t.autofit = False
        tblPr = t._tbl.tblPr
        bordes = OxmlElement("w:tblBorders")
        for b in ("top", "left", "bottom", "right", "insideH", "insideV"):
            e = OxmlElement(f"w:{b}")
            e.set(self._qn("w:val"), "single")
            e.set(self._qn("w:sz"), "4")
            e.set(self._qn("w:space"), "0")
            e.set(self._qn("w:color"), "C3C2B7")
            bordes.append(e)
        tblPr.append(bordes)
        mar = OxmlElement("w:tblCellMar")
        for lado, v in (("top", 20), ("bottom", 20), ("left", 70), ("right", 70)):
            e = OxmlElement(f"w:{lado}")
            e.set(self._qn("w:w"), str(v))
            e.set(self._qn("w:type"), "dxa")
            mar.append(e)
        tblPr.append(mar)
        alinear = alinear or ["l"] + ["c"] * (len(encabezado) - 1)
        for i, fila in enumerate([encabezado] + filas):
            tr = t.rows[i]._tr
            trPr = tr.get_or_add_trPr()
            cs = OxmlElement("w:cantSplit")
            trPr.append(cs)
            if i == 0:
                th = OxmlElement("w:tblHeader")
                trPr.append(th)
            for j, texto in enumerate(fila):
                celda = t.rows[i].cells[j]
                celda.width = self.Cm(anchos[j])
                par = celda.paragraphs[0]
                par.paragraph_format.space_after = self.Pt(0)
                par.paragraph_format.line_spacing = 1.0
                par.alignment = {"l": self.ALIGN.LEFT, "c": self.ALIGN.CENTER, "r": self.ALIGN.RIGHT}[alinear[j]]
                self._runs(par, str(texto), tam=tam)
                if i == 0 or (i - 1) in destacar:
                    for r in par.runs:
                        r.bold = True
                    fondo = "F0EFEC" if i == 0 else "E3EEFB"
                    tcPr = celda._tc.get_or_add_tcPr()
                    shd = OxmlElement("w:shd")
                    shd.set(self._qn("w:val"), "clear")
                    shd.set(self._qn("w:color"), "auto")
                    shd.set(self._qn("w:fill"), fondo)
                    tcPr.append(shd)
        if len(filas) <= 14:  # tablas cortas: todas las filas "con la siguiente", así Word no las parte
            for fila in t.rows[:-1]:
                for celda in fila.cells:
                    for par in celda.paragraphs:
                        par.paragraph_format.keep_with_next = True
        self.doc.add_paragraph().paragraph_format.space_after = self.Pt(2)
        return self.ntab

    def salto(self):
        from docx.enum.text import WD_BREAK

        self.doc.add_paragraph().add_run().add_break(WD_BREAK.PAGE)

    def guardar(self, ruta):
        self.doc.save(str(ruta))


# ---------------------------------------------------------------------------
# Contenido
# ---------------------------------------------------------------------------

def f1(x, nd=1):
    """Número con coma decimal (es-CL)."""
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "—"
    return f"{x:,.{nd}f}".replace(",", "X").replace(".", ",").replace("X", ".")


def pct(x, nd=1):
    if x is None or not math.isfinite(x):
        return "—"
    if x == 0:
        return "0 %"
    return f"{f1(100 * x, nd)} %"





def hora(t0: datetime, h: float) -> str:
    dias = ["lun", "mar", "mié", "jue", "vie", "sáb", "dom"]
    t = t0 + timedelta(hours=h)
    return f"{dias[t.weekday()]} {t.day} {t:%H:%M}"


def git_commit() -> str:
    try:
        sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=RAIZ, capture_output=True, text=True).stdout.strip()
        sucio = subprocess.run(["git", "status", "--porcelain", "--", "src", "config", "scripts"], cwd=RAIZ,
                               capture_output=True, text=True).stdout.strip()
        return sha + (" + cambios sin commit" if sucio else "")
    except OSError:
        return "—"


ETIQUETAS = {
    "unica": ("Vertida única sin protección (moldaje 72 h)", "18", "—", "—"),
    "unica-M10": ("Única + manta", "18", "10 cm", "—"),
    "unica-M10-T10": ("Única + manta + preenfriado", "10", "10 cm", "—"),
    "unica-M10-A25": ("Única + manta + tubos en todo", "18", "10 cm", "V1–V5, 25 °C, 14 d"),
    "unica-M5-T10-A25": ("Única + manta 5 cm + preenfriado + tubos en todo", "10", "5 cm", "V1–V5, 25 °C, 14 d"),
    "unica-M10-T10-A25": ("Única + manta + preenfriado + tubos en todo", "10", "10 cm", "V1–V5, 25 °C, 14 d"),
    "unica-M10-T8-A25": ("Ídem, colocación a 8 °C", "8", "10 cm", "V1–V5, 25 °C, 14 d"),
    "unica-M10-T10-A25-S75": ("Ídem, tubos cada 0,75 m", "10", "10 cm", "V1–V5, 25 °C, 14 d"),
    "unica-M10-T8-A25-S75": ("Ídem, 8 °C y tubos cada 0,75 m", "8", "10 cm", "V1–V5, 25 °C, 14 d"),
    "unica-M10-T10-N25": ("Única + manta + preenfriado + tubos en el núcleo (14 d)", "10", "10 cm", "V3–V5, 25 °C, 14 d"),
    "unica-M10-T10-N20": ("Ídem con agua a 20 °C", "10", "10 cm", "V3–V5, 20 °C, 14 d"),
    "unica-M10-T8-N20": ("Ídem con agua a 20 °C y 8 °C", "8", "10 cm", "V3–V5, 20 °C, 14 d"),
    "unica-M10-T8-N25": ("Tubos en el núcleo (14 d), colocación a 8 °C", "8", "10 cm", "V3–V5, 25 °C, 14 d"),
    "unica-M10-T12-N25": ("Tubos en el núcleo (14 d), colocación a 12 °C", "12", "10 cm", "V3–V5, 25 °C, 14 d"),
    "unica-M10-T14-N25": ("Tubos en el núcleo (14 d), colocación a 14 °C", "14", "10 cm", "V3–V5, 25 °C, 14 d"),
    "unica-M10-T16-N25": ("Tubos en el núcleo (14 d), colocación a 16 °C", "16", "10 cm", "V3–V5, 25 °C, 14 d"),
    "unica-M10-T10-N25-R28": ("Propuesta anterior (10 °C, retiro el día 28)", "10", "10 cm", "V3–V5, programa"),
    "5V-M10-T10-A25-e12": ("5 vertidas, esperas de 12 h (más lenta)", "10", "10 cm", "V1–V5, 25 °C, 14 d"),
    "5V-M10-T10-A25-e24": ("5 vertidas, esperas de 24 h (más lenta)", "10", "10 cm", "V1–V5, 25 °C, 14 d"),
}
RONDAS_RETIRO = ("retiro_tubos_era5", "retiro_tubos2_era5", "retiro_tubos3_era5")
DIAS = ["lunes", "martes", "miércoles", "jueves", "viernes", "sábado", "domingo"]
MESES = ["enero", "febrero", "marzo", "abril", "mayo", "junio", "julio", "agosto", "septiembre", "octubre",
         "noviembre", "diciembre"]
AUTORES = "Trinidad Díaz · Pablo Arancibia · Eduardo Vergara"
CURSO = "Métodos Computacionales en Obras Civiles (IOC4201) · Universidad de los Andes"
PROFESOR = "José Antonio Abell Mena"  # de Canvas (IOC4201-202620-834); se puede cambiar con --profesor
try:
    import ifcopenshell as _ifc

    IFC_VERSION = "v" + _ifc.version.split(".post")[0]
except ImportError:
    IFC_VERSION = ""


def fecha_larga(t: datetime) -> str:
    return f"{DIAS[t.weekday()]} {t.day} de {MESES[t.month - 1]} de {t.year}, {t:%H:%M}"


def contar_tests() -> int:
    r = subprocess.run([sys.executable, "-m", "pytest", "--collect-only", "-q"], cwd=RAIZ, capture_output=True, text=True)
    m = re.search(r"(\d+) tests? collected", r.stdout)
    return int(m.group(1)) if m else -1


def plan_de(cfg: dict, nombre: str) -> dict:
    return next(p for p in cfg["planes"] if p.get("nombre") == nombre)


def programa_texto(tubo: dict) -> str:
    """'25 °C hasta el día 7, 20 °C hasta el 14, 15 °C hasta el 20,5' (horas desde que parte el circuito)."""
    if tubo.get("programa"):
        partes = [f"{f1(t['T_agua_C'], 0)} °C hasta el día {f1(t['hasta_h'] / 24, 1).replace(',0', '')}"
                  for t in tubo["programa"]]
        return ", ".join(partes)
    return f"{f1(tubo.get('T_agua_C', 15), 0)} °C durante {f1(tubo.get('durante_h', 168) / 24, 0)} días"


def retiro_de(plan: dict) -> float:
    """Horas desde el fin de la colocación hasta el retiro de la aislación (o del moldaje si no hay aislación)."""
    cub = plan.get("cubiertas", {})
    for clave in ("aislacion", "moldaje"):
        r = cub.get(clave, {}).get("retiro_h")
        if isinstance(r, (int, float)):
            return float(r)
    return float("nan")


def antes_despues(carpeta: Path, filas: list[dict], t_ret: float, delta: float) -> dict:
    """M_Δ antes y después del retiro (con la corrección de malla), desde las series de cada escenario."""
    antes, despues, peor = [], [], None
    for f in filas:
        s = serie(carpeta, f["plan"], f["curva"], f["clima"])
        if s is None:
            continue
        t, dT = s["t_h"], s["dT"]
        a = float(np.nanmax(dT[t < t_ret])) + delta
        d = float(np.nanmax(dT[t >= t_ret])) + delta if np.any(t >= t_ret) else float("nan")
        antes.append(a)
        despues.append(d)
        if peor is None or d > peor[0]:
            peor = (d, f)
    return {"antes_max": max(antes, default=float("nan")), "despues_max": max(despues, default=float("nan")),
            "antes_min": min(antes, default=float("nan")), "despues_min": min(despues, default=float("nan")),
            "peor_despues": peor[1] if peor else None, "n": len(antes)}


def resumen_pesos(filas: list[dict], solo: str | None = None) -> dict:
    """P_T, P_Δ y P_f con los pesos de la enumeración; con `solo` = "SS", sólo esa campaña (pesos renormalizados)."""
    fs = [f for f in filas if solo is None or f["curva"].endswith("_" + solo)]
    w = np.array([f["w"] for f in fs])
    if not len(fs) or w.sum() <= 0:
        return {}
    w = w / w.sum()
    return {"n": len(fs), "P_T": float(w @ [f["I_T"] for f in fs]), "P_D": float(w @ [f["I_D"] for f in fs]),
            "P_F": float(w @ [f["I_F"] for f in fs]), "MD_max": max(f["M_D"] for f in fs),
            "MT_max": max(f["M_T"] for f in fs)}


def criterio_retiro(carpeta: Path, filas: list[dict], t_ret: float) -> dict:
    """Lo que medirían los sensores al momento del retiro en los escenarios que después cumplen."""
    dT, nucleo_aire, nucleo = [], [], []
    for f in filas:
        if f["M_D"] > LIM_D:
            continue
        s = serie(carpeta, f["plan"], f["curva"], f["clima"])
        if s is None:
            continue
        t = s["t_h"]
        i, j = int(np.searchsorted(t, t_ret)), int(np.searchsorted(t, t_ret + 72))
        if i >= len(t):
            continue
        dT.append(float(s["dT"][i]))
        nucleo.append(float(s["Tmax"][i]))
        nucleo_aire.append(float(s["Tmax"][i] - np.min(s["T_amb"][i:j])))
    return {"dT_max": max(dT, default=float("nan")), "nucleo_aire_max": max(nucleo_aire, default=float("nan")),
            "nucleo_max": max(nucleo, default=float("nan")), "n": len(dT)}


def preparar() -> dict:
    """Lee resultados, calcula las métricas y dibuja las figuras (una sola vez)."""
    IMG.mkdir(parents=True, exist_ok=True)
    X = {}
    cfg_fin = cargar_estudio(RAIZ / "config" / f"{CONF_PROP}.yaml")
    X["prop"] = cfg_fin["planes"][0]["nombre"]
    X["prop_cfg"] = cfg_fin["planes"][0]
    X["t0"] = datetime.fromisoformat(cfg_fin["t0"])
    dirs = {k: RES / v for k, v in (("fam", "estudio_era5"), ("conf", CONF_PROP), ("ant", CONF_ANT), ("est", "estres_era5"),
                                     ("var", "sensibilidad_era5"),
                                     ("ver", "verificacion_ensemble"))}
    X["dirs"] = dirs
    filas = {k: filas_de(d) for k, d in dirs.items()}
    # Exploración del retiro en rondas (cada una declarada antes de ver sus resultados): plan -> carpeta.
    X["rondas"] = [r for r in RONDAS_RETIRO if (RAIZ / "config" / f"{r}.yaml").exists()]
    X["carpeta_ret"], filas["ret"] = {}, []
    for r in X["rondas"]:
        for f in filas_de(RES / r):
            X["carpeta_ret"][f["plan"]] = RES / r
            filas["ret"].append(f)
    prop = X["prop"]

    # Corrección de malla: la mayor diferencia (0,11 m - 0,275 m) de M_Δ en los escenarios de la propuesta verificados
    # con malla fina, medida sobre todo el horizonte (incluye el retiro).
    X["sens"] = sensibilidades([dirs["conf"]], {dirs["conf"]: filas["conf"]})
    X["prop_ant"] = prop_ant = cargar_estudio(RAIZ / "config" / f"{CONF_ANT}.yaml")["planes"][0]["nombre"]
    X["sens_ant"] = [s_ for s_ in sensibilidades([dirs["ant"]], {dirs["ant"]: filas["ant"]}) if s_["plan"] == prop_ant]
    # δ es la mayor diferencia medida en TODAS las corridas con malla fina (propuesta actual y anterior): es un error
    # del método numérico, no del plan, y así la regla no depende de qué escenarios se refinaron.
    deltas = [s_["delta"] for s_ in X["sens"] if s_["plan"] == prop] + [s_["delta"] for s_ in X["sens_ant"]]
    X["delta"] = d = max(deltas) if deltas else 0.0
    X["delta_ant_max"] = max((s_["delta"] for s_ in X["sens_ant"]), default=float("nan"))
    if X["sens_ant"] and prop_ant != prop:
        pa = max(X["sens_ant"], key=lambda s_: s_["delta"])
        X["delta_ant_txt"] = (f"{pa['curva']} {pa['clima'].split('_')[-1]}, el mayor de sus "
                              f"{len(X['sens_ant'])} escenarios con malla fina")
    if prop_ant != prop and filas["ant"]:
        finos_ant = {(s_["plan"], s_["curva"], s_["clima"]): s_["MD_fina"] for s_ in X["sens_ant"]}
        X["C_ant"] = confiabilidad(corregir(filas["ant"], d, finos_ant), prop_ant)
    else:
        X["C_ant"] = {}
    finos = {(s_["plan"], s_["curva"], s_["clima"]): s_["MD_fina"] for s_ in X["sens"]}
    c = {k: corregir(v, d, finos if k == "conf" else None) for k, v in filas.items()}
    X["tabla_delta"] = tabla_delta(filas["conf"], prop, d, finos)
    X["filas_c"] = c
    X["fam"] = {k: v for k, v in resumir(c["fam"]).items()
                if v["n"] == max((u["n"] for u in resumir(c["fam"]).values()), default=0)}
    X["fam_gruesa"] = {k: v for k, v in resumir(filas["fam"]).items() if k in X["fam"]}
    X["MT_familia"] = max((f["M_T"] for f in filas["fam"]), default=float("nan"))
    X["C"] = confiabilidad(c["conf"], prop)
    _cg = confiabilidad(corregir(filas["conf"], 0.0), prop)
    X["pf_gruesa"] = _cg["P_F"] if _cg else None
    X["solo_ss"] = resumen_pesos([f for f in c["conf"] if f["plan"] == prop], "SS")
    # Un estudio a medio correr daría un P_f parcial: sólo se informan los planes con todos sus escenarios.
    def completos(res):
        n_max = max((v["n"] for v in res.values()), default=0)
        return {k: v for k, v in res.items() if v["n"] == n_max}

    X["est"] = completos(resumir(c["est"])).get(prop, {})
    X["var"] = completos(resumir(c["var"]))
    X["ver"] = completos(resumir(c["ver"])).get(prop, {})
    X["ret"] = resumir(c["ret"])

    print("Geometría, cronograma y clima...")
    X["vox"], X["geo"] = geometria()
    X["inicios"], X["fines"], X["D"], _ = cronograma_vertidas(RAIZ / "config" / f"{CONF_PROP}.yaml", prop)
    X["retiro_h"] = retiro_de(X["prop_cfg"])
    X["t_ret"] = X["D"] + X["retiro_h"]
    X["ad"] = antes_despues(dirs["conf"], [f for f in filas["conf"] if f["plan"] == prop], X["t_ret"], d)
    X["crit"] = criterio_retiro(dirs["conf"], [f for f in c["conf"] if f["plan"] == prop], X["t_ret"])
    X["ret_cfg"], X["ronda_de"] = {}, {}
    for n_r, r in enumerate(X["rondas"], start=1):
        for p in cargar_estudio(RAIZ / "config" / f"{r}.yaml")["planes"]:
            X["ret_cfg"][p["nombre"]] = p
            X["ronda_de"][p["nombre"]] = n_r
    X["ret_ad"] = {}
    for nombre, p in X["ret_cfg"].items():
        fs = [f for f in filas["ret"] if f["plan"] == nombre]
        if fs:
            X["ret_ad"][nombre] = antes_despues(X["carpeta_ret"][nombre], fs, X["D"] + retiro_de(p), d)
    # Familia: núcleo al retirar la manta y M_Δ antes/después del retiro, en los planes cuyos tubos (si los hay) se
    # apagan antes del retiro. Sirve para afirmar, con datos, qué pasa sin seguir enfriando hasta el retiro.
    X["fam_cfg"] = {p["nombre"]: p for p in cargar_estudio(RAIZ / "config" / "estudio_era5.yaml")["planes"]
                    if "nombre" in p}
    X["fam_retiro"] = {}
    for nombre, r in X["fam_gruesa"].items():
        p = X["fam_cfg"].get(nombre)
        if nombre == prop or p is None or not math.isfinite(retiro_de(p)):
            continue
        if not isinstance(p.get("cubiertas", {}).get("aislacion", {}).get("retiro_h"), (int, float)):
            continue  # sin manta no hay retiro de la aislación que evaluar (el desmolde temprano es otra cosa)
        if any(t.get("programa") for t in p.get("tubos", [])):
            continue
        t_r = r["D"] + retiro_de(p)
        fs = [f for f in filas["fam"] if f["plan"] == nombre]
        nuc = []
        for f in fs:
            s = serie(dirs["fam"], nombre, f["curva"], f["clima"])
            if s is not None and s["t_h"][-1] >= t_r:
                nuc.append(float(np.interp(t_r, s["t_h"], s["Tmax"])))
        X["fam_retiro"][nombre] = {"dia": t_r / 24, "nucleo_max": max(nuc, default=float("nan")),
                                   "nucleo_min": min(nuc, default=float("nan")),
                                   **antes_despues(dirs["fam"], fs, t_r, d), "P_F": X["fam"][nombre]["P_F"]}
    historias, X["meta_era5"] = historias_reales(cfg_fin)
    X["n_temporadas"] = len(historias)
    ens = historias_ensemble_de(cfg_fin)
    X["meta_ens"] = ens[1]

    from mcplan.curvas import cargar_fuentes

    fuentes = cargar_fuentes(RAIZ / "data" / "hidratacion" / "source_terms.txt", desfase_h=2.0, cola="exponencial")
    X["fuentes_E"] = {n: f.energia(336.0) / 1e6 for n, f in fuentes.items()}
    X["mat"] = margen_material(c["conf"], c["est"], fuentes, prop, mapa=representantes(dirs["est"]))
    X["terciles"] = representantes(dirs["fam"])

    print("Figuras...")
    figs = X["figs"] = {"geo": fig_geometria(X["vox"], X["geo"], IMG / "f_geometria.png"),
                        "curvas": fig_curvas(fuentes, IMG / "f_curvas.png"),
                        "clima": fig_clima(historias, ens[0], IMG / "f_clima.png")}
    figs["crono"] = fig_cronograma(X["inicios"], X["fines"], X["D"], IMG / "f_cronograma.png",
                                   retiro_d=X["t_ret"] / 24, tubos=X["prop_cfg"]["tubos"],
                                   monitoreo_d=X["t_ret"] / 24 + 3)
    C = X["C"]
    X["sensores"], X["zonas"] = [], Counter()
    if C:
        peor = C["peor"]
        s_peor = serie(dirs["conf"], prop, peor["curva"], peor["clima"])
        series_anios = {a: serie(dirs["conf"], prop, peor["curva"], a) for a in C["anios"]}
        series_anios = {a: s for a, s in series_anios.items() if s is not None}
        rot = f"{peor['curva']}, temporada {peor['clima'].split('_')[-1]}"
        if s_peor is not None:
            figs["temp"] = fig_temperaturas(s_peor, rot, series_anios, peor["clima"], IMG / "f_temperaturas.png",
                                            t_ret=X["t_ret"], delta=d)
        figs["conf"] = fig_confiabilidad(C, IMG / "f_confiabilidad.png")
    if X["fam"]:
        figs["dpf"] = fig_D_pf(X["fam"], prop, IMG / "f_D_pf.png", {p: v[0] for p, v in ETIQUETAS.items()})
        X["sensores"], cal, fri, X["zonas"] = plan_sensores(C)
        figs["sens"] = fig_sensores(X["vox"], X["geo"], cal, fri, X["sensores"], IMG / "f_sensores.png")
    paneles = [(f"Ronda {n_r}", {}) for n_r in range(1, len(X["rondas"]) + 1)]
    for nombre in X["ret_cfg"]:
        if nombre in X["ret"]:
            pr = X["ret"][nombre]["peor"]
            s = serie(X["carpeta_ret"][nombre], nombre, pr["curva"], pr["clima"])
            if s is not None:
                paneles[X["ronda_de"][nombre] - 1][1][nombre] = (s, (X["D"] + retiro_de(X["ret_cfg"][nombre])) / 24)
    if any(series for _, series in paneles):
        figs["retiro"] = fig_retiro(paneles, IMG / "f_retiro.png", d)
    detalle = sorted(dirs["conf"].glob(f"detalle_{prop}_*"))
    d0 = detalle[0] if detalle else None
    X["corte"] = d0 / "cortes_y0.png" if d0 and (d0 / "cortes_y0.png").exists() else None
    X["vista"] = d0 / "vista_3d.png" if d0 and (d0 / "vista_3d.png").exists() else None
    X["tubos"] = json.loads((d0 / "tubos.json").read_text(encoding="utf-8")) if d0 and (d0 / "tubos.json").exists() else {}
    res = [abs(f["residuo_energia"]) for f in filas["conf"] if f["plan"] == prop
           and math.isfinite(f.get("residuo_energia", float("nan")))]
    if res and max(res) > 0:
        e = math.floor(math.log10(max(res)))
        X["residuo_txt"] = f"{f1(max(res) / 10 ** e, 1)}·10^{e}^".replace("^-", "^−")
    else:
        X["residuo_txt"] = "0" if res else "—"
    X["n_tests"] = contar_tests()
    X["segundos"] = float(np.median([f["segundos"] for f in filas["conf"]])) if filas["conf"] else float("nan")
    X["n_corridas"] = sum(len(v) for v in filas.values())
    # Tiempos medidos: suma de los segundos de cada corrida (un proceso cada una) y de las corridas con malla fina.
    X["horas_cpu"] = (sum(float(f["segundos"]) for v in filas.values() for f in v)
                      + sum(x["segundos_fina"] for x in X["sens"] + X["sens_ant"])) / 3600
    finas_s = [x["segundos_fina"] for x in X["sens"] if abs(x["dx_fina"] - 0.11) < 1e-6]
    X["segundos_fina"] = float(np.median(finas_s)) if finas_s else float("nan")
    X["hw"] = hardware()
    X["curvas_ctl"] = control_curvas(RAIZ / "data" / "hidratacion" / "source_terms.txt",
                                     float(cfg_fin["curvas"].get("desfase_h", 2.0)))
    X["eq1"] = ecuacion(r"$\rho\,c_p\,\dfrac{\partial T}{\partial t}=\nabla\cdot(k\,\nabla T)+\dot q_h(a)-\dot q_w,"
                        r"\qquad -k\,\dfrac{\partial T}{\partial n}=U\,(T-T_{eq}),"
                        r"\qquad U=\left(\dfrac{1}{h_{ext}}+R_{capas}+\dfrac{\Delta x}{2k}\right)^{-1}$",
                        IMG / "eq_modelo.png", tam=11)
    X["eq2"] = ecuacion(r"$M_\Delta=\max_{t\leq H}\left[\max_{\Omega_c}T-\min_{S_{25}}T\right],\qquad "
                        r"P_f=\sum_{m,k}p_m\,q_k\,\mathbf{1}\left[M_T>70\,^\circ\mathrm{C}\ \cup\ M_\Delta>20\,^\circ\mathrm{C}\right]$",
                        IMG / "eq_criterio.png", tam=11)
    X["eq3"] = ecuacion(r"$C_i\,\dfrac{T_i^{n+1}-T_i^n}{\Delta t}=\sum_j G_{ij}\,(T_j^n-T_i^n)+\sum_{caras}U A\,(T_{eq}-T_i^n)"
                        r"+\dot q_h V_i-\dot m\,c_w\,(T_{w,i+1}-T_{w,i}),\qquad \Delta t\leq 0{,}9\,\min_i\dfrac{C_i}{\sum G_i}$",
                        IMG / "eq_discreta.png", tam=10.5)
    X["eq4"] = ecuacion(r"$T_{w,i+1}=T_{c,i}+(T_{w,i}-T_{c,i})\,e^{-U_L L_i/(\dot m c_w)},\qquad "
                        r"U_L=\left[\dfrac{1}{h_w\pi D_i}+\dfrac{\ln(D_o/D_i)}{2\pi k_t}+\dfrac{\ln(r_{eq}/r_o)}{2\pi k}\right]^{-1}$",
                        IMG / "eq_tubos.png", tam=10.5)
    return X


def _frase_tcol(fr: dict) -> str:
    """Efecto de la temperatura de colocación antes y después del retiro (planes unica-M10-T{T}-N25 de la familia)."""
    pares = []
    for nombre, v in fr.items():
        m = re.fullmatch(r"unica-M10-T(\d+)-N25", nombre)
        if m and math.isfinite(v["antes_max"]) and math.isfinite(v["despues_max"]):
            pares.append((int(m.group(1)), v["antes_max"], v["despues_max"]))
    if len(pares) < 2:
        return ""
    pares.sort()
    Ts = "/".join(str(p[0]) for p in pares)
    return (f"Con tubos en el núcleo durante 14 días, colocar a {Ts} °C da un M~Δ~ máximo de "
            f"{'/'.join(f1(p[1]) for p in pares)} °C con la manta puesta y de {'/'.join(f1(p[2]) for p in pares)} °C "
            f"después del retiro: la temperatura de colocación mueve el pico de los primeros días "
            f"({f1((pares[-1][1] - pares[0][1]) / (pares[-1][0] - pares[0][0]), 2)} °C por °C), mientras el salto del "
            f"retiro cambia {f1((pares[-1][2] - pares[0][2]) / (pares[-1][0] - pares[0][0]), 2)} °C por °C.")


def _tabla_planes_familia(X):
    fam = X["fam"]
    orden = [p for p in ETIQUETAS if p in fam] + [p for p in fam if p not in ETIQUETAS]
    filas = []
    for p in orden:
        r = fam[p]
        Tp = X["prop_cfg"].get("T_colocacion") if p == X["prop"] else None
        et, Tc, manta, tub = ETIQUETAS.get(p, ("**Propuesta final** (tubos hasta el retiro)" if p == X["prop"] else p,
                                               f1(float(Tp), 0) if Tp is not None else "—", "10 cm",
                                               "V3–V5, programa"))
        filas.append([et, f1(r["D"]), Tc, manta, tub, pct(r["P_T"], 0), pct(r["P_D"], 0), pct(r["P_F"], 0),
                      "sí" if r["P_F"] <= 0.05 else "no"])
    return orden, filas


def _etiqueta(X, p):
    return ETIQUETAS.get(p, ("**Propuesta final** (tubos hasta el retiro)" if p == X["prop"] else p,))[0]


def _tabla_maximos_familia(X, orden):
    """Por plan: M_T y M_Δ máximos con su instante, ubicación y escenario."""
    filas = []
    for p in orden:
        r = X["fam"][p]
        pt, pd_ = r["peor_T"], r["peor"]
        filas.append([_etiqueta(X, p), f1(pt["M_T"]), f1(pt["t_MT_h"], 0), _zona(pt["z_MT"]),
                      f"{pt['curva']} {pt['clima']}", f1(pd_["M_D"]), f1(pd_["t_MD_h"], 0), _zona(pd_["z_Tmin"]),
                      f"{pd_['curva']} {pd_['clima']}"])
    return filas


def escribir(W: Informe, X: dict, R: dict):
    ref = lambda k: str(R.get(k, "?"))  # noqa: E731
    C, fam, figs, d = X["C"], X["fam"], X["figs"], X["delta"]
    geo, inicios, fines, D, t0 = X["geo"], X["inicios"], X["fines"], X["D"], X["t0"]
    prop, pc = X["prop"], X["prop_cfg"]
    t_ret, dia_ret = X["t_ret"], X["t_ret"] / 24
    ad, sens, mat, est = X["ad"], X["sens"], X["mat"], X["est"]
    tubo3 = next(t for t in pc["tubos"] if t["vertida"] == 3)
    H_geo = max(g["z1"] for g in geo)  # altura total de la geometría real (IFC)
    Tb = f1(float(pc.get("T_colocacion", 10)), 0)  # máximo a la descarga = temperatura de colocación del plan
    T_obj = float(pc.get("T_colocacion", 10)) - 2  # la planta apunta 2 °C bajo el máximo
    # Hielo por m³: ρc·(T_mezcla − T_obj) / (L + c_agua·T_obj), con la mezcla sin hielo a 16–18 °C y ρc del modelo.
    hielo = tuple(RHO_CP / 1000 * (Ts - T_obj) / (334 + 4.19 * T_obj) for Ts in (16.0, 18.0))
    V_tot = sum(g["V"] for g in geo)
    prog = programa_texto(tubo3)
    cumple = bool(C) and C["P_F"] <= 0.05
    pf = C["P_F"] if C else float("nan")
    n_sens = sum(s["n"] for s in X["sensores"]) + 6
    # Resultado de la propuesta en una frase (de los datos, sin suponer que cumple).
    if not C:
        frase_prop = "está en evaluación"
    elif C["k_F"] == 0:
        frase_prop = (f"cumple en los {C['n']} escenarios (P~f~ = 0 %), con M~Δ~ máximo de {f1(C['MD_max'])} °C y M~T~ "
                      f"máximo de {f1(C['MT_max'])} °C")
    elif cumple:
        frase_prop = (f"cumple el criterio con P~f~ = {pct(pf)} ({int(C['k_F'])} de {C['n']} escenarios superan 20 °C, "
                      f"con M~Δ~ máximo de {f1(C['MD_max'])} °C) y M~T~ máximo de {f1(C['MT_max'])} °C")
    else:
        frase_prop = (f"no cumple: P~f~ = {pct(pf)} ({int(C['k_F'])} de {C['n']} escenarios), sobre el 5 % admitido"
                      + (f"; con la malla del estudio, sin corrección, P~f~ sería {pct(X['pf_gruesa'])}"
                         if X.get("pf_gruesa") is not None else ""))
    # Planes de la familia que apagan los tubos antes del retiro: ¿alguno cumple?
    fr = X.get("fam_retiro", {})
    fr_ok = [n for n, v in fr.items() if v["P_F"] <= 0.05]
    # Planes cuyo peor M_Δ (con δ) aparece después de retirar la manta y supera 20 °C: en ellos manda el retiro.
    fr_gob = {n: v for n, v in fr.items() if v["despues_max"] > LIM_D and v["despues_max"] >= v["antes_max"]}
    nuc = [v[k] for v in fr_gob.values() for k in ("nucleo_min", "nucleo_max") if math.isfinite(v[k])]
    nuc_txt = f"{f1(min(nuc), 0)}–{f1(max(nuc), 0)} °C" if nuc else "—"
    if fr and not fr_ok and fr_gob:
        frase_retiro = ("El retiro de las protecciones, evaluado dentro del horizonte, gobierna el diseño: ningún plan "
                        f"que apague los tubos antes del retiro cumple, y en {len(fr_gob)} de {len(fr)} el peor "
                        f"diferencial aparece después de retirar la manta, con el núcleo aún a {nuc_txt}.")
    elif fr and not fr_ok:
        frase_retiro = ("Ningún plan que apague los tubos antes del retiro cumple; en todos ellos el peor diferencial "
                        "ocurre durante el pico de la primera semana, antes del retiro.")
    elif fr_gob:
        frase_retiro = ("El retiro de las protecciones, evaluado dentro del horizonte, es una fase crítica del diseño: en "
                        f"{len(fr_gob)} de {len(fr)} planes sin tubos al retirar, el peor diferencial aparece después del "
                        f"retiro (núcleo a {nuc_txt}).")
    else:
        frase_retiro = "El retiro de las protecciones se evalúa dentro del horizonte."
    X["_fr_ok"], X["_nuc_txt"] = fr_ok, nuc_txt
    dias_r = sorted({round((D + retiro_de(p)) / 24) for p in X["ret_cfg"].values() if math.isfinite(retiro_de(p))})
    dias_ret_txt = (", ".join(str(x) for x in dias_r[:-1]) + f" o {dias_r[-1]}") if len(dias_r) > 1 else \
        (str(dias_r[0]) if dias_r else "—")
    aguas = sorted({t["T_agua_C"] for p in X["ret_cfg"].values() for tb in p.get("tubos", [])
                    for t in tb.get("programa", [])[1:]}, reverse=True)
    aguas_txt = (", ".join(f1(a, 0) for a in aguas[:-1]) + f" o {f1(aguas[-1], 0)}") if len(aguas) > 1 else \
        (f1(aguas[0], 0) if aguas else "25")

    # ---------------- Título y resumen
    W.p("**Planificación probabilística de la vertida de una fundación de hormigón masivo con un solver térmico "
        "propio, BIM y 30 temporadas de clima real**", tam=16, alinear="centro", color="#184f95", despues=4)
    W.p(AUTORES, tam=11, alinear="centro", despues=1)
    W.p(f"{CURSO}" + (f" · Profesor: {PROFESOR}" if PROFESOR else ""), tam=9.5, alinear="centro", color=TINTA2,
        despues=1)
    W.p("Grupo 4 — sitio Litueche · Proyecto 2: Planificación de vertidas de hormigón masivo mediante simulación · "
        f"{MESES[datetime.now().month - 1]} de {datetime.now().year}", tam=9.5, alinear="centro", color=TINTA2, despues=1)
    W.p(f"Repositorio github.com/tbdiaz/Proyecto-Hormigon-Masivo · commit {git_commit()}", tam=9, alinear="centro",
        color=TINTA2, despues=10)
    resumen = (
        f"Se planifica la construcción de la fundación de un aerogenerador (636,8 m³ en cinco vertidas definidas en IFC) "
        f"para que la temperatura máxima no supere 70 °C y el diferencial entre el núcleo y la superficie, medido a "
        f"25 mm bajo las caras expuestas, no supere 20 °C. Un solver de volúmenes finitos propio, alimentado por el "
        f"IFC, resuelve la conducción con calor de hidratación por edad, colocación por capas, moldaje y mantas con su "
        f"retiro, tubos de enfriamiento y clima horario. La incertidumbre se representa con 8 curvas de hidratación "
        f"medidas y 30 temporadas reales de reanálisis ERA5 alineadas al inicio del {fecha_larga(t0).split(',')[0]}; "
        f"la probabilidad de falla se obtiene por enumeración ponderada de escenarios. Con malla fina en los escenarios "
        f"críticos (+{f1(d, 1)} °C en el resto), una vertida única continua de {f1(D)} h, con hormigón a {Tb} °C, manta de 10 cm y tubos en el "
        f"núcleo que siguen enfriando hasta el retiro de la manta el día {f1(dia_ret, 0)} {frase_prop}. "
        + frase_retiro + " Se entrega el plan de control con sensores, umbrales y contingencias.")
    W.p("**Resumen.** " + resumen, alinear="justificado")
    W.p("**Palabras clave:** hormigón masivo; diferencial térmico; volúmenes finitos; confiabilidad; tubos de "
        "enfriamiento; BIM/IFC.", despues=8)
    X["_resumen"] = resumen

    # ---------------- 1. Introducción
    W.h1("1. Introducción")
    W.p("En una fundación de hormigón masivo, el calor de hidratación del cemento se acumula en el núcleo, que se "
        "calienta durante días mientras la superficie se enfría hacia el aire. La diferencia de temperatura entre "
        "ambos produce tensiones de tracción en la superficie y es la causa principal de fisuración térmica temprana. "
        "Ju y Lei (2019) midieron en una losa de 18 m × 18 m y 3,3 m de espesor, sin control térmico, un máximo de "
        "≈73 °C en el núcleo al cuarto día y un diferencial núcleo–superficie de ≈35 °C. La Especificación para "
        "Hormigón Masivo Estructural del ICH (2019) exige para estos elementos un Plan de Control de Temperatura, "
        "propone el método de Schmidt (diferencias finitas) para estimar las temperaturas y limita el diferencial entre "
        "20 y 30 °C según la resistencia; ACI 301 fija 70 °C como máximo después de la colocación, y ACI 207.4R "
        "describe los controles disponibles: preenfriado, aislación y tubos de enfriamiento. Smolana et al. (2021) "
        "muestran, en losas de fundación reales, cómo la secuencia de construcción y la aislación reducen el riesgo de "
        "fisuras tempranas.", alinear="justificado")
    W.p("El proyecto adopta límites fijos y conservadores: T~máx~ ≤ 70 °C y ΔT~máx~ ≤ 20 °C durante todo el "
        "horizonte, con T~mín~ medida a 25 mm bajo las superficies exteriores. La pregunta de ingeniería es cómo "
        "ejecutar las vertidas desde el jueves 19 de noviembre de 2026 para mantener el riesgo térmico dentro de esos "
        "límites, considerando la incertidumbre del hormigón, de la ejecución y del clima. El objetivo es minimizar la "
        "duración de la colocación D con una probabilidad de falla P~f~ ≤ 5 %.", alinear="justificado")
    W.p("Este trabajo aporta: (i) un solver tridimensional propio, verificado, que lee la geometría y la secuencia "
        "desde IFC; (ii) un modelo de incertidumbre sin distribuciones inventadas, con curvas de hidratación medidas y "
        "temporadas climáticas reales; (iii) una evaluación que incluye el retiro de las protecciones dentro del "
        "horizonte, que resultó ser la fase crítica, y (iv) un plan de control derivado del modelo. La sección 2 "
        "describe el caso y los datos; las secciones 3 a 5, el modelo, su discretización y su verificación; la 6, el "
        "tratamiento de la incertidumbre; la 7, los resultados; la 8, la discusión y el plan de control, y la 9, las "
        "conclusiones.", alinear="justificado")

    # ---------------- 2. Caso y datos
    W.h1("2. Caso y datos")
    W.h2("2.1 Geometría y modelo BIM")
    W.p("La geometría es la oficial del curso, modelada en Bonsai/Blender y exportada a IFC: una losa octogonal de 18 m "
        "en dos vertidas de 0,55 m (V1, V2), dos troncos con caras inclinadas a ≈22° (V3, V4) y un pedestal cilíndrico "
        "de 7 m de diámetro y 1 m de alto (V5), con 4,30 m de altura total (Figura " + ref("f:geo") + "). Cada vertida "
        "es un IfcFooting con el property set UANDES_HormigonMasivo (VertidaID, Inicio_h, TemperaturaColocacion_C, "
        "TipoHormigon, Moldaje, Aislacion, RetiroMoldaje_h). El programa lee la geometría y las propiedades con "
        f"IfcOpenShell {IFC_VERSION} (IfcOpenShell, {IFC_VERSION}) y voxeliza cada vertida: una celda pertenece a la "
        "vertida que contiene su centro. Con "
        "Δx = 0,275 m el error de volumen total es +0,38 % (+0,06 % con 0,11 m).", alinear="justificado")
    W.figura(figs["geo"], "Las cinco vertidas del IFC y su volumen. Cada vertida conserva su color en todo el artículo.",
             ref="geo")
    W.h2("2.2 Sitio y fecha")
    pm = X["meta_era5"].get("punto_modelo", {})
    W.p(f"El sitio asignado es Litueche, Región de O'Higgins (−34,105457; −71,649539); el punto de grilla de ERA5 está "
        f"en ({f1(pm.get('latitud'), 2)}; {f1(pm.get('longitud'), 2)}; {f1(pm.get('elevacion_m'), 0)} m). La primera "
        f"vertida comienza el **{fecha_larga(t0)}** (hora local, America/Santiago, UTC−3), común a todos los planes, "
        f"como fija la pauta.", alinear="justificado")
    W.h2("2.3 Parámetros")
    filas = [
        ["Hormigón: k, ρ, c~p~", "2,222 W/mK; 2400 kg/m³; 970 J/kgK (ρc~p~ = 2,33 MJ/m³K)", "bibliográfico (Ju y Lei 2019)"],
        ["Probetas de las curvas", "350,77 kg/m³ de cemento; ρ = 2476 kg/m³; c~p~ = 1023 J/kgK", "medido (encabezado)"],
        ["Suelo: k, ρ, c~p~", "1,389 W/mK; 1800 kg/m³; 1273 J/kgK; 3 m bajo la base y 2 m de margen", "bibliográfico (Ju y Lei 2019)"],
        ["Suelo: T inicial y del fondo", "climatología ERA5 de las 4 capas de suelo a la hora de inicio", "medido (reanálisis)"],
        ["Moldaje", "terciado 18 mm, k = 0,12 W/mK (R = 0,15 m²K/W), caras verticales de V1, V2 y V5", "supuesto"],
        ["Manta", "10 cm, k = 0,04 W/mK (R = 2,5 m²K/W); α = 0,4; ε = 0,9", "supuesto (ficha típica)"],
        ["Superficie del hormigón", "absortividad 0,65; emisividad 0,9", "supuesto (rango típico)"],
        ["Convección", "h~c~ = 5,6 + 4v (v: viento a 10 m)", "bibliográfico (Jürges, 1924)"],
        ["Tubos", "HDPE Ø32 mm, pared 2 mm, k = 0,45 W/mK; 20 L/min por circuito", "supuesto (catálogo)"],
        ["Tasa de colocación", "27 m³/h, una bomba y una cuadrilla", "ajustado (Ju y Lei: ≈637 m³ en ≈24 h)"],
        ["Calor de hidratación", "8 curvas Q(t) por edad (sección 2.4)", "medido (probetas del curso)"],
    ]
    W.tabla("Parámetros del modelo y su origen.", ["Parámetro", "Valor", "Origen"], filas, [3.7, 9.0, 4.5],
            alinear=["l", "l", "l"], tam=8, ref="param")
    W.p("Las propiedades térmicas del hormigón y del suelo vienen del caso de Ju y Lei (2019), que tiene la misma "
        "geometría y un hormigón de dosis parecida. El encabezado de las curvas informa para las probetas una densidad "
        "y un calor específico mayores, con ρc~p~ = 2,53 MJ/m³K frente a 2,33 MJ/m³K del modelo. Con el mismo calor, un "
        "ρc~p~ menor da más aumento de temperatura, por lo que el valor del modelo es conservador; la sección 7.4 "
        "repite la propuesta con el ρc~p~ de las probetas. Los valores marcados como supuestos (moldaje, manta, "
        "superficie, tubos) son de catálogo o de rango típico y no se calibraron; la manta se especifica en el plan "
        "por su resistencia térmica (R ≥ 2,5 m²K/W) y no por su espesor, y la sección 7.4 prueba R = 2,0 y 3,0 m²K/W.",
        alinear="justificado")
    W.h2("2.4 Calor de hidratación")
    W.p("El calor proviene de las 8 curvas Q(t) entregadas por el profesor (Equipo docente, 2026; versión corregida "
        "del 14 de noviembre): cuatro probetas (M1–M4) en dos campañas, SS (240 h) y 22A (163 h, marcada de menor confianza), "
        "obtenidas por inversión con volúmenes finitos y regularización de Tikhonov con una conductancia de borde "
        "β = 2,5 W/m²K. Las curvas liberan 61–100 MJ/m³ en 14 días, equivalentes a 26–43 °C de aumento adiabático con "
        "el ρc~p~ del modelo (Figura " + ref("f:curvas") + "); las SS siguen liberando calor hasta ≈200 h y las 22A se "
        "apagan cerca de las 85 h. Cada celda recibe Q(a) según su propia edad a = t − t~colocación~, sin madurez; "
        "después del registro se continúa con una cola exponencial.", alinear="justificado")
    cc = X["curvas_ctl"]
    if cc:
        W.p(f"**Origen temporal y control.** Según el encabezado, el tiempo 0 del archivo es el inicio de la ventana "
            f"conservada, {f1(cc['desfase_h'], 0)} h después del inicio del registro. Se supone, como supuesto "
            f"declarado, que el registro empezó con el hormigón recién mezclado: la edad es a = t~archivo~ + "
            f"{f1(cc['desfase_h'], 0)} h y entre 0 y {f1(cc['desfase_h'], 0)} h se usa una rampa lineal de 0 a "
            f"Q({f1(cc['desfase_h'], 0)} h), que aporta a lo más {f1(cc['rampa_max'], 2)} MJ/m³ "
            f"({pct(cc['rampa_frac'], 1)} de la energía de la curva). Si el registro hubiera empezado más tarde, la "
            f"fuente debería desplazarse más. La integral de cada curva reproduce la energía que informa el "
            f"encabezado para β = 2,5 W/m²K con un error máximo de {pct(cc['err_max'], 2)} en las "
            f"{cc['n_tabla']} curvas.", alinear="justificado")
    W.p("**Por qué equiprobables.** El encabezado informa una sola dosificación (350,77 kg/m³ de cemento) y la misma "
        "probeta y método de inversión para las cuatro probetas, por lo que se tratan como réplicas de una misma "
        "mezcla ensayadas en dos campañas: su dispersión representa la variabilidad entre probetas y entre campañas, "
        "no mezclas que se puedan elegir. Sin información para preferir una sobre otra, cada curva pesa 1/8. La "
        "campaña 22A es de menor confianza (su único sensor de ambiente no reproduce el ciclo diario); por eso se "
        "informa además el resultado con sólo las 4 curvas SS (sección 7.4).", alinear="justificado")
    W.figura(figs["curvas"], "Las 8 curvas de hidratación (arriba) y el aumento adiabático que implican con "
                             "ρc~p~ = 2,33 MJ/m³K (abajo).", ref="curvas")
    W.h2("2.5 Clima")
    me = X["meta_ens"]
    W.p(f"El clima de los escenarios son las **{X['n_temporadas']} temporadas reales 1996–2025** del reanálisis ERA5 "
        f"(Hersbach et al., 2020; Open-Meteo Historical Weather API), alineadas al 19 de noviembre de cada año: "
        f"temperatura del aire, viento a 10 m, radiación de onda corta y nubosidad, cada hora. Son bloques históricos "
        f"comparables, que conservan el ciclo diario, los frentes y la coherencia entre variables. Para la familia de "
        f"planes se agrupan en terciles de un índice de frío (media de las mínimas de los primeros 7 días) y se "
        f"representa cada tercil por su temporada mediana, con peso 1/3; la propuesta final se evalúa con las 30. El "
        f"pronóstico más reciente, el ensemble ECMWF IFS de 51 miembros (ECMWF, 2026) descargado el "
        f"{str(me.get('descargado_utc', '—'))[:10]}, sirve de verificación para los primeros días (Figura "
        f"{ref('f:clima')}).", alinear="justificado")
    W.figura(figs["clima"], "Temperatura del aire y viento de las 30 temporadas ERA5 alineadas al inicio, y el "
                            "pronóstico ensemble más reciente para 2026.", ref="clima")

    # ---------------- 3. Modelo continuo
    W.h1("3. Modelo continuo")
    W.p("La temperatura T(x, t) en el hormigón ya colocado Ω~c~(t) y en el suelo Ω~s~ cumple el balance de energía de "
        "un sólido rígido con fuente volumétrica: la variación de energía interna de un volumen es igual al flujo de "
        "calor que entra por conducción (ley de Fourier), más el calor de hidratación, menos lo que extraen los tubos:",
        alinear="justificado")
    W.ecuacion(X["eq1"])
    W.p("donde q̇~h~(a) es la curva del escenario evaluada en la edad de cada punto y q̇~w~ el sumidero de los tubos. "
        "Cada capa de hormigón aparece a su temperatura de colocación con edad cero; el suelo parte del perfil ERA5. "
        "En las caras expuestas rige una condición de Robin en serie con las capas de protección; T~eq~ es la "
        "temperatura sol-aire que reúne convección, onda larga hacia el cielo y radiación solar absorbida. Las "
        "condiciones cambian por etapa, como resume la Tabla " + ref("t:bordes") + ".", alinear="justificado")
    filas = [
        ["Cara superior de la última capa colocada", "Robin sol-aire; con manta, R = 2,5 m²K/W en serie", "colocación → retiro"],
        ["Caras verticales de V1, V2, V5", "Robin + moldaje (R = 0,15) + manta (R = 2,5)", "colocación → retiro"],
        ["Caras inclinadas de V3, V4", "Robin + manta", "colocación → retiro"],
        ["Cualquier cara, después del retiro", "Robin sol-aire sin capas", f"día {f1(dia_ret, 0)} → fin"],
        ["Capa sobre capa; hormigón–suelo", "continuidad (k armónico)", "siempre"],
        ["Superficie del suelo", "Robin sol-aire", "siempre"],
        ["Fondo del suelo (3 m)", "Dirichlet: T profunda de ERA5", "siempre"],
        ["Costados del suelo", "flujo nulo", "siempre"],
    ]
    W.tabla("Condiciones de borde por superficie y etapa (normal saliente).", ["Superficie", "Condición", "Vigencia"],
            filas, [5.6, 7.6, 4.0], alinear=["l", "l", "l"], tam=8, ref="bordes")
    W.p("**Ambiente.** Cada mecanismo entra una sola vez. La convección usa h~c~ = 5,6 + 4v con el viento horario; la "
        "onda larga se linealiza como h~r~ = 4εσT~m~³ hacia una temperatura radiante que mezcla el cielo y el entorno "
        "según la orientación de la cara (fracción de cielo 1 en caras superiores y 0,5 en laterales); la temperatura "
        "del cielo sigue la fórmula de Swinbank (1963) para cielo despejado, corregida por la nubosidad horaria; la radiación "
        "solar αG entra sólo en las caras horizontales hacia arriba, porque ERA5 entrega la radiación sobre un plano "
        "horizontal. Así, h~ext~ = h~c~ + h~r~ y T~eq~ = (h~c~ T~aire~ + h~r~ T~rad~ + αG)/h~ext~. En noches despejadas "
        "T~eq~ queda varios grados bajo el aire, y ese es el mecanismo que enfría la superficie después del retiro de la "
        "manta.", alinear="justificado")
    W.p("**Condiciones iniciales.** El suelo parte de un perfil que interpola, en profundidad, las cuatro capas de "
        "suelo de ERA5 a la hora de inicio (climatología de las temporadas), y su fondo, a 3 m bajo la base, se fija a "
        "la temperatura profunda. El hormigón no existe hasta que se coloca: cada celda se activa cuando la capa que la "
        "contiene llega a su cota, con la temperatura de colocación del plan y edad cero, y desde ese instante recibe "
        "su propio calor de hidratación.", alinear="justificado")
    W.p("**Tubos.** Un tubo de 32 mm es mucho más delgado que una celda, por lo que se modela como un sumidero lineal: "
        "el agua entra a T~w,0~ y se calienta a lo largo del recorrido, intercambiando con cada celda que atraviesa "
        "según una resistencia agua–pared–hormigón por metro (convección interior de Dittus y Boelter, 1930; conducción en la "
        "pared, y el radio equivalente de Peaceman, 1978: r~eq~ = 0,198 Δx):", alinear="justificado")
    W.ecuacion(X["eq4"])
    W.p("La temperatura de entrada puede seguir un **programa por tramos**, por ejemplo 25 °C la primera semana y "
        "luego más fría, para seguir enfriando el núcleo antes de retirar la aislación. Las métricas de cada escenario "
        "son la temperatura máxima M~T~ y el diferencial máximo M~Δ~ en todo el horizonte H, que llega hasta 7 días "
        "después del último retiro de protecciones, de modo que el diferencial que induce el retiro forma parte de la "
        "evaluación:", alinear="justificado")
    W.ecuacion(X["eq2"])
    W.p("S~25~ son los puntos a 25 mm bajo cada cara expuesta, incluidas las cubiertas por moldaje o manta, en caras "
        "planas; las aristas de la escalera de vóxeles se informan aparte. P~f~ suma, sin contar dos veces, los pesos "
        "de los escenarios que exceden cualquiera de los dos límites.", alinear="justificado")

    # ---------------- 4. Discretización
    W.h1("4. Discretización y algoritmos")
    W.p("El dominio se divide en vóxeles cúbicos de lado Δx y cada celda i cumple un balance de volúmenes finitos con "
        "Euler explícito:", alinear="justificado")
    W.ecuacion(X["eq3"])
    W.p("Entre celdas vecinas la conductancia es G~ij~ = k~arm~ A/Δx; en las caras expuestas U A, con U de la sección "
        "3, y en el fondo del suelo una conductancia hacia la temperatura fija. El paso se recalcula en cada evento "
        "(capa colocada, retiro, cambio del agua) y nunca supera 600 s por precisión. Los eventos se ordenan y la "
        "integración avanza de uno a otro con la geometría fija; dentro de cada tramo las conductancias no cambian, lo "
        "que permite vectorizar el paso completo con NumPy (Harris et al., 2020). El algoritmo es:",
        alinear="justificado")
    for linea in ("1. Leer IFC → voxelizar → asignar suelo, aire y vertidas.",
                  "2. Cronograma: cada etapa sube por capas a 27 m³/h; t_act de cada celda.",
                  "3. Eventos = colocaciones ∪ retiros ∪ tramos del agua ∪ guardados.",
                  "4. Para cada tramo entre eventos:",
                  "     a. armar G, U·A por categoría de cara, sumideros de los tubos y Δt ≤ 0,9 Δt_crít;",
                  "     b. repetir: bordes + conducción + hidratación + agua a lo largo del tubo → T^(n+1);",
                  "     c. registrar Tmax, Tmin a 25 mm, balance de energía.",
                  "5. Métricas M_T, M_Δ, ubicación de los extremos y residuo de energía."):
        par = W.p(linea, tam=8.5, despues=0)
        for r in par.runs:
            r.font.name = "Consolas"
    n_finas = len(sens) + len(X["sens_ant"])
    W.p(f"**Costo.** Las corridas se hicieron en computadores portátiles (el principal: {X['hw']}). Con Δx = 0,275 m el dominio tiene ≈100 mil "
        f"celdas y un escenario con retiro (≈{f1(t_ret / 24 + 7, 0)} días) toma una mediana de "
        f"{f1(X['segundos'], 0)} s en un proceso; con Δx = 0,11 m, {f1(X['segundos_fina'] / 60, 0)} min. Los "
        f"escenarios son independientes y se corren 6 en paralelo. Las {X['n_corridas']} simulaciones de los estudios "
        f"y las {n_finas} con malla fina suman {f1(X['horas_cpu'], 0)} h de proceso (medidas por el programa en cada "
        f"corrida), es decir, del orden de {f1(X['horas_cpu'] / 6, 0)} h de reloj con 6 procesos.",
        alinear="justificado", despues=6)
    W.p("**Paso de tiempo.** Euler explícito es estable si Δt ≤ C~i~/ΣG~i~ en todas las celdas. En el interior del "
        "hormigón, con Δx = 0,275 m, ese límite es Δx²/(6α) ≈ 3,7 h, y las celdas de borde con manta lo bajan poco, "
        "porque la manta agrega resistencia; por eso el paso lo fija la precisión (600 s) y no la estabilidad. Con "
        "0,11 m el límite interior baja a ≈35 min, todavía sobre 600 s, y el costo crece con el número de celdas, "
        "(0,275/0,11)³ ≈ 15,6 veces. El paso de 600 s se verificó contra 300 s (sección 5).", alinear="justificado")
    W.p("**Estudio.** Un segundo algoritmo recorre planes y escenarios. Cada simulación tiene una firma con el plan, "
        "la curva, el clima y los parámetros que la afectan; si ya está en la caché, no se repite:",
        alinear="justificado")
    for linea in ("para cada plan u declarado en el archivo de configuración:",
                  "     para cada curva m (p_m) y cada clima k (q_k):",
                  "          firma ← hash(u, m, k, parámetros que afectan); simular si no está en caché",
                  "          guardar M_T, M_Δ, extremos, serie horaria y residuo de energía",
                  "     P_f(u) ← Σ p_m q_k · 1[M_T > 70 ∪ M_Δ + δ > 20]",
                  "elegir u* = el de menor D(u) con P_f(u) ≤ 5 %; empate → mayor margen 20 − máx M_Δ"):
        par = W.p(linea, tam=8.5, despues=0)
        for r in par.runs:
            r.font.name = "Consolas"
    W.p("", despues=4)

    # ---------------- 5. Verificación
    W.h1("5. Verificación")
    W.p("**Verificación, calibración y validación.** Verificar es comprobar que el código resuelve bien las ecuaciones "
        "del modelo; calibrar, ajustar parámetros para reproducir mediciones; validar, comparar las predicciones con "
        "mediciones independientes del sistema real. Este trabajo hace verificación: tests contra soluciones exactas, "
        "balance de energía y sensibilidad numérica. No calibra: los parámetros son medidos, bibliográficos o "
        "supuestos declarados (Tabla " + ref("t:param") + "), y las curvas vienen de una inversión sobre probetas, no "
        "de un ajuste del modelo de la fundación. Y **no hay validación integral**, porque no existen mediciones "
        "térmicas de esta fundación: la comparación con Ju y Lei (2019) es sólo de orden de magnitud, y la validación "
        "queda a cargo del monitoreo en obra (sección 8.2), que compara lo medido con lo previsto.",
        alinear="justificado")
    W.p(f"El código tiene {X['n_tests']} tests automáticos. Los principales comparan con soluciones exactas: conducción "
        f"1D, 2D y 3D (orden 2 en espacio y 1 en tiempo), una caja con moldaje y desmolde contra un solver de "
        f"referencia, un bloque adiabático cuya temperatura debe seguir la integral de la curva, el agua a lo largo del "
        f"tubo contra la solución exponencial, y el programa del agua por tramos (un tramo da lo mismo que agua "
        f"constante; con agua a la temperatura del hormigón no hay intercambio). Cada simulación cierra su balance de "
        f"energía con un residuo relativo del orden de 10^−15^.", alinear="justificado")
    finas = sorted((x for x in sens if x["plan"] == prop), key=lambda x: -x["MD_fina"])
    detalle = sorted((x for x in finas if x["otras"] or x["intermedias"]), key=lambda x: x["t_MD_fina"])
    criticos = [x for x in finas if x not in detalle]
    W.p(f"El residuo de energía compara el cambio de energía interna con el calor liberado menos lo que salió por los "
        f"bordes y el agua; en los {C['n'] if C else 0} escenarios de la propuesta es a lo más {X['residuo_txt']}, de "
        "modo que las diferencias entre planes son físicas y no numéricas. Como la métrica de diseño depende de un "
        "gradiente fuerte en los primeros centímetros bajo la cara, se mide además la sensibilidad a la malla, al paso "
        "de tiempo y a la profundidad del punto frío "
        + ("en los escenarios que gobiernan las dos fases (el pico con la manta puesta y el salto después del retiro)"
           if len({x["t_MD_fina"] < t_ret for x in detalle}) == 2 else
           "en el escenario de mayor M~Δ~ ("
           + ("el pico con la manta puesta" if detalle and detalle[0]["t_MD_fina"] < t_ret else "el salto del retiro")
           + ")")
        + (f" y, sólo con malla fina, en los {len(criticos)} escenarios siguientes más cercanos al límite"
           + (", entre ellos el de mayor M~Δ~ después del retiro"
              if any(x["t_MD_fina"] >= t_ret for x in criticos) and all(x["t_MD_fina"] < t_ret for x in detalle)
              else "") if criticos else "")
        + ".", alinear="justificado")
    if finas:
        def fase(x):
            return "antes del retiro" if x["t_MD_fina"] < t_ret else "después del retiro"

        def dx_txt(v):
            return f1(v, 4) if round(v * 1e3, 6) % 1 else f1(v, 3)

        filas = []
        for x in detalle:
            esc = f"{x['curva']} {x['clima'].split('_')[-1]}"
            filas.append([esc, fase(x), "0,275", "600", "25", f1(x["MT_gruesa"], 2), f1(x["MD_gruesa"], 2), ""])
            for m in sorted(x["intermedias"], key=lambda m: -m["dx"]):
                filas.append(["", "", dx_txt(m["dx"]), "600", "25", f1(m["M_T"], 2), f1(m["M_D"], 2),
                              f"{'+' if m['delta'] >= 0 else '−'}{f1(abs(m['delta']), 2)}"])
            filas.append(["", "", dx_txt(x["dx_fina"]), "600", "25", f1(x["MT_fina"], 2), f1(x["MD_fina"], 2),
                          f"+{f1(x['delta'], 2)}"])
            filas += [["", "", dx_txt(float(r["dx"])), f1(float(r["dt_s"]), 0), f1(1000 * float(r["profundidad_m"]), 0),
                       f1(float(r["M_T"]), 2), f1(float(r["M_D"]), 2), ""] for r in x["otras"]]
        if detalle:
            W.tabla("Sensibilidad numérica en "
                    + ("los escenarios de la propuesta que gobiernan cada fase"
                       if len({x["t_MD_fina"] < t_ret for x in detalle}) == 2 else
                       "el escenario de la propuesta de mayor M~Δ~")
                    + " (todo el "
                    "horizonte, incluido el retiro): malla, paso de tiempo y profundidad de T~mín~, por separado. La "
                    "fase es la del instante de M~Δ~; δ es la diferencia con la malla de 0,275 m.",
                    ["Escenario", "Fase", "Δx [m]", "Δt [s]", "Prof. [mm]", "M~T~ [°C]", "M~Δ~ [°C]", "δ [°C]"], filas,
                    [2.4, 2.6, 1.6, 1.4, 1.6, 1.9, 1.9, 1.6], alinear=["l", "l"] + ["c"] * 6, tam=7.5, ref="malla")
        dts = [abs(float(r["M_D"]) - x["MD_gruesa"]) for x in detalle for r in x["otras"] if float(r["dt_s"]) != 600]
        prs = [abs(float(r["M_D"]) - x["MD_gruesa"]) for x in detalle for r in x["otras"]
               if abs(float(r["profundidad_m"]) - 0.025) > 1e-6]
        tendencia = []
        for x in detalle:
            ms = sorted(x["intermedias"], key=lambda m: -m["dx"])
            if not ms:
                continue
            pasos = [ms[0]["delta"]] + [b["delta"] - a["delta"] for a, b in zip(ms, ms[1:])] \
                + [x["delta"] - ms[-1]["delta"]]
            dxs = [0.275] + [m["dx"] for m in ms] + [x["dx_fina"]]
            txt = (f"en {x['curva']} {x['clima'].split('_')[-1]} ({fase(x)}) δ vale "
                   + ", ".join(f"{f1(m['delta'], 2)} °C con {dx_txt(m['dx'])} m" for m in ms)
                   + f" y {f1(x['delta'], 2)} °C con {dx_txt(x['dx_fina'])} m")
            if all(pq > 0 for pq in pasos):
                txt += (", con incrementos de " + " y ".join(f"{f1(pq, 2)} °C" for pq in pasos)
                        + " en pasos de malla de " + " y ".join(f"{f1(a - b, 3)} m" for a, b in zip(dxs, dxs[1:]))
                        + (" (el último incremento es menor: δ tiende a estabilizarse)" if pasos[-1] < pasos[-2]
                           else " (el último incremento no es menor: δ todavía no se estabiliza)"))
            else:
                txt += " (la secuencia no es monótona, por lo que no se le asigna un orden de convergencia)"
            tendencia.append(txt)
        if detalle:
            W.p("La malla gruesa subestima el enfriamiento junto a la superficie: refinar a "
                f"{f1(detalle[0]['dx_fina'], 2)} m sube M~Δ~ en "
                + " y ".join(f"{f1(x['delta'], 2)} °C ({fase(x)})" for x in detalle)
                + ", mientras M~T~ casi no cambia"
                + (f"; el paso de 300 s cambia M~Δ~ en a lo más {f1(max(dts), 2)} °C" if dts else "")
                + (f" y medir a 50 mm en vez de 25 mm, en {f1(max(prs), 2)} °C" if prs else "") + ". "
                + ("Con tres mallas, " + "; ".join(tendencia) + ". Tres resoluciones bastan para ver la tendencia, no "
                   "para afirmar un orden de convergencia." if tendencia else ""),
                alinear="justificado")
        if criticos:
            dxc = f1(criticos[0]["dx_fina"], 2)
            filas = [[f"{x['curva']} {x['clima'].split('_')[-1]}", fase(x), f1(x["MD_gruesa"], 2),
                      f1(x["MD_gruesa"] + d, 2), f1(x["MD_fina"], 2), f"+{f1(x['delta'], 2)}",
                      "sí" if x["MD_fina"] > LIM_D else "no"] for x in criticos]
            W.tabla(f"Malla de {dxc} m en los escenarios de la propuesta más cercanos al límite, además de los de "
                    f"la Tabla {ref('t:malla')}: M~Δ~ con la malla del estudio, con la corrección uniforme "
                    f"(+{f1(d, 2)} °C) y con la malla fina.",
                    ["Escenario", "Fase", "M~Δ~ 0,275 m", "+ δ", f"M~Δ~ {dxc} m", "δ [°C]", "¿> 20 °C?"], filas,
                    [2.6, 2.8, 2.0, 1.6, 2.0, 1.6, 1.8], alinear=["l", "l"] + ["c"] * 5, tam=7.5, ref="malla2")
        ks = sum(1 for x in finas if x["MD_fina"] > LIM_D)
        ds = [x["delta"] for x in finas]
        W.p(f"En los {len(finas)} escenarios corridos con malla fina, δ va de {f1(min(ds), 2)} a {f1(max(ds), 2)} °C y "
            + ("ninguno supera 20 °C" if ks == 0 else f"{ks} superan 20 °C") + ". Por eso **en esos escenarios se usa "
            f"directamente M~Δ~ de la malla fina, y a los demás se les suma δ = {f1(d, 2)} °C** antes de compararlos "
            "con 20 °C: la mayor diferencia medida en todas las corridas con malla fina"
            + (f", que se dio en la propuesta anterior ({X['delta_ant_txt']}). Se usa la mayor porque es un error del "
               "método numérico y no del plan" if X.get("delta_ant_txt") and X["delta_ant_max"] > max(ds) + 1e-9
               else "") + ".", alinear="justificado")
        td = X.get("tabla_delta") or {}
        if td:
            filas = [[f1(r["delta"], 2) + (" (medida)" if abs(r["delta"] - round(d, 2)) < 1e-9 else ""),
                      f"{pct(r['pf_unif'])} ({r['k_unif']})", f"{pct(r['pf_mixto'])} ({r['k_mixto']})",
                      "sí" if r["pf_mixto"] <= 0.05 else "no"] for r in td["filas"]]
            W.tabla(f"¿Puede el error de malla cambiar la clasificación? P~f~ de la propuesta (entre paréntesis, "
                    f"escenarios que superan los límites, de {td['n']}) según la corrección δ: sumada a los {td['n']} "
                    f"escenarios, o con malla fina directa en {td['n_finos']} y δ en el resto (lo que usa el "
                    f"artículo). ¿Cumple? se refiere a la segunda.",
                    ["δ [°C]", f"P~f~, δ en los {td['n']}", f"P~f~, malla fina en {td['n_finos']} + δ", "¿Cumple?"],
                    filas, [3.4, 4.2, 5.0, 2.0], alinear=["c"] * 4, tam=8, ref="delta")
            cu, cm = td["crit_unif"], td["crit_mixto"]
            X["_delta_crit"] = cm
            cola = (f" Una malla más fina que {f1(finas[-1]['dx_fina'], 2)} m no se corrió: con dos resoluciones"
                    + (" (tres en los escenarios que gobiernan)" if tendencia else "")
                    + " se evalúa la sensibilidad, y esta tabla acota cuánto podría cambiar la decisión.")
            if cumple:
                W.p("La tabla responde la pregunta que importa para la decisión: "
                    + (f"con lo que usa el artículo, P~f~ supera 5 % sólo si la corrección de los escenarios sin malla "
                       f"fina llega a **δ~crít~ = {f1(cm, 2)} °C**, {f1(cm / d, 1)} veces la mayor medida" if cm else
                       "con lo que usa el artículo, ninguna corrección de hasta 8 °C lleva P~f~ sobre 5 %")
                    + (f"; si en cambio se sumara δ a los {td['n']} escenarios, bastaría δ~crít~ = {f1(cu, 2)} °C"
                       if cu else "")
                    + ". Esa es la **condición de validez** de la conclusión: la propuesta cumple mientras el error de "
                    "malla que queda en los escenarios no refinados sea menor que δ~crít~." + cola,
                    alinear="justificado")
            else:
                W.p(f"La tabla responde la pregunta que importa para la decisión, y la respuesta es que **sí cambia la "
                    f"clasificación**: con la corrección medida (δ = {f1(d, 2)} °C) la propuesta **no cumple** "
                    f"(P~f~ = {pct(pf)}). "
                    + (f"Cumpliría sólo si el error de malla de los escenarios sin malla fina fuese menor que "
                       f"δ~crít~ = {f1(cm, 2)} °C, por debajo de la mediana de las diferencias medidas "
                       f"({f1(float(np.median(ds)), 2)} °C)" if cm is not None and cm < float(np.median(ds)) else
                       f"Cumpliría sólo si el error de malla de los escenarios sin malla fina fuese menor que "
                       f"δ~crít~ = {f1(cm, 2)} °C" if cm is not None else "")
                    + (f". De las {len(ds)} diferencias medidas, {sum(1 for x in ds if x < cm)} son menores que "
                       f"δ~crít~; las que lo superan son de la temporada "
                       f"{', '.join(sorted({x['clima'].split('_')[-1] for x in finas if x['delta'] >= cm}))}, y nada "
                       "indica qué otras temporadas no refinadas se comportarían igual, así que se mantiene la regla "
                       "conservadora declarada antes de ver estos resultados."
                       if cm is not None else ".") + cola, alinear="justificado")
    else:
        W.p("La corrida con malla fina no está disponible: las M~Δ~ se informan sin corrección de malla.",
            alinear="justificado")

    # ---------------- 6. Incertidumbre
    W.h1("6. Incertidumbre y optimización")
    W.p("**Decisión y objetivo.** Un plan u fija la partición en etapas, las esperas, la temperatura de colocación, "
        "las protecciones y su retiro, y los tubos con su programa. Con una bomba, la duración D(u) es el tiempo hasta "
        "la última colocación; se busca el menor D(u) con P~f~(u) ≤ 5 %, con los mismos escenarios y pesos para todos "
        "los planes.", alinear="justificado")
    W.p("**Escenarios y pesos.** El material se representa por las 8 curvas medidas, equiprobables (p~m~ = 1/8; sección 2.4): cada "
        "escenario aplica una curva completa a todo el hormigón. El clima, por las temporadas ERA5: 3 terciles con "
        "q~k~ = 1/3 para comparar planes y las 30 temporadas con q~k~ = 1/30 para la propuesta. Material y clima se "
        "suponen independientes, w~mk~ = p~m~ q~k~. La ejecución (temperatura de colocación, interrupciones, fallas de "
        "equipos) se trata con variantes del plan y contingencias, no con distribuciones.", alinear="justificado")
    W.p("**Probabilidad.** P~f~ es la enumeración ponderada de esos escenarios: es la probabilidad exacta del modelo "
        "discreto declarado, no una estimación muestral, por lo que no se le asignan intervalos binomiales ni se "
        "remuestrean las curvas. Su resolución es el peso de un escenario: 1/24 = 4,2 % en la familia y 1/240 = 0,42 % "
        "en la propuesta. La robustez frente a lo que el modelo no representa se examina con sensibilidades "
        "deterministas: sólo la campaña SS (22A tiene menor confianza), el ρc~p~ medido en las probetas, curvas con "
        "+20 % de calor (cota de la inversión con β = 3,0 W/m²K), la temperatura de colocación y la resistencia de la "
        "manta.", alinear="justificado")
    W.p("**Búsqueda.** Las familias de planes se declararon en los archivos de configuración antes de ver sus "
        "resultados: (a) niveles de protección de la vertida única y la variante más lenta de 5 vertidas con esperas "
        f"de 12 y 24 h; (b) temperatura de colocación de 8 a 16 °C en pasos de 2 °C; (c) en {len(X['rondas'])} rondas, "
        f"el retiro de la manta los días {dias_ret_txt}, con los tubos del núcleo funcionando hasta el retiro y el agua "
        f"bajando por tramos desde 25 °C hasta {aguas_txt} °C. La vertida única es el límite inferior de D con una "
        f"bomba (D = {f1(D)} h): todo plan sin esperas es la misma colocación, por lo que, entre los planes que "
        "cumplen, la decisión se desempata por el margen frente al límite.", alinear="justificado")
    W.p("**Por qué enumeración.** Con 8 curvas y 30 temporadas el espacio de escenarios es finito y se recorre "
        "completo, sin error de muestreo. Un Monte Carlo sobre distribuciones ajustadas a 8 probetas agregaría supuestos "
        "(forma de la distribución, correlaciones) que los datos no permiten verificar, y su cola, que es lo que "
        "importa para P~f~ ≤ 5 %, sería sobre todo una extrapolación. Lo que queda fuera de ese conjunto se trata con "
        "sensibilidades.", alinear="justificado")
    W.p("**Margen frente al material.** Un hormigón más exotérmico que las 8 curvas no está representado. Para medir "
        "el margen se repite la propuesta con las curvas escaladas en +20 % (la diferencia entre β = 2,5 y 3,0 W/m²K) "
        "y se ajusta, para las curvas SS, el peor M~Δ~ contra la energía a 14 días: la energía con la que esa recta "
        "llega a 20 °C es E~crít~, una cota determinista que se controla con los sensores de las primeras 48 h.",
        alinear="justificado")

    # ---------------- 7. Resultados
    W.h1("7. Resultados")
    W.h2("7.1 Familia de planes")
    if fam:
        orden, filas = _tabla_planes_familia(X)
        W.tabla("Familia de planes con 8 curvas × 3 terciles climáticos (24 escenarios por plan), evaluada hasta 7 días "
                f"después del retiro de las protecciones, con M~Δ~ corregido por malla (+{f1(d, 2)} °C): "
                "probabilidades marginales (P~T~, P~Δ~) y conjunta (P~f~).",
                ["Plan", "D [h]", "T col. [°C]", "Manta", "Tubos", "P~T~", "P~Δ~", "P~f~", "¿Cumple?"], filas,
                [4.8, 1.1, 1.2, 1.2, 2.9, 1.2, 1.2, 1.2, 1.7], alinear=["l"] + ["c"] * 8,
                destacar=tuple(i for i, p in enumerate(orden) if p == prop), tam=7.5, ref="familia")
        W.tabla("Máximos de cada plan de la familia: valor, instante (horas desde t~0~), ubicación y escenario (curva y "
                "tercil climático). La ubicación de M~Δ~ es la del punto frío a 25 mm bajo la cara.",
                ["Plan", "M~T~ [°C]", "t [h]", "Ubicación", "Escenario", "M~Δ~ [°C]", "t [h]", "Punto frío",
                 "Escenario"], _tabla_maximos_familia(X, orden),
                [4.4, 1.2, 0.9, 1.9, 2.0, 1.2, 0.9, 1.9, 2.0], alinear=["l"] + ["c"] * 8,
                destacar=tuple(i for i, p in enumerate(orden) if p == prop), tam=7, ref="maximos")
        if "dpf" in figs:
            W.figura(figs["dpf"], "Duración de la colocación D frente a P~f~ en la familia de planes, con la línea del "
                                  "5 %. La vertida única (D = " + f1(D) + " h) es la colocación más rápida posible con "
                                  "una bomba; los planes con la misma D se separan levemente para que se vean.",
                     ref="dpf")
        n_ok = [p for p in fam if fam[p]["P_F"] <= 0.05]
        pt_max = max(fam[p]["P_T"] for p in fam)
        W.p(f"Con las curvas medidas, la temperatura máxima nunca es el problema "
            f"({'P~T~ = 0 en todos' if pt_max == 0 else f'P~T~ ≤ {pct(pt_max)}'} y M~T~ ≤ {f1(X['MT_familia'])} °C); "
            f"gobierna el diferencial (Tablas {ref('t:familia')} y {ref('t:maximos')}, Figura {ref('f:dpf')}). Al "
            f"incluir el retiro de las protecciones en la evaluación, "
            + (f"sólo {len(n_ok)} plan{'es' if len(n_ok) != 1 else ''} cumple{'n' if len(n_ok) != 1 else ''}: "
               + ", ".join("la propuesta final" if p == prop else p for p in n_ok) + ". " if n_ok else "ningún plan cumple. ")
            + (f"De los {len(fr)} planes que apagan los tubos a los 14 días (o no los tienen) y retiran la manta el día "
               f"{f1(min(v['dia'] for v in fr.values()), 0)}, el mayor M~Δ~ con la manta puesta es "
               f"{f1(min(v['antes_max'] for v in fr.values()))}–{f1(max(v['antes_max'] for v in fr.values()))} °C, y "
               f"después del retiro, {f1(min(v['despues_max'] for v in fr.values()))}–"
               f"{f1(max(v['despues_max'] for v in fr.values()))} °C. " if fr else "")
            + (f"En {len(fr_gob)} de ellos el peor diferencial aparece después del retiro: "
               f"el núcleo sigue a {nuc_txt} al retirar la manta y la superficie cae hacia el aire nocturno. "
               if fr_gob else "")
            + (f"En los otros {len(fr) - len(fr_gob)} el máximo ocurre antes, durante el pico de la primera semana. "
               if fr and len(fr_gob) < len(fr) else "")
            + "Las variantes más lentas son peores: mientras se espera, las capas inferiores se enfrían hacia el suelo y "
            "la siguiente se calienta encima, lo que suma un diferencial entre capas al del núcleo con la superficie. "
            + _frase_tcol(fr), alinear="justificado")
    W.h2("7.2 Retiro de las protecciones")
    filas = []
    for nombre, p in X["ret_cfg"].items():
        if nombre not in X["ret"]:
            continue
        r, a = X["ret"][nombre], X["ret_ad"].get(nombre, {})
        tub = next((t for t in p["tubos"] if t["vertida"] == 3), {})
        filas.append([f"R{X['ronda_de'][nombre]} · {nombre}", f"día {f1((D + retiro_de(p)) / 24, 0)}",
                      " → ".join(f1(t["T_agua_C"], 0) for t in tub.get("programa", [])) or f1(tub.get("T_agua_C", 25), 0)
                      + (" + V1–V2" if any(t["vertida"] in (1, 2) for t in p["tubos"]) else ""),
                      f1(a.get("antes_max")), f1(a.get("despues_max")), pct(r["P_F"], 0), f1(r["peor"]["dT_agua_max"], 1)])
    if filas:
        escs = []
        for n_r, r in enumerate(X["rondas"], start=1):
            cr = cargar_estudio(RAIZ / "config" / f"{r}.yaml")
            cur = cr["curvas"].get("usar")
            cl = cr["clima"]
            n_cur = len(cur) if isinstance(cur, list) else 8
            if cl.get("temporadas"):
                escs.append(f"R{n_r}: {n_cur} curvas × las temporadas {', '.join(map(str, cl['temporadas']))}")
            else:
                escs.append(f"R{n_r}: {n_cur} curvas × {cl.get('historias', 3)} terciles")
        W.tabla(f"Exploración del retiro en {len(X['rondas'])} rondas, cada una declarada antes de ver sus resultados "
                f"({'; '.join(escs)}; 12 escenarios por variante, con el mismo peso). M~Δ~ máximo antes y después del "
                "retiro y P~f~ de ese filtro, con la corrección de malla, y diferencia máxima hormigón–agua [°C].",
                ["Variante", "Retiro", "Agua V3–V5 [°C]", "M~Δ~ antes", "M~Δ~ después", "P~f~", "ΔT agua"], filas,
                [4.6, 1.6, 3.4, 1.7, 1.9, 1.3, 1.6], alinear=["l", "c", "c", "c", "c", "c", "c"], tam=7.5, ref="retiro")
        if "retiro" in figs:
            W.figura(figs["retiro"], "ΔT(t) en el peor escenario de cada variante de retiro, un panel por ronda (malla "
                                     "del estudio); las líneas verticales marcan cada retiro.", ref="retiro")
    W.p("Lo que controla el salto es la diferencia entre el núcleo y la mínima nocturna del aire al momento del "
        "retiro: bajo 10 cm de manta el núcleo se enfría muy lento, y mantener los tubos del núcleo funcionando hasta "
        "el retiro baja esa diferencia, siempre que el agua no enfríe de más el hormigón junto a las caras.",
        alinear="justificado")
    if X["ret"]:
        ok = [n for n in X["ret_cfg"] if n in X["ret"] and X["ret"][n]["P_F"] <= 0.05]
        # Variante con el mayor M_Δ antes del retiro entre las que enfrían sólo el núcleo (V3-V5).
        solo_nucleo = [n for n in X["ret_ad"] if all(t["vertida"] >= 3 for t in X["ret_cfg"][n].get("tubos", []))]
        a_max = max(solo_nucleo, key=lambda n: X["ret_ad"][n]["antes_max"], default=None)
        if a_max:  # el mecanismo descrito abajo vale si esa variante enfría el agua antes del día 14
            prog0 = next((t.get("programa") for t in X["ret_cfg"][a_max]["tubos"] if t.get("programa")), None)
            if not prog0 or prog0[0]["hasta_h"] >= 336:
                a_max = None
        txt = (f"De las {len(X['ret'])} variantes, "
               + (f"{'cumple' if len(ok) == 1 else 'cumplen'} con la corrección de malla: {', '.join(ok)}. " if ok else
                  "ninguna cumple con la corrección de malla. "))
        if a_max:
            txt += (f"Bajar el agua temprano tiene un costo: la variante {a_max} sube el M~Δ~ antes del retiro a "
                    f"{f1(X['ret_ad'][a_max]['antes_max'])} °C, porque el agua fría enfría la zona junto a la cara "
                    "inclinada de V4 mientras el núcleo sigue caliente. ")
        txt += ("La segunda ronda mantuvo por eso el agua a 25 °C hasta los días 10 a 14, pasado el pico, y recién "
                "entonces la bajó por tramos.")
        Ca = X.get("C_ant") or {}
        if Ca and X["prop_ant"] != prop:
            txt += (f" Su mejor variante, evaluada con las 8 curvas y las {len(Ca['anios'])} temporadas, cumplía con la "
                    f"malla del estudio, pero con la malla fina el salto del retiro subió hasta "
                    f"{f1(X['delta_ant_max'], 2)} °C más y P~f~ quedó en {pct(Ca['P_F'])}: no cumple.")
        W.p(txt, alinear="justificado")
        r3 = [n for n in X["ret_cfg"] if X["ronda_de"][n] == len(X["rondas"]) and n in X["ret"]]
        if Ca and len(r3) > 1:
            def _ultimo_agua(n):
                pr = next((t.get("programa") for t in X["ret_cfg"][n]["tubos"] if t.get("programa")), None) or [{}]
                return pr[-1].get("T_agua_C")

            def _par(clave):  # variantes iguales salvo en `clave`: efecto de esa palanca sobre el peor M_Δ
                otras = {"T": lambda n: X["ret_cfg"][n].get("T_colocacion"), "agua": _ultimo_agua,
                         "retiro": lambda n: retiro_de(X["ret_cfg"][n])}
                for a in r3:
                    for b in r3:
                        if a < b and otras[clave](a) != otras[clave](b) and all(
                                f(a) == f(b) for k, f in otras.items() if k != clave):
                            return (a, b) if otras[clave](a) > otras[clave](b) else (b, a)
                return None

            partes7 = []
            for clave, txt_ in (("T", "colocar más frío"), ("agua", "bajar el agua del último tramo"),
                                ("retiro", "retirar más tarde")):
                par = _par(clave)
                if par:
                    alto, bajo = par
                    m_a, m_b = X["ret"][alto]["MD_max"], X["ret"][bajo]["MD_max"]
                    if clave == "retiro":
                        alto, bajo, m_a, m_b = bajo, alto, m_b, m_a  # el "cambio" es retirar más tarde
                    partes7.append(f"{txt_} {'baja' if m_b < m_a else 'sube'} el peor M~Δ~ de {f1(m_a)} a {f1(m_b)} °C")
            mejor = min(r3, key=lambda n: X["ret"][n]["MD_max"])
            W.p("La tercera ronda se declaró después de ver ese resultado y se filtró en los escenarios que fallaban. "
                + ("En ella, " + "; ".join(partes7) + ". " if partes7 else "")
                + f"La de menor M~Δ~ es {mejor} ({X['ret_cfg'][mejor].get('descripcion', '').replace(' -> ', ' → ')}), "
                + ("que es la propuesta final" if mejor == prop else f"y la propuesta final es {prop}")
                + ", evaluada luego con las 8 curvas y las 30 temporadas.", alinear="justificado")
    W.h2("7.3 Confiabilidad de la propuesta final")
    if C:
        pe = C["peor"]
        filas = [["P(T~máx~ > 70 °C)", f"{int(C['k_T'])} de {C['n']}", pct(C["P_T"])],
                 ["P(ΔT~máx~ > 20 °C)", f"{int(C['k_D'])} de {C['n']}", pct(C["P_D"])],
                 ["P~f~ (evento conjunto)", f"{int(C['k_F'])} de {C['n']}", pct(C["P_F"])],
                 ["M~Δ~ antes del retiro (máx.)", "", f"{f1(ad['antes_max'])} °C"],
                 ["M~Δ~ después del retiro (máx.)", "", f"{f1(ad['despues_max'])} °C"],
                 ["M~T~ (P50 / P95 / máx.)", "", f"{f1(C['MT_p50'])} / {f1(C['MT_p95'])} / {f1(C['MT_max'])} °C"],
                 ["M~Δ~ (P50 / P95 / máx.)", "", f"{f1(C['MD_p50'])} / {f1(C['MD_p95'])} / {f1(C['MD_max'])} °C"]]
        W.tabla(f"Propuesta final con 8 curvas × {len(C['anios'])} temporadas = {C['n']} escenarios (pesos 1/240), "
                "M~Δ~ corregido por malla. Los percentiles son de la distribución ponderada de escenarios.",
                ["Indicador", "Escenarios que exceden", "Valor"], filas, [6.0, 4.4, 6.0], alinear=["l", "c", "c"],
                ref="conf")
        W.p(f"El peor escenario ({pe['curva']}, temporada {pe['clima'].split('_')[-1]}) alcanza M~Δ~ = {f1(pe['M_D'])} °C "
            f"a las {f1(pe['t_MD_h'], 0)} h, con el punto frío en ({f1(pe['x_Tmin'], 2)}; {f1(pe['y_Tmin'], 2)}; "
            f"{f1(pe['z_Tmin'], 2)}) m, en {_zona(pe['z_Tmin'])} (Figuras {ref('f:temp')} y {ref('f:conf')})"
            + (f". La coordenada es la del modelo discretizado: con Δx = 0,275 m la cara superior del pedestal queda "
               f"en z = {f1(pe['z_Tmin'] + 0.025, 2)} m por la escalera de vóxeles, y no en {f1(H_geo, 2)} m; sobre la "
               f"geometría real, el punto está 25 mm bajo esa cara, en z = {f1(H_geo - 0.025, 3)} m"
               if pe["z_Tmin"] > H_geo else "")
            + ". Entre "
            "temporadas, el pico de los primeros días cambia poco para una misma curva: con la manta, la superficie "
            "está desacoplada del aire y el diferencial lo fijan la curva y los tubos.", alinear="justificado")
        md22 = max((f["M_D"] for f in C["filas"] if f["curva"].endswith("_22A")), default=float("nan"))
        i_pe = C["curvas"].index(pe["curva"])
        fila_pe = C["MD"][i_pe]
        orden_t = np.argsort(fila_pe)[::-1]
        peores_t = ", ".join(C["anios"][j].split("_")[-1] for j in orden_t[:3])
        W.p((f"Después del retiro, en cambio, la temporada sí importa. Con la manta puesta el mayor diferencial es "
             f"{f1(ad['antes_max'])} °C y después del retiro llega a {f1(ad['despues_max'])} °C: el salto lo decide la "
             f"noche más fría de los días que siguen al retiro. " if ad["despues_max"] > ad["antes_max"] else
             f"Con la manta puesta el mayor diferencial es {f1(ad['antes_max'])} °C y después del retiro, "
             f"{f1(ad['despues_max'])} °C: retirando el día {f1(dia_ret, 0)}, con el núcleo ya frío, el retiro deja de "
             "gobernar y manda el pico de los primeros días. ")
            + f"Con la curva {pe['curva']}, M~Δ~ va de "
            f"{f1(np.nanmin(fila_pe))} a {f1(np.nanmax(fila_pe))} °C según la temporada, y las más exigentes son "
            f"{peores_t}. El plan no fija una fecha rígida: el retiro se condiciona a lo que miden los sensores "
            "y al pronóstico de las mínimas (sección 8.3). "
            + (f"Las curvas 22A, que liberan menos calor y lo hacen antes, llegan a lo más a {f1(md22)} °C; el riesgo se "
               f"concentra en las SS. " if math.isfinite(md22) else "")
            + (f"{pe['curva']} es además la curva de mayor energía a 14 días ({f1(X['fuentes_E'][pe['curva']], 0)} "
               "MJ/m³)." if max(X["fuentes_E"], key=X["fuentes_E"].get) == pe["curva"] else ""),
            alinear="justificado")
        if "temp" in figs:
            W.figura(figs["temp"], "Izquierda: núcleo, superficie a 25 mm y aire en el peor escenario. Derecha: ΔT(t) con "
                                   "la misma curva en las 30 temporadas (malla del estudio; la línea punteada vertical es "
                                   "el retiro).", ref="temp")
        W.figura(figs["conf"], f"M~Δ~ (con la corrección de malla) y M~T~ de los {C['n']} escenarios de la propuesta: cada "
                               "punto es una temporada.", ref="conf")
        if X["corte"] is not None:
            W.figura(X["corte"], "Cortes por el eje del peor escenario en la primera semana después de la colocación "
                                 "(temperatura en °C).",
                     ancho_cm=14.5, ref="corte")
    W.h2("7.4 Sensibilidades")
    filas = []
    ss = X["solo_ss"]
    if ss:
        filas.append(["Sólo campaña SS (pesos renormalizados)", pct(ss["P_F"]), f1(ss["MD_max"])])
    if prop in fam:
        filas.append(["Propuesta con 3 terciles (referencia de las filas siguientes)", pct(fam[prop]["P_F"]),
                      f1(fam[prop]["MD_max"])])
    T_base = float(X["prop_cfg"].get("T_colocacion", 10))
    var_T = sorted(((int(m.group(1)), n) for n in X["var"] for m in [re.fullmatch(re.escape(prop) + r"-T(\d+)", n)]
                    if m), key=lambda x: x[0])
    for nombre, etiq in ((f"{prop}-rhocp", "ρc~p~ de las probetas (2,53 MJ/m³K)"),
                         *((n, f"Colocación a {T} °C") for T, n in var_T),
                         (f"{prop}-mantaR20", "Manta con R = 2,0 m²K/W (húmeda o de menor rendimiento)"),
                         (f"{prop}-mantaR30", "Manta con R = 3,0 m²K/W")):
        if nombre in X["var"]:
            r = X["var"][nombre]
            filas.append([etiq + " (3 terciles)", pct(r["P_F"]), f1(r["MD_max"])])
    if est:
        filas.append(["Curvas con +20 % de calor (3 terciles)", pct(est["P_F"]), f1(est["MD_max"])])
    if X["ver"]:
        filas.append(["Clima del ensemble ECMWF (3 terciles de 51 miembros)", pct(X["ver"]["P_F"]), f1(X["ver"]["MD_max"])])
    if filas:
        W.tabla("Sensibilidades de la propuesta final (8 curvas salvo indicación; M~Δ~ corregido por malla).",
                ["Caso", "P~f~", "M~Δ~ máx. [°C]"], filas, [10.0, 3.0, 3.6], alinear=["l", "c", "c"], ref="sensib")
    partes = []
    if ss:
        partes.append(f"Si se descarta la campaña 22A por su menor confianza y se renormalizan los pesos a las 4 curvas "
                      f"SS (120 escenarios), P~f~ pasa a {pct(ss['P_F'])} con M~Δ~ máximo de {f1(ss['MD_max'])} °C: "
                      + (f"la propuesta cumple también con ese criterio más exigente (con las 8 curvas, "
                         f"P~f~ = {pct(pf)})." if ss["P_F"] <= 0.05 else
                         "con ese criterio más exigente la propuesta no cumple, y el resultado depende de aceptar la "
                         "campaña 22A."))
    rc = X["var"].get(f"{prop}-rhocp")
    if rc:
        partes.append(f"Con el ρc~p~ de las probetas, el mismo calor sube menos la temperatura y M~Δ~ máximo queda en "
                      f"{f1(rc['MD_max'])} °C (P~f~ = {pct(rc['P_F'])}), lo que confirma que el valor del modelo es el "
                      "lado conservador.")
    if var_T:
        vs = [(T, X["var"][n]) for T, n in var_T]
        falla = next((T for T, r in vs if T > T_base and r["P_F"] > 0.05), None)
        txts = [f"a {T} °C da M~Δ~ máximo de {f1(r['MD_max'])} °C (P~f~ = {pct(r['P_F'])})" for T, r in vs]
        partes.append("Colocar " + (", ".join(txts[:-1]) + " y " + txts[-1] if len(txts) > 1 else txts[0]) + ": "
                      + (f"con hormigón a {falla} °C el plan ya no cumple, de modo que el límite de {f1(T_base, 0)} °C "
                         "a la descarga es obligatorio y se controla camión a camión." if falla is not None else
                         "un camión algo más tibio que lo especificado no cambia la conclusión, pero sí reduce el margen; "
                         f"por eso el límite de {f1(T_base, 0)} °C a la descarga se controla camión a camión."))
    r20, r30 = X["var"].get(f"{prop}-mantaR20"), X["var"].get(f"{prop}-mantaR30")
    if r20 and r30:
        b = fam.get(prop)
        base = (f" (con R = 2,5, en los mismos 24 escenarios, P~f~ = {pct(b['P_F'])} y {f1(b['MD_max'])} °C)"
                if b else "")
        poco = b is not None and abs(r20["MD_max"] - b["MD_max"]) < 0.5 and abs(r30["MD_max"] - b["MD_max"]) < 0.5
        partes.append(f"Con una manta de R = 2,0 m²K/W, húmeda o de menor rendimiento, P~f~ = {pct(r20['P_F'])} y "
                      f"M~Δ~ máximo de {f1(r20['MD_max'])} °C; con R = 3,0, P~f~ = {pct(r30['P_F'])} y "
                      f"{f1(r30['MD_max'])} °C{base}: "
                      + ("la resistencia de la manta casi no cambia el resultado"
                         + (", porque el peor diferencial ocurre después de retirarla."
                            if b["peor"]["t_MD_h"] > t_ret else ".") if poco else
                         "perder un 20 % de la resistencia no cambia la conclusión, pero reduce el margen; la R se "
                         "controla en obra." if r20["P_F"] <= 0.05 else
                         "una manta un 20 % menos aislante ya no cumple, de modo que R ≥ 2,5 m²K/W es un requisito "
                         "del plan: la manta se mantiene seca bajo el film y se revisa en cada turno."))
    if est:
        partes.append(f"Con las curvas escaladas en +20 %, P~f~ = {pct(est['P_F'])} y M~Δ~ máximo de "
                      f"{f1(est['MD_max'])} °C: " + ("el plan resiste ese exceso de calor." if est["P_F"] <= 0.05 else
                                                     "el plan no resiste ese exceso de calor, por lo que el monitoreo "
                                                     "temprano es la defensa frente a un hormigón más exotérmico."))
    if X["ver"]:
        partes.append(f"Con el pronóstico ensemble más reciente (ECMWF, 51 miembros, continuado con la climatología "
                      f"ERA5), P~f~ = {pct(X['ver']['P_F'])} y M~Δ~ máximo de {f1(X['ver']['MD_max'])} °C, "
                      + ("dentro del rango de las 30 temporadas históricas." if C and X["ver"]["MD_max"] <= C["MD_max"]
                         else "sobre el máximo de las 30 temporadas históricas.")
                      + " El pronóstico cubre sólo las primeras dos semanas: el retiro queda en la continuación "
                      "climatológica, más suave que una temporada real, por lo que esta verificación vale para el pico "
                      "de los primeros días y no para el retiro.")
    if partes:
        W.p(" ".join(partes), alinear="justificado")
    if mat:
        W.p(f"Ajustando M~Δ~ contra la energía liberada a 14 días por las curvas SS, con y sin el +20 % (recta con "
            f"R² = {f1(mat['r2'], 3)}), el diferencial llegaría a 20 °C con **E~crít~ ≈ {f1(mat['E_crit'], 0)} MJ/m³**, "
            + (f"un {pct(mat['exceso'], 0)} más que la probeta más exotérmica ({mat['curva_max']}, "
               f"{f1(mat['E_max'], 0)} MJ/m³). Es un margen determinista, no una probabilidad: indica cuánto más "
               "exotérmico que las probetas puede ser el hormigón de la obra antes de que el plan deje de cumplir."
               if mat["exceso"] >= 0 else
               f"un {pct(-mat['exceso'], 0)} menos que la probeta más exotérmica ({mat['curva_max']}, "
               f"{f1(mat['E_max'], 0)} MJ/m³): con la corrección de malla, el plan no tiene margen ni siquiera frente "
               "a las probetas medidas."), alinear="justificado")

    # ---------------- 8. Discusión y plan de control
    W.h1("8. Discusión y plan de control")
    # Con manta = la columna "Manta" de la tabla de la familia no es "—" (5 o 10 cm; la propuesta lleva 10 cm).
    con_manta = [p for p in fam if ETIQUETAS.get(p, ("", "", "10 cm"))[2] != "—"]
    espesores = sorted({ETIQUETAS.get(p, ("", "", "10 cm"))[2] for p in con_manta})
    n25 = sum(1 for p in con_manta if fam[p]["MD_max"] <= 25.0 and fam[p]["MT_max"] <= LIM_T)
    W.p("**Mecanismo.** Ninguna medida basta sola. La manta desacopla la superficie del aire; el preenfriado baja el "
        "pico del núcleo; los tubos con agua tibia sólo en el núcleo (V3–V5) recortan el pico sin enfriar el anillo "
        "exterior de la losa, que es la zona más fría. El punto frío del período protegido aparece donde el serpentín de "
        "V4 se acerca a la cara inclinada, y el del retiro, en las caras que quedan expuestas a la noche. "
        + ("La fase crítica es el retiro: por eso los tubos siguen funcionando hasta que la manta se retira."
           if ad and ad["despues_max"] > ad["antes_max"] else
           f"El retiro, que gobernaba los planes anteriores, deja de ser crítico cuando se hace el día "
           f"{f1(dia_ret, 0)} con los tubos funcionando hasta entonces."), alinear="justificado")
    W.p("**Comparación.** Ju y Lei (2019) midieron, sin control térmico, un diferencial de ≈35 °C en una losa de "
        f"espesor parecido; el plan propuesto lo mantiene bajo 20 °C en "
        f"{'todos los' if C and C['k_F'] == 0 else 'casi todos los' if cumple else f'el {pct(1 - pf, 0)} de los'} "
        "escenarios, con controles que ACI 207.4R describe como habituales: preenfriado, aislación y tubos. El límite de "
        "20 °C es el extremo exigente del rango de 20–30 °C de la especificación del ICH"
        + (f"; con un límite de 25 °C cumplirían {n25} de los {len(con_manta)} planes de la Tabla {ref('t:familia')} "
           f"que llevan manta ({' o '.join(espesores)}), en todos sus escenarios. " if con_manta else ". ")
        + ("En la familia, la manta sola (con o sin preenfriado) no cumple el límite de 20 °C, y " if
           all(p in fam and fam[p]["P_F"] > 0.05 for p in ("unica-M10", "unica-M10-T10")) else "En la familia, ")
        + "el retiro de las protecciones aparece como un momento de riesgo que debe planificarse, igual que en la "
        "práctica de obra.", alinear="justificado")
    W.p("**Factibilidad.** Todo el plan usa equipos corrientes: hielo en escamas en la planta, mantas de obra, tubos "
        "HDPE como los de riego, un estanque y un enfriador del tamaño de los que se arriendan para obras. El costo "
        "principal es la operación del circuito de agua durante "
        f"{f1(max(t['hasta_h'] for t in tubo3['programa']) / 24, 0)} días, con turnos de revisión y energía de "
        "respaldo, y el de mantener la manta en buen estado, sellada en los bordes, todo ese tiempo. A cambio, la "
        "colocación se hace en un solo día, sin juntas frías, que es la opción más rápida posible con una bomba.",
        alinear="justificado")
    W.p("**Limitaciones.** La fuente depende sólo de la edad (sin madurez); el ρc~p~ del modelo es menor que el de las "
        "probetas, lo que es conservador; las propiedades del hormigón y del suelo son bibliográficas y no hay "
        "validación con mediciones de la fundación; la malla fina se corrió sólo en "
        f"{len([x for x in sens if x['plan'] == prop])} escenarios y en los demás se usa la corrección δ; el ΔT es un "
        "indicador de riesgo, no un análisis de tensiones; y con 8 curvas el "
        "modelo de material es tan bueno como esas probetas: un hormigón más exotérmico que E~crít~ no está cubierto, por "
        "lo que el monitoreo de las primeras 48 h es parte del plan.", alinear="justificado")
    W.h2("8.1 Procedimiento de vertida y protecciones")
    filas = [[f"V{v}", f"{f1(geo[v - 1]['z0'], 2)}–{f1(geo[v - 1]['z1'], 2)}", f1(geo[v - 1]["V"]),
              f"{f1(inicios[v - 1])}–{f1(fines[v - 1])}", f"{hora(t0, inicios[v - 1])} → {hora(t0, fines[v - 1])}"]
             for v in range(1, 6)]
    filas.append(["Total", "0–4,30", f1(sum(g["V"] for g in geo)), f"0–{f1(D)}", f"{hora(t0, 0)} → {hora(t0, D)}"])
    W.tabla("Colocación continua por capas a 27 m³/h (una bomba y una cuadrilla, sin juntas frías).",
            ["Vertida", "Cota z [m]", "Volumen [m³]", "Horas desde t~0~", f"Reloj ({MESES[t0.month - 1][:3]} {t0.year})"],
            filas, [2.0, 2.6, 2.6, 3.2, 5.6], destacar=(5,), ref="secuencia")
    W.figura(figs["crono"], "Cronograma: colocación (izquierda, en horas) y protecciones, programa del agua y monitoreo "
                            "(derecha, en días).", ref="crono")
    tub_txt = ""
    if X["tubos"]:
        Pmax = sum(s_["potencia_max_kW"] for s_ in X["tubos"].values())
        Tsal = max((s_["T_salida_max_C"] or 0) for s_ in X["tubos"].values())
        tub_txt = (f" En el peor escenario la potencia extraída llega a {f1(Pmax)} kW y el agua vuelve a hasta "
                   f"{f1(Tsal)} °C: ese es el tamaño del enfriador del circuito.")
    W.vinetas([
        "**Colocación:** continua, de día y de noche, con dos turnos, iluminación y una bomba de respaldo; capas de "
        "≤ 0,5 m vibradas. Una interrupción larga equivale a las variantes con esperas, que no cumplen.",
        f"**Temperatura de colocación:** máximo {Tb} °C a la descarga, medida en cada camión; la planta apunta a "
        f"{f1(T_obj, 0)} °C con hielo en escamas (≈{f1(hielo[0], 0)}–{f1(hielo[1], 0)} kg/m³ según un balance de "
        f"energía con los materiales a 16–18 °C, ≈{f1(hielo[0] * V_tot / 1000, 0)}–{f1(hielo[1] * V_tot / 1000, 0)} t "
        "en total), agua fría y áridos a la sombra.",
        "**Moldaje y manta:** terciado de 18 mm en las caras verticales de V1, V2 y V5; manta de 10 cm con film de "
        "polietileno (curado y protección contra la lluvia) en todas las caras, a medida que sube el hormigón; se "
        f"retiran el **día {f1(dia_ret, 0)}**, de día, sólo si se cumple el criterio de la sección 8.3.",
        f"**Tubos:** un serpentín HDPE Ø32 mm a media altura de V3, V4 y V5, filas cada 1 m, 20 L/min por circuito, con "
        f"agua a {prog}; circuito cerrado con estanque, bomba, calefactor y enfriador, y respaldo de energía.{tub_txt} "
        "Al terminar se inyectan con lechada.",
    ])
    W.p(f"El plan queda también en el modelo BIM: geometry/{X['ifc'].name} es el IFC del curso con el property set de "
        "cada vertida actualizado (Anexo B).", alinear="justificado")
    W.h2("8.2 Monitoreo")
    if X["sensores"]:
        W.p(f"Se instalan {n_sens} sensores de temperatura (termocuplas tipo T o digitales, ±0,5 °C) donde el modelo "
            f"predice el punto más caliente y los puntos fríos que definen el diferencial (Figura {ref('f:sens')}, Tabla "
            f"{ref('t:sensores')}), con redundancia en los de más peso, más el aire, el suelo y el agua de los tubos. "
            "Registro cada 15 min hasta las 72 h, cada 30 min hasta el día 14 y cada hora hasta 3 días después del "
            "retiro, con transmisión inalámbrica como en el sistema de Wang et al. (2025), alarmas y comparación con "
            "el modelo dos veces al día.", alinear="justificado")
        W.p("**Responsables.** El jefe de terreno responde por la ejecución del plan y autoriza el retiro de la manta. "
            "Un ingeniero de control térmico, con un reemplazo designado para cada turno, revisa los registros dos veces "
            "al día, los compara con la curva prevista, activa las contingencias de la Tabla "
            f"{ref('t:contingencias')} y firma el criterio de retiro. El laboratorista mide la temperatura de cada camión "
            f"a la descarga y lo rechaza si supera {Tb} °C. El operador del circuito de agua registra caudales y "
            "temperaturas en cada cambio de tramo, y el instalador de los sensores verifica su lectura antes de la "
            "primera colocación.", alinear="justificado")
        if "sens" in figs:
            W.figura(figs["sens"], "Sensores sobre la sección y la planta, junto a los puntos más calientes (naranja) y "
                                   "más fríos a 25 mm (azul) de los escenarios de la propuesta.", ref="sens")
        filas = [[s_["id"], s_["tipo"], str(s_["n"]), f"({f1(s_['x'], 2)}; {f1(s_['y'], 2)}; {f1(s_['z'], 2)})",
                  s_["donde"]] for s_ in X["sensores"]]
        filas.append(["W1–W6", "agua", "6", "entrada y salida de V3–V5", "temperatura del agua y caudal de cada circuito"])
        W.tabla("Sensores: identificador, tipo, cantidad, posición (x; y; z) en metros y propósito.",
                ["ID", "Tipo", "N.º", "Posición [m]", "Propósito"], filas, [1.3, 1.5, 0.9, 3.8, 9.7],
                alinear=["c", "l", "c", "c", "l"], tam=7.5, ref="sensores")
    W.h2("8.3 Criterio de retiro y contingencias")
    cr = X["crit"]
    if cr and math.isfinite(cr.get("dT_max", float("nan"))):
        W.p(f"La manta y el moldaje se retiran el día {f1(dia_ret, 0)} sólo si los sensores muestran lo que el modelo "
            f"tenía en esa fecha en los escenarios que cumplen: núcleo ≤ {math.floor(cr['nucleo_max'])} °C, "
            f"ΔT medido ≤ {math.floor(cr['dT_max'] + d)} °C y núcleo − mínima del aire prevista para los 3 días "
            f"siguientes ≤ {math.floor(cr['nucleo_aire_max'])} °C. Si no se cumplen, el agua sigue en su último tramo y "
            "se evalúa de nuevo cada 2 días.", alinear="justificado")
    mas_tibios = [(T, X["var"][n]) for T, n in var_T if T > T_base]
    if mas_tibios:
        T_f, r_f = next(((T, r) for T, r in mas_tibios if r["P_F"] > 0.05), mas_tibios[-1])
        txt12 = (f"a {T_f} °C, P~f~ = {pct(r_f['P_F'])}" if r_f["P_F"] > 0.05 else
                 f"colocar a {T_f} °C reduce el margen")
    else:
        txt12 = "la temperatura de colocación mueve el pico"
    filas = [
        [f"Hormigón a la descarga > {Tb} °C", "Se rechaza el camión o se ajusta el hielo en planta.",
         f"{txt12} (Tabla {ref('t:sensib')})"],
        ["Falla de la bomba de hormigón o de la planta", "Bomba de respaldo; la interrupción no debe pasar de 2 h; si pasa, "
         "se trata la junta y se reevalúa el plan.", "las esperas de 12–24 h no cumplen"],
        ["Falla del circuito de agua", "Bomba y generador de respaldo; reparar en < 12 h; no retirar la manta mientras "
         "tanto.", "sin tubos después del día 14 el retiro no cumple"],
        ["ΔT medido sobre la curva prevista o ≥ 18 °C", "Revisar mantas y sellos; segunda manta en la cara del sensor más "
         "frío; recalibrar el modelo con lo medido.", f"peor escenario previsto: {f1(C['MD_max']) if C else '—'} °C"],
        ["ΔT medido ≥ 19 °C", "Segunda manta en todo el sector; subir el agua del circuito vecino al punto frío.",
         "el punto frío está junto al serpentín de V4"],
        ["Núcleo ≥ 60 °C o agua–hormigón > 20 °C", "Revisar caudal y temperatura del agua.", f"M~T~ previsto ≤ "
         f"{f1(C['MT_max']) if C else '—'} °C"],
        ["Helada o frente frío previsto en la fecha de retiro", "Postergar el retiro y mantener el agua.",
         "el salto depende del aire nocturno"],
        ["Sensor dañado", "Usar el redundante; reemplazar si es accesible.", "redundancia"],
    ]
    W.tabla("Umbrales y acciones de contingencia.", ["Disparador", "Acción", "Respaldo en el modelo"], filas,
            [4.6, 8.4, 4.2], alinear=["l", "l", "l"], tam=7.5, ref="contingencias")

    # ---------------- 9. Conclusiones
    W.h1("9. Conclusiones")
    n_finas = len([x for x in sens if x["plan"] == prop])
    dcrit = X.get("_delta_crit")
    W.p("**Lo demostrado**, dentro del modelo declarado (8 curvas medidas, 30 temporadas observadas y los parámetros "
        f"de la Tabla {ref('t:param')}):", despues=2)
    W.vinetas([
        f"Con las 8 curvas medidas, la temperatura máxima nunca limita (M~T~ ≤ {f1(X['MT_familia'])} °C en todos los "
        "planes); el diseño lo gobierna el diferencial de 20 °C.",
        frase_retiro,
        f"La propuesta final —vertida única de {f1(D)} h, hormigón a ≤ {Tb} °C, manta de 10 cm y tubos en V3–V5 con agua "
        f"a {prog}, retiro el día {f1(dia_ret, 0)}—, evaluada con malla fina en los {n_finas} escenarios que "
        f"gobiernan o están más cerca del límite y la corrección δ = {f1(d, 2)} °C en el resto, {frase_prop}.",
        ("**Condición de validez:** la conclusión se mantiene mientras el error de malla de los escenarios no refinados "
         f"sea menor que δ~crít~ = {f1(dcrit, 2)} °C, {f1(dcrit / d, 1)} veces la mayor diferencia medida."
         if cumple and dcrit and d > 0 else
         "El error de malla no cambia la clasificación en el rango explorado (Tabla " + ref("t:delta") + ")."
         if cumple else
         f"**El error de malla cambia la clasificación:** con la malla del estudio la propuesta cumple, pero con la "
         f"malla fina y la mayor corrección medida (δ = {f1(d, 2)} °C) no; cumpliría sólo si ese error fuese menor "
         f"que δ~crít~ = {f1(dcrit, 2)} °C (Tabla {ref('t:delta')})." if dcrit is not None else
         f"**El error de malla cambia la clasificación** (Tabla {ref('t:delta')})."),
        "P~f~ es una enumeración exacta sobre las curvas medidas y las temporadas observadas, sin distribuciones "
        "supuestas" + (f"; con sólo la campaña SS resulta {pct(X['solo_ss']['P_F'])}" if X["solo_ss"] else "") + ".",
    ])
    r20, r30 = X["var"].get(f"{prop}-mantaR20"), X["var"].get(f"{prop}-mantaR30")
    W.p("**Lo pendiente o fuera del modelo:**", despues=2)
    W.vinetas([
        "Validación: no hay mediciones de esta fundación; el monitoreo (sección 8.2) compara lo medido con lo previsto "
        "desde la primera colocación.",
        (f"Material: con +20 % de calor el plan no cumple (P~f~ = {pct(est['P_F'])}) y el margen frente a un hormigón "
         f"más exotérmico es acotado (E~crít~ ≈ {f1(mat['E_crit'], 0)} MJ/m³, {pct(abs(mat['exceso']), 0)} "
         f"{'sobre' if mat['exceso'] >= 0 else 'bajo'} la probeta más exotérmica): los sensores de las primeras 48 h y "
         "las contingencias son parte del plan."
         if est and mat and est["P_F"] > 0.05 else
         "Material: un hormigón más exotérmico que las probetas se controla con los sensores de las primeras 48 h."),
        "Datos por confirmar: el origen temporal de las curvas (desfase de "
        f"{f1(X['curvas_ctl'].get('desfase_h', 2.0), 0)} h) y que la dosificación de la obra sea la de las probetas.",
        *([("Manta: entre R = 2,0 y 3,0 m²K/W el resultado casi no cambia (M~Δ~ máximo "
            f"{f1(r20['MD_max'])} °C con R = 2,0" + (f" y {f1(r30['MD_max'])} °C con R = 3,0" if r30 else "")
            + f", frente a {f1(fam[prop]['MD_max'])} °C con R = 2,5); igual se especifica R ≥ 2,5 m²K/W.")
           if prop in fam and abs(r20["MD_max"] - fam[prop]["MD_max"]) < 0.5 else
           ("Manta: con R = 2,0 m²K/W la propuesta no cumple; R ≥ 2,5 m²K/W es un requisito que se controla en obra.")
           if r20["P_F"] > 0.05 else
           "Manta: con R = 2,0 m²K/W la propuesta sigue cumpliendo, con menos margen."] if r20 else []),
        f"Malla: no se corrió una resolución más fina que {f1(min(x['dx_fina'] for x in sens), 2)} m; la Tabla "
        f"{ref('t:delta')} acota su efecto sobre la decisión." if sens else "Malla: no se corrió una malla fina.",
    ])

    # ---------------- Referencias y anexos
    W.h1("Referencias")
    W.vinetas([
        "ACI Committee 207 (2005). *Cooling and Insulating Systems for Mass Concrete* (ACI 207.4R-05). American Concrete "
        "Institute.",
        "ACI Committee 301 (2020). *Specifications for Concrete Construction* (ACI 301-20). American Concrete Institute.",
        "Dittus, F. W. y Boelter, L. M. K. (1930). Heat transfer in automobile radiators of the tubular type. "
        "*University of California Publications in Engineering*, 2(13), 443–461.",
        "ECMWF (2026). IFS Ensemble (ENS), 51 miembros, 0,25°. European Centre for Medium-Range Weather Forecasts; "
        "obtenido mediante Open-Meteo Ensemble API.",
        "Equipo docente del curso (2026). Curvas de calor de hidratación Q(t) de las probetas M1–M4, campañas SS y 22A "
        "(source_terms.txt, versión corregida del 14 de noviembre de 2026). " + CURSO + ".",
        "Harris, C. R. et al. (2020). Array programming with NumPy. *Nature*, 585, 357–362.",
        "Hersbach, H. et al. (2020). The ERA5 global reanalysis. *Quarterly Journal of the Royal Meteorological "
        "Society*, 146, 1999–2049.",
        f"IfcOpenShell ({IFC_VERSION}). IfcOpenShell: the open source IFC toolkit and geometry engine. "
        "https://ifcopenshell.org.",
        "Instituto del Cemento y del Hormigón de Chile (2019). *Especificación para Hormigón Masivo Estructural*. Santiago.",
        "Ju, Y. y Lei, H. (2019). Actual temperature evolution of thick raft concrete foundations and cracking risk "
        "analysis. *Advances in Materials Science and Engineering*, 7029671. doi:10.1155/2019/7029671.",
        "Jürges, W. (1924). Der Wärmeübergang an einer ebenen Wand. *Beihefte zum Gesundheits-Ingenieur*, Reihe 1, "
        "Heft 19.",
        "Open-Meteo (2026). Historical Weather API y Ensemble API. https://open-meteo.com.",
        "Peaceman, D. W. (1978). Interpretation of well-block pressures in numerical reservoir simulation. *Society of "
        "Petroleum Engineers Journal*, 18(3), 183–194.",
        "Smolana, A., Klemczak, B., Azenha, M. y Schlicke, D. (2021). Experiences and analysis of the construction "
        "process of mass foundation slabs aimed at reducing the risk of early age cracks. *Journal of Building "
        "Engineering*, 44, 102947. doi:10.1016/j.jobe.2021.102947.",
        "Swinbank, W. C. (1963). Long-wave radiation from clear skies. *Quarterly Journal of the Royal Meteorological "
        "Society*, 89(381), 339–348.",
        "Wang, T., Li, D., Zhou, J. y Zhang, J. (2025). Temperature monitoring of mass concrete structure using "
        "wireless sensing system. *Structural Control and Health Monitoring*, 2025(1), 7847074. "
        "doi:10.1155/stc/7847074.",
    ])
    def _g(v):
        return f"{v:g}"

    finas_txt = []
    for x in sorted((x for x in sens if x["plan"] == prop), key=lambda x: (not x["otras"], -x["MD_fina"])):
        mallas = [0.275] + sorted((m["dx"] for m in x["intermedias"]), reverse=True) + [x["dx_fina"]]
        opciones = "" if mallas == [0.275, 0.11] and x["otras"] else " --dx-variantes " + " ".join(map(_g, mallas))
        if not x["otras"]:
            opciones += " --dt-variantes 600 --profundidades 0.025"
        finas_txt.append(f"{x['curva']} {x['clima']}{opciones}")
    if not finas_txt and C:
        finas_txt = [f"{C['peor']['curva']} {C['peor']['clima']}"]
    W.h1("Anexo A. Reproducción")
    W.p("Todo se regenera desde el repositorio (Python 3.12: pip install -r requirements.txt; pytest):", despues=2)
    for linea in ("mcplan clima era5 --fecha-inicio 2026-11-19T08:00 --anios 1996 2025 --dias 40 --margen 7",
                  "mcplan clima ensemble",
                  "mcplan estudio config/estudio_era5.yaml            # familia de planes",
                  "mcplan estudio config/retiro_tubos_era5.yaml       # retiro: primera ronda",
                  "mcplan estudio config/retiro_tubos2_era5.yaml      # retiro: segunda ronda",
                  f"mcplan estudio config/{CONF_ANT}.yaml      # propuesta anterior: 8 curvas x 30 temporadas",
                  "mcplan estudio config/retiro_tubos3_era5.yaml      # retiro: tercera ronda (temporadas críticas)",
                  f"mcplan estudio config/{CONF_PROP}.yaml --planes {prop}   # propuesta final: 8 x 30",
                  "mcplan estudio config/estres_era5.yaml             # +20 % de calor",
                  "mcplan estudio config/sensibilidad_era5.yaml       # rho cp, colocación, R de la manta",
                  "mcplan estudio config/verificacion_ensemble.yaml   # clima del ensemble",
                  "# malla fina de la propuesta (una línea por escenario; la de la anterior, en el LEEME):",
                  *(f"mcplan estudio config/{CONF_PROP}.yaml --procesos 1 --sensibilidad {prop} {e}"
                    for e in finas_txt),
                  "python scripts/articulo_f01.py --pdf               # este artículo"):
        par = W.p(linea, tam=7.5, despues=0)
        for r in par.runs:
            r.font.name = "Consolas"
    W.p(f"Clima: snapshots inmutables con SHA-256 en data/clima/snapshots (ERA5 descargado el "
        f"{str(X['meta_era5'].get('descargado_utc', '—'))[:10]}). Curvas: data/hidratacion/source_terms.txt, sin "
        "modificar. Cada escenario es una fila de escenarios.csv con su peso, métricas, ubicación de los extremos y "
        "residuo de energía; la corrida de detalle deja los campos 3D en VTK para ParaView.", tam=9.5)
    W.h1("Anexo B. Property set de las vertidas en el IFC del plan")
    filas = [[f"V{v}", f1(X["ifc_valores"][f"Vertida_0{v}"]["Inicio_h"], 2), "10",
              X["ifc_valores"][f"Vertida_0{v}"]["Moldaje"], f1(X["ifc_valores"][f"Vertida_0{v}"]["RetiroMoldaje_h"], 1)]
             for v in range(1, 6)]
    W.tabla("UANDES_HormigonMasivo en geometry/" + X["ifc"].name + " (Aislacion: manta de 10 cm con film en todas las "
            "caras; RetiroMoldaje_h desde el inicio de cada vertida).",
            ["Vertida", "Inicio_h", "T col. [°C]", "Moldaje", "RetiroMoldaje_h"], filas, [2.0, 2.4, 2.4, 6.0, 3.4],
            ref="ifc")
    W.h1("Anexo C. Declaración de uso de herramientas de IA")
    W.p(f"Curso: {CURSO}. Actividad: F01, artículo del Proyecto 2. Estudiantes: {AUTORES.replace(' · ', ', ')}. "
        "Sí se utilizaron herramientas de IA generativa.", alinear="justificado")
    W.tabla("Herramientas de IA utilizadas (una fila por herramienta y versión).",
            ["Herramienta y versión", "Quién la usó", "Propósito", "Cómo se integró"],
            [["Claude Opus 5.5 (Anthropic), mediante Claude Code", "Eduardo Vergara",
              "programación del solver y de los scripts, análisis de resultados y redacción del artículo",
              "el código queda en el repositorio, verificado con los tests de la sección 5; los números del texto se "
              "calculan desde los resultados de las simulaciones (Anexo A); el grupo revisó el texto"]],
            [4.0, 2.6, 4.6, 6.0], alinear=["l"] * 4, tam=8, ref="ia")
    W.p("Los integrantes declaran que revisaron críticamente el contenido, que lo entienden y pueden defenderlo, y que "
        "asumen la responsabilidad por su calidad, veracidad y originalidad.", alinear="justificado")
    W.p("Firmas:  Trinidad Díaz ____________    Pablo Arancibia ____________    Eduardo Vergara ____________",
        despues=2)


def ifc_del_plan(X: dict) -> Path:
    """Copia del IFC oficial con el plan final en el pset UANDES_HormigonMasivo de cada vertida."""
    from mcplan.ifc_io import agregar_propiedades_curso

    valores = {}
    for v in range(1, 6):
        ini = X["inicios"][v - 1]
        valores[f"Vertida_0{v}"] = {
            "VertidaID": v, "Inicio_h": round(ini, 2), "TemperaturaColocacion_C": float(X["prop_cfg"].get("T_colocacion", 10)), "TipoHormigon": "M1",
            "Moldaje": "Terciado fenólico 18 mm" if v in (1, 2, 5) else "Sin moldaje (cara inclinada)",
            "Aislacion": "Manta 10 cm (R = 2,5 m2K/W) con film, todas las caras expuestas",
            "RetiroMoldaje_h": round(X["t_ret"] - ini, 1)}
    salida = RAIZ / "geometry" / "foundation_plan_final.ifc"
    agregar_propiedades_curso(RAIZ / "geometry" / "foundation.ifc", salida, valores)
    X["ifc_valores"] = valores
    return salida


def contar_palabras(docx: Path) -> dict:
    """Palabras del resumen y del cuerpo (secciones 1 a 9, sin tablas, leyendas, referencias ni anexos)."""
    from docx import Document

    doc = Document(str(docx))
    cuerpo, dentro = 0, False
    for p in doc.paragraphs:
        texto = p.text.strip()
        if p.style.name.startswith("Heading 1"):
            dentro = bool(re.match(r"^[1-9]\. ", texto))
            continue
        if dentro and p.style.name != "Caption" and texto:
            cuerpo += len(texto.split())
    return {"cuerpo": cuerpo}


def armar(args) -> Path:
    X = preparar()
    X["ifc"] = ifc_del_plan(X)
    print(f"IFC del plan: {X['ifc']}")
    print("Documento...")
    W = Informe()
    escribir(W, X, {})  # primera pasada: numera figuras y tablas
    W2 = Informe()
    escribir(W2, X, W.refs)
    SALIDA.mkdir(parents=True, exist_ok=True)
    W2.guardar(DOCX)
    n = contar_palabras(DOCX)
    n["resumen"] = len(X["_resumen"].split())
    print(f"Artículo: {DOCX}  (resumen {n['resumen']} palabras; cuerpo {n['cuerpo']} palabras)")
    return DOCX


def exportar_pdf(docx: Path) -> tuple[Path, int]:
    """Exporta con Word (COM, Windows) y devuelve la ruta del PDF y su número de páginas."""
    pdf = docx.with_suffix(".pdf")
    ps = (f"$w = New-Object -ComObject Word.Application; $w.Visible = $false; "
          f"try {{ $d = $w.Documents.Open('{docx}', $false, $true); $d.Fields.Update() | Out-Null; "
          f"$d.ExportAsFixedFormat('{pdf}', 17); $n = $d.ComputeStatistics(2); $d.Close($false); Write-Output $n }} "
          f"finally {{ $w.Quit() }}")
    r = subprocess.run(["powershell", "-NoProfile", "-Command", ps], capture_output=True, text=True)
    paginas = int(re.findall(r"\d+", r.stdout)[-1]) if re.findall(r"\d+", r.stdout) else -1
    if r.returncode != 0 or paginas < 0:
        print(r.stdout, r.stderr)
    return pdf, paginas


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdf", action="store_true", help="exportar a PDF con Word y contar las páginas")
    ap.add_argument("--profesor", help="nombre del profesor para la portada (por defecto, el de Canvas)")
    args = ap.parse_args()
    if args.profesor is not None:
        global PROFESOR
        PROFESOR = args.profesor.strip()
    docx = armar(args)
    if args.pdf:
        pdf, n = exportar_pdf(docx)
        print(f"PDF: {pdf} ({n} páginas)")


if __name__ == "__main__":
    main()
