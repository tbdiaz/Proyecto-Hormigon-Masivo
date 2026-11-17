"""Plan de Control de Temperatura (informe final, Word), armado desde los resultados del estudio.

Uso, desde la raíz del repositorio y después de correr los estudios (ver docs/informe_final/LEEME.md):

    python scripts/informe_final.py          # docs/informe_final/Plan_Control_Temperatura_G4.docx + figuras
    python scripts/informe_final.py --pdf    # además lo exporta a PDF con Word (Windows) y cuenta las páginas

Todo número del informe sale de los resultados (resultados/estudio_2025, confiabilidad_2025, estres_2025,
retiro_2025): si se vuelven a correr los estudios (p. ej. con otra fecha), basta con volver a correr este script.
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
SALIDA = RAIZ / "docs" / "informe_final"
IMG = SALIDA / "img"
DOCX = SALIDA / "Plan_Control_Temperatura_G4.docx"
PROPUESTA = "unica-M10-T10-N25"
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
            fina = next((r for r in filas if float(r["dx"]) < 0.2 and float(r["dt_s"]) == 600
                         and abs(float(r["profundidad_m"]) - 0.025) < 1e-6), None)
            gruesa = next((r for r in filas if abs(float(r["dx"]) - 0.275) < 1e-6 and float(r["dt_s"]) == 600
                           and abs(float(r["profundidad_m"]) - 0.025) < 1e-6), None)
            esc = next((f for f in filas_por_carpeta.get(carpeta, []) if (f["plan"], f["curva"], f["clima"]) == (plan, curva, clima)), None)
            MD_g = float(gruesa["M_D"]) if gruesa else (esc["M_D"] if esc else None)
            MT_g = float(gruesa["M_T"]) if gruesa else (esc["M_T"] if esc else None)
            if fina is None or MD_g is None:
                continue
            out.append({"plan": plan, "curva": curva, "clima": clima, "dx_fina": float(fina["dx"]),
                        "MD_gruesa": MD_g, "MD_fina": float(fina["M_D"]), "MT_gruesa": MT_g, "MT_fina": float(fina["M_T"]),
                        "delta": float(fina["M_D"]) - MD_g, "otras": [r for r in filas if r is not fina and r is not gruesa],
                        "t_MD_fina": float(fina["t_MD_h"])})
    return out


def _brecha_clima(X: dict) -> float:
    finas = [x["MD_fina"] for x in X.get("sens", []) if x["plan"] == PROPUESTA]
    return (max(finas) - min(finas)) if len(finas) > 1 else float("nan")


def corregir(filas: list[dict], delta: float) -> list[dict]:
    """M_Δ con la corrección de malla (+delta) y los indicadores de falla recalculados con ella."""
    out = []
    for f in filas:
        g = dict(f)
        g["M_D_gruesa"] = f["M_D"]
        g["M_D"] = f["M_D"] + delta
        g["M_D_14d"] = f["M_D_14d"] + delta
        g["I_D"] = float(g["M_D"] > LIM_D)
        g["I_F"] = float(g["I_D"] or f["I_T"])
        out.append(g)
    return out


def serie(carpeta: Path, plan, curva, clima):
    archivos = sorted(glob.glob(str(carpeta / "series" / f"{plan}__{curva}__{clima}__*.npz")))
    return dict(np.load(archivos[-1])) if archivos else None


def cp_superior(k: int, n: int, conf=0.95) -> float:
    """Cota superior de Clopper-Pearson (unilateral) para k fallas en n ensayos, sin scipy (bisección)."""
    if n == 0:
        return float("nan")
    if k == 0:
        return 1 - (1 - conf) ** (1 / n)
    acum = lambda p: sum(math.comb(n, i) * p**i * (1 - p) ** (n - i) for i in range(k + 1))  # noqa: E731
    lo, hi = k / n, 1.0
    for _ in range(80):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if acum(mid) > 1 - conf else (lo, mid)
    return hi


def confiabilidad(filas: list[dict]) -> dict:
    """Estadística del plan propuesto con las 30 temporadas reales: enumeración exacta, cotas y margen."""
    fs = [f for f in filas if f["plan"] == PROPUESTA]
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
    n, k_T, k_D = len(fs), sum(f["I_T"] for f in fs), sum(f["I_D"] for f in fs)
    k_F = sum(f["I_F"] for f in fs)
    peor = max(fs, key=lambda f: f["M_D"])
    peor_T = max(fs, key=lambda f: f["M_T"])
    real = [f for f in fs if f["clima"].endswith("2025")]
    return {
        "n": n, "curvas": curvas, "anios": anios, "MD": MD, "MT": MT, "k_T": k_T, "k_D": k_D, "k_F": k_F,
        "P_T": float(w @ [f["I_T"] for f in fs]), "P_D": float(w @ [f["I_D"] for f in fs]),
        "P_F": float(w @ [f["I_F"] for f in fs]),
        "cp_total": cp_superior(int(k_F), n), "cp_anios": cp_superior(0, len(anios)), "cp_curvas": cp_superior(0, len(curvas)),
        "MD_p50": float(np.nanmedian(MD)), "MD_p95": float(np.nanquantile(MD, 0.95)), "MD_max": float(np.nanmax(MD)),
        "MT_p50": float(np.nanmedian(MT)), "MT_p95": float(np.nanquantile(MT, 0.95)), "MT_max": float(np.nanmax(MT)),
        "peor": peor, "peor_T": peor_T, "real": real, "filas": fs,
        "MD_real_max": max((f["M_D"] for f in real), default=float("nan")),
        "MT_real_max": max((f["M_T"] for f in real), default=float("nan")),
    }


def _p_normal_mayor(x, mu, sd):
    return 0.5 * math.erfc((x - mu) / (sd * math.sqrt(2))) if sd > 0 else float(x < mu)


def margen_material(filas_base: list[dict], filas_estres: list[dict], fuentes: dict, plan=PROPUESTA, escala=1.2,
                    n_boot=4000, seed=11) -> dict:
    """¿Cuánto más calor tendría que liberar el hormigón para que M_Δ llegue a 20 °C?

    El clima casi no mueve M_Δ con la manta (las 30 temporadas, ±0,05 °C por curva): la incertidumbre que domina
    es la del material. Se ajusta M_Δ (peor temporada o tercil) contra la energía a 14 días de cada curva SS, con y
    sin el +20 % de la prueba de estrés, y se despeja la energía crítica E*. Luego se estima P(E > E*) con la
    energía de las probetas (normal ajustada a sus 8 valores, y a las 4 SS; bootstrap para la cota).
    """
    E = {n: fuentes[n].energia(336.0) / 1e6 for n in fuentes}  # MJ/m3
    base, estres = {}, {}
    for f in filas_base:
        if f["plan"] == plan:
            base[f["curva"]] = max(base.get(f["curva"], -np.inf), f["M_D"])
    for f in filas_estres:
        if f["plan"] == plan:
            estres[f["curva"]] = max(estres.get(f["curva"], -np.inf), f["M_D"])
    ss = [c for c in base if c.endswith("_SS") and c in estres]
    if len(ss) < 2:
        return {}
    x = np.array([E[c] for c in ss] + [escala * E[c] for c in ss])
    y = np.array([base[c] for c in ss] + [estres[c] for c in ss])
    b, a = np.polyfit(x, y, 1)
    r2 = 1 - np.sum((y - (a + b * x)) ** 2) / np.sum((y - y.mean()) ** 2)
    E_crit = (LIM_D - a) / b
    todas = np.array(list(E.values()))
    solo_ss = np.array([E[c] for c in E if c.endswith("_SS")])
    rng = np.random.default_rng(seed)
    boot = [_p_normal_mayor(E_crit, m.mean(), m.std(ddof=1))
            for m in (rng.choice(todas, len(todas)) for _ in range(n_boot))]
    e_max = max(E.values())
    return {"E": E, "base": base, "estres": estres, "a": float(a), "b": float(b), "r2": float(r2),
            "E_crit": float(E_crit), "E_max": float(e_max), "curva_max": max(E, key=E.get),
            "exceso": float(E_crit / e_max - 1), "mu": float(todas.mean()), "sd": float(todas.std(ddof=1)),
            "p": _p_normal_mayor(E_crit, todas.mean(), todas.std(ddof=1)), "p_95": float(np.quantile(boot, 0.95)),
            "p_ss": _p_normal_mayor(E_crit, solo_ss.mean(), solo_ss.std(ddof=1))}


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
    from mcplan.clima import historias_era5, leer_snapshot

    era5 = leer_snapshot(Path(cfg["_ruta"]).parent / cfg["clima"]["climatologia"])
    return historias_era5(era5, cfg["t0"], 24 * 30), era5["metadata"]


# ---------------------------------------------------------------------------
# Figuras
# ---------------------------------------------------------------------------

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


def fig_clima(historias, ruta, dias=14):
    nombres = [n for n, _ in historias]
    T = np.array([c.T_aire for _, c in historias])[:, : 24 * dias + 1]
    v = np.array([c.viento for _, c in historias])[:, : 24 * dias + 1]
    t = np.arange(T.shape[1]) / 24
    real = nombres.index("ERA5 2025")
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(16.5 * CM, 8.4 * CM), sharex=True,
                                 gridspec_kw={"height_ratios": [1.6, 1]})
    for ax, X, unidad in ((a1, T, "T del aire [°C]"), (a2, v, "viento a 10 m [m/s]")):
        lo, mid, hi = np.percentile(X, [10, 50, 90], axis=0)
        ax.fill_between(t, lo, hi, color=GRIS_CLARO, alpha=0.55, lw=0, label="1996–2025: 10–90 %")
        ax.plot(t, mid, color=GRIS, lw=1.0, label="1996–2025: mediana")
        ax.plot(t, X[real], color=AZUL, lw=1.4, label="2025 (el clima que ocurrió)")
        ax.set_ylabel(unidad)
    a1.legend(ncol=3, loc="upper left", bbox_to_anchor=(0, 1.22))
    a2.set_xlabel("días desde el inicio de la vertida (mié 19 nov 2025, 08:00)")
    a2.set_xlim(0, dias)
    a2.set_xticks(range(0, dias + 1, 2))
    fig.subplots_adjust(hspace=0.12)
    return _guardar(fig, ruta)


def fig_cronograma(inicios, fines, D, ruta, retiro_d=26, tubos_d=14, monitoreo_d=29):
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
    filas = [("Monitoreo", 0, monitoreo_d, GRIS), ("Moldaje + manta 10 cm", 0, retiro_d, TINTA2),
             ("Tubos V5 (agua 25 °C)", inicios[4] / 24, inicios[4] / 24 + tubos_d, COLOR_V[5]),
             ("Tubos V4 (agua 25 °C)", inicios[3] / 24, inicios[3] / 24 + tubos_d, COLOR_V[4]),
             ("Tubos V3 (agua 25 °C)", inicios[2] / 24, inicios[2] / 24 + tubos_d, COLOR_V[3]),
             ("Horizonte de evaluación", 0, 14, EJE)]
    for k, (nombre, a, b, c) in enumerate(filas):
        a2.barh(k, b - a, left=a, height=0.6, color=c, lw=0)
        a2.text(b + 0.3, k, nombre, va="center", fontsize=7.5, color=TINTA)
    a2.axvline(retiro_d, color=TINTA, lw=0.8)
    a2.text(retiro_d, len(filas) - 0.35, f" retiro día {retiro_d}", fontsize=7.5, color=TINTA, va="bottom")
    a2.set_yticks([])
    a2.set_xlim(0, 40)
    a2.set_xticks(range(0, 31, 5))
    a2.set_xlabel("días desde t₀")
    a2.set_title("Protecciones, tubos y monitoreo")
    a2.grid(axis="y", visible=False)
    a2.spines["left"].set_visible(False)
    fig.subplots_adjust(wspace=0.22)
    return _guardar(fig, ruta)


def fig_temperaturas(s_peor, rotulo_peor, series_anios, anio_peor, ruta):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(16.5 * CM, 6.4 * CM))
    t = s_peor["t_h"] / 24
    a1.plot(t, s_peor["Tmax"], color=NARANJA, lw=1.6, label="T máx. (núcleo)")
    a1.plot(t, s_peor["Tmin"], color=AZUL, lw=1.6, label="T mín. a 25 mm")
    a1.plot(t, s_peor["T_amb"], color=GRIS, lw=0.8, label="aire")
    a1.set_xlim(0, 14)
    a1.set_xticks(range(0, 15, 2))
    a1.set_xlabel("días desde t₀")
    a1.set_ylabel("°C")
    a1.set_title(f"Peor escenario: {rotulo_peor}")
    a1.legend(loc="upper right", ncol=1)
    for nombre, s in series_anios.items():
        if nombre in (anio_peor, "ERA5_2025"):
            continue
        a2.plot(s["t_h"] / 24, s["dT"], color=GRIS_CLARO, lw=0.7)
    etiquetas = ((anio_peor, NARANJA, f"peor temporada ({anio_peor.split('_')[-1]})"),
                 ("ERA5_2025", AZUL, "2025 (clima real)"))
    for nombre, color, et in etiquetas:
        if nombre in series_anios and not (nombre == "ERA5_2025" and anio_peor == "ERA5_2025"):
            s = series_anios[nombre]
            a2.plot(s["t_h"] / 24, s["dT"], color=color, lw=1.6, label=et)
    a2.axhline(LIM_D, color=TINTA, lw=1.0)
    a2.text(0.2, LIM_D + 0.4, "límite 20 °C", fontsize=7.5, color=TINTA)
    a2.set_xlim(0, 14)
    a2.set_xticks(range(0, 15, 2))
    a2.set_ylim(0, 22)
    a2.set_xlabel("días desde t₀")
    a2.set_ylabel("ΔT = T máx − T mín 25 mm [°C]")
    curva = rotulo_peor.split(",")[0]
    a2.set_title(f"ΔT(t) con {curva}: las 30 temporadas reales")
    a2.legend(loc="lower right")
    fig.subplots_adjust(wspace=0.25)
    return _guardar(fig, ruta)


def fig_confiabilidad(C, ruta):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(16.5 * CM, 6.6 * CM), sharey=True)
    curvas, anios = C["curvas"], C["anios"]
    rng = np.random.default_rng(3)
    for i, c in enumerate(curvas):
        y = len(curvas) - 1 - i
        jit = rng.uniform(-0.22, 0.22, len(anios))
        for ax, M in ((a1, C["MD"]), (a2, C["MT"])):
            otros = [j for j, a in enumerate(anios) if not a.endswith("2025")]
            ax.scatter(M[i, otros], y + jit[otros], s=9, color=AZUL, alpha=0.55, lw=0)
            j25 = [j for j, a in enumerate(anios) if a.endswith("2025")]
            ax.scatter(M[i, j25], [y] * len(j25), s=26, color=NARANJA, edgecolor="white", lw=0.8, zorder=3)
    a1.axvline(LIM_D, color=TINTA, lw=1.0)
    a1.text(LIM_D - 0.15, -0.45, "límite 20 °C", ha="right", fontsize=7.5, color=TINTA)
    a2.axvline(LIM_T, color=TINTA, lw=1.0)
    a2.text(LIM_T - 0.5, -0.45, "límite 70 °C", ha="right", fontsize=7.5, color=TINTA)
    a1.set_yticks(range(len(curvas)), list(reversed(curvas)))
    a1.set_xlabel(r"M$_\Delta$: ΔT máximo en 14 días [°C]")
    a2.set_xlabel(r"M$_T$: temperatura máxima [°C]")
    a1.set_ylim(-0.7, len(curvas) - 0.4)
    a1.set_xlim(min(10, np.nanmin(C["MD"]) - 1), 21)
    a2.set_xlim(min(30, np.nanmin(C["MT"]) - 2), 72)
    a1.set_title("Diferencial (8 curvas × 30 temporadas)")
    a2.set_title("Temperatura máxima")
    a1.scatter([], [], s=9, color=AZUL, alpha=0.55, label="temporadas 1996–2024")
    a1.scatter([], [], s=26, color=NARANJA, label="2025 (clima real)")
    a1.legend(loc="lower left", bbox_to_anchor=(0, 1.08), ncol=2)
    a1.grid(axis="y", visible=False)
    a2.grid(axis="y", visible=False)
    fig.subplots_adjust(wspace=0.08)
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


def fig_retiro(series_planes, ruta, retiros_d):
    fig, ax = plt.subplots(figsize=(16.5 * CM, 5.6 * CM))
    colores = [AZUL, NARANJA, AQUA, AMARILLO]
    for (nombre, s), c, d in zip(series_planes.items(), colores, retiros_d):
        ax.plot(s["t_h"] / 24, s["dT"], color=c, lw=1.5, label=nombre)
        ax.axvline(d, color=c, lw=0.8, alpha=0.7)
    ax.axhline(LIM_D, color=TINTA, lw=1.0)
    ax.text(37.8, LIM_D + 0.5, "límite 20 °C", fontsize=7.5, color=TINTA, ha="right")
    ax.set_xlim(0, 38)
    ax.set_ylim(0, 24)
    ax.set_xlabel("días desde t₀ (cada línea vertical marca un retiro)")
    ax.set_ylabel("ΔT [°C]")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=4, fontsize=7.5,
              title="retiro de moldaje y mantas el:", title_fontsize=7.5)
    ax.set_ylim(0, 28)
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
        r = p.add_run("Plan de Control de Temperatura · Fundación de aerogenerador · Grupo 4, Litueche")
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


def chica(p) -> str:
    """Probabilidad muy chica en notación científica legible: 3,2·10^−6^."""
    if p is None or not math.isfinite(p):
        return "—"
    if p < 1e-12:
        return "< 10^−12^"
    e = math.floor(math.log10(p))
    return f"{f1(p / 10**e, 1)}·10^{str(e).replace('-', '−')}^"


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
    "unica-M10-A25": ("Única + manta + tubos en todo", "18", "10 cm", "V1–V5, 25 °C"),
    "unica-M5-T10-A25": ("Única + manta 5 cm + preenfriado + tubos en todo", "10", "5 cm", "V1–V5, 25 °C"),
    "unica-M10-T10-A25": ("Única + manta + preenfriado + tubos en todo", "10", "10 cm", "V1–V5, 25 °C"),
    "unica-M10-T8-A25": ("Ídem, colocación a 8 °C", "8", "10 cm", "V1–V5, 25 °C"),
    "unica-M10-T10-A25-S75": ("Ídem, tubos cada 0,75 m", "10", "10 cm", "V1–V5, 25 °C, 0,75 m"),
    "unica-M10-T8-A25-S75": ("Ídem, 8 °C y tubos cada 0,75 m", "8", "10 cm", "V1–V5, 25 °C, 0,75 m"),
    "unica-M10-T10-N25": ("**Propuesta:** única + manta + preenfriado + tubos en el núcleo", "10", "10 cm", "V3–V5, 25 °C"),
    "unica-M10-T10-N20": ("Propuesta con agua a 20 °C", "10", "10 cm", "V3–V5, 20 °C"),
    "unica-M10-T8-N20": ("Propuesta con agua a 20 °C y 8 °C", "8", "10 cm", "V3–V5, 20 °C"),
    "unica-M10-T8-N25": ("Propuesta colocando a 8 °C (meta de la planta)", "8", "10 cm", "V3–V5, 25 °C"),
    "unica-M10-T12-N25": ("Propuesta colocando a 12 °C", "12", "10 cm", "V3–V5, 25 °C"),
    "unica-M10-T14-N25": ("Propuesta colocando a 14 °C", "14", "10 cm", "V3–V5, 25 °C"),
    "unica-M10-T16-N25": ("Propuesta colocando a 16 °C", "16", "10 cm", "V3–V5, 25 °C"),
    "5V-M10-T10-A25-e12": ("5 vertidas, esperas de 12 h (más lenta)", "10", "10 cm", "V1–V5, 25 °C"),
    "5V-M10-T10-A25-e24": ("5 vertidas, esperas de 24 h (más lenta)", "10", "10 cm", "V1–V5, 25 °C"),
}


RETIROS = [("propuesta-retiro-d14", 14, 312.0), ("propuesta-retiro-d21", 21, 480.0), ("propuesta-retiro-d26", 26, 600.0),
           ("propuesta-retiro-d33", 33, 768.0)]
CARPETAS_RETIRO = ("retiro_2025", "retiro_tardio_2025")


def serie_retiro(plan, curva, clima):
    for c in CARPETAS_RETIRO:
        s_ = serie(RES / c, plan, curva, clima)
        if s_ is not None:
            return s_
    return None
DIAS = ["lunes", "martes", "miércoles", "jueves", "viernes", "sábado", "domingo"]
MESES = ["enero", "febrero", "marzo", "abril", "mayo", "junio", "julio", "agosto", "septiembre", "octubre",
         "noviembre", "diciembre"]


def fecha_larga(t: datetime) -> str:
    return f"{DIAS[t.weekday()]} {t.day} de {MESES[t.month - 1]} de {t.year}, {t:%H:%M}"


def contar_tests() -> int:
    r = subprocess.run([sys.executable, "-m", "pytest", "--collect-only", "-q"], cwd=RAIZ, capture_output=True, text=True)
    m = re.search(r"(\d+) tests? collected", r.stdout)
    return int(m.group(1)) if m else -1


def criterio_retiro(ret_dir: Path, plan: str, retiro_h: float, D: float, filas: list[dict]) -> dict:
    """Lo que medirían los sensores al momento del retiro, sólo en los escenarios que después cumplen (con la
    corrección de malla): esas condiciones son las que el criterio de obra exige antes de retirar."""
    dT, nucleo_aire = [], []
    for f in filas:
        if f["plan"] != plan or f["M_D"] > LIM_D:
            continue
        s = serie_retiro(plan, f["curva"], f["clima"])
        if s is None:
            continue
        t = s["t_h"]
        i = int(np.searchsorted(t, retiro_h + D))
        j = int(np.searchsorted(t, retiro_h + D + 72))
        if i >= len(t):
            continue
        dT.append(float(s["dT"][i]))
        nucleo_aire.append(float(s["Tmax"][i] - np.min(s["T_amb"][i:j])))
    return {"dT_max": max(dT, default=float("nan")), "nucleo_aire_max": max(nucleo_aire, default=float("nan")),
            "n": len(dT)}


def preparar() -> dict:
    """Lee resultados, calcula estadísticas y dibuja las figuras (una sola vez)."""
    IMG.mkdir(parents=True, exist_ok=True)
    X = {}
    cfg_fam = RAIZ / "config" / "estudio_2025.yaml"
    X["fam_dir"], X["conf_dir"] = RES / "estudio_2025", RES / "confiabilidad_2025"
    X["est_dir"], X["ret_dir"] = RES / "estres_2025", RES / "retiro_2025"
    filas_fam, filas_conf = filas_de(X["fam_dir"]), filas_de(X["conf_dir"])
    filas_est = filas_de(X["est_dir"])
    filas_ret = [f for c in CARPETAS_RETIRO for f in filas_de(RES / c)]
    # Corrección de malla: cuánto sube M_Δ al pasar de Δx = 0,275 a 0,11 m en los escenarios verificados con malla
    # fina; se usa la mayor de las diferencias de la propuesta (conservador) en todas las corridas del estudio.
    X["sens"] = sensibilidades([X["conf_dir"], X["fam_dir"]], {X["conf_dir"]: filas_conf, X["fam_dir"]: filas_fam})
    deltas = [s_["delta"] for s_ in X["sens"] if s_["plan"] == PROPUESTA]
    X["delta"] = max(deltas) if deltas else 0.0
    d = X["delta"]
    X["filas_fam_c"] = corregir(filas_fam, d)
    X["fam"] = resumir(X["filas_fam_c"])
    X["fam_gruesa"] = resumir(filas_fam)
    X["MT_familia"] = max((f["M_T"] for f in filas_fam), default=float("nan"))
    X["C"] = confiabilidad(corregir(filas_conf, d))
    X["C_gruesa"] = confiabilidad(filas_conf)
    filas_est = corregir(filas_est, d)
    X["est"] = resumir(filas_est).get(PROPUESTA, {})
    X["est8"] = resumir(filas_est).get("unica-M10-T8-N25", {})
    filas_ret = corregir(filas_ret, d)
    X["ret"] = resumir(filas_ret)
    ruta_def = X["fam_dir"] / "definicion.json"
    X["definicion"] = json.loads(ruta_def.read_text(encoding="utf-8")) if ruta_def.exists() else {}

    print("Geometría, cronograma y clima...")
    X["vox"], X["geo"] = geometria()
    X["inicios"], X["fines"], X["D"], cfg = cronograma_vertidas(cfg_fam)
    X["t0"] = datetime.fromisoformat(cfg["t0"])
    historias, X["meta_era5"] = historias_reales(cfg)
    c25 = dict(historias)["ERA5 2025"]
    X["T25"], X["v25"] = np.array(c25.T_aire[: 24 * 14 + 1]), np.array(c25.viento[: 24 * 14 + 1])
    from mcplan.clima import indice_frio

    indices = {n: indice_frio(c) for n, c in historias}
    X["rango_frio_2025"] = sorted(indices, key=indices.get).index("ERA5 2025") + 1
    X["n_temporadas"] = len(historias)

    # Retiro: el primer día que cumple en todos los escenarios con 1 °C de margen (si no, el primero que cumple).
    ret = X["ret"]
    evaluados = [(p, d, h) for p, d, h in RETIROS if p in ret]
    ok = [(p, d, h) for p, d, h in evaluados if ret[p]["P_D"] == 0]
    con_margen = [(p, d, h) for p, d, h in ok if ret[p]["MD_max"] <= LIM_D - 1.0]
    # Si ninguna fecha cumple en todos los escenarios, la de menos fallas (y menor M_Δ): retirar más tarde no
    # siempre ayuda, porque bajo la manta el núcleo sigue caliente por semanas.
    menos_fallas = sorted(evaluados, key=lambda e: (ret[e[0]]["P_D"], ret[e[0]]["MD_max"]))[:1]
    X["retiro"] = (con_margen or ok or menos_fallas or [RETIROS[-1]])[0]
    X["retiro_ok"] = bool(ok)
    X["crit_retiro"] = criterio_retiro(X["ret_dir"], X["retiro"][0], X["retiro"][2], X["D"], filas_ret) \
        if X["retiro"][0] in ret else {}

    # T máxima de colocación: la más alta que cumple en los 24 escenarios, con todas las menores cumpliendo también.
    fam = X["fam"]
    X["tmax_col"] = None
    for p, T in ((PROPUESTA, 10), ("unica-M10-T12-N25", 12), ("unica-M10-T14-N25", 14), ("unica-M10-T16-N25", 16)):
        if p in fam and fam[p]["P_F"] == 0:
            X["tmax_col"] = T
        else:
            break

    from mcplan.curvas import cargar_fuentes

    fuentes = cargar_fuentes(RAIZ / "data" / "hidratacion" / "source_terms.txt", desfase_h=2.0, cola="exponencial")
    X["mat"] = margen_material(X["C"]["filas"] if X["C"] else [], filas_est, fuentes, PROPUESTA)
    X["mat8"] = margen_material(X["filas_fam_c"], filas_est, fuentes, "unica-M10-T8-N25")
    print("Figuras...")
    X["figs"] = figs = {"geo": fig_geometria(X["vox"], X["geo"], IMG / "f_geometria.png"),
                        "curvas": fig_curvas(fuentes, IMG / "f_curvas.png"),
                        "clima": fig_clima(historias, IMG / "f_clima.png"),
                        "crono": fig_cronograma(X["inicios"], X["fines"], X["D"], IMG / "f_cronograma.png",
                                                retiro_d=X["retiro"][1])}
    C = X["C"]
    X["sensores"], X["zonas"] = [], Counter()
    if C:
        peor = C["peor"]
        s_peor = serie(X["conf_dir"], PROPUESTA, peor["curva"], peor["clima"])
        series_anios = {a: serie(X["conf_dir"], PROPUESTA, peor["curva"], a) for a in C["anios"]}
        series_anios = {a: s for a, s in series_anios.items() if s is not None}
        rot = f"{peor['curva']}, temporada {peor['clima'].split('_')[-1]}"
        if s_peor is not None:
            figs["temp"] = fig_temperaturas(s_peor, rot, series_anios, peor["clima"], IMG / "f_temperaturas.png")
        figs["conf"] = fig_confiabilidad(C, IMG / "f_confiabilidad.png")
        X["sensores"], cal, fri, X["zonas"] = plan_sensores(C)
        figs["sens"] = fig_sensores(X["vox"], X["geo"], cal, fri, X["sensores"], IMG / "f_sensores.png")
    hechos = [(p, d, h) for p, d, h in RETIROS if p in ret]
    if hechos:
        sr, xs = {}, []
        for p, d, h in hechos:
            pr = ret[p]["peor"]
            s = serie_retiro(p, pr["curva"], pr["clima"])
            if s is not None:
                sr[f"día {d} ({pr['curva']}, {pr['clima']})"] = s
                xs.append((h + X["D"]) / 24)
        if sr:
            figs["retiro"] = fig_retiro(sr, IMG / "f_retiro.png", xs)
    detalle = sorted(X["conf_dir"].glob(f"detalle_{PROPUESTA}_*"))
    d0 = detalle[0] if detalle else None
    X["corte"] = d0 / "cortes_y0.png" if d0 and (d0 / "cortes_y0.png").exists() else None
    X["vista"] = d0 / "vista_3d.png" if d0 and (d0 / "vista_3d.png").exists() else None
    X["tubos"] = json.loads((d0 / "tubos.json").read_text(encoding="utf-8")) if d0 and (d0 / "tubos.json").exists() else {}
    X["n_tests"] = contar_tests()
    X["eq1"] = ecuacion(r"$\rho\,c_p\,\dfrac{\partial T}{\partial t}=\nabla\cdot(k\,\nabla T)+\dot q_h(a),"
                        r"\qquad -k\,\dfrac{\partial T}{\partial n}=U\,(T-T_{eq}),"
                        r"\qquad U=\left(\dfrac{1}{h_{ext}}+R_{capas}+\dfrac{\Delta x}{2k}\right)^{-1}$",
                        IMG / "eq_modelo.png", tam=11)
    X["eq2"] = ecuacion(r"$M_\Delta=\max_{t\leq 14\,d}\left[\max_{\Omega_c}T-\min_{S_{25}}T\right],\qquad "
                        r"P_f=\sum_{m,k}p_m\,q_k\,\mathbf{1}\left[M_T>70\,^\circ\mathrm{C}\ \cup\ M_\Delta>20\,^\circ\mathrm{C}\right]$",
                        IMG / "eq_criterio.png", tam=11)
    return X


def escribir(W: Informe, X: dict, R: dict):
    """Todo el texto del informe. R: números de figuras y tablas por etiqueta (de una primera pasada)."""
    ref = lambda k: str(R.get(k, "?"))  # noqa: E731
    C, fam, est, ret, figs = X["C"], X["fam"], X["est"], X["ret"], X["figs"]
    geo, inicios, fines, D, t0 = X["geo"], X["inicios"], X["fines"], X["D"], X["t0"]
    sensores, tmax_col = X["sensores"], X["tmax_col"]
    tcol = f"{tmax_col} °C" if tmax_col else "—"
    plan_ret, dia_ret, h_ret = X["retiro"]
    prop = fam.get(PROPUESTA, {})
    n_sens = sum(s["n"] for s in sensores) + 6

    # ---------------- Portada
    W.p("")
    W.p("**PLAN DE CONTROL DE TEMPERATURA**", tam=22, alinear="centro", color="#184f95", despues=2)
    W.p("Fundación de hormigón masivo de un aerogenerador", tam=14, alinear="centro", despues=2)
    W.p("Proyecto 2 · Planificación de vertidas de hormigón masivo mediante simulación", tam=11, alinear="centro",
        color=TINTA2, despues=12)
    if X["vista"] is not None:
        W.imagen(X["vista"], 13.0)
        W.p("Campo de temperatura simulado de la propuesta en su peor escenario.", tam=8.5, alinear="centro",
            color=TINTA2, despues=12)
    W.p("**Grupo 4 — Sitio Litueche** (34,105° S; 71,650° O)", tam=12, alinear="centro", despues=10)
    for linea in ("Integrantes: ____________________________________________",
                  "Curso: ____________________________   Profesor: ____________________________",
                  "Fecha: 16 de noviembre de 2026"):
        W.p(linea, tam=11, alinear="centro", despues=6)
    W.p(f"Repositorio: github.com/tbdiaz/Proyecto-Hormigon-Masivo · commit {git_commit()}", tam=9,
        alinear="centro", color=TINTA2)
    W.salto()

    # ---------------- Resumen
    W.h1("Resumen ejecutivo")
    W.p(f"Este Plan de Control de Temperatura (PCT), inspirado en la Especificación para Hormigón Masivo "
        f"Estructural del ICH (2019) como pide el curso, define cómo construir la fundación del aerogenerador (636,8 m³ en las cinco "
        f"vertidas del IFC) para que **T~máx~ ≤ 70 °C** y **ΔT~máx~ ≤ 20 °C** durante los 14 días siguientes a la "
        f"primera vertida. La construcción se analiza el **{fecha_larga(t0)}**: la fecha de la presentación, un año "
        f"antes, para que todo el horizonte tenga **clima real medido** (reanálisis ERA5) y el plan pueda verificarse "
        f"contra lo que efectivamente ocurrió.", alinear="justificado")
    filas = [
        ["Secuencia", f"Una sola colocación continua V1 → V5, por capas desde abajo, a 27 m³/h: termina a las "
                      f"{f1(D)} h ({hora(t0, D)})."],
        ["Temperatura de colocación", (f"**Máximo 10 °C a la descarga**; la planta apunta a 8 °C (hormigón "
                                       f"preenfriado con hielo, ≈ 50–60 kg/m³)." if tmax_col == 10 else
                                       f"Objetivo **10 °C** (hormigón preenfriado con hielo, ≈ 50 kg/m³); máximo "
                                       f"admisible **{tcol}** a la descarga.")],
        ["Moldaje y aislación", "Terciado de 18 mm en V1, V2 y V5 y **manta de 10 cm (R = 2,5 m²K/W) en todas las caras "
                                f"expuestas**, colocadas a medida que sube el hormigón; se retiran **desde el día {dia_ret}**, según los sensores."],
        ["Tubos de enfriamiento", "Serpentín HDPE Ø32 mm a media altura de **V3, V4 y V5** (filas cada 1 m), 20 L/min "
                                  "por circuito, **agua a 25 ± 2 °C** durante 14 días. V1 y V2 sin tubos."],
        ["Monitoreo", f"{n_sens} sensores de temperatura (núcleo, 25 mm bajo la superficie, suelo, aire y agua) en los puntos críticos "
                      "que entrega el modelo; registro cada 15 min hasta las 72 h y luego cada 30 min."],
    ]
    W.tabla("Decisiones del plan propuesto.", ["Decisión", "Valor"], filas, [4.0, 13.2], alinear=["l", "l"])
    if C:
        filas = [
            ["Temperatura máxima prevista, M~T~", f"P50 {f1(C['MT_p50'])} °C · P95 {f1(C['MT_p95'])} °C · máx. {f1(C['MT_max'])} °C"],
            ["Diferencial máximo previsto, M~Δ~ (con corrección de malla)",
             f"P50 {f1(C['MD_p50'])} °C · P95 {f1(C['MD_p95'])} °C · máx. {f1(C['MD_max'])} °C"],
            ["P(T~máx~ > 70 °C)", f"{pct(C['P_T'])}  ({int(C['k_T'])} de {C['n']} escenarios)"],
            ["P(ΔT~máx~ > 20 °C)", f"{pct(C['P_D'])}  ({int(C['k_D'])} de {C['n']})"
                                   + (f"; por el margen del material ≈ {pct(X['mat8']['p'], 1)} con la planta a 8 °C "
                                      f"y {pct(X['mat']['p'], 1)} si todo llegara a 10 °C (sección 9)"
                                      if X.get("mat") and X.get("mat8") else "")],
            ["Con el clima real de 2025", (f"**cumple**: con las 8 curvas, M~Δ~ ≤ {f1(C['MD_real_max'])} °C (con la "
                                           f"corrección conservadora de malla) y M~T~ ≤ {f1(C['MT_real_max'])} °C"
                                           if C["MD_real_max"] <= LIM_D else
                                           f"**no cumple**: M~Δ~ hasta {f1(C['MD_real_max'])} °C")
             + next((f"; calculado con malla fina, {x['curva']} da {f1(x['MD_fina'])} °C" for x in X["sens"]
                     if x["plan"] == PROPUESTA and x["clima"].endswith("2025")), "")],
        ]
        W.tabla(f"Resultados de la propuesta: 8 curvas de hidratación × {len(C['anios'])} temporadas reales "
                f"(1996–2025) = {C['n']} simulaciones 3D, con el hormigón en el máximo de 10 °C.",
                ["Indicador", "Valor"], filas, [5.2, 12.0],
                alinear=["l", "l"])
    W.p(f"**Por qué este plan.** Con las curvas de hidratación medidas, la temperatura máxima nunca es el problema "
        f"(M~T~ ≤ {f1(X['MT_familia'])} °C en todos los planes evaluados); lo que gobierna es el diferencial entre el "
        f"núcleo y la superficie. La vertida única es la más rápida posible con una bomba y evita que las capas "
        f"inferiores se enfríen mientras se espera. Ninguna medida basta sola: la manta protege la superficie, el "
        f"preenfriado baja el núcleo y los tubos con agua tibia, sólo en el núcleo (V3–V5), recortan el pico sin "
        f"enfriar el anillo exterior de la losa, que ya es la zona más fría (Tabla {ref('t:familia')}).",
        alinear="justificado")

    # ---------------- 1. Alcance y criterios
    W.h1("1. Objetivo, alcance y criterios")
    W.p("El PCT responde a la pregunta del proyecto: **¿cómo ejecutar las vertidas para mantener el riesgo térmico "
        "dentro de los límites, considerando la incertidumbre del hormigón, de la ejecución y del clima?** Los criterios "
        "del curso son T~máx~ ≤ 70 °C (igual al límite general de ACI 301) y ΔT~máx~ ≤ 20 °C fijo durante los 14 días, "
        "más conservador que la Tabla 4 del ICH (20, 25 o 30 °C según la resistencia). Cumplir el ΔT se trata como un "
        "control del riesgo de agrietamiento térmico, no como una demostración mecánica.", alinear="justificado")
    W.ecuacion(X["eq2"])
    W.p("T~mín~ se mide en las celdas a **25 mm bajo las caras expuestas** (también las cubiertas por moldaje o manta), "
        "en caras planas; las aristas de la escalera de vóxeles se informan aparte. M~T~ es el máximo de T en el "
        "hormigón. P~f~ suma los pesos de los escenarios que fallan en cualquiera de los dos criterios, sin contar dos "
        "veces una corrida. La meta del plan es **P~f~ ≤ 5 %**.", alinear="justificado")

    # ---------------- 2. Fundación, sitio y fecha
    W.h1("2. Fundación, sitio y fecha de construcción")
    W.p("La geometría es la oficial del curso (Bonsai/IFC), basada en la fundación de Ju y Lei (2019): losa octogonal "
        "de 18 m en dos vertidas de 0,55 m, dos troncos con caras inclinadas (≈22°) y un pedestal cilíndrico de 7 m de "
        "diámetro. El solver lee el IFC con IfcOpenShell (geometría y property sets de cada vertida) y lo voxeliza; con "
        "Δx = 0,275 m el error de volumen es +0,38 % (+0,06 % con 0,11 m).", alinear="justificado")
    W.figura(figs["geo"], "Las cinco vertidas del IFC y su volumen. Cada vertida conserva su color en todo el informe.",
             ref="geo")
    pm = X["meta_era5"].get("punto_modelo", {})
    W.p(f"**Sitio:** Litueche (−34,105457; −71,649539); el punto de grilla de ERA5 está en "
        f"({f1(pm.get('latitud'), 2)}; {f1(pm.get('longitud'), 2)}; {f1(pm.get('elevacion_m'), 0)} m), con "
        f"hora local America/Santiago. **Fecha:** t~0~ = {fecha_larga(t0)}. El curso fija como fecha hipotética la de "
        f"la defensa; el grupo la trasladó al mismo día de 2025 para trabajar con clima real en todo el horizonte, en "
        f"vez de un pronóstico que hoy cubre sólo 11 de los 14 días. El flujo con pronóstico (ensemble ECMWF de 51 "
        f"miembros más climatología) está implementado y se usa igual para la fecha final.", alinear="justificado")

    # ---------------- 3. Materiales
    W.h1("3. Mezcla, materiales y parámetros térmicos")
    su = X["definicion"].get("suelo", {})
    filas = [
        ["Hormigón", "k = 2,222 W/mK · ρ = 2400 kg/m³ · c~p~ = 970 J/kgK (ρc~p~ = 2,33 MJ/m³K)", "Ju y Lei (2019)"],
        ["Mezcla de las probetas", "350,77 kg/m³ de cemento; ρ = 2476 kg/m³; c~p~ = 1023 J/kgK", "encabezado de las curvas"],
        ["Suelo", "k = 1,389 W/mK · ρ = 1800 kg/m³ · c~p~ = 1273 J/kgK; 3 m bajo la base y 2 m de margen", "Ju y Lei (2019)"],
        ["Suelo: estado inicial",
         (f"real (ERA5) el 19 nov 2025: {', '.join(f1(x) for x in su['T_capas_C'])} °C a "
          f"{', '.join(f1(x, 2) for x in su['profundidades_m'])} m; fondo a 3 m: {f1(su['T_fondo_C'])} °C")
         if su else "perfil real de ERA5", "ERA5"],
        ["Moldaje", "terciado fenólico de 18 mm, k = 0,12 W/mK (R = 0,15 m²K/W), caras verticales de V1, V2 y V5", "supuesto"],
        ["Manta de aislación", "10 cm, k = 0,04 W/mK (R = 2,5 m²K/W); α = 0,4; ε = 0,9", "valor típico, por confirmar"],
        ["Superficie del hormigón", "absortividad solar 0,65; emisividad 0,9", "rango típico"],
        ["Tubos", "HDPE Ø32 mm; agua a 20 L/min por circuito", "docs/tubos.md"],
    ]
    W.tabla("Parámetros térmicos del modelo.", ["Material", "Valor", "Fuente"], filas, [3.6, 10.2, 3.4],
            alinear=["l", "l", "l"])
    W.p("**Temperatura de colocación.** Máximo 10 °C a la descarga, apuntando a 8 °C en planta (sección 9.1). Con "
        "ρc~p~ = 2,33 MJ/m³K, bajar el hormigón de ≈18 °C a 8–10 °C exige retirar 18,6–23,3 MJ/m³. Con hielo en "
        "escamas en reemplazo de parte del agua de amasado (334 kJ/kg al fundirse y ≈40 kJ/kg al entibiarse) son "
        "**≈ 50–62 kg de hielo por m³ (≈ 32–40 t en total)**, junto con agua fría y áridos a la sombra. Es un balance "
        "de energía que la planta debe confirmar con su dosificación.", alinear="justificado")

    # ---------------- 4. Curvas
    W.h1("4. Calor de hidratación y aumento adiabático")
    W.p("El calor viene de las **8 curvas Q(t) entregadas por el profesor** (versión corregida del 14 nov 2026): cuatro "
        "probetas (M1–M4) en dos campañas, SS (240 h) y 22A (163 h, marcada de menor confianza), obtenidas por "
        "inversión con volúmenes finitos y regularización de Tikhonov con β = 2,5 W/m²K. Cada celda recibe Q(a) según "
        "su propia edad a = t − t~colocación~, sin madurez ni Arrhenius (limitación declarada: la fuente no se acelera "
        "en el núcleo caliente). Después del registro sigue una cola exponencial. El aumento adiabático es "
        "ΔT~ad~(a) = ∫Q da / ρc~p~.", alinear="justificado")
    W.figura(figs["curvas"], "Las 8 curvas de hidratación (arriba) y el aumento adiabático que implican con "
                             "ρc~p~ = 2,33 MJ/m³K (abajo), con su valor a los 14 días.", ref="curvas")
    W.p("Las curvas liberan 61–100 MJ/m³ (26–43 °C adiabáticos), menos que la de Ju y Lei (60 °C). Las SS siguen "
        "liberando calor hasta ≈200 h; las 22A se apagan cerca de las 85 h. Se usan como **8 escenarios de material "
        "equiprobables**: cada escenario aplica una curva completa a todo el hormigón.", alinear="justificado")

    # ---------------- 5. Clima y bordes
    W.h1("5. Clima real y condiciones de borde")
    T25, v25 = X["T25"], X["v25"]
    W.p(f"El clima es el **reanálisis ERA5** (Hersbach et al., 2020; Open-Meteo Historical Weather API): temperatura "
        f"del aire a 2 m, viento a 10 m, radiación de onda corta y nubosidad, cada hora, en las **{X['n_temporadas']} "
        f"temporadas 1996–2025**, desde el 19 de noviembre de cada año. Cada temporada es una historia completa y "
        f"real, sin promediar, que conserva el ciclo diario, los frentes fríos y la coherencia entre variables; la de "
        f"2025 es la que ocurrió. En sus primeros 14 días el aire varió entre {f1(T25.min())} y {f1(T25.max())} °C "
        f"(media {f1(T25.mean())} °C) y el viento promedió {f1(v25.mean())} m/s. Por su índice de frío (media de las "
        f"mínimas de los primeros 7 días), 2025 fue la temporada n.º {X['rango_frio_2025']} más fría de "
        f"{X['n_temporadas']}.", alinear="justificado")
    W.figura(figs["clima"], f"Clima real del horizonte: 2025 frente a la variabilidad de las {X['n_temporadas']} "
                            "temporadas 1996–2025.", ref="clima")
    W.p("**Bordes** (normal saliente). Caras expuestas: Robin con convección por viento h~c~ = 5,6 + 4v (Jürges), "
        "onda larga hacia el cielo según la nubosidad y radiación solar absorbida (α = 0,65), en serie con el moldaje "
        "y la manta mientras están puestos. Hormigón–suelo y capa sobre capa: continuidad (media armónica de k). "
        "Superficie del suelo: el mismo balance con el clima. Fondo del suelo a 3 m: Dirichlet con la temperatura "
        "profunda de ERA5; costados: flujo nulo. Al colocar una capa, la cara superior de la anterior deja de ser borde.",
        alinear="justificado")

    # ---------------- 6. Modelo numérico
    W.h1("6. Modelo numérico y verificación")
    W.ecuacion(X["eq1"])
    W.p("Solver propio en Python: volúmenes finitos sobre vóxeles cúbicos con Euler explícito y "
        "Δt ≤ min(0,9 Δt~crít~, 600 s), recalculado en cada evento (capa colocada, retiro de una protección). Las celdas "
        "aparecen al colocarse, a su temperatura de colocación y con edad cero. Los tubos son sumideros lineales con la "
        "resistencia agua–pared–hormigón y el calentamiento del agua a lo largo del circuito. Verificación:",
        alinear="justificado")
    vin = [f"{X['n_tests']} tests automáticos: soluciones exactas 1D, 2D y 3D (orden 2 en espacio y 1 en tiempo), caja "
           f"con moldaje y desmolde, bloque adiabático igual a la curva del ensayo, voxelización y clima;",
           "balance de energía en cada simulación: residuo relativo del orden de 10^−15^ (columna residuo_energia de "
           "cada escenario);"]
    sens, d = X["sens"], X["delta"]
    prop_s = [x for x in sens if x["plan"] == PROPUESTA]
    if prop_s:
        pr = max(prop_s, key=lambda x: x["MD_fina"])
        otras_clima = [x for x in prop_s if x is not pr]
        otras_num = [r for x in prop_s for r in x["otras"]]
        dif = max([abs(float(r["M_D"]) - pr["MD_gruesa"]) for r in otras_num], default=0.0)
        txt = (f"malla: refinar de Δx = 0,275 a {f1(pr['dx_fina'], 2)} m sube M~Δ~ del peor escenario "
               f"({pr['curva']}, {pr['clima'].replace('_', ' ')}) de {f1(pr['MD_gruesa'], 2)} a {f1(pr['MD_fina'], 2)} °C")
        if otras_clima:
            txt += " (" + "; ".join(f"{x['clima'].replace('_', ' ')}: +{f1(x['delta'], 2)} °C" for x in otras_clima) + ")"
        txt += (f", mientras M~T~ cambia {f1(max(abs(x['MT_fina'] - x['MT_gruesa']) for x in prop_s), 1)} °C o menos"
                + (f"; Δt = 300 s y medir T~mín~ a 50 mm lo cambian en {f1(dif, 2)} °C o menos" if otras_num else "")
                + f". La malla gruesa subestima el enfriamiento junto a la superficie, así que **todas las M~Δ~ del "
                f"estudio se corrigen sumando {f1(d, 2)} °C** (la mayor diferencia observada) antes de compararlas con "
                f"20 °C. La verificación con Δx = 0,055 m queda para el cluster.")
        vin.append(txt)
    W.vinetas(vin)

    # ---------------- 7. Plan de vertidas
    W.h1("7. Procedimiento de vertidas y protecciones")
    W.h2("7.1 Secuencia y horario")
    filas = [[f"V{v}", f"{f1(geo[v - 1]['z0'], 2)}–{f1(geo[v - 1]['z1'], 2)}", f1(geo[v - 1]["V"]),
              f"{f1(inicios[v - 1])}–{f1(fines[v - 1])}", f"{hora(t0, inicios[v - 1])} → {hora(t0, fines[v - 1])}"]
             for v in range(1, 6)]
    filas.append(["Total", "0–4,30", f1(sum(g["V"] for g in geo)), f"0–{f1(D)}", f"{hora(t0, 0)} → {hora(t0, D)}"])
    W.tabla("Colocación continua por capas a 27 m³/h (una bomba y una cuadrilla, sin juntas frías).",
            ["Vertida", "Cota z [m]", "Volumen [m³]", "Horas desde t~0~", f"Reloj ({MESES[t0.month - 1][:3]} {t0.year})"],
            filas, [2.0, 2.6, 2.6, 3.2, 5.6], destacar=(5,), ref="secuencia")
    W.figura(figs["crono"], "Cronograma: colocación (izquierda, en horas) y protecciones, tubos y monitoreo "
                            "(derecha, en días).", ref="crono")
    W.p(f"El plan queda también en el modelo BIM: geometry/{X['ifc'].name} es el IFC del curso con el property set "
        f"UANDES_HormigonMasivo de cada vertida actualizado (Inicio_h de la tabla, TemperaturaColocacion_C = 10, "
        f"Moldaje, Aislacion y RetiroMoldaje_h hasta el día {dia_ret}); se revisa con mcplan ifc-info.",
        alinear="justificado")
    W.vinetas([
        "La colocación es continua, de día y de noche: dos turnos, iluminación y una **bomba de respaldo en obra**. "
        f"Una interrupción larga equivale a las variantes con esperas, que no cumplen (Tabla {ref('t:familia')}).",
        "Capas de ≤ 0,5 m, vibradas; la superficie final se termina y se cubre con la manta apenas se puede caminar "
        "sobre ella.",
        f"El hormigón se recibe sólo con **T ≤ {tcol} a la descarga**; se mide en cada camión, junto con el cono.",
    ])
    W.h2("7.2 Moldaje y aislación")
    W.p(f"Moldaje de terciado de 18 mm en las caras verticales de V1, V2 y V5 (V3 y V4 tienen caras inclinadas y no lo "
        f"necesitan). La **manta de 10 cm** cubre todas las caras expuestas: la cara superior de la losa y de los "
        f"troncos a medida que quedan terminadas, y las caras moldeadas por fuera del moldaje. Bajo la manta va un film "
        f"de polietileno, que retiene la humedad (curado) y evita que la manta se moje, porque mojada pierde gran parte "
        f"de su aislación; encima, lastre contra el viento. Moldaje y mantas se retiran **desde el día {dia_ret}**, en "
        f"forma gradual y según lo que midan los sensores (sección 10.2).", alinear="justificado")
    W.h2("7.3 Tubos de enfriamiento")
    tubos = X["tubos"]
    texto_tubos = ""
    if tubos:
        largos = ", ".join(f"V{n.split('V')[-1]}: {f1(s['largo_m'], 0)} m" for n, s in tubos.items())
        Pmax = sum(s["potencia_max_kW"] for s in tubos.values())
        Tsal = max((s["T_salida_max_C"] or 0) for s in tubos.values())
        E = sum(s["energia_MJ"] for s in tubos.values())
        texto_tubos = (f" En el peor escenario los circuitos ({largos}) extraen en total {f1(E, 0)} MJ, con una potencia "
                       f"máxima conjunta de {f1(Pmax)} kW y agua de retorno de hasta {f1(Tsal)} °C: ese es el tamaño del "
                       f"disipador.")
    W.p("Un serpentín de HDPE Ø32 mm a media altura de cada una de **V3, V4 y V5**, con filas cada 1,0 m y 20 L/min "
        "por circuito, desde que empieza cada vertida y durante 14 días. El agua entra a **25 ± 2 °C**: al comienzo "
        "entibia el hormigón colocado a 10 °C y, cuando el núcleo supera los 25 °C, extrae calor. Se necesita un "
        f"circuito cerrado con estanque, bomba, calefactor y un disipador (aerotermo).{texto_tubos} No se ponen tubos "
        f"en V1 ni V2: enfriarían el anillo exterior de la losa, que ya es la zona más fría (Tabla {ref('t:familia')}). "
        "Al terminar, los tubos se inyectan con lechada.", alinear="justificado")

    # ---------------- 8. Temperaturas previstas
    W.h1("8. Temperaturas previstas")
    if C:
        pe = C["peor"]
        W.p(f"La Figura {ref('f:temp')} muestra el peor escenario ({pe['curva']}, temporada {pe['clima'].split('_')[-1]}): "
            f"el núcleo llega a {f1(pe['M_T'])} °C a las {f1(pe['t_MT_h'], 0)} h y el diferencial máximo, "
            f"{f1(pe['M_D_gruesa'])} °C con la malla del estudio ({f1(pe['M_D'])} °C con la corrección de malla), "
            f"ocurre a las {f1(pe['t_MD_h'], 0)} h, con el punto frío en "
            f"({f1(pe['x_Tmin'], 2)}; {f1(pe['y_Tmin'], 2)}; {f1(pe['z_Tmin'], 2)}) m, en {_zona(pe['z_Tmin'])}. "
            f"El diferencial sube durante los primeros días y después se estabiliza o baja; en ninguna de las "
            f"{len(C['anios'])} temporadas llega a 20 °C. En los cortes (Figura {ref('f:corte')}) los serpentines se "
            f"ven como franjas más frías a media altura de V3, V4 y V5: el punto frío aparece donde el serpentín de V4 "
            f"llega cerca de la cara inclinada, porque allí el agua a 25 °C enfría más de lo que la manta protege. Por "
            f"eso el clima influye poco y por eso los sensores S1 y S2 van en esa cara.", alinear="justificado")
        if "temp" in figs:
            W.figura(figs["temp"], "Izquierda: núcleo, superficie (25 mm) y aire en el peor escenario. Derecha: ΔT(t) "
                                   "con la misma curva en las 30 temporadas reales. Malla del estudio (0,275 m): con "
                                   f"malla fina el máximo sube {f1(X['delta'], 1)} °C.", ref="temp")
        if X["corte"] is not None:
            W.figura(X["corte"], "Cortes por el eje en cuatro instantes (fin de la colocación, +1, +3 y +6 días), "
                                 "en el peor escenario.", ancho_cm=15.0, ref="corte")

    # ---------------- 9. Confiabilidad
    W.h1("9. Análisis de confiabilidad")
    W.p("**Modelo de incertidumbre.** *Material:* las 8 curvas, equiprobables (p~m~ = 1/8). *Clima:* las temporadas "
        "reales, equiprobables (q~k~ = 1/30), independientes del material: cada escenario pesa 1/240. *Ejecución:* la "
        f"temperatura de colocación y las fallas de equipos se tratan como variantes del plan (Tabla {ref('t:familia')}) "
        f"y como contingencias (Tabla {ref('t:contingencias')}). No se usa ninguna distribución normal arbitraria: "
        "cada escenario es una curva medida y un clima que ocurrió. Cada simulación es independiente (el estudio "
        "las reparte en procesos paralelos; en el cluster, en un arreglo de trabajos).", alinear="justificado")
    if C:
        W.figura(figs["conf"], f"M~Δ~ (con la corrección de malla, +{f1(X['delta'], 2)} °C) y M~T~ de las 240 "
                               "simulaciones de la propuesta: cada punto es una temporada; en naranja, el clima real "
                               "de 2025.", ref="conf")
        filas = [
            ["P(T~máx~ > 70 °C)", f"{int(C['k_T'])}/{C['n']}", pct(C["P_T"]), f"≤ {pct(cp_superior(int(C['k_T']), C['n']))}"],
            ["P(ΔT~máx~ > 20 °C)", f"{int(C['k_D'])}/{C['n']}", pct(C["P_D"]), f"≤ {pct(cp_superior(int(C['k_D']), C['n']))}"],
            ["P~f~ (cualquiera de los dos)", f"{int(C['k_F'])}/{C['n']}", pct(C["P_F"]), f"≤ {pct(C['cp_total'])}"],
        ]
        W.tabla("Probabilidades de la propuesta (enumeración de los 240 escenarios, M~Δ~ con la corrección de malla) "
                "y cota superior de Clopper–Pearson al 95 % si los 240 se tratan como ensayos independientes.",
                ["Evento", "Fallas", "Probabilidad", "Cota 95 %"], filas, [5.6, 2.6, 3.2, 3.6], ref="prob")
        mrg = LIM_D - C["MD_max"]
        W.p(f"**Incertidumbre estadística.** La enumeración es la probabilidad exacta del modelo discreto, pero sus "
            f"escenarios son una muestra finita. Sin usar el margen, 0 fallas en 240 sólo acotan P~f~ ≤ "
            f"{pct(C['cp_total'])} (Clopper–Pearson, 95 %), y la cota sube a {pct(C['cp_anios'])} si se reconocen sólo "
            f"las 30 temporadas como independientes, o a {pct(C['cp_curvas'])} con sólo las 8 curvas. Lo que acota de "
            f"verdad es el **margen**, y hay que mirarlo por fuente. Con la manta, el **clima mueve poco** el "
            f"diferencial: con la malla del estudio, las 30 temporadas cambian M~Δ~ en ±"
            f"{f1(float(np.nanmax(np.nanstd(X['C_gruesa']['MD'], axis=1))), 2)} °C por curva (desviación típica), y con malla "
            f"fina la temporada fría queda {f1(_brecha_clima(X), 1)} °C sobre las demás verificadas; por eso la corrección "
            f"de malla usa el valor de la temporada fría para todas. Lo que manda es **el material**: M~Δ~ va de "
            f"{f1(np.nanmin(C['MD']))} a {f1(C['MD_max'])} °C según la curva.", alinear="justificado")
        mat = X.get("mat")
        if mat:
            W.p(f"Ajustando M~Δ~ (corregido por malla) contra la energía liberada a 14 días por las curvas SS, con y sin "
                f"el +20 % de la prueba de estrés (recta con R² = {f1(mat['r2'], 3)}), el diferencial llegaría a 20 °C con "
                f"**E~crít~ = {f1(mat['E_crit'], 0)} MJ/m³**, un {pct(mat['exceso'], 0)} más que la probeta más exotérmica "
                f"({mat['curva_max']}, {f1(mat['E_max'], 0)} MJ/m³). Con la energía de las 8 probetas (media "
                f"{f1(mat['mu'], 0)} MJ/m³, desviación {f1(mat['sd'], 0)} MJ/m³, ajustada a una normal sólo para "
                f"este cálculo), P(E > E~crít~) ≈ **{pct(mat['p'], 1)}**, con cota 95 % de {pct(mat['p_95'], 1)} al "
                f"remuestrear las probetas (bootstrap de 4000 réplicas); con sólo las 4 SS, {pct(mat['p_ss'], 1)}. "
                f"Ese es el valor con el hormigón justo en el máximo admisible de 10 °C.", alinear="justificado")
            m8 = X.get("mat8") or {}
            ref_m = m8 or mat
            texto8 = (f"A la meta de la planta, 8 °C, la recta se desplaza: E~crít~ sube a **{f1(m8['E_crit'], 0)} MJ/m³** "
                      f"({pct(m8['exceso'], 0)} sobre la probeta más exotérmica) y P(E > E~crít~) baja a "
                      f"**{pct(m8['p'], 1)}** (cota 95 %: {pct(m8['p_95'], 1)}; sólo SS: {pct(m8['p_ss'], 1)}). "
                      if m8 else "")
            W.p(texto8 + "Es la estimación más honesta de P~f~ con la información disponible: "
                + ("queda bajo la meta de 5 % incluso en su cota." if ref_m["p_95"] <= 0.05 else
                   "su valor central queda bajo la meta de 5 %, pero su cota no; el monitoreo, con umbrales que "
                   "detectan un hormigón más exotérmico en las primeras 48 h, y las contingencias cubren esa "
                   "diferencia." if ref_m["p"] <= 0.05 else
                   "supera la meta de 5 %: el plan depende del monitoreo y de las contingencias."),
                alinear="justificado")
    if est:
        W.p(f"**Prueba de estrés del material.** El encabezado de las curvas da su energía con β = 3,0 W/m²K: +20 %. "
            f"Con las 8 curvas ×1,2 (24 escenarios, 3 terciles climáticos), la propuesta llega a M~Δ~ = "
            f"{f1(est['MD_max'])} °C y M~T~ = {f1(est['MT_max'])} °C: P~f~ = {pct(est['P_F'])}"
            + (f" (a 8 °C, P~f~ = {pct(X['est8']['P_F'])})" if X.get("est8") else "") + ". "
            + ("El plan resiste también ese extremo." if est["P_F"] == 0 else
               "Ese extremo ya no se cubre por completo: es la principal razón para monitorear y tener contingencias."),
            alinear="justificado")
    W.h2("9.1 Comparación con otros planes")
    if fam:
        orden = [p for p in ETIQUETAS if p in fam]
        filas = []
        for p in orden:
            r = fam[p]
            et, Tc, manta, tub = ETIQUETAS[p]
            filas.append([et, f1(r["D"]), Tc, manta, tub, pct(r["P_D"], 0), pct(r["P_F"], 0),
                          f"{f1(r['MD_min'])}–{f1(r['MD_max'])}", "sí" if r["P_F"] <= 0.05 else "no"])
        W.tabla("Familia de planes con 3 historias climáticas (terciles de las 30 temporadas, peso 1/3 cada uno) y las "
                f"8 curvas: 24 escenarios por plan, con M~Δ~ corregido por malla (+{f1(X['delta'], 2)} °C). P~T~ = 0 "
                f"en todos (M~T~ ≤ {f1(X['MT_familia'])} °C).",
                ["Plan", "D [h]", "T col. [°C]", "Manta", "Tubos", "P~Δ~", "P~f~", "M~Δ~ [°C]", "¿Cumple?"], filas,
                [5.4, 1.2, 1.3, 1.3, 2.6, 1.2, 1.2, 1.9, 1.3], alinear=["l"] + ["c"] * 8,
                destacar=tuple(i for i, p in enumerate(orden) if p == PROPUESTA), tam=8, ref="familia")
        if prop:
            temps = "; ".join(f"a {ETIQUETAS[p][1]} °C, P~f~ = {pct(fam[p]['P_F'], 0)} (M~Δ~ hasta {f1(fam[p]['MD_max'])} °C)"
                              for p in ("unica-M10-T12-N25", "unica-M10-T14-N25", "unica-M10-T16-N25") if p in fam)
            t8 = fam.get("unica-M10-T8-N25")
            lentas = [p for p in fam if p.startswith("5V")]
            cumplen = [q for q in fam if fam[q]["P_F"] <= 0.05]
            mejor = min(cumplen, key=lambda q: fam[q]["MD_max"]) if cumplen else None
            margen = ("la de mayor margen entre las que cumplen" if mejor == PROPUESTA else
                      f"cumple con M~Δ~ hasta {f1(prop['MD_max'])} °C")
            W.p(f"La propuesta es la **más rápida posible** (D = {f1(prop['D'])} h: con una bomba, todo plan sin esperas "
                f"es esta misma colocación) y {margen}. "
                + (f"Las **variantes más lentas** (5 vertidas con esperas de 12 o 24 h) tienen P~f~ de "
                   f"{' y '.join(pct(fam[p]['P_F'], 0) for p in lentas)}: mientras se espera, las capas inferiores se "
                   f"enfrían hacia el suelo y la siguiente se calienta encima. " if lentas else "")
                + (f"Colocar más caliente reduce el margen: {temps}"
                   + next((f" (con malla fina, colocando a 12 °C el peor escenario llega a {f1(x['MD_fina'], 2)} °C)"
                           for x in X["sens"] if x["plan"] == "unica-M10-T12-N25"), "")
                   + f". Por eso la **temperatura máxima de colocación es {tcol}**, y la planta apunta a 8 °C para "
                   f"cumplirla con holgura"
                   + (f" (a 8 °C, M~Δ~ queda en {f1(t8['MD_min'])}–{f1(t8['MD_max'])} °C, con "
                      f"{f1(LIM_D - t8['MD_max'])} °C de margen)" if t8 else "")
                   + "; el análisis de confiabilidad de la sección 9 se hizo con el hormigón justo en el máximo."
                   if temps else ""), alinear="justificado")

    # ---------------- 10. Monitoreo, retiro y contingencias
    W.h1("10. Monitoreo, retiro de protecciones y contingencias")
    W.h2("10.1 Sensores y registro")
    if sensores:
        W.p(f"Los sensores van donde el modelo predice el punto más caliente y los puntos fríos que definen el ΔT "
            f"(Figura {ref('f:sens')}), con redundancia en los de más peso, más el aire, el suelo y el agua de los "
            f"tubos. Se usan termocuplas tipo T o sensores digitales (±0,5 °C) fijados a barras antes de hormigonar; "
            f"los de 25 mm se montan con un separador a esa distancia de la cara.", alinear="justificado")
        if "sens" in figs:
            W.figura(figs["sens"], "Sensores sobre la sección por el eje y la planta, junto a los puntos más calientes "
                                   "(naranja) y más fríos a 25 mm (azul) de los 240 escenarios.", ref="sens")
        filas = [[s["id"], s["tipo"], str(s["n"]), f"({f1(s['x'], 2)}; {f1(s['y'], 2)}; {f1(s['z'], 2)})", s["donde"]]
                 for s in sensores]
        filas.append(["W1–W6", "agua", "6", "entrada y salida de V3–V5", "temperatura del agua y caudal de cada circuito"])
        W.tabla("Sensores: identificador, tipo, cantidad (redundancia), posición (x; y; z) en metros y propósito.",
                ["ID", "Tipo", "N.º", "Posición [m]", "Propósito"], filas, [1.3, 1.5, 0.9, 3.8, 9.7],
                alinear=["c", "l", "c", "c", "l"], tam=8, ref="sensores")
    W.p("**Frecuencia de registro:** cada 15 min desde la primera vertida hasta las 72 h (cubre la colocación y la "
        "subida inicial; el ICH pide 30 min o menos en las primeras 48 h), cada 30 min hasta el día 14 y cada hora hasta "
        "3 días después del retiro de las protecciones. Registro continuo con respaldo de energía, transmisión remota "
        "con alarmas y comparación con el modelo dos veces al día.", alinear="justificado")
    W.h2("10.2 Desmolde y retiro de las protecciones")
    hechos = [(p, d, h) for p, d, h in RETIROS if p in ret]
    if hechos:
        filas = [[f"día {d}", pct(ret[p]["P_D"], 0), f1(ret[p]["MD14_max"]), f1(ret[p]["MD_max"]),
                  "sí" if ret[p]["P_D"] == 0 else "no"] for p, d, _ in hechos]
        W.tabla("Retiro del moldaje y las mantas en distintos días (24 escenarios por fecha, simulados 30 días; 37 "
                f"para el día 33), con M~Δ~ corregido por malla (+{f1(X['delta'], 2)} °C).",
                ["Retiro", "P(ΔT > 20 °C)", "M~Δ~ en 14 días [°C]", "M~Δ~ hasta el final [°C]", "¿Admisible?"], filas,
                [3.0, 3.2, 3.4, 3.6, 2.6], ref="retiro")
        if "retiro" in figs:
            W.figura(figs["retiro"], "ΔT(t) en el peor escenario de cada fecha de retiro: al sacar las mantas, la "
                                     "superficie se enfría y el diferencial salta.", ref="retiro")
    cr = X["crit_retiro"]
    if cr and math.isfinite(cr.get("dT_max", float("nan"))):
        cond = (f"ΔT medido (núcleo − sensor más frío) ≤ {math.floor(cr['dT_max'] + X['delta'])} °C y núcleo − mínima del aire "
                f"prevista para los 3 días siguientes ≤ {math.floor(cr['nucleo_aire_max'])} °C")
    else:
        cond = "ΔT medido y diferencia núcleo–aire dentro de lo que el modelo muestra seguro"
    aviso = ("" if X.get("retiro_ok") else
             f"Con la corrección de malla, ningún retiro de una sola vez cumple en todos los escenarios: bajo la manta el "
             f"núcleo sigue caliente por semanas, así que retirar más tarde no basta"
             + (" (el día 33 falla en más escenarios que el 26)" if "propuesta-retiro-d33" in ret and "propuesta-retiro-d26" in ret
                and ret["propuesta-retiro-d33"]["P_D"] > ret["propuesta-retiro-d26"]["P_D"] else "")
             + f". El día {dia_ret} es el de menos fallas ({pct(ret[plan_ret]['P_D'], 1)}), y el retiro debe "
             f"ser **gradual y condicionado a lo que midan los sensores**. ")
    W.p(f"**Criterio:** {aviso}Moldaje y mantas se retiran el **día {dia_ret}** como muy pronto, de día y sin pronóstico "
        f"de helada, y sólo si los sensores confirman las condiciones que el modelo tenía en esa fecha en los escenarios "
        f"que después cumplen: {cond}. Si no se cumplen, se espera y se evalúa de nuevo cada 2 días. Se retira por "
        f"etapas: primero se reemplaza la manta de 10 cm por una de 5 cm (o se retiran sólo las caras laterales) y, "
        f"cuando el ΔT medido vuelve a estabilizarse, el resto; el retiro gradual se debe verificar con el modelo antes "
        f"de la obra. El registro horario sigue 3 días después del último retiro.",
        alinear="justificado")
    W.h2("10.3 Umbrales de alerta y acciones de contingencia")
    mt_prev = C["MT_max"] if C else float("nan")
    filas = [
        [f"Hormigón a la descarga > {tcol}", "Se rechaza el camión o se ajusta el hielo en planta.",
         f"colocar más caliente baja el margen (Tabla {ref('t:familia')})"],
        ["Falla de la bomba de hormigón o de la planta", "Bomba de respaldo en obra; la interrupción no debe pasar de 2 h. "
         "Si pasa, se trata la junta y se reevalúa el plan con el modelo.", "las esperas de 12–24 h no cumplen"],
        ["Falla del circuito de agua (bomba o calefactor)", "Bomba y generador de respaldo; reparar en menos de 12 h. "
         "Sin calefactor el agua llega fría: se reduce el caudal hasta repararlo.", "sin tubos, o con agua a 20 °C, no cumple"],
        [f"ΔT medido sobre la curva del peor escenario (Figura {ref('f:temp')}) o ≥ 18 °C (alerta)", "Revisar mantas y "
         "sellos; segunda manta en la cara del sensor más frío; recalibrar el modelo con lo medido.",
         f"peor escenario previsto: {f1(C['MD_max']) if C else '—'} °C"],
        ["ΔT medido ≥ 19 °C (acción)", "Segunda manta en todas las caras del sector; subir a 27–28 °C el agua del "
         "circuito vecino al punto frío (V4) o reducir su caudal.", "el punto frío está junto al serpentín de V4"],
        ["T del núcleo ≥ 60 °C, o agua–hormigón > 20 °C", "Revisar caudal y temperatura del agua; aumentar el caudal "
         "del circuito del núcleo.", f"M~T~ previsto ≤ {f1(mt_prev)} °C"],
        ["Frente frío, lluvia o viento fuerte en el pronóstico", "Verificar lastre y film de las mantas; no retirar "
         "protecciones; adelantar las lecturas.", "las temporadas frías ya están en el análisis"],
        ["Sensor dañado o lecturas fuera de rango", "Usar el redundante; reemplazar si es accesible.", "redundancia"],
        ["T medidas a 24–48 h difieren > 3 °C del modelo", "Recalibrar la curva con lo medido y reevaluar el plan (el "
         "flujo lo permite en minutos).", "modelo reproducible"],
    ]
    W.tabla("Umbrales y acciones: los umbrales dejan margen antes de los límites del proyecto.",
            ["Disparador", "Acción", "Respaldo en el modelo"], filas, [4.6, 8.4, 4.2], alinear=["l", "l", "l"], tam=8,
            ref="contingencias")

    # ---------------- 11. Limitaciones y 12. reproducibilidad
    W.h1("11. Supuestos y limitaciones")
    W.vinetas([
        "La fuente de calor depende sólo de la edad (sin madurez): no se acelera en el núcleo caliente ni se frena con "
        "frío. El ρc~p~ del modelo (2,33 MJ/m³K) es menor que el de las probetas (2,53): el aumento adiabático sale ≈9 % "
        "mayor, del lado conservador.",
        "Las curvas son de probetas con 351 kg/m³ de cemento; si la dosificación de la obra es otra, hay que reescalarlas.",
        "Las propiedades del hormigón y del suelo son las de Ju y Lei (2019), no medidas en el sitio. El suelo bajo la "
        "losa es el que más enfría la base: conviene medir su temperatura antes de hormigonar.",
        "El ΔT es un indicador de riesgo, no un análisis de tensiones.",
        "La corrección de malla (+" + f1(X["delta"], 1) + " °C) se midió en el pico de los primeros días; después del "
        "retiro de la manta el gradiente junto a la superficie es más fuerte y la malla fina podría dar más. El retiro "
        "gradual no se simuló: queda para la verificación en el cluster, junto con Δx = 0,055 m.",
        "El margen frente a un hormigón más exotérmico que las probetas es estrecho"
        + (f" (≈ {pct(X['mat8']['exceso'], 0)} de energía a 8 °C)" if X.get("mat8") else "")
        + ": comparar la temperatura medida del núcleo con la prevista en las primeras 48 h es la forma de detectarlo "
        "a tiempo y activar las contingencias.",
        "La fecha usa el clima real de 2025; para la fecha final se repite el estudio con el ensemble más reciente "
        "(config/estudio.yaml), sin cambiar el programa.",
    ])
    W.h1("12. Reproducibilidad")
    W.p("Todo se regenera desde el repositorio (Python 3.12: pip install -r requirements.txt; pytest):", despues=2)
    for linea in ("mcplan clima era5 --fecha-inicio 2025-11-19T08:00 --anios 1996 2025 --dias 40 --margen 7",
                  "mcplan estudio config/estudio_2025.yaml            # familia: 17 planes x 24 escenarios",
                  "mcplan estudio config/confiabilidad_2025.yaml      # propuesta: 8 curvas x 30 temporadas",
                  "mcplan estudio config/estres_2025.yaml             # +20 % de calor",
                  "mcplan estudio config/retiro_2025.yaml             # retiro de protecciones (días 14, 21, 26)",
                  "mcplan estudio config/retiro_tardio_2025.yaml      # retiro el día 33",
                  "mcplan estudio config/confiabilidad_2025.yaml --sensibilidad unica-M10-T10-N25 M4_SS ERA5_2013",
                  "mcplan estudio config/confiabilidad_2025.yaml --detalle unica-M10-T10-N25 M4_SS ERA5_2013",
                  "python scripts/informe_final.py --pdf              # este informe"):
        par = W.p(linea, tam=8, despues=0)
        for r in par.runs:
            r.font.name = "Consolas"
    me = X["meta_era5"]
    W.p(f"Clima: {me.get('fuente', 'ERA5')}, descargado el {me.get('descargado_utc', '—')} (snapshot inmutable con "
        f"SHA-256, en data/clima/). Curvas: data/hidratacion/source_terms.txt, sin modificar. Cada escenario es una "
        f"fila de escenarios.csv (en la carpeta de cada estudio, bajo resultados/) con su peso, métricas, ubicación de "
        f"los extremos y residuo de energía. La corrida de detalle (--detalle) deja además los campos 3D en VTK para "
        f"ParaView.", tam=9.5)
    W.h1("Referencias")
    W.vinetas([
        "Ju, Y. y Lei, H. (2019). Actual temperature evolution of thick raft concrete foundations and cracking risk "
        "analysis. *Advances in Materials Science and Engineering*, 7029671. doi:10.1155/2019/7029671.",
        "Instituto del Cemento y del Hormigón de Chile (2019). *Especificación para Hormigón Masivo Estructural*. Santiago.",
        "ACI Committee 301 (2020). *Specifications for Concrete Construction* (ACI 301-20).",
        "ACI Committee 207 (2005). *Cooling and Insulating Systems for Mass Concrete* (ACI 207.4R-05).",
        "Hersbach, H. et al. (2020). The ERA5 global reanalysis. *Quarterly Journal of the Royal Meteorological "
        "Society*, 146, 1999–2049. Datos vía Open-Meteo Historical Weather API.",
        "Curvas de hidratación Q(t) de las probetas M1–M4, campañas SS y 22A (profesor del curso, versión del 14 nov 2026).",
    ])


def ifc_del_plan(X: dict) -> Path:
    """Copia del IFC oficial con el plan final en el pset UANDES_HormigonMasivo de cada vertida."""
    from mcplan.ifc_io import agregar_propiedades_curso

    _, dia, h = X["retiro"]
    retiro_abs = h + X["D"]  # el retiro se cuenta desde el fin de la colocación
    valores = {}
    for v in range(1, 6):
        ini = X["inicios"][v - 1]
        valores[f"Vertida_0{v}"] = {
            "VertidaID": v, "Inicio_h": round(ini, 2), "TemperaturaColocacion_C": 10.0, "TipoHormigon": "M1",
            "Moldaje": "Terciado fenólico 18 mm" if v in (1, 2, 5) else "Sin moldaje (cara inclinada)",
            "Aislacion": "Manta 10 cm (R = 2,5 m2K/W) con film, todas las caras expuestas",
            "RetiroMoldaje_h": round(retiro_abs - ini, 1)}
    salida = RAIZ / "geometry" / "foundation_plan_final.ifc"
    agregar_propiedades_curso(RAIZ / "geometry" / "foundation.ifc", salida, valores)
    return salida


def armar(args) -> Path:
    X = preparar()
    X["ifc"] = ifc_del_plan(X)
    print(f"IFC del plan: {X['ifc']}")
    print("Documento...")
    W = Informe()
    escribir(W, X, {})  # primera pasada: numera figuras y tablas
    W2 = Informe()
    escribir(W2, X, W.refs)  # segunda: las referencias del texto ya apuntan al número correcto
    SALIDA.mkdir(parents=True, exist_ok=True)
    W2.guardar(DOCX)
    print(f"Informe: {DOCX}")
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
    args = ap.parse_args()
    docx = armar(args)
    if args.pdf:
        pdf, n = exportar_pdf(docx)
        print(f"PDF: {pdf} ({n} páginas)")


if __name__ == "__main__":
    main()
