"""Figuras y tablas del estudio de decisión (para la presentación y el artículo)."""

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

COLORES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]


def fig_curvas(fuentes, rho_cp, ruta):
    """Las 8 fuentes Q(a) y su energía acumulada como aumento adiabático equivalente E/(ρ c_p)."""
    from .hidratacion import HidratacionJuLei2019

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.8))
    probetas = sorted({f.info["probeta"] for f in fuentes.values()})
    for n, f in fuentes.items():
        color = COLORES[probetas.index(f.info["probeta"]) % 8]
        ls = "-" if f.info["campana"] == "SS" else "--"
        a = f.edad_h
        m = a <= 336
        a1.plot(a[m], f.Q_tabla[m], ls, color=color, lw=1.4, label=n)
        a2.plot(a[m], f.energia(a[m]) / rho_cp, ls, color=color, lw=1.4)
        fin = f.info["fin_registro_h"]
        a1.axvline(fin, color=color, lw=0.5, ls=":")
    jl = HidratacionJuLei2019()
    a = np.linspace(0, 336, 600)
    a1.plot(a, rho_cp * jl.dtheta_dte(a) / 3600, color="0.3", lw=2, alpha=0.6, label="Ju y Lei (2019), te = t")
    a2.plot(a, jl.theta(a), color="0.3", lw=2, alpha=0.6, label="Ju y Lei (2019)")
    a1.set_xlabel("edad a [h]")
    a1.set_ylabel("Q [W/m³]")
    a1.set_title("Fuente de calor de cada curva (línea: SS; trazos: 22A; punteado: fin del registro)", fontsize=9.5)
    a1.legend(fontsize=7.5, ncol=2)
    a2.set_xlabel("edad a [h]")
    a2.set_ylabel("E(a) / (ρ c_p) [°C]")
    a2.set_title(f"Calor liberado como aumento adiabático equivalente (ρ c_p = {rho_cp / 1e6:.2f} MJ/m³K)",
                 fontsize=9.5)
    a2.legend(fontsize=8)
    for ax in (a1, a2):
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(ruta, dpi=150)
    plt.close(fig)


def fig_clima(todas, elegidas, ruta, dias=14):
    """Todas las historias candidatas (gris) y las representantes, con su peso."""
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(12, 7), gridspec_kw={"height_ratios": [2, 1]})
    n = dias * 24 + 1
    for _, c in todas:
        a1.plot(c.t_h[:n], c.T_aire[:n], color="0.75", lw=0.5)
    colores = {"fría": "#2a78d6", "media": "#1baf7a", "cálida": "#e34948"}
    for h in elegidas:
        c = h["clima"]
        a1.plot(c.t_h[:n], c.T_aire[:n], lw=1.6, color=colores.get(h["etiqueta"], "k"),
                label=f"{h['nombre']}: peso {h['peso']:.3f}")
    a1.set_ylabel("T aire [°C]")
    a1.set_xlabel("t desde t0 [h]")
    a1.legend(fontsize=8)
    a1.grid(alpha=0.3)
    a1.set_title(f"{len(todas)} historias candidatas (gris) y las {len(elegidas)} representantes", fontsize=10)
    from .clima import indice_frio

    idx = sorted((indice_frio(c), nombre) for nombre, c in todas)
    grupo = {m: h["etiqueta"] for h in elegidas for m in h["miembros"]}
    a2.bar(range(len(idx)), [v for v, _ in idx], color=[colores.get(grupo[nm], "0.5") for _, nm in idx])
    a2.set_ylim(min(v for v, _ in idx) - 0.5, max(v for v, _ in idx) + 0.5)
    a2.set_ylabel("índice de frío [°C]")
    a2.set_xlabel("historias ordenadas (media de las mínimas diarias de los primeros 7 días)")
    a2.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(ruta, dpi=150)
    plt.close(fig)


def fig_D_Pf(resumen, ruta, dec=None):
    fig, ax = plt.subplots(figsize=(8, 5))
    for r in resumen:
        color = "#1b7f4b" if r["cumple"] else "#b3261e"
        ax.plot(r["D_h"], 100 * r["P_f"], "o", ms=8, color=color, mfc=color if r["cumple"] else "white")
        ax.annotate(r["plan"], (r["D_h"], 100 * r["P_f"]), xytext=(5, 5), textcoords="offset points", fontsize=8)
    ax.axhline(5, color="0.3", ls="--", lw=1, label="P_f = 5 %")
    ax.set_xlabel("D: duración de la colocación [h desde t0]")
    ax.set_ylabel("P_f conjunta [%]")
    ax.set_ylim(-3, 103)
    titulo = "Duración vs. probabilidad de falla conjunta (lleno: cumple)"
    if dec and dec.get("propuesto"):
        titulo += f" · propuesto: {dec['propuesto']}"
    elif dec:
        titulo += " · ningún plan evaluado cumple"
    ax.set_title(titulo, fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(ruta, dpi=150)
    plt.close(fig)


def fig_cronogramas(cronos, ruta, planes=None):
    planes = planes or list(cronos)
    fig, ax = plt.subplots(figsize=(11, 0.55 * len(planes) + 1.5))
    for y, p in enumerate(planes):
        for k, e in enumerate(cronos[p]["etapas"]):
            ax.barh(y, e["duracion_h"], left=e["inicio_h"], color=COLORES[k % 8], edgecolor="white", height=0.6)
        ax.text(cronos[p]["D_h"] + 2, y, f"D = {cronos[p]['D_h']:.1f} h", va="center", fontsize=8)
    ax.set_yticks(range(len(planes)), planes)
    ax.invert_yaxis()
    ax.set_xlabel("horas desde t0 (cada color es una etapa; los huecos son esperas)")
    ax.set_title("Cronogramas de colocación", fontsize=10)
    ax.grid(alpha=0.3, axis="x")
    fig.tight_layout()
    fig.savefig(ruta, dpi=150)
    plt.close(fig)


def fig_escenarios(filas, planes, ruta):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 0.45 * len(planes) + 2), sharey=True)
    climas = sorted({f["clima"] for f in filas})
    marcas = dict(zip(climas, "o^s"))
    for y, p in enumerate(planes):
        fs = [f for f in filas if f["plan"] == p]
        for f in fs:
            m = marcas.get(f["clima"], "o")
            a1.plot(f["M_T"], y, m, color="#c0392b", alpha=0.6, ms=5)
            a2.plot(f["M_D"], y, m, color="#2a78d6", alpha=0.6, ms=5)
    a1.axvline(70, color="0.3", ls="--")
    a2.axvline(20, color="0.3", ls="--")
    a1.set_yticks(range(len(planes)), planes)
    a1.invert_yaxis()
    a1.set_xlabel("M_T [°C] (límite 70)")
    a2.set_xlabel("M_Δ [°C] (límite 20)")
    for c, m in marcas.items():
        a2.plot([], [], m, color="0.4", label=f"clima {c}")
    a2.legend(fontsize=8)
    for ax in (a1, a2):
        ax.grid(alpha=0.3)
    fig.suptitle("Métricas por escenario (cada punto: una curva × una historia climática)", fontsize=10)
    fig.tight_layout()
    fig.savefig(ruta, dpi=150)
    plt.close(fig)


def fig_nucleo(serie, titulo, ruta, D_h=None):
    t = serie["t_h"]
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 6.5), sharex=True, gridspec_kw={"height_ratios": [1.6, 1]})
    a1.plot(t, serie["Tmax"], color="#c0392b", lw=2, label="núcleo: Tmax del hormigón colocado")
    a1.plot(t, serie["Tmin"], color="#2a78d6", lw=1.6, label="superficie: Tmin a 25 mm (caras planas)")
    a1.plot(t, serie["T_amb"], color="0.6", lw=0.8, label="T aire")
    a1.axhline(70, color="#c0392b", ls="--", lw=1, label="70 °C")
    a2.plot(t, serie["dT"], color="#b3261e", lw=1.6, label="ΔT = Tmax − Tmin")
    a2.axhline(20, color="#b3261e", ls="--", lw=1, label="20 °C")
    if D_h is not None:
        for ax in (a1, a2):
            ax.axvline(D_h, color="0.3", lw=1, ls=":")
        a1.text(D_h, a1.get_ylim()[1], " fin de la colocación", va="top", fontsize=8)
    a1.set_ylabel("T [°C]")
    a2.set_ylabel("ΔT [°C]")
    a2.set_xlabel("t desde t0 [h]")
    a1.set_title(titulo, fontsize=10)
    for ax in (a1, a2):
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(ruta, dpi=150)
    plt.close(fig)


def _pct(x):
    return f"{100 * x:.2f}"


def escribir_resumen(salida, esc, resumen, filas, dec, sensibilidades):
    """resumen.md (tablas para las diapositivas), resumen_planes.csv y definicion.json."""
    salida = Path(salida)
    with open(salida / "resumen_planes.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, list(resumen[0]))
        w.writeheader()
        w.writerows(resumen)
    definicion = {
        "t0": esc["t0"], "zona_horaria": "America/Santiago", "tasa_colocacion_m3_h": esc["tasa"],
        "curvas": {n: {**f.info, "peso": esc["p_curva"][n]} for n, f in esc["fuentes"].items()},
        "clima": {**esc["origen_clima"], "historias": [
            {k: h[k] for k in ("clave", "nombre", "representante", "peso", "indice_C", "rango_indice_C", "miembros")}
            for h in esc["historias"]]},
        "suelo": esc["info_suelo"],
        "planes": {p: c for p, c in esc["cronos"].items()},
        "decision": dec,
    }
    (salida / "definicion.json").write_text(json.dumps(definicion, indent=1, ensure_ascii=False, default=str),
                                            encoding="utf-8")
    L = [f"# Estudio de decisión — t0 = {esc['t0']} (America/Santiago)", ""]
    L += ["## Planes", "", "| Plan | Descripción | D [h] | P_T [%] | P_Δ [%] | P_f conjunta [%] | ¿Cumple? |",
          "|---|---|---|---|---|---|---|"]
    for r in sorted(resumen, key=lambda r: r["D_h"]):
        estado = "sí" if r["cumple"] else ("no" if r["completo"] else "incompleto")
        L.append(f"| {r['plan']} | {r['descripcion']} | {r['D_h']:.1f} | {_pct(r['P_T'])} | {_pct(r['P_D'])} | "
                 f"{_pct(r['P_f'])} | {estado} |")
    L += ["", "**Decisión:** " + (f"plan propuesto **{dec['propuesto']}** (menor D con P_f ≤ 5 % entre los evaluados)"
                                  if dec.get("propuesto") else
                                  "ningún plan evaluado cumple P_f ≤ 5 %; menor riesgo: " + str(dec.get("menor_riesgo"))), ""]
    L += ["## Escenarios y pesos", "", "| Curva | p | | Historia | q | representa |", "|---|---|---|---|---|---|"]
    curvas = list(esc["p_curva"].items())
    hist = esc["historias"]
    for i in range(max(len(curvas), len(hist))):
        nombre_c, peso_c = curvas[i] if i < len(curvas) else ("", None)
        h = hist[i] if i < len(hist) else None
        texto_p = "" if peso_c is None else f"{peso_c:.4f}"
        texto_h = ["", "", ""]
        if h:
            lo, hi = h["rango_indice_C"]
            texto_h = [h["nombre"], f"{h['peso']:.4f}", f"{len(h['miembros'])} historias, índice {lo:.1f}–{hi:.1f} °C"]
        L.append(f"| {nombre_c} | {texto_p} | | {texto_h[0]} | {texto_h[1]} | {texto_h[2]} |")
    L += ["", "## Sensibilidad a los pesos (mismas corridas)", "", "| Variante | " +
          " | ".join(r["plan"] for r in sorted(resumen, key=lambda r: r["D_h"])) + " |",
          "|---|" + "---|" * len(resumen)]
    for nombre, res_v in sensibilidades.items():
        por = {r["plan"]: r for r in res_v}
        L.append(f"| {nombre} | " + " | ".join(f"{_pct(por[r['plan']]['P_f'])} %" for r in
                                                sorted(resumen, key=lambda r: r["D_h"])) + " |")
    L += ["", "## Métricas por escenario", "", "| Plan | Curva | Clima | w | M_T [°C] | t [h] | M_Δ [°C] | t [h] | I_T | I_Δ | I_F |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    for f in sorted(filas, key=lambda f: (f["D_h"], f["plan"], f["curva"], f["clima"])):
        L.append(f"| {f['plan']} | {f['curva']} | {f['clima']} | {f['w']:.4f} | {f['M_T']:.1f} | {f['t_MT_h']:.0f} | "
                 f"{f['M_D']:.1f} | {f['t_MD_h']:.0f} | {int(f['I_T'])} | {int(f['I_D'])} | {int(f['I_F'])} |")
    residuo = max(f["residuo_energia"] for f in filas)
    L += ["", f"Residuo relativo máximo del balance de energía en todas las corridas: {residuo:.1e}.",
          "Picos al final del horizonte: " + (", ".join(r["plan"] for r in resumen if r["pico_al_final"]) or "ninguno") + "."]
    (salida / "resumen.md").write_text("\n".join(L) + "\n", encoding="utf-8")
