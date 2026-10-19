"""Gráficos 2D con Matplotlib (curvas, perfiles y mapa z-t)."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # sólo archivos PNG, no abre ventanas
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from .metricas import MetricasTermicas  # noqa: E402
from .solver1d import Resultado  # noqa: E402


def _limites_hormigon(res: Resultado):
    idx = np.flatnonzero(res.malla.es_hormigon)
    return res.malla.z_caras[idx[0]], res.malla.z_caras[idx[-1] + 1]


def graficar_series(res: Resultado, met: MetricasTermicas, ruta):
    t = res.t_h
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    ax1.plot(t, met.Tmax_t, label="Tmax en el hormigón (núcleo)", color="C3")
    ax1.plot(t, met.T_25mm_sup_t, label="25 mm bajo superficie superior", color="C0")
    ax1.plot(t, met.T_25mm_inf_t, label="25 mm sobre base (contacto suelo)", color="C2")
    if np.isfinite(res.T_exterior_sup).any():
        ax1.plot(t, res.T_exterior_sup, label="Ambiente", color="0.6", lw=0.8)
    ax1.axhline(met.limite_Tmax, color="C3", ls="--", lw=1, label=f"Límite {met.limite_Tmax:.0f} °C")
    ax1.set_ylabel("Temperatura [°C]")
    ax1.legend(fontsize=8, loc="upper right")
    ax1.grid(alpha=0.3)

    ax2.plot(t, met.dT_t, color="C1", label="ΔT = Tmax − Tmin(25 mm)")
    ax2.axhline(met.limite_dT, color="C1", ls="--", lw=1, label=f"Límite {met.limite_dT:.0f} °C")
    ax2.plot(met.t_dTmax_h, met.dTmax, "o", color="C1")
    ax2.set_xlabel("Tiempo desde la colocación [h]")
    ax2.set_ylabel("ΔT [°C]")
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(ruta, dpi=150)
    plt.close(fig)


def graficar_perfiles(res: Resultado, ruta, dias=(1, 3, 4, 9, 14)):
    """Perfiles verticales T(z), comparables con la Fig. 7 de Ju y Lei (2019)."""
    m = res.malla
    z_base, z_tope = _limites_hormigon(res)
    altura = np.concatenate(([m.z_caras[0]], m.z, [m.z_caras[-1]])) - z_base
    fig, ax = plt.subplots(figsize=(6, 7))
    ax.axhspan(0.0, z_tope - z_base, color="0.9", label="Hormigón")
    for d in dias:
        i = int(np.argmin(np.abs(res.t_h - 24 * d)))
        T_ext = np.concatenate(([res.T_cara_inf[i]], res.T[i], [res.T_cara_sup[i]]))
        ax.plot(T_ext, altura, label=f"{res.t_h[i] / 24:.1f} d")
    ax.set_xlabel("Temperatura [°C]")
    ax.set_ylabel("Altura sobre la base del hormigón [m]  (negativo = suelo)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(ruta, dpi=150)
    plt.close(fig)


def graficar_mapa(res: Resultado, ruta):
    fig, ax = plt.subplots(figsize=(9, 5))
    z_base, z_tope = _limites_hormigon(res)
    malla = ax.pcolormesh(res.t_h, res.malla.z - z_base, res.T.T, shading="nearest", cmap="inferno")
    for z in (0.0, z_tope - z_base):
        ax.axhline(z, color="w", lw=0.8, ls="--")
    fig.colorbar(malla, ax=ax, label="Temperatura [°C]")
    ax.set_xlabel("Tiempo [h]")
    ax.set_ylabel("Altura sobre la base del hormigón [m]")
    ax.set_title("Campo de temperatura T(z, t)")
    fig.tight_layout()
    fig.savefig(ruta, dpi=150)
    plt.close(fig)


def graficar_todo(res: Resultado, met: MetricasTermicas, carpeta) -> list[Path]:
    carpeta = Path(carpeta)
    carpeta.mkdir(parents=True, exist_ok=True)
    rutas = [carpeta / "series.png", carpeta / "perfiles.png", carpeta / "mapa_zt.png"]
    graficar_series(res, met, rutas[0])
    graficar_perfiles(res, rutas[1])
    graficar_mapa(res, rutas[2])
    return rutas
