"""Gráficos con Matplotlib: columna 1D (curvas, perfiles, mapa z-t) y rectángulo 2D (campos, GIF)."""

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


def _nombre_punto(p) -> str:
    return f"({p[0]:g}, {p[1]:g})"


def _dibujar_campo(ax, res, i):
    """Mapa de T en el cuadro i (cada cuadrado es una celda) con los puntos marcados."""
    imagen = ax.imshow(res.T[i].T, origin="lower", extent=res.extension, cmap="inferno",
                       vmin=min(res.Tb, res.T.min()), vmax=res.T.max(), interpolation="nearest")
    for j, p in enumerate(res.puntos):
        ax.plot(*p, "o", color=f"C{j}", mec="w", ms=7)
        ax.annotate(_nombre_punto(p), p, xytext=(6, 6), textcoords="offset points", color="w", fontsize=8,
                    bbox={"boxstyle": "round,pad=0.2", "fc": "k", "alpha": 0.5, "lw": 0})
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    return imagen


def graficar_puntos_rectangulo(res, ruta, limite_Tmax=70.0, limite_dT=20.0):
    """T(t) en los puntos pedidos: numérica (línea) y exacta (círculos, sólo sin hidratación).

    Con hidratación agrega Tmax del campo y un segundo panel con ΔT, como la
    serie del modelo 1D.
    """
    from .rectangulo import unidad_tiempo

    f, unidad = unidad_tiempo(res.t[-1])
    t = res.t / f
    cada = max(1, len(t) // 25)
    if res.hidratacion:
        fig, (ax, ax_dT) = plt.subplots(2, 1, figsize=(9, 7.5), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
    else:
        fig, ax = plt.subplots(figsize=(9, 5))
    ax.axhspan(res.Tb - res.tolerancia, res.Tb + res.tolerancia, color="0.9", label=f"Tb ± {res.tolerancia:g} °C")
    ax.axhline(res.Tb, ls="--", color="0.4", lw=1, label=f"Borde Tb = {res.Tb:g} °C")
    for j, p in enumerate(res.puntos):
        ax.plot(t, res.T_puntos[:, j], color=f"C{j}", label=f"T{_nombre_punto(p)} numérica")
        if res.T_puntos_exacta is not None:
            ax.plot(t[::cada], res.T_puntos_exacta[::cada, j], "o", mfc="none", color=f"C{j}", ms=5,
                    label=f"T{_nombre_punto(p)} exacta")
    if res.hidratacion:
        ax.plot(t, res.Tmax_t, color="C3", lw=1, ls="-.", label="Tmax del campo")
        ax.plot(t, res.Tmin_superficie, color="C2", lw=1, ls="-.",
                label=f"Tmin a {res.profundidad * 1000:g} mm del borde")
        ax.axhline(limite_Tmax, color="C3", ls="--", lw=1, label=f"Límite {limite_Tmax:g} °C")
    if res.t_estabilizacion is not None:
        ax.axvline(res.t_estabilizacion / f, ls=":", color="0.3",
                   label=f"|T − Tb| < {res.tolerancia:g} °C en todo: {res.t_estabilizacion / f:.1f} {unidad}")
    ax.set_ylabel("Temperatura [°C]")
    ax.set_title(f"Rectángulo {res.Lx:g} × {res.Ly:g} m, Ti = {res.Ti:g} °C, borde a {res.Tb:g} °C"
                 f"{', con hidratación' if res.hidratacion else ''} (α = {res.alpha:.3g} m²/s, dx = {res.dx:g} m)",
                 fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    if res.hidratacion:
        ax_dT.plot(t, res.dT_t, color="C1", label=f"ΔT = Tmax − Tmin({res.profundidad * 1000:g} mm)")
        ax_dT.axhline(limite_dT, color="C1", ls="--", lw=1, label=f"Límite {limite_dT:g} °C")
        ax_dT.set_ylabel("ΔT [°C]")
        ax_dT.legend(fontsize=8)
        ax_dT.grid(alpha=0.3)
        ax = ax_dT
    ax.set_xlabel(f"Tiempo [{unidad}]")
    fig.tight_layout()
    fig.savefig(ruta, dpi=150)
    plt.close(fig)


def graficar_campos_rectangulo(res, ruta):
    """Campo T(x, y) en cuatro instantes.

    Sin hidratación: 0; 0,1; 0,5 y 2 veces la constante de tiempo τ.
    Con hidratación: 0, el momento de Tmax, y τ y 3τ después.
    """
    from .rectangulo import unidad_tiempo

    f, unidad = unidad_tiempo(res.t[-1])
    if res.hidratacion:
        t_max = res.t[res.posicion_Tmax()[0]]
        instantes = (0.0, t_max, t_max + res.tau, t_max + 3 * res.tau)
    else:
        instantes = tuple(fr * res.tau for fr in (0.0, 0.1, 0.5, 2.0))
    fig, ejes = plt.subplots(2, 2, figsize=(11, 6.5), sharex=True, sharey=True)
    for ax, instante in zip(ejes.flat, instantes):
        i = int(np.argmin(np.abs(res.t - instante)))
        imagen = _dibujar_campo(ax, res, i)
        ax.set_title(f"t = {res.t[i] / f:.2f} {unidad}  ({res.t[i] / res.tau:.2f} τ)"
                     + (f"   Tmax = {res.T[i].max():.1f} °C" if res.hidratacion else ""))
        ax.label_outer()
    fig.colorbar(imagen, ax=ejes, label="Temperatura [°C]", shrink=0.9)
    fig.savefig(ruta, dpi=150, bbox_inches="tight")
    plt.close(fig)


def animar_rectangulo(res, ruta, fps=15):
    """GIF de la evolución: campo a la izquierda y T(t) de los puntos a la derecha."""
    from matplotlib.animation import FuncAnimation, PillowWriter

    from .rectangulo import unidad_tiempo

    f, unidad = unidad_tiempo(res.t[-1])
    t = res.t / f
    fig, (ax_campo, ax_t) = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios": [1.6, 1]})
    imagen = _dibujar_campo(ax_campo, res, 0)
    fig.colorbar(imagen, ax=ax_campo, label="T [°C]", shrink=0.85)

    ax_t.axhline(res.Tb, ls="--", color="0.4", lw=1)
    lineas, marcas = [], []
    for j, p in enumerate(res.puntos):
        # De fondo, la curva completa (exacta si existe): se ve hacia dónde va.
        fondo = res.T_puntos[:, j] if res.T_puntos_exacta is None else res.T_puntos_exacta[:, j]
        ax_t.plot(t, fondo, color=f"C{j}", lw=0.8, alpha=0.35)
        lineas.append(ax_t.plot([], [], color=f"C{j}", label=f"T{_nombre_punto(p)}")[0])
        marcas.append(ax_t.plot([], [], "o", color=f"C{j}")[0])
    ax_t.set_xlim(0, t[-1])
    ax_t.set_ylim(min(res.Tb, res.T.min()) - 0.5, res.T.max() + 0.5)
    ax_t.set_xlabel(f"Tiempo [{unidad}]")
    ax_t.set_ylabel("T [°C]")
    ax_t.legend(fontsize=8, loc="upper right")
    ax_t.grid(alpha=0.3)
    fig.tight_layout()

    def cuadro(i):
        imagen.set_data(res.T[i].T)
        ax_campo.set_title(f"t = {t[i]:.2f} {unidad}" + (f"   Tmax = {res.T[i].max():.1f} °C" if res.hidratacion else ""))
        for j, (linea, marca) in enumerate(zip(lineas, marcas)):
            linea.set_data(t[:i + 1], res.T_puntos[:i + 1, j])
            marca.set_data([t[i]], [res.T_puntos[i, j]])
        return [imagen, *lineas, *marcas]

    animacion = FuncAnimation(fig, cuadro, frames=res.indices_animacion(), blit=False)
    animacion.save(ruta, writer=PillowWriter(fps=fps), dpi=80)
    plt.close(fig)


def graficar_todo(res: Resultado, met: MetricasTermicas, carpeta) -> list[Path]:
    carpeta = Path(carpeta)
    carpeta.mkdir(parents=True, exist_ok=True)
    rutas = [carpeta / "series.png", carpeta / "perfiles.png", carpeta / "mapa_zt.png"]
    graficar_series(res, met, rutas[0])
    graficar_perfiles(res, rutas[1])
    graficar_mapa(res, rutas[2])
    return rutas
