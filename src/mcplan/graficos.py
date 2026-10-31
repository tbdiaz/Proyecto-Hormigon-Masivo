"""Gráficos con Matplotlib: columna 1D (curvas, perfiles, mapa z-t), rectángulo 2D y bloque 3D (cortes, GIF)."""

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


def _guardar_gif(fig, actualizar, cuadros, ruta, fps, dpi=80):
    """Dibuja cada cuadro con Agg y arma el GIF con Pillow.

    Es unas 3 veces más rápido que matplotlib.animation.PillowWriter (no pasa
    por savefig en cada cuadro).
    """
    from PIL import Image

    fig.set_dpi(dpi)
    imagenes = []
    for i in cuadros:
        actualizar(i)
        fig.canvas.draw()
        rgb = np.asarray(fig.canvas.buffer_rgba())[..., :3]
        imagenes.append(Image.fromarray(rgb).quantize(256, method=Image.Quantize.FASTOCTREE))
    imagenes[0].save(ruta, save_all=True, append_images=imagenes[1:], duration=round(1000 / fps), loop=0)
    plt.close(fig)


def animar_rectangulo(res, ruta, fps=15):
    """GIF de la evolución: campo a la izquierda y T(t) de los puntos a la derecha."""
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

    _guardar_gif(fig, cuadro, res.indices_animacion(), ruta, fps)


# ---------------------------------------------------------------------------
# Bloque 3D (bloque.py)
# ---------------------------------------------------------------------------

_COLOR_CASO = {"dirichlet": "C0", "neumann": "C3", "robin": "C2"}


def _color_caso(res, i):
    """Color por caso: Dirichlet en tonos azules (uno por temperatura), Neumann rojo, Robin verde."""
    if res.caso.tipo == "dirichlet":
        return ("C0", "C9", "C4", "C7")[i % 4]
    return _COLOR_CASO[res.caso.tipo]


def _nombre_corto(res) -> str:
    c = res.caso
    return {"dirichlet": f"Dirichlet {c.T_borde:g} °C", "neumann": "Neumann (adiabático)",
            "robin": f"Robin β = {c.h / res.k:.3g} 1/m"}[c.tipo]


def _instantes_bloque(res):
    """T/2 (máximo del pulso), el momento de Tmax, T (fin del pulso) y 2T, lo que quepa en la simulación."""
    T_s = res.pulso.T_h * 3600.0
    instantes = sorted({T_s / 2, res.t_Tmax, T_s, 2 * T_s})
    return [i for i in instantes if i <= res.t[-1]]


def graficar_puntos_bloque(res, ruta, limite_Tmax=70.0, limite_dT=20.0):
    """T(t) en el centro y en los puntos (numérica en línea, exacta en círculos), Tmax, Tmin y ΔT."""
    from .rectangulo import unidad_tiempo

    f, unidad = unidad_tiempo(res.t[-1])
    t = res.t / f
    cada = max(1, len(t) // 30)
    fig, (ax, ax_dT) = plt.subplots(2, 1, figsize=(9.5, 7.5), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
    series = [("centro (0, 0, 0)" if res.origen == "centro" else "centro", res.T_centro, res.T_centro_exacta)]
    for j, p in enumerate(res.puntos):
        series.append((f"({p[0]:g}, {p[1]:g}, {p[2]:g})", res.T_puntos[:, j],
                       None if res.T_puntos_exacta is None else res.T_puntos_exacta[:, j]))
    for j, (nombre, num, exacta) in enumerate(series):
        ax.plot(t, num, color=f"C{j}", label=f"T{nombre} numérica" if j else f"T {nombre} numérica")
        if exacta is not None:
            ax.plot(t[::cada], exacta[::cada], "o", mfc="none", color=f"C{j}", ms=5, label="exacta" if j == 0 else None)
    ax.plot(t, res.Tmax_t, color="C3", lw=1, ls="-.", label="Tmax del bloque")
    ax.plot(t, res.Tmin_cascara_t, color="C2", lw=1, ls="-.", label=f"Tmin a {res.profundidad * 1000:g} mm de las caras")
    ax.axhline(res.T_eq, color="0.4", ls="--", lw=1, label=f"Equilibrio {res.T_eq:.1f} °C")
    ax.axhline(limite_Tmax, color="C3", ls="--", lw=1, label=f"Límite {limite_Tmax:g} °C")
    if res.t_estabilizacion is not None:
        ax.axvline(res.t_estabilizacion / f, ls=":", color="0.3",
                   label=f"|T − T_eq| < {res.tolerancia:g} °C en todo: {res.t_estabilizacion / f:.2f} {unidad}")
    ax.set_ylabel("Temperatura [°C]")
    ax.set_title(f"Bloque {' × '.join(f'{Li:g}' for Li in res.L)} m con el pulso de calor — {res.caso.etiqueta(res.k)}",
                 fontsize=10)
    ax.legend(fontsize=7, ncol=2)
    ax.grid(alpha=0.3)
    ax_dT.plot(t, res.dT_t, color="C1", label=f"ΔT = Tmax − Tmin({res.profundidad * 1000:g} mm)")
    ax_dT.axhline(limite_dT, color="C1", ls="--", lw=1, label=f"Límite {limite_dT:g} °C")
    ax_dT.set_ylabel("ΔT [°C]")
    ax_dT.set_xlabel(f"Tiempo [{unidad}]")
    ax_dT.legend(fontsize=8)
    ax_dT.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(ruta, dpi=150)
    plt.close(fig)


def _dibujar_cortes(ax_xy, ax_xz, res, i, vmin, vmax):
    """Corte horizontal (z = centro) y vertical (y = centro) del instante i."""
    kw = {"origin": "lower", "cmap": "inferno", "vmin": vmin, "vmax": vmax, "interpolation": "nearest"}
    im_xy = ax_xy.imshow(res.corte_xy[i].T, extent=res.extension((0, 1)), **kw)
    im_xz = ax_xz.imshow(res.corte_xz[i].T, extent=res.extension((0, 2)), **kw)
    ax_xy.set_ylabel("y [m]")
    ax_xz.set_ylabel("z [m]")
    ax_xz.set_xlabel("x [m]")
    return im_xy, im_xz


def graficar_cortes_bloque(res, ruta):
    """Cortes centrales del campo en T/2, t(Tmax), T y 2T: arriba el plano horizontal, abajo el vertical."""
    instantes = _instantes_bloque(res)
    vmin, vmax = min(res.corte_xy.min(), res.corte_xz.min()), max(res.corte_xy.max(), res.corte_xz.max())
    fig, ejes = plt.subplots(2, len(instantes), figsize=(3.4 * len(instantes) + 1.5, 4.6), squeeze=False,
                             gridspec_kw={"height_ratios": [res.L[1], 2.5 * res.L[2]], "hspace": 0.35,
                                          "wspace": 0.12})
    for c, instante in enumerate(instantes):
        i = int(np.argmin(np.abs(res.t - instante)))
        im, _ = _dibujar_cortes(ejes[0, c], ejes[1, c], res, i, vmin, vmax)
        nota = "  (Tmax)" if abs(res.t[i] - res.t_Tmax) < 1.0 else ""
        ejes[0, c].set_title(f"t = {res.t[i] / 3600:.1f} h{nota}\nTmax = {res.Tmax_t[i]:.1f} °C", fontsize=10)
        if c:
            for ax in ejes[:, c]:
                ax.set_ylabel("")
                ax.tick_params(labelleft=False)
    fig.suptitle(f"{res.caso.etiqueta(res.k)} — arriba: plano horizontal central; abajo: plano vertical central",
                 fontsize=11, y=1.06)
    fig.colorbar(im, ax=ejes, label="Temperatura [°C]", shrink=0.9)
    fig.savefig(ruta, dpi=150, bbox_inches="tight")
    plt.close(fig)


def graficar_comparacion_bloque(resultados, ruta, limite_dT=20.0):
    """Los casos juntos: T en el centro (numérica y exacta), Tmin a 25 mm de las caras y ΔT."""
    from .rectangulo import unidad_tiempo

    f, unidad = unidad_tiempo(max(r.t[-1] for r in resultados))
    fig, (ax_c, ax_s, ax_d) = plt.subplots(3, 1, figsize=(9.5, 10), sharex=True)
    for n, r in enumerate(resultados):
        color, t = _color_caso(r, n), r.t / f
        cada = max(1, len(t) // 30)
        ax_c.plot(t, r.T_centro, color=color, label=_nombre_corto(r))
        if r.T_centro_exacta is not None:
            ax_c.plot(t[::cada], r.T_centro_exacta[::cada], "o", mfc="none", ms=4, color=color)
        ax_s.plot(t, r.Tmin_cascara_t, color=color, label=_nombre_corto(r))
        ax_d.plot(t, r.dT_t, color=color, label=_nombre_corto(r))
    ax_c.set_ylabel("T en el centro [°C]")
    ax_c.set_title("Mismo bloque y mismo calor, distinta condición de borde (líneas: numérica; círculos: exacta)",
                   fontsize=10)
    ax_s.set_ylabel(f"Tmin a {resultados[0].profundidad * 1000:g} mm de las caras [°C]")
    ax_d.axhline(limite_dT, color="0.3", ls="--", lw=1, label=f"Límite {limite_dT:g} °C")
    ax_d.set_ylabel("ΔT [°C]")
    ax_d.set_xlabel(f"Tiempo [{unidad}]")
    for ax in (ax_c, ax_s, ax_d):
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(ruta, dpi=150)
    plt.close(fig)


def animar_bloque(resultados, ruta, fps=12, cuadros=100):
    """GIF con los casos lado a lado: cortes horizontal y vertical, y T en el centro."""
    t = resultados[0].t
    vmin = min(min(r.corte_xy.min(), r.corte_xz.min()) for r in resultados)
    vmax = max(max(r.corte_xy.max(), r.corte_xz.max()) for r in resultados)
    n = len(resultados)
    fig = plt.figure(figsize=(3.2 * n + 1.2, 7.4))
    grilla = fig.add_gridspec(3, n, height_ratios=[3, 1, 2.2], hspace=0.45)
    cortes, titulos = [], []
    for c, r in enumerate(resultados):
        ax_xy, ax_xz = fig.add_subplot(grilla[0, c]), fig.add_subplot(grilla[1, c])
        cortes.append(_dibujar_cortes(ax_xy, ax_xz, r, 0, vmin, vmax))
        titulos.append(ax_xy)
        if c:
            for ax in (ax_xy, ax_xz):
                ax.set_ylabel("")
                ax.tick_params(labelleft=False)
    fig.colorbar(cortes[-1][0], ax=[a for a in fig.axes], label="T [°C]", shrink=0.6)
    ax_t = fig.add_subplot(grilla[2, :])
    lineas, marcas = [], []
    for c, r in enumerate(resultados):
        color = _color_caso(r, c)
        ax_t.plot(r.t / 3600, r.T_centro, color=color, lw=0.8, alpha=0.3)
        lineas.append(ax_t.plot([], [], color=color, label=_nombre_corto(r))[0])
        marcas.append(ax_t.plot([], [], "o", color=color)[0])
    ax_t.set_xlim(0, t[-1] / 3600)
    ax_t.set_ylim(vmin - 1, vmax + 1)
    ax_t.set_xlabel("Tiempo [h]")
    ax_t.set_ylabel("T en el centro [°C]")
    ax_t.legend(fontsize=8, loc="upper right")
    ax_t.grid(alpha=0.3)

    def cuadro(i):
        fig.suptitle(f"Bloque {' × '.join(f'{Li:g}' for Li in resultados[0].L)} m con el pulso de calor — "
                     f"t = {t[i] / 3600:.0f} h", fontsize=11)
        for c, r in enumerate(resultados):
            cortes[c][0].set_data(r.corte_xy[i].T)
            cortes[c][1].set_data(r.corte_xz[i].T)
            titulos[c].set_title(f"{_nombre_corto(r)}\nTmax = {r.Tmax_t[i]:.1f} °C", fontsize=9)
            lineas[c].set_data(r.t[:i + 1] / 3600, r.T_centro[:i + 1])
            marcas[c].set_data([r.t[i] / 3600], [r.T_centro[i]])

    indices = np.unique(np.linspace(0, len(t) - 1, cuadros + 1).round().astype(int))
    _guardar_gif(fig, cuadro, indices, ruta, fps)


# ---------------------------------------------------------------------------
# Voxelización (voxelizacion.py)
# ---------------------------------------------------------------------------

def graficar_voxelizacion(vox, ruta):
    """Cortes de la grilla coloreados por vertida: planta a media altura de cada vertida y corte vertical central."""
    from matplotlib.colors import ListedColormap

    colores = ["#ffffff", "#d9984d", "#4d8ccc", "#73b866", "#d16b73", "#9478c7"][:vox.n_vertidas + 1]
    mapa = ListedColormap(colores)
    x, y, z = (vox.centros(e) for e in range(3))
    ext_xz = (x[0] - vox.dx / 2, x[-1] + vox.dx / 2, z[0] - vox.dx / 2, z[-1] + vox.dx / 2)
    ext_xy = (x[0] - vox.dx / 2, x[-1] + vox.dx / 2, y[0] - vox.dx / 2, y[-1] + vox.dx / 2)
    fig = plt.figure(figsize=(13, 7.2))
    grilla = fig.add_gridspec(2, 3, height_ratios=[1, 1.15])
    ax = fig.add_subplot(grilla[0, :])
    j = int(np.argmin(np.abs(y)))
    ax.imshow(vox.vertida_id[:, j, :].T, origin="lower", extent=ext_xz, cmap=mapa, vmin=0,
              vmax=vox.n_vertidas, interpolation="nearest", aspect="auto")
    ax.set_title(f"Corte vertical y = 0 (dx = {vox.dx:g} m): cada color es una vertida; los bordes inclinados quedan en escalera",
                 fontsize=10)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("z [m]")
    for v in range(1, vox.n_vertidas + 1):
        ax.plot([], [], "s", color=colores[v], label=f"{vox.nombres[v - 1]} (id {v})")
    ax.legend(fontsize=8, ncol=2, loc="upper left")
    for c, v in enumerate((1, 3, 5)):
        if v > vox.n_vertidas:
            break
        k = int(np.round(np.mean(np.nonzero(np.any(vox.vertida_id == v, axis=(0, 1)))[0])))
        axp = fig.add_subplot(grilla[1, c])
        axp.imshow(vox.vertida_id[:, :, k].T, origin="lower", extent=ext_xy, cmap=mapa, vmin=0,
                   vmax=vox.n_vertidas, interpolation="nearest")
        axp.set_title(f"Planta z = {z[k]:.2f} m ({vox.nombres[v - 1]})", fontsize=10)
        axp.set_xlabel("x [m]")
        if c == 0:
            axp.set_ylabel("y [m]")
    fig.tight_layout()
    fig.savefig(ruta, dpi=150)
    plt.close(fig)


def graficar_error_volumen(voxelizaciones, ruta):
    """Error de volumen (total y por vertida) en función de dx."""
    dxs = [v.dx for v in voxelizaciones]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.loglog(dxs, [abs(v.error_volumen_total) * 100 for v in voxelizaciones], "ko-", lw=2, label="total")
    for i in range(1, voxelizaciones[0].n_vertidas + 1):
        ax.loglog(dxs, [max(abs(v.error_volumen(i)) * 100, 1e-4) for v in voxelizaciones], ".--", alpha=0.7,
                  label=voxelizaciones[0].nombres[i - 1])
    ax.set_xlabel("dx [m]")
    ax.set_ylabel("|error de volumen| [%]")
    ax.set_title("Volumen voxelizado vs. IFC (con dx = 0,05 m las juntas caen en caras de celda)", fontsize=10)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8)
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
