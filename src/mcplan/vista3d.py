"""Vista 3D de la fundación voxelizada: caras exteriores de las celdas coloreadas por temperatura.

- `graficar_vista_3d`: figura PNG (Matplotlib) con un cuarto de la fundación cortado
  (x > 0, y < 0) para ver el interior.
- `exportar_html_3d`: página HTML interactiva (Plotly desde CDN, no requiere instalar
  nada en Python): una traza por vertida, que se oculta o muestra con un clic en la
  leyenda, y un deslizador con los instantes guardados.

Para que la vista sea liviana, la grilla se submuestrea a celdas de ~dx_vista.
"""

import json
from pathlib import Path

import numpy as np

_DIRECCIONES = [  # (eje, signo): cara de la celda en esa dirección
    (0, -1), (0, 1), (1, -1), (1, 1), (2, -1), (2, 1),
]


def _submuestreo(dx: float, dx_vista: float) -> int:
    return max(1, int(round(dx_vista / dx)))


def caras_visibles(mostrar: np.ndarray, grupo: np.ndarray | None = None):
    """Caras de las celdas de `mostrar` cuya vecina no se muestra (o es de otro `grupo`).

    Devuelve (indices (n, 3) de la celda, eje, signo) de cada cara.
    """
    grupo = np.zeros(mostrar.shape, dtype=np.int8) if grupo is None else grupo
    celdas, ejes, signos = [], [], []
    for eje, s in _DIRECCIONES:
        vecina_mostrar = np.zeros_like(mostrar)
        vecina_grupo = np.zeros_like(grupo)
        a = [slice(None)] * 3
        b = [slice(None)] * 3
        if s > 0:
            a[eje], b[eje] = slice(None, -1), slice(1, None)
        else:
            a[eje], b[eje] = slice(1, None), slice(None, -1)
        vecina_mostrar[tuple(a)] = mostrar[tuple(b)]
        vecina_grupo[tuple(a)] = grupo[tuple(b)]
        visible = mostrar & ~(vecina_mostrar & (vecina_grupo == grupo))
        idx = np.argwhere(visible)
        celdas.append(idx)
        ejes.append(np.full(len(idx), eje))
        signos.append(np.full(len(idx), s))
    return np.concatenate(celdas), np.concatenate(ejes), np.concatenate(signos)


def _cuadrilateros(idx, ejes, signos, origen, d):
    """Vértices (n, 4, 3) de cada cara, en metros."""
    base = origen + idx * d  # esquina inferior de la celda
    q = np.zeros((len(idx), 4, 3))
    for n_eje in range(3):
        sel = ejes == n_eje
        u, v = [e for e in range(3) if e != n_eje]
        esquina = base[sel].copy()
        esquina[:, n_eje] += np.where(signos[sel] > 0, d, 0.0)
        for c, (du, dv) in enumerate(((0, 0), (1, 0), (1, 1), (0, 1))):
            p = esquina.copy()
            p[:, u] += du * d
            p[:, v] += dv * d
            q[sel, c] = p
    return q


def _vista(res, t_h, dx_vista, cortar=True):
    """Campo submuestreado del instante t_h (o el guardado más cercano) y máscara de lo que se muestra."""
    from .simulacion3d import HORMIGON

    dom = res.dominio
    t = min(res.campos, key=lambda x: abs(x - t_h))
    s = _submuestreo(dom.dx, dx_vista)
    T = res.campos[t][::s, ::s, ::s]
    vid = dom.vertida_id[::s, ::s, ::s]
    mostrar = (dom.material[::s, ::s, ::s] == HORMIGON) & np.isfinite(T)
    d = dom.dx * s
    origen = dom.origen + (np.array([0.5, 0.5, 0.5]) * dom.dx - 0.5 * d)  # celdas submuestreadas centradas
    if cortar:
        xc = origen[0] + (np.arange(T.shape[0]) + 0.5) * d
        yc = origen[1] + (np.arange(T.shape[1]) + 0.5) * d
        mostrar &= ~((xc[:, None, None] > 0) & (yc[None, :, None] < 0))
    return t, T, vid, mostrar, origen, d


def graficar_vista_3d(res, instantes_h, ruta, dx_vista=0.5, vmin=None, vmax=None):
    """Vista 3D con un cuarto cortado en uno o más instantes."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import colormaps
    from matplotlib.colors import Normalize
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    vmin = 15.0 if vmin is None else vmin
    vmax = float(np.nanmax(res.Tmax)) if vmax is None else vmax
    norma, mapa = Normalize(vmin, vmax), colormaps["inferno"]
    fig = plt.figure(figsize=(6.2 * len(instantes_h), 5.6))
    for n, t_h in enumerate(instantes_h, start=1):
        t, T, _, mostrar, origen, d = _vista(res, t_h, dx_vista)
        idx, ejes, signos = caras_visibles(mostrar)
        quads = _cuadrilateros(idx, ejes, signos, origen, d)
        valores = T[tuple(idx.T)]
        sombra = np.select([ejes == 2, ejes == 0], [1.0, 0.85], 0.7)  # luz simple por orientación
        colores = mapa(norma(valores))
        colores[:, :3] *= sombra[:, None]
        ax = fig.add_subplot(1, len(instantes_h), n, projection="3d")
        ax.add_collection3d(Poly3DCollection(quads, facecolors=colores, edgecolors="none"))
        ax.set_xlim(-9.5, 9.5)
        ax.set_ylim(-9.5, 9.5)
        ax.set_zlim(0, 4.5)
        ax.set_box_aspect((19, 19, 4.5 * 1.8))
        ax.view_init(elev=24, azim=-58)
        activas = res.plan.vertidas_activas(t)
        ax.set_title(f"t = {t:g} h  (vertidas {activas[0]}–{activas[-1]})" if activas else f"t = {t:g} h", fontsize=10)
        ax.set_xlabel("x [m]")
        ax.set_ylabel("y [m]")
        ax.set_zlabel("z [m]")
    barra = fig.colorbar(plt.cm.ScalarMappable(norma, mapa), ax=fig.axes, shrink=0.6, pad=0.06)
    barra.set_label("T [°C]")
    fig.suptitle("Fundación 3D con un cuarto cortado (x > 0, y < 0): temperatura en las caras visibles", fontsize=11)
    fig.savefig(ruta, dpi=140, bbox_inches="tight")
    plt.close(fig)


def exportar_html_3d(res, ruta, dx_vista=0.5):
    """HTML interactivo: una traza por vertida (clic en la leyenda = ocultar/mostrar) y un deslizador de tiempo."""
    dom = res.dominio
    nombres = dom.vox.nombres
    vmax = float(np.nanmax(res.Tmax))
    cuadros = []
    for t in sorted(res.campos):
        _, T, vid, mostrar, origen, d = _vista(res, t, dx_vista)
        trazas = []
        for v in range(1, len(nombres) + 1):
            idx, ejes, signos = caras_visibles(mostrar & (vid == v))
            q = _cuadrilateros(idx, ejes, signos, origen, d).reshape(-1, 3)
            val = np.repeat(T[tuple(idx.T)], 4) if len(idx) else np.zeros(0)
            base = np.arange(len(idx)) * 4
            trazas.append({
                "type": "mesh3d", "name": nombres[v - 1], "showlegend": True, "showscale": v == 1,
                "x": np.round(q[:, 0], 3).tolist(), "y": np.round(q[:, 1], 3).tolist(),
                "z": np.round(q[:, 2], 3).tolist(),
                "i": np.concatenate([base, base]).tolist(), "j": np.concatenate([base + 1, base + 2]).tolist(),
                "k": np.concatenate([base + 2, base + 3]).tolist(),
                "intensity": np.round(val, 2).tolist(), "colorscale": "Inferno", "cmin": 15, "cmax": vmax,
                "colorbar": {"title": {"text": "T [°C]"}}, "flatshading": True,
                "hovertemplate": nombres[v - 1] + "<br>T = %{intensity:.1f} °C<extra></extra>",
            })
        cuadros.append({"name": f"{t:g} h", "data": trazas})
    figura = {
        "data": cuadros[0]["data"],
        "frames": cuadros,
        "layout": {
            "title": {"text": "Fundación: temperatura con un cuarto cortado "
                              "(clic en la leyenda para ocultar/mostrar una vertida)"},
            "scene": {"aspectmode": "manual", "aspectratio": {"x": 1, "y": 1, "z": 0.4},
                      "xaxis": {"range": [-9.5, 9.5]}, "yaxis": {"range": [-9.5, 9.5]}, "zaxis": {"range": [0, 4.5]}},
            "legend": {"x": 0, "y": 1},
            "sliders": [{"currentvalue": {"prefix": "t = "}, "steps": [
                {"label": c["name"], "method": "animate",
                 "args": [[c["name"]], {"mode": "immediate", "frame": {"duration": 0, "redraw": True}}]}
                for c in cuadros]}],
        },
    }
    html = f"""<!DOCTYPE html>
<html lang="es"><head><meta charset="utf-8"><title>Fundación 3D</title>
<script src="https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"></script>
<style>html,body{{margin:0;height:100%;background:#fff;font-family:sans-serif}}#g{{width:100%;height:100%}}</style>
</head><body><div id="g"></div><script>
const fig = {json.dumps(figura)};
Plotly.newPlot("g", fig.data, fig.layout).then(() => Plotly.addFrames("g", fig.frames));
</script></body></html>"""
    Path(ruta).write_text(html, encoding="utf-8")
