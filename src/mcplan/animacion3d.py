"""Animación 3D (artística) del ejercicio del bloque: render con Blender de los cuatro casos.

Flujo de `mcplan bloque-3d`:
1. Se simulan los casos con el solver 3D (bloque.py) y, en cada instante
   guardado, se lee la temperatura en las superficies visibles de cada bloque
   con un cuarto cortado (x > 0, y < 0): caras exteriores y los dos planos de
   corte, que muestran el interior. Las temperaturas pasan a colores (inferno,
   escala común) y se guardan en un .npz.
2. Blender (scripts/render_bloque_3d.py, sin interfaz) arma la escena, la
   renderiza con EEVEE cuadro a cuadro y junta un MP4.
3. Con Pillow se arma también un GIF liviano.

Los colores son la temperatura que calculó nuestro solver: no hay nada inventado,
sólo la presentación es artística.
"""

import glob
import os
import shutil
from pathlib import Path

import numpy as np

from .bloque import Interpolador, casos_pizarra, simular_bloque, t_final_comun

# Rectángulos visibles del bloque sin el cuarto x > 0, y < 0 (coordenadas desde el centro):
# (eje fijo, signo del valor: -1, 0 o +1 medio lado, eje u, rango u, eje v, rango v, normal hacia afuera)
# Los rangos están en fracciones del medio lado: -1 .. 1.
_RECTANGULOS = [
    (2, +1, 0, (-1, 0), 1, (-1, 1), (0, 0, 1)),  # techo, mitad x < 0
    (2, +1, 0, (0, 1), 1, (0, 1), (0, 0, 1)),  # techo, cuarto x > 0, y > 0
    (0, -1, 1, (-1, 1), 2, (-1, 1), (-1, 0, 0)),  # cara x = -lx
    (1, +1, 0, (-1, 1), 2, (-1, 1), (0, 1, 0)),  # cara y = +ly
    (0, +1, 1, (0, 1), 2, (-1, 1), (1, 0, 0)),  # cara x = +lx, sólo y > 0
    (1, -1, 0, (-1, 0), 2, (-1, 1), (0, -1, 0)),  # cara y = -ly, sólo x < 0
    (0, 0, 1, (-1, 0), 2, (-1, 1), (1, 0, 0)),  # corte x = 0 (mira hacia +x)
    (1, 0, 0, (0, 1), 2, (-1, 1), (0, -1, 0)),  # corte y = 0 (mira hacia -y)
]


def superficie_cortada(L, paso):
    """Malla (vértices, cuadriláteros) de las caras visibles del bloque con un cuarto cortado.

    Coordenadas desde el centro del bloque; cada cuadrilátero queda ordenado con
    la normal hacia afuera, para que Blender lo ilumine del lado correcto.
    """
    medio = np.asarray(L, dtype=float) / 2
    vertices, caras = [], []
    for fijo, signo, eu, (u0, u1), ev, (v0, v1), normal in _RECTANGULOS:
        nu = max(2, round((u1 - u0) * medio[eu] / paso) + 1)
        nv = max(2, round((v1 - v0) * medio[ev] / paso) + 1)
        U, V = np.meshgrid(np.linspace(u0, u1, nu) * medio[eu], np.linspace(v0, v1, nv) * medio[ev], indexing="ij")
        p = np.zeros((nu * nv, 3))
        p[:, fijo] = signo * medio[fijo]
        p[:, eu], p[:, ev] = U.ravel(), V.ravel()
        base = sum(len(v) for v in vertices)
        i, j = np.meshgrid(np.arange(nu - 1), np.arange(nv - 1), indexing="ij")
        a = (i * nv + j).ravel() + base
        q = np.stack([a, a + nv, a + nv + 1, a + 1], axis=1)
        # Orientación: si (b - a) x (d - a) no apunta hacia afuera, se invierte el orden.
        pa, pb, pd = p[q[0, 0] - base], p[q[0, 1] - base], p[q[0, 3] - base]
        if np.dot(np.cross(pb - pa, pd - pa), normal) < 0:
            q = q[:, ::-1]
        vertices.append(p)
        caras.append(q)
    return np.vstack(vertices), np.vstack(caras)


def _srgb_a_lineal(c):
    """Blender guarda los colores de atributo en espacio lineal."""
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def generar_datos(ruta_npz, L=(3.0, 3.0, 0.6), dx=0.05, beta=None, t_final=None, paso_malla=0.05,
                  k=2.222, rho=2400.0, cp=970.0, Ti=20.0, pulso=None, mapa="inferno"):
    """Simula los cuatro casos y guarda, por cuadro, los colores de la superficie de cada bloque."""
    from matplotlib import colormaps

    from .hidratacion import CalorPulsoCoseno

    pulso = CalorPulsoCoseno() if pulso is None else pulso
    casos = casos_pizarra(k, beta)
    t_final = t_final or t_final_comun(casos, L, k, rho, cp, Ti, pulso)
    vertices, caras = superficie_cortada(L, paso_malla)
    n = tuple(round(Li / dx) for Li in L)
    lector = Interpolador(vertices + np.asarray(L) / 2, L, n)  # el interpolador usa coordenadas desde la esquina

    temperaturas, tmax, etiquetas = [], [], []
    for caso in casos:
        muestras = []
        r = simular_bloque(caso, L, dx, k, rho, cp, Ti, pulso, ((0.0, 0.0, 0.0),), "centro", t_final,
                           con_exacta=False, observador=lambda t, P: muestras.append(lector(P)))
        temperaturas.append(np.array(muestras, dtype=np.float32))
        tmax.append(r.Tmax_t)
        etiquetas.append({"dirichlet": f"Dirichlet {caso.T_borde:g} °C", "neumann": "Neumann (adiabático)",
                          "robin": f"Robin  β = {caso.h / k:.1f} 1/m".replace(".", ",")}[caso.tipo])
        t = r.t
    temperaturas = np.array(temperaturas)  # (casos, cuadros, vértices)
    T_min, T_max = float(min(Ti, temperaturas.min())), float(temperaturas.max())
    rgb = colormaps[mapa]((temperaturas - T_min) / (T_max - T_min))[..., :3]
    barra = colormaps[mapa](np.linspace(0, 1, 64))[:, :3]
    np.savez_compressed(
        ruta_npz, vertices=vertices.astype(np.float32), caras=caras.astype(np.int32),
        colores=_srgb_a_lineal(rgb).astype(np.float16), barra=_srgb_a_lineal(barra).astype(np.float32),
        t=t, Q=pulso.Q(t / 3600.0), tmax=np.array(tmax), etiquetas=np.array(etiquetas),
        escala=np.array([T_min, T_max]), L=np.asarray(L), Q_max=pulso.Q_max, T_h=pulso.T_h,
    )
    return ruta_npz, len(t), (T_min, T_max)


def buscar_blender(indicado=None) -> str:
    """Ruta a blender: la indicada, la del PATH o la instalación más nueva en Program Files."""
    if indicado:
        return indicado
    en_path = shutil.which("blender")
    if en_path:
        return en_path
    candidatos = sorted(glob.glob(r"C:\Program Files\Blender Foundation\Blender *\blender.exe"))
    if candidatos:
        return candidatos[-1]
    raise ValueError("No se encontró Blender: instálelo o indique la ruta con --blender")


def script_render() -> Path:
    return Path(__file__).resolve().parents[2] / "scripts" / "render_bloque_3d.py"


def crear_gif(carpeta_cuadros, ruta_gif, ancho=800, fps=15, cada=1):
    """GIF liviano a partir de los PNG renderizados (reescalados a `ancho`)."""
    from PIL import Image

    archivos = sorted(glob.glob(os.path.join(carpeta_cuadros, "cuadro_*.png")))[::cada]
    if not archivos:
        raise ValueError(f"No hay cuadros renderizados en {carpeta_cuadros}")
    imagenes = []
    for archivo in archivos:
        with Image.open(archivo) as im:
            alto = round(im.height * ancho / im.width)
            imagenes.append(im.convert("RGB").resize((ancho, alto), Image.LANCZOS)
                            .quantize(256, method=Image.Quantize.FASTOCTREE))
    # La última imagen se queda un poco más: se lee el estado final.
    duraciones = [round(1000 * cada / fps)] * (len(imagenes) - 1) + [2500]
    imagenes[0].save(ruta_gif, save_all=True, append_images=imagenes[1:], duration=duraciones, loop=0)
    return ruta_gif
