"""Exportación a VTK (ParaView / PyVista) y vista con PyVista (E11).

Los campos 3D se escriben como VTK ImageData (.vti, XML con datos binarios comprimidos
con zlib), un archivo por instante, y un .pvd que los junta en una serie de tiempo:
ParaView abre el .pvd y deja recorrer los instantes. El formato se escribe aquí mismo
(sin depender de VTK); PyVista sólo se usa para mirar los resultados.

Cada celda lleva: T [°C] (NaN en el aire y en vertidas aún no colocadas), te [h],
vertida_id, material (0 aire, 1 suelo, 2 hormigón) y superficie (1 = hormigón con
una cara al aire en ese instante).
"""

import re
import struct
import zlib
from pathlib import Path

import numpy as np

_TIPOS = {np.dtype("float32"): "Float32", np.dtype("float64"): "Float64", np.dtype("int8"): "Int8",
          np.dtype("uint8"): "UInt8", np.dtype("int32"): "Int32"}
_BLOQUE = 1 << 16


def _comprimir(datos: bytes) -> bytes:
    """Bloques zlib con el encabezado de VTK (UInt64): [n_bloques, tamaño, último, comprimidos...]."""
    bloques = [datos[i:i + _BLOQUE] for i in range(0, len(datos), _BLOQUE)] or [b""]
    comprimidos = [zlib.compress(b, 6) for b in bloques]
    cabeza = struct.pack(f"<{3 + len(comprimidos)}Q", len(bloques), _BLOQUE, len(bloques[-1]), *map(len, comprimidos))
    return cabeza + b"".join(comprimidos)


def escribir_vti(ruta, origen, dx, celdas: dict) -> Path:
    """celdas: nombre -> arreglo (nx, ny, nz). VTK recorre x primero: se aplanan en orden Fortran."""
    forma = next(iter(celdas.values())).shape
    nx, ny, nz = forma
    xml, binario, desplazamiento = [], [], 0
    for nombre, a in celdas.items():
        a = np.asarray(a)
        if a.shape != forma:
            raise ValueError(f"{nombre}: forma {a.shape}, se esperaba {forma}")
        crudo = a.astype(a.dtype.newbyteorder("<")).ravel(order="F").tobytes()
        bloque = _comprimir(crudo)
        xml.append(f'        <DataArray type="{_TIPOS[a.dtype]}" Name="{nombre}" format="appended" '
                   f'offset="{desplazamiento}"/>')
        binario.append(bloque)
        desplazamiento += len(bloque)
    ox, oy, oz = (float(c) for c in origen)
    cabeza = (
        '<?xml version="1.0"?>\n'
        '<VTKFile type="ImageData" version="1.0" byte_order="LittleEndian" header_type="UInt64" '
        'compressor="vtkZLibDataCompressor">\n'
        f'  <ImageData WholeExtent="0 {nx} 0 {ny} 0 {nz}" Origin="{ox} {oy} {oz}" Spacing="{dx} {dx} {dx}">\n'
        f'    <Piece Extent="0 {nx} 0 {ny} 0 {nz}">\n'
        '      <CellData Scalars="T">\n' + "\n".join(xml) + "\n      </CellData>\n"
        "    </Piece>\n  </ImageData>\n"
        '  <AppendedData encoding="raw">\n   _'
    )
    ruta = Path(ruta)
    with open(ruta, "wb") as f:
        f.write(cabeza.encode("ascii"))
        for b in binario:
            f.write(b)
        f.write(b"\n  </AppendedData>\n</VTKFile>\n")
    return ruta


def escribir_pvd(ruta, entradas) -> Path:
    """entradas: [(t_h, nombre de archivo .vti)]."""
    filas = "\n".join(f'    <DataSet timestep="{t:g}" part="0" file="{n}"/>' for t, n in entradas)
    texto = ('<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n'
             f"  <Collection>\n{filas}\n  </Collection>\n</VTKFile>\n")
    Path(ruta).write_text(texto, encoding="ascii")
    return Path(ruta)


def superficie_expuesta(material, existe):
    """Hormigón que existe y tiene alguna vecina inexistente (aire o vertida aún no colocada)."""
    hormigon = (material == 2) & existe
    sup = np.zeros_like(hormigon)
    for eje in range(3):
        for paso in (1, -1):
            vecina = np.roll(existe, paso, axis=eje)
            borde = [slice(None)] * 3
            borde[eje] = 0 if paso == 1 else -1
            vecina[tuple(borde)] = False
            sup |= hormigon & ~vecina
    return sup.astype(np.uint8)


def exportar_vtk(res, carpeta) -> Path:
    """Un .vti por instante guardado en res.campos (instantes_3d_h y el de la Tmax) y un .pvd."""
    dom = res.dominio
    carpeta = Path(carpeta)
    carpeta.mkdir(parents=True, exist_ok=True)
    entradas = []
    for t, campo in res.campos.items():
        existe = np.isfinite(campo)
        nombre = f"T_{int(round(t)):03d}h.vti"
        escribir_vti(carpeta / nombre, dom.origen, dom.dx, {
            "T": campo.astype(np.float32), "te": res.edades[t].astype(np.float32),
            "vertida_id": dom.vertida_id.astype(np.int8), "material": dom.material.astype(np.int8),
            "superficie": superficie_expuesta(dom.material, existe),
        })
        entradas.append((t, nombre))
    return escribir_pvd(carpeta / "fundacion.pvd", entradas)


def instantes_vtk(carpeta) -> list[tuple[float, Path]]:
    carpeta = Path(carpeta)
    texto = (carpeta / "fundacion.pvd").read_text(encoding="ascii")
    return [(float(t), carpeta / n) for t, n in re.findall(r'timestep="([^"]+)" part="0" file="([^"]+)"', texto)]


def nucleo_y_superficie(malla) -> dict:
    """Celda más caliente del hormigón (núcleo) y celda de superficie más fría, en una malla de PyVista."""
    T = np.asarray(malla.cell_data["T"])
    material = np.asarray(malla.cell_data["material"])
    centros = malla.cell_centers().points
    es_h = (material == 2) & np.isfinite(T)
    i_n = int(np.flatnonzero(es_h)[np.nanargmax(T[es_h])])
    sup = es_h & (np.asarray(malla.cell_data["superficie"]) == 1)
    i_s = int(np.flatnonzero(sup)[np.nanargmin(T[sup])]) if sup.any() else None
    return {"nucleo": (centros[i_n].tolist(), float(T[i_n])),
            "superficie": None if i_s is None else (centros[i_s].tolist(), float(T[i_s]))}


def _para_dibujar(malla):
    """Superficie con un juego de puntos por cara y T en los puntos.

    VTK pinta los datos por celda con un "texture buffer" que algunos drivers (p. ej. Intel
    Iris de 2020) leen mal: todo sale oscuro y la barra de color se repite. Con los puntos
    sin compartir, T por punto da el mismo color exacto por celda sin usar ese camino.
    """
    return malla.extract_surface(algorithm="dataset_surface").shrink(1.0).cell_data_to_point_data()


def _barra_de_color(p, vmin, vmax, a0, a1, z, fijo, a_lo_largo_de_y=False):
    """Barra de color dibujada como una franja con T en los puntos (la barra de VTK sufre el mismo problema).

    La franja va a lo largo de x (o de y) a la altura z, en la coordenada `fijo` del otro eje.
    """
    import pyvista as pv

    s = np.linspace(a0, a1, 129)
    S_, Z = np.meshgrid(s, [z, z + 0.5], indexing="ij")
    F = np.full(S_.shape, fijo)
    X, Y = (F, S_) if a_lo_largo_de_y else (S_, F)
    barra = pv.StructuredGrid(X[:, None, :], Y[:, None, :], Z[:, None, :])
    barra.point_data["T"] = (vmin + (S_ - a0) / (a1 - a0) * (vmax - vmin))[:, None, :].ravel(order="F")
    p.add_mesh(barra, scalars="T", cmap="inferno", clim=(vmin, vmax), show_scalar_bar=False, lighting=False)
    marcas = np.linspace(vmin, vmax, 5)
    pos = [a0 + (m - vmin) / (vmax - vmin) * (a1 - a0) for m in marcas]
    puntos = [[fijo, q, z - 0.25] if a_lo_largo_de_y else [q, fijo, z - 0.25] for q in pos]
    p.add_point_labels(np.array(puntos), [f"{m:.0f} °C" for m in marcas], font_size=13, shape_opacity=0.0,
                       always_visible=True, show_points=False)


def ver_pyvista(archivo, corte=None, recorte=None, posicion=0.0, sin_suelo=False, captura=None, interactivo=False,
                escala_z=1.0):
    """Abre un .vti con PyVista: slice (corte) o clipping (recorte) por un plano x, y o z, y marca núcleo y superficie.

    escala_z exagera la vertical sólo en la imagen (las coordenadas informadas son las reales).
    """
    import pyvista as pv

    malla = pv.read(str(archivo))
    solido = malla.threshold(1.5 if sin_suelo else 0.5, scalars="material")  # saca el aire (y el suelo)
    info = nucleo_y_superficie(solido)
    normal = {"x": (1, 0, 0), "y": (0, 1, 0), "z": (0, 0, 1)}
    plano = corte or recorte
    if plano:
        origen = [0.0, 0.0, 2.0]
        origen["xyz".index(plano)] = posicion
        vista = solido.slice(normal=normal[plano], origin=origen) if corte else solido.clip(
            normal=normal[plano], origin=origen)
    else:
        vista = solido
    if captura is None and not interactivo:
        return info
    ez = np.array([1.0, 1.0, escala_z])
    p = pv.Plotter(off_screen=not interactivo, window_size=(1400, 900))
    T = np.asarray(solido.cell_data["T"])
    vmin, vmax = 10.0, float(max(40.0, np.ceil(np.nanmax(T))))
    p.add_mesh(_para_dibujar(vista).scale(ez), scalars="T", cmap="inferno", clim=(vmin, vmax),
               show_scalar_bar=False, ambient=0.35, diffuse=0.65, specular=0.0)
    if plano:
        p.add_mesh(solido.outline().scale(ez), color="gray")
    x0, x1, y0, y1, z0, _ = solido.bounds
    if plano == "x":
        _barra_de_color(p, vmin, vmax, y0 + 2, y1 - 2, z0 * escala_z - 1.6, x1, a_lo_largo_de_y=True)
    else:
        _barra_de_color(p, vmin, vmax, x0 + 2, x1 - 2, z0 * escala_z - 1.6, y1 if plano == "y" else y0)
    puntos, etiquetas = [np.array(info["nucleo"][0]) * ez], [f"núcleo {info['nucleo'][1]:.1f} °C"]
    if info["superficie"]:
        puntos.append(np.array(info["superficie"][0]) * ez)
        etiquetas.append(f"superficie {info['superficie'][1]:.1f} °C")
    p.add_point_labels(np.array(puntos), etiquetas, point_size=14, font_size=16, point_color="cyan",
                       render_points_as_spheres=True, always_visible=True, shape_opacity=0.6)
    p.add_axes()
    # Vista elegida para que se vea la cara del corte (recorte en y = 0: se mira desde +y).
    vista_desde = {"x": (1.0, -0.35, 0.45), "y": (0.35, 1.0, 0.45), "z": (0.5, -0.7, 1.2)}.get(plano, (0.8, -0.9, 0.7))
    p.camera_position = [tuple(30 * c for c in vista_desde), (0, 0, 0.5 * escala_z), (0, 0, 1)]
    p.reset_camera()
    p.camera.zoom(1.15)
    p.show(screenshot=None if captura is None else str(captura), auto_close=True)
    return info
