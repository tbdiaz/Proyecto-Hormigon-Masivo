"""Render 3D (artístico) del ejercicio del bloque con Blender + EEVEE, sin interfaz.

Lo llama `mcplan bloque-3d`; también se puede usar directo:

    blender -b --factory-startup --python scripts/render_bloque_3d.py -- datos.npz carpeta
            [--res 1280 720] [--muestras 32] [--cada 1] [--fps 15] [--mp4 ruta.mp4]

datos.npz lo genera mcplan.animacion3d.generar_datos: la malla de la superficie
del bloque con un cuarto cortado y, por cuadro, el color de cada vértice (la
temperatura que calculó el solver). Aquí sólo se arma la escena y se renderiza.
"""

import math
import os
import sys

import bpy
import numpy as np

args = sys.argv[sys.argv.index("--") + 1:]
# Rutas absolutas: Blender no resuelve las relativas contra la carpeta actual.
RUTA_NPZ, CARPETA = os.path.abspath(args[0]), os.path.abspath(args[1])


def opcion(nombre, n=1, defecto=None, tipo=float):
    if nombre in args:
        i = args.index(nombre)
        valores = [tipo(v) for v in args[i + 1:i + 1 + n]]
        return valores if n > 1 else valores[0]
    return defecto


ANCHO, ALTO = opcion("--res", 2, [1280, 720], int)
MUESTRAS = opcion("--muestras", 1, 32, int)
CADA = opcion("--cada", 1, 1, int)
FPS = opcion("--fps", 1, 15, int)
RUTA_MP4 = opcion("--mp4", 1, None, str)

d = np.load(RUTA_NPZ)
VERTICES, CARAS, COLORES = d["vertices"], d["caras"], d["colores"]
T, Q, TMAX, ETIQUETAS = d["t"], d["Q"], d["tmax"], [str(e) for e in d["etiquetas"]]
T_MIN, T_MAX = (float(v) for v in d["escala"])
LX, LY, LZ = (float(v) for v in d["L"])
N_CASOS = COLORES.shape[0]
os.makedirs(CARPETA, exist_ok=True)

# ---------------------------------------------------------------------------
# Escena
# ---------------------------------------------------------------------------
bpy.ops.wm.read_factory_settings(use_empty=True)
escena = bpy.context.scene
escena.render.engine = "BLENDER_EEVEE"
escena.eevee.taa_render_samples = MUESTRAS
if hasattr(escena.eevee, "use_raytracing"):
    escena.eevee.use_raytracing = True  # reflejos en la mesa
escena.render.resolution_x, escena.render.resolution_y = ANCHO, ALTO
escena.render.image_settings.file_format = "PNG"
escena.view_settings.view_transform = "Standard"  # que los colores de la escala se vean tal cual

mundo = bpy.data.worlds.new("Fondo")
mundo.use_nodes = True
mundo.node_tree.nodes["Background"].inputs["Color"].default_value = (0.004, 0.005, 0.009, 1)
escena.world = mundo

fuente = None
for candidata in (r"C:\Windows\Fonts\segoeui.ttf", r"C:\Windows\Fonts\arial.ttf"):
    if os.path.exists(candidata):
        fuente = bpy.data.fonts.load(candidata)
        break


def material_termico():
    """El color de cada vértice es la temperatura; además brilla un poco (hormigón 'caliente')."""
    m = bpy.data.materials.new("Hormigon termico")
    m.use_nodes = True
    nodos = m.node_tree.nodes
    bsdf = nodos["Principled BSDF"]
    atributo = nodos.new("ShaderNodeAttribute")
    atributo.attribute_name = "T_color"
    m.node_tree.links.new(atributo.outputs["Color"], bsdf.inputs["Base Color"])
    m.node_tree.links.new(atributo.outputs["Color"], bsdf.inputs["Emission Color"])
    bsdf.inputs["Emission Strength"].default_value = 0.35
    bsdf.inputs["Roughness"].default_value = 0.5
    return m


def material_emision(color=(1, 1, 1), fuerza=1.0, atributo=None):
    m = bpy.data.materials.new("Emision")
    m.use_nodes = True
    nodos = m.node_tree.nodes
    nodos.remove(nodos["Principled BSDF"])
    emision = nodos.new("ShaderNodeEmission")
    emision.inputs["Strength"].default_value = fuerza
    emision.inputs["Color"].default_value = (*color, 1)
    if atributo:
        nodo = nodos.new("ShaderNodeAttribute")
        nodo.attribute_name = atributo
        m.node_tree.links.new(nodo.outputs["Color"], emision.inputs["Color"])
    m.node_tree.links.new(emision.outputs["Emission"], nodos["Material Output"].inputs["Surface"])
    return m


def texto(cuerpo, tamano, material, alinear="CENTER"):
    curva = bpy.data.curves.new("Texto", type="FONT")
    curva.body, curva.size = cuerpo, tamano
    curva.align_x, curva.align_y = alinear, "CENTER"
    if fuente:
        curva.font = fuente
    obj = bpy.data.objects.new("Texto", curva)
    escena.collection.objects.link(obj)
    obj.data.materials.append(material)
    return obj


def malla(nombre, vertices, caras, material, colores=None):
    me = bpy.data.meshes.new(nombre)
    me.from_pydata(vertices.tolist(), [], caras.tolist())
    me.update()
    atributo = me.color_attributes.new("T_color", "FLOAT_COLOR", "POINT")
    if colores is not None:
        atributo.data.foreach_set("color", colores.ravel())
    obj = bpy.data.objects.new(nombre, me)
    escena.collection.objects.link(obj)
    me.materials.append(material)
    return obj, me, atributo


def rgba(rgb):
    return np.concatenate([rgb.astype(np.float32), np.ones((len(rgb), 1), np.float32)], axis=1)


blanco = material_emision((0.95, 0.95, 0.95), 1.0)
gris = material_emision((0.55, 0.57, 0.62), 1.0)
termico = material_termico()

# Cuatro bloques en 2 x 2, con el cuarto cortado mirando a la cámara.
separacion = LX + 1.6
posiciones = [(-separacion / 2, separacion / 2), (separacion / 2, separacion / 2),
              (-separacion / 2, -separacion / 2), (separacion / 2, -separacion / 2)]
centro = bpy.data.objects.new("Centro", None)
escena.collection.objects.link(centro)
centro.location = (0, 0.4, 0.6)  # la cámara apunta un poco por encima: la escena baja en el cuadro

bloques, textos_tmax = [], []
for c in range(N_CASOS):
    x, y = posiciones[c]
    obj, me, atributo = malla(f"Bloque {ETIQUETAS[c]}", VERTICES, CARAS, termico, rgba(COLORES[c, 0]))
    obj.location = (x, y, 0)
    bloques.append((me, atributo))
    for cuerpo, z, tam, mat in ((ETIQUETAS[c], 2.0, 0.36, blanco), ("", 1.5, 0.32, gris)):
        etiqueta = texto(cuerpo, tam, mat)
        etiqueta.location = (x, y, LZ / 2 + z)
        restriccion = etiqueta.constraints.new("TRACK_TO")
        restriccion.target, restriccion.track_axis, restriccion.up_axis = centro, "TRACK_NEGATIVE_Z", "UP_Y"
        if not cuerpo:
            textos_tmax.append(etiqueta)

# Mesa oscura y brillante (refleja los bloques) y luces.
bpy.ops.mesh.primitive_plane_add(size=80, location=(0, 0, -LZ / 2 - 0.003))
mesa = bpy.context.active_object
mat_mesa = bpy.data.materials.new("Mesa")
mat_mesa.use_nodes = True
bsdf = mat_mesa.node_tree.nodes["Principled BSDF"]
bsdf.inputs["Base Color"].default_value = (0.010, 0.011, 0.015, 1)
bsdf.inputs["Roughness"].default_value = 0.32
mesa.data.materials.append(mat_mesa)

for ubicacion, potencia, tamano in (((-2, 3, 12), 500, 10), ((10, -9, 6), 900, 8), ((-9, 7, 5), 300, 6)):
    luz = bpy.data.lights.new("Luz", type="AREA")
    luz.energy, luz.size = potencia, tamano
    objeto_luz = bpy.data.objects.new("Luz", luz)
    escena.collection.objects.link(objeto_luz)
    objeto_luz.location = ubicacion
    restriccion = objeto_luz.constraints.new("TRACK_TO")
    restriccion.target, restriccion.track_axis, restriccion.up_axis = centro, "TRACK_NEGATIVE_Z", "UP_Y"

# Cámara con paneo lento desde el frente derecho (se ve el corte de cada bloque).
datos_camara = bpy.data.cameras.new("Camara")
datos_camara.lens = 40
camara = bpy.data.objects.new("Camara", datos_camara)
escena.collection.objects.link(camara)
escena.camera = camara
restriccion = camara.constraints.new("TRACK_TO")
restriccion.target, restriccion.track_axis, restriccion.up_axis = centro, "TRACK_NEGATIVE_Z", "UP_Y"
# Las etiquetas de los bloques miran a la cámara.
for obj in escena.objects:
    if obj.type == "FONT":
        obj.constraints[0].target = camara
        obj.constraints[0].track_axis = "TRACK_Z"


def ubicar_camara(fraccion):
    angulo = math.radians(-78 + 56 * fraccion)  # de -78° a -22°: siempre del lado del corte
    radio, altura = 15.5, 9.5
    camara.location = (radio * math.cos(angulo), radio * math.sin(angulo), altura)


# HUD: textos y barra de colores pegados a la cámara (a 1 m delante del lente).
medio_ancho = 18.0 / datos_camara.lens  # sensor de 36 mm
medio_alto = medio_ancho * ALTO / ANCHO


def hud(obj, x, y):
    obj.parent = camara
    obj.location = (x * medio_ancho, y * medio_alto, -1.0)
    obj.constraints.clear()
    return obj


hud(texto("Mismo calor, distinto borde", 0.032, blanco, "LEFT"), -0.95, 0.88)
medidas = f"{LX:g} × {LY:g} × {LZ:g}".replace(".", ",")
hud(texto(f"Bloque {medidas} m de hormigón · pulso Q(t) = Qmax/2·(1 − cos 2πt/T), "
          f"Qmax = {float(d['Q_max']):g} W/m³, T = {float(d['T_h']):g} h · cuarto cortado para ver el interior",
          0.0125, gris, "LEFT"), -0.95, 0.78)
reloj = hud(texto("", 0.036, blanco, "LEFT"), -0.95, -0.80)
texto_Q = hud(texto("", 0.016, gris, "LEFT"), -0.95, -0.92)
hud(texto("Solver propio (volúmenes finitos 3D, Euler explícito) · Grupo 4", 0.011, gris, "RIGHT"), 0.95, -0.94)

# Barra de colores vertical a la derecha (medidas en metros a 1 m del lente).
barra = d["barra"]
n = len(barra)
alto_barra, ancho_barra = 1.1 * medio_alto, 0.012
vertices_barra, caras_barra, colores_barra = [], [], []
for i in range(n):
    yy = -alto_barra / 2 + alto_barra * i / (n - 1)
    vertices_barra += [(0, yy, 0), (ancho_barra, yy, 0)]
    colores_barra += [barra[i], barra[i]]
    if i:
        a = 2 * (i - 1)
        caras_barra.append((a, a + 1, a + 3, a + 2))
obj_barra, _, _ = malla("Escala", np.array(vertices_barra), np.array(caras_barra),
                        material_emision(atributo="T_color"), rgba(np.array(colores_barra)))
obj_barra.parent = camara
obj_barra.location = (0.86 * medio_ancho, 0.0, -1.0)
for valor in np.linspace(T_MIN, T_MAX, 5):
    yy = -alto_barra / 2 + alto_barra * (valor - T_MIN) / (T_MAX - T_MIN)
    marca = texto(f"{valor:.0f} °C", 0.014, gris, "LEFT")
    marca.parent = camara
    marca.location = (0.86 * medio_ancho + ancho_barra + 0.008, yy, -1.0)
    marca.constraints.clear()

# ---------------------------------------------------------------------------
# Render cuadro a cuadro
# ---------------------------------------------------------------------------
indices = list(range(0, len(T), CADA))
for k, i in enumerate(indices):
    for c, (me, atributo) in enumerate(bloques):
        atributo.data.foreach_set("color", rgba(COLORES[c, i]).ravel())
        me.update()
        textos_tmax[c].data.body = f"Tmax {TMAX[c, i]:.1f} °C".replace(".", ",")
    reloj.data.body = f"t = {T[i] / 3600:.0f} h"
    texto_Q.data.body = (f"Q = {Q[i]:.0f} W/m³" + ("  (pulso apagado)" if Q[i] == 0 and T[i] > 0 else ""))
    ubicar_camara(k / max(1, len(indices) - 1))
    escena.render.filepath = os.path.join(CARPETA, f"cuadro_{k:04d}.png")
    bpy.ops.render.render(write_still=True)
    if not os.path.exists(escena.render.filepath):
        raise RuntimeError(f"Blender no guardó {escena.render.filepath}")
    print(f"cuadro {k + 1}/{len(indices)}", flush=True)

# ---------------------------------------------------------------------------
# MP4 con el editor de video (no vuelve a renderizar: junta los PNG)
# ---------------------------------------------------------------------------
if RUTA_MP4:
    video = bpy.data.scenes.new("Video")
    editor = video.sequence_editor_create()
    archivos = [f"cuadro_{k:04d}.png" for k in range(len(indices))]
    archivos += [archivos[-1]] * (2 * FPS)  # el último cuadro se queda 2 s
    franja = editor.strips.new_image("Cuadros", os.path.join(CARPETA, archivos[0]), 1, 1)
    for archivo in archivos[1:]:
        franja.elements.append(archivo)
    video.frame_start, video.frame_end = 1, len(archivos)
    video.render.fps = FPS
    video.render.resolution_x, video.render.resolution_y = ANCHO, ALTO
    if hasattr(video.render.image_settings, "media_type"):  # Blender 5: primero se declara que es video
        video.render.image_settings.media_type = "VIDEO"
    video.render.image_settings.file_format = "FFMPEG"
    video.render.ffmpeg.format = "MPEG4"
    video.render.ffmpeg.codec = "H264"
    video.render.ffmpeg.constant_rate_factor = "HIGH"
    video.render.filepath = os.path.abspath(RUTA_MP4)
    bpy.ops.render.render(animation=True, scene=video.name)
    print("MP4:", os.path.abspath(RUTA_MP4), flush=True)
