"""Capturas de evidencia del modelo BIM (guía del curso, sección 32): 3D, superior, lateral y explotada.

Uso, desde la raíz del repositorio:

    blender -b --python scripts/capturas_fundacion_bonsai.py -- geometry/foundation.ifc docs/img

Abre el IFC con Bonsai (File > Open IFC Project) y renderiza con el motor
Workbench, con un color fijo por vertida para distinguirlas. Los colores son
sólo de visualización (propiedad del objeto en Blender): no modifica el IFC.
"""

import os
import sys
from math import radians

import bpy
from mathutils import Vector

ruta_ifc, carpeta = sys.argv[sys.argv.index("--") + 1:][:2]
os.makedirs(carpeta, exist_ok=True)

COLORES = [(0.85, 0.60, 0.30, 1), (0.30, 0.55, 0.80, 1), (0.45, 0.72, 0.40, 1),
           (0.82, 0.42, 0.45, 1), (0.58, 0.47, 0.78, 1)]

bpy.ops.bim.load_project(filepath=os.path.abspath(ruta_ifc))
vertidas = sorted((o for o in bpy.data.objects if "Vertida_" in o.name), key=lambda o: o.name)
for obj, color in zip(vertidas, COLORES):
    obj.color = color

escena = bpy.context.scene
escena.render.engine = "BLENDER_WORKBENCH"
escena.render.resolution_x, escena.render.resolution_y = 1400, 900
sombreado = escena.display.shading
sombreado.light = "STUDIO"
sombreado.color_type = "OBJECT"
sombreado.show_object_outline = True
sombreado.show_cavity = True

camara_datos = bpy.data.cameras.new("Camara")
camara = bpy.data.objects.new("Camara", camara_datos)
escena.collection.objects.link(camara)
escena.camera = camara


def mirar(desde, hacia):
    camara.location = Vector(desde)
    camara.rotation_euler = (Vector(hacia) - Vector(desde)).to_track_quat("-Z", "Y").to_euler()


def capturar(nombre):
    escena.render.filepath = os.path.abspath(os.path.join(carpeta, nombre))
    bpy.ops.render.render(write_still=True)
    print("Captura:", escena.render.filepath)


# 3D en perspectiva
camara_datos.type = "PERSP"
camara_datos.lens = 35
mirar((24.0, -26.0, 17.0), (0.0, 0.0, 1.2))
capturar("fundacion_3d.png")

# Vista superior (Numpad 7) y lateral (Numpad 1), ortográficas.
# ortho_scale es el ancho visible: la superior necesita ~22 m de alto visible.
camara_datos.type = "ORTHO"
camara_datos.ortho_scale = 22.0 * 1400 / 900
camara.location, camara.rotation_euler = (0.0, 0.0, 40.0), (0.0, 0.0, 0.0)
capturar("fundacion_superior.png")
camara_datos.ortho_scale = 21.0
camara.location, camara.rotation_euler = (0.0, -40.0, 2.15), (radians(90.0), 0.0, 0.0)
capturar("fundacion_lateral.png")

# Vista explotada con etiquetas: cada vertida se separa verticalmente
# (sólo para la imagen; nada se guarda).
camara_datos.type = "PERSP"
mirar((26.0, -28.0, 16.0), (0.0, 0.0, 3.5))
for i, obj in enumerate(vertidas):
    obj.location.z += 1.2 * i
bpy.context.view_layer.update()  # recalcula matrix_world con las nuevas posiciones
for i, obj in enumerate(vertidas):
    esquinas = [obj.matrix_world @ Vector(c) for c in obj.bound_box]
    z_medio = sum(c.z for c in esquinas) / 8.0
    x_borde = max(c.x for c in esquinas)
    texto = bpy.data.objects.new(f"Etiqueta {obj.name}", bpy.data.curves.new(obj.name, type="FONT"))
    texto.data.body = obj.name.split("/")[-1]
    texto.data.size = 0.9
    texto.data.align_y = "CENTER"
    texto.location = (x_borde + 0.8, 0.0, z_medio)
    texto.rotation_euler = camara.rotation_euler
    texto.color = COLORES[i]
    escena.collection.objects.link(texto)
capturar("fundacion_explotada.png")
