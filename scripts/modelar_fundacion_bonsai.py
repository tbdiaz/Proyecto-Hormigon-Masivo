"""Modela la fundación en Blender + Bonsai siguiendo la guía del curso (secciones 8 a 17).

Uso, desde la raíz del repositorio (no abre la interfaz gráfica):

    blender -b --python scripts/modelar_fundacion_bonsai.py -- geometry/course_geometry.json geometry/foundation.ifc

Cada paso es el mismo que se hace a mano en la interfaz:

    File > New IFC Project > New Metric (m) Project       -> bim.new_project(preset="metric_m")
    Add > Mesh > Cylinder / Cone  (F9: Vertices, Radius, Depth)
    N > Item: Location Z = z0 + h/2, Rotation Z           -> Blender ubica la primitiva por su centro
    Object Information > IFC Class = IfcFooting > Assign   -> bim.assign_class
    Attributes > Name = Vertida_XX                         -> atributo IFC Name
    File > Save IFC Project                                -> bim.save_project

Las propiedades del curso (UANDES_HormigonMasivo) se agregan después, fuera de
Blender, con:  mcplan ifc-propiedades geometry/foundation.ifc
"""

import json
import os
import sys
from math import radians

import bonsai.tool as tool
import bpy
import ifcopenshell.api.attribute

ruta_geometria, ruta_salida = sys.argv[sys.argv.index("--") + 1:][:2]
with open(ruta_geometria, encoding="utf-8") as f:
    geometria = json.load(f)

bpy.ops.bim.new_project(preset="metric_m")

for v in geometria["vertidas"]:
    alto = v["z1"] - v["z0"]
    centro = (0.0, 0.0, v["z0"] + alto / 2.0)
    rotacion = (0.0, 0.0, radians(v["rotacion_z_deg"]))
    if v["primitiva"] == "cilindro":
        bpy.ops.mesh.primitive_cylinder_add(
            vertices=v["lados"], radius=v["radio1"], depth=alto, location=centro, rotation=rotacion)
    else:
        bpy.ops.mesh.primitive_cone_add(
            vertices=v["lados"], radius1=v["radio1"], radius2=v["radio2"], depth=alto,
            location=centro, rotation=rotacion)

    obj = bpy.context.active_object
    obj.name = v["nombre"]
    bpy.ops.bim.assign_class(obj=obj.name, ifc_class="IfcFooting")

    # Nombre IFC explícito (no basta con el nombre del objeto en el Outliner).
    elemento = tool.Ifc.get_entity(obj)
    ifcopenshell.api.attribute.edit_attributes(tool.Ifc.get(), product=elemento, attributes={"Name": v["nombre"]})
    print(f"OK {v['nombre']}: {elemento.is_a()} {elemento.GlobalId}")

bpy.ops.bim.save_project(filepath=os.path.abspath(ruta_salida), should_save_as=True)
print(f"Guardado: {ruta_salida}")
