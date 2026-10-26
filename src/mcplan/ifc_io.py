"""Lectura y creación del modelo IFC con IfcOpenShell.

Flujo del proyecto:  Bonsai/Blender -> IFC -> IfcOpenShell -> geometría + propiedades -> solver

El solver nunca lee Blender: sólo lee el archivo .ifc (formato abierto). Por eso
cualquier cambio hecho en Bonsai y guardado en el IFC se refleja aquí.

Convenciones del curso (guía "Modelar la fundación en Blender + Bonsai"):
- cada vertida es un elemento IFC (IfcFooting) cuyo atributo Name es "Vertida_XX";
- sus datos constructivos van en el property set UANDES_HormigonMasivo. No se usa
  el prefijo "Pset_" porque buildingSMART lo reserva para los property sets oficiales.
"""

import json
from dataclasses import dataclass, field
from math import cos, pi, sin
from pathlib import Path

import numpy as np

PSET_VERTIDA = "UANDES_HormigonMasivo"
PREFIJO_VERTIDA = "Vertida_"
PROPIEDADES_VERTIDA = (
    "VertidaID",
    "Inicio_h",
    "TemperaturaColocacion_C",
    "TipoHormigon",
    "Moldaje",
    "Aislacion",
    "RetiroMoldaje_h",
)

# Valores iniciales de la guía del curso (sección 24). Son el punto de partida del
# plan, no una recomendación: el grupo los cambiará al estudiar alternativas.
PROPIEDADES_INICIALES = {
    f"Vertida_0{i}": {"VertidaID": i, "Inicio_h": 6.0 * (i - 1), "TemperaturaColocacion_C": 18.0,
                      "TipoHormigon": "M1", "Moldaje": "TBD", "Aislacion": "TBD", "RetiroMoldaje_h": 72.0}
    for i in range(1, 6)
}


# --------------------------------------------------------------------------- #
# Propiedades del curso
# --------------------------------------------------------------------------- #
def _asignar_propiedades(f, valores: dict) -> list[str]:
    """Crea o actualiza UANDES_HormigonMasivo en cada vertida nombrada en `valores`."""
    import ifcopenshell.api.pset
    import ifcopenshell.util.element

    productos = {p.Name: p for p in f.by_type("IfcProduct") if p.Name}
    asignadas = []
    for nombre, props in valores.items():
        producto = productos.get(nombre)
        if producto is None:
            continue
        existente = ifcopenshell.util.element.get_psets(producto).get(PSET_VERTIDA)
        if existente:
            pset = f.by_id(existente["id"])
        else:
            pset = ifcopenshell.api.pset.add_pset(f, product=producto, name=PSET_VERTIDA)
        ifcopenshell.api.pset.edit_pset(f, pset=pset, properties=props)
        asignadas.append(nombre)
    return asignadas


def agregar_propiedades_curso(entrada, salida=None, valores=PROPIEDADES_INICIALES) -> list[str]:
    """Agrega (o actualiza) las propiedades del curso a las vertidas de un IFC.

    Equivale al script add_course_properties.py de la guía. Si no se da `salida`,
    modifica el archivo de entrada.
    """
    import ifcopenshell

    f = ifcopenshell.open(str(entrada))
    asignadas = _asignar_propiedades(f, valores)
    f.write(str(salida or entrada))
    return asignadas


# --------------------------------------------------------------------------- #
# Creación del modelo con IfcOpenShell (sin Blender)
# --------------------------------------------------------------------------- #
def cargar_geometria(ruta) -> dict:
    return json.loads(Path(ruta).read_text(encoding="utf-8"))


def malla_prisma(lados, radio1, radio2, z0, z1, rotacion_z_deg=0.0):
    """Vértices y caras de un prisma (radio1 = radio2) o tronco (radio1 != radio2) de n lados.

    Caras orientadas hacia afuera: base, tapa y n caras laterales.
    """
    rot = rotacion_z_deg * pi / 180.0
    anillo = lambda r, z: [(r * cos(2 * pi * k / lados + rot), r * sin(2 * pi * k / lados + rot), z)  # noqa: E731
                           for k in range(lados)]
    vertices = anillo(radio1, z0) + anillo(radio2, z1)
    caras = [list(range(lados - 1, -1, -1)), list(range(lados, 2 * lados))]
    caras += [[k, (k + 1) % lados, lados + (k + 1) % lados, lados + k] for k in range(lados)]
    return vertices, caras


def crear_modelo_base(ruta, ruta_geometria, valores=PROPIEDADES_INICIALES):
    """Escribe un IFC4 de la fundación oficial directamente con IfcOpenShell.

    Produce la misma geometría que scripts/modelar_fundacion_bonsai.py pero sin
    Blender. Lo usan los tests, para no depender del archivo que el grupo edita.
    """
    from ifcopenshell import api
    from ifcopenshell.api import aggregate, context, geometry, project, root, spatial, unit  # noqa: F401

    f = api.project.create_file(version="IFC4")
    proyecto = api.root.create_entity(f, ifc_class="IfcProject", name="Fundacion aerogenerador")
    api.unit.assign_unit(f, units=[api.unit.add_si_unit(f, unit_type=u)
                                   for u in ("LENGTHUNIT", "AREAUNIT", "VOLUMEUNIT")])
    modelo = api.context.add_context(f, context_type="Model")
    cuerpo = api.context.add_context(
        f, context_type="Model", context_identifier="Body", target_view="MODEL_VIEW", parent=modelo)
    site = api.root.create_entity(f, ifc_class="IfcSite", name="Site")
    edificio = api.root.create_entity(f, ifc_class="IfcBuilding", name="Building")
    nivel = api.root.create_entity(f, ifc_class="IfcBuildingStorey", name="Storey")
    api.aggregate.assign_object(f, products=[site], relating_object=proyecto)
    api.aggregate.assign_object(f, products=[edificio], relating_object=site)
    api.aggregate.assign_object(f, products=[nivel], relating_object=edificio)

    for v in cargar_geometria(ruta_geometria)["vertidas"]:
        vertices, caras = malla_prisma(v["lados"], v["radio1"], v["radio2"], v["z0"], v["z1"], v["rotacion_z_deg"])
        rep = api.geometry.add_mesh_representation(f, context=cuerpo, vertices=[vertices], faces=[caras])
        elem = api.root.create_entity(f, ifc_class="IfcFooting", name=v["nombre"])
        api.geometry.edit_object_placement(f, product=elem)  # placement en el origen
        api.geometry.assign_representation(f, product=elem, representation=rep)
        api.spatial.assign_container(f, products=[elem], relating_structure=nivel)

    _asignar_propiedades(f, valores)
    ruta = Path(ruta)
    ruta.parent.mkdir(parents=True, exist_ok=True)
    f.write(str(ruta))
    return ruta


# --------------------------------------------------------------------------- #
# Lectura
# --------------------------------------------------------------------------- #
def _ajustes_geometria():
    import ifcopenshell.geom

    ajustes = ifcopenshell.geom.settings()
    ajustes.set("use-world-coords", True)  # coordenadas globales (aplica el placement)
    return ajustes


def resumir_ifc(ruta) -> dict:
    """Lee el IFC y devuelve, por elemento: clase, nombre, property sets, bbox y volumen.

    El resultado es un diccionario serializable a JSON, indexado por GlobalId
    (identificador estable: no cambia cuando se edita el elemento en Bonsai).
    """
    import ifcopenshell
    import ifcopenshell.geom
    import ifcopenshell.util.element
    import ifcopenshell.util.shape

    f = ifcopenshell.open(str(ruta))
    ajustes = _ajustes_geometria()
    elementos = {}
    for e in f.by_type("IfcElement"):
        if e.is_a("IfcFeatureElement") or not e.Representation:
            continue
        psets = {
            nombre: {k: v for k, v in props.items() if k != "id"}
            for nombre, props in ifcopenshell.util.element.get_psets(e).items()
        }
        forma = ifcopenshell.geom.create_shape(ajustes, e)
        vertices = ifcopenshell.util.shape.get_vertices(forma.geometry)
        elementos[e.GlobalId] = {
            "clase": e.is_a(),
            "nombre": e.Name,
            "psets": psets,
            "bbox_min": [round(float(x), 6) for x in vertices.min(axis=0)],
            "bbox_max": [round(float(x), 6) for x in vertices.max(axis=0)],
            "volumen_m3": round(float(ifcopenshell.util.shape.get_volume(forma.geometry)), 6),
        }
    return {"archivo": str(ruta), "schema": f.schema, "elementos": elementos}


@dataclass
class VertidaIFC:
    """Lo que el resto del programa necesita de una vertida, sin detalles del esquema IFC."""

    guid: str
    nombre: str
    vertices: np.ndarray  # (n, 3) en coordenadas globales [m]
    triangulos: np.ndarray  # (m, 3) índices de vértices
    propiedades: dict = field(default_factory=dict)  # contenido de UANDES_HormigonMasivo


def leer_vertidas(ruta) -> list[VertidaIFC]:
    """Vertidas del IFC (elementos con Name "Vertida_XX"), trianguladas, ordenadas por nombre."""
    import ifcopenshell
    import ifcopenshell.geom
    import ifcopenshell.util.element
    import ifcopenshell.util.shape

    f = ifcopenshell.open(str(ruta))
    ajustes = _ajustes_geometria()
    vertidas = []
    for p in f.by_type("IfcProduct"):
        if not (p.Name or "").startswith(PREFIJO_VERTIDA) or not p.Representation:
            continue
        forma = ifcopenshell.geom.create_shape(ajustes, p)
        props = ifcopenshell.util.element.get_psets(p).get(PSET_VERTIDA, {})
        vertidas.append(VertidaIFC(
            guid=p.GlobalId,
            nombre=p.Name,
            vertices=ifcopenshell.util.shape.get_vertices(forma.geometry),
            triangulos=ifcopenshell.util.shape.get_faces(forma.geometry),
            propiedades={k: v for k, v in props.items() if k != "id"},
        ))
    return sorted(vertidas, key=lambda v: v.nombre)


def extension_vertical(vertida: VertidaIFC, x: float, y: float, tol: float = 1e-9):
    """Tramo [z_min, z_max] que ocupa el sólido sobre la vertical (x, y), o None.

    Lanza un rayo vertical y busca en qué triángulos de la superficie cae (x, y)
    usando coordenadas baricéntricas en planta. Es exacto para sólidos convexos,
    como todas las vertidas de la geometría del curso.
    """
    tri = vertida.vertices[vertida.triangulos]  # (m, 3 vértices, 3 coordenadas)
    a, b, c = tri[:, 0], tri[:, 1], tri[:, 2]
    det = (b[:, 1] - c[:, 1]) * (a[:, 0] - c[:, 0]) + (c[:, 0] - b[:, 0]) * (a[:, 1] - c[:, 1])
    validos = np.abs(det) > 1e-12  # descarta caras verticales (su proyección no tiene área)
    det = np.where(validos, det, 1.0)
    l1 = ((b[:, 1] - c[:, 1]) * (x - c[:, 0]) + (c[:, 0] - b[:, 0]) * (y - c[:, 1])) / det
    l2 = ((c[:, 1] - a[:, 1]) * (x - c[:, 0]) + (a[:, 0] - c[:, 0]) * (y - c[:, 1])) / det
    l3 = 1.0 - l1 - l2
    dentro = validos & (l1 >= -tol) & (l2 >= -tol) & (l3 >= -tol)
    if not dentro.any():
        return None
    z = l1[dentro] * a[dentro, 2] + l2[dentro] * b[dentro, 2] + l3[dentro] * c[dentro, 2]
    return float(z.min()), float(z.max())


def columna_desde_ifc(ruta, x: float = 0.0, y: float = 0.0, tol: float = 1e-6) -> list[dict]:
    """Vertidas que atraviesa la vertical (x, y), ordenadas de abajo hacia arriba.

    Es una "voxelización" de una sola columna. Devuelve
    [{"nombre", "guid", "z0", "z1", "T_colocacion", "inicio_h"}, ...].
    """
    columna = []
    for v in leer_vertidas(ruta):
        tramo = extension_vertical(v, x, y)
        if tramo is None or tramo[1] - tramo[0] < tol:
            continue
        if "TemperaturaColocacion_C" not in v.propiedades:
            raise ValueError(f"{v.nombre} no tiene {PSET_VERTIDA}.TemperaturaColocacion_C. "
                             "Agrega las propiedades con: mcplan ifc-propiedades <archivo.ifc>")
        columna.append({
            "nombre": v.nombre,
            "guid": v.guid,
            "z0": tramo[0],
            "z1": tramo[1],
            "T_colocacion": float(v.propiedades["TemperaturaColocacion_C"]),
            "inicio_h": float(v.propiedades.get("Inicio_h", 0.0)),
        })
    columna.sort(key=lambda c: c["z0"])
    if not columna:
        raise ValueError(f"Ninguna vertida pasa por (x={x}, y={y})")
    for abajo, arriba in zip(columna, columna[1:]):
        if abs(arriba["z0"] - abajo["z1"]) > tol:
            raise ValueError(f"Hueco o traslape entre {abajo['nombre']} y {arriba['nombre']} "
                             f"(z = {abajo['z1']:.3f} vs {arriba['z0']:.3f})")
    return columna


# --------------------------------------------------------------------------- #
# Comparación de lecturas (para detectar cambios hechos en Bonsai)
# --------------------------------------------------------------------------- #
def guardar_resumen(resumen: dict, ruta) -> None:
    Path(ruta).parent.mkdir(parents=True, exist_ok=True)
    Path(ruta).write_text(json.dumps(resumen, indent=2, ensure_ascii=False), encoding="utf-8")


def cargar_resumen(ruta) -> dict:
    return json.loads(Path(ruta).read_text(encoding="utf-8"))


def comparar_resumenes(antes: dict, despues: dict) -> list[str]:
    """Lista legible de diferencias entre dos lecturas del IFC (vacía si no hay cambios)."""
    cambios = []
    ea, ed = antes["elementos"], despues["elementos"]
    for gid in sorted(ea.keys() - ed.keys()):
        cambios.append(f"ELIMINADO  {ea[gid]['nombre']} ({gid})")
    for gid in sorted(ed.keys() - ea.keys()):
        cambios.append(f"NUEVO      {ed[gid]['nombre']} ({gid})")
    for gid in sorted(ea.keys() & ed.keys()):
        a, d = ea[gid], ed[gid]
        etiqueta = d["nombre"] or gid
        for campo in ("clase", "nombre", "bbox_min", "bbox_max", "volumen_m3"):
            if a[campo] != d[campo]:
                cambios.append(f"MODIFICADO {etiqueta}: {campo} {a[campo]} -> {d[campo]}")
        for pset in sorted(a["psets"].keys() | d["psets"].keys()):
            pa, pd = a["psets"].get(pset, {}), d["psets"].get(pset, {})
            for prop in sorted(pa.keys() | pd.keys()):
                if pa.get(prop) != pd.get(prop):
                    cambios.append(f"MODIFICADO {etiqueta}: {pset}.{prop} {pa.get(prop)} -> {pd.get(prop)}")
    return cambios


def formatear_resumen(resumen: dict) -> str:
    elementos = resumen["elementos"]
    lineas = [f"Archivo: {resumen['archivo']}  (schema {resumen['schema']})",
              f"Elementos con geometría: {len(elementos)}", ""]
    for gid, e in sorted(elementos.items(), key=lambda kv: kv[1]["nombre"] or ""):
        dims = np.subtract(e["bbox_max"], e["bbox_min"])
        lineas.append(f"- {e['nombre']}  [{e['clase']}]  GlobalId={gid}")
        lineas.append(f"    bbox: {dims[0]:.2f} x {dims[1]:.2f} x {dims[2]:.2f} m, "
                      f"z = {e['bbox_min'][2]:.2f} .. {e['bbox_max'][2]:.2f} m, volumen = {e['volumen_m3']:.2f} m3")
        for pset, props in e["psets"].items():
            lineas.append(f"    {pset}: " + ", ".join(f"{k}={v}" for k, v in props.items()))
    if elementos:
        lo = np.min([e["bbox_min"] for e in elementos.values()], axis=0)
        hi = np.max([e["bbox_max"] for e in elementos.values()], axis=0)
        n_vertidas = sum((e["nombre"] or "").startswith(PREFIJO_VERTIDA) for e in elementos.values())
        lineas += ["", f"GLOBAL  bbox: {hi[0] - lo[0]:.2f} x {hi[1] - lo[1]:.2f} x {hi[2] - lo[2]:.2f} m, "
                       f"volumen total = {sum(e['volumen_m3'] for e in elementos.values()):.2f} m3",
                   f"Número de vertidas: {n_vertidas}"]
    return "\n".join(lineas)
