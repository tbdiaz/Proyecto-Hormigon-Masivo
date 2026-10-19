"""Lectura y creación del modelo IFC con IfcOpenShell.

Flujo del proyecto:  Bonsai/Blender -> IFC -> IfcOpenShell -> geometría + propiedades -> solver

El solver nunca lee Blender: sólo lee el archivo .ifc (formato abierto). Por eso
cualquier cambio hecho en Bonsai y guardado en el IFC se refleja aquí.

Cada vertida es un elemento IFC con el property set MCPlan_Vertida. No usamos el
prefijo "Pset_" porque buildingSMART lo reserva para los property sets oficiales.
"""

import json
from pathlib import Path

import numpy as np

PSET_VERTIDA = "MCPlan_Vertida"
PROPIEDADES_VERTIDA = (
    "VertidaID",
    "TipoHormigon",
    "InicioPlanificado",
    "TemperaturaColocacion_C",
    "Moldaje",
    "Aislacion",
    "RetiroMoldaje_h",
)

# Coordenada oficial del Grupo 4 (Sitios de construcción por grupo).
SITIO_G4 = {"nombre": "Litueche", "latitud": -34.105457, "longitud": -71.649539}

# Geometría base. Fuente: Ju y Lei (2019): planta 18 m x 18 m, espesor 3,30 m,
# parte superior cilíndrica de 7 m de diámetro.
# SUPUESTO PENDIENTE DE VALIDAR con la geometría oficial del curso: el espesor
# total en el eje es 3,30 m y se divide en 3 vertidas de 1,10 m (el paper vierte
# en tres capas y ubica sensores cada 1,10 m). La tercera vertida es el pedestal.
GEOMETRIA_BASE = {"lado_losa": 18.0, "diametro_pedestal": 7.0, "espesor_total": 3.30}

# Valores iniciales de la planificación: son decisiones del plan, no datos.
VERTIDAS_BASE = [
    {"nombre": "V1 - Losa inferior", "forma": "losa", "z0": 0.00, "z1": 1.10,
     "props": {"VertidaID": "V1", "InicioPlanificado": "TBD+0h"}},
    {"nombre": "V2 - Losa superior", "forma": "losa", "z0": 1.10, "z1": 2.20,
     "props": {"VertidaID": "V2", "InicioPlanificado": "TBD+48h"}},
    {"nombre": "V3 - Pedestal", "forma": "pedestal", "z0": 2.20, "z1": 3.30,
     "props": {"VertidaID": "V3", "InicioPlanificado": "TBD+96h"}},
]
PROPS_COMUNES = {
    "TipoHormigon": "C35 (Ju y Lei 2019)",
    "TemperaturaColocacion_C": 20.0,
    "Moldaje": "Metalico",
    "Aislacion": "Ninguna",
    "RetiroMoldaje_h": 72.0,
}


# --------------------------------------------------------------------------- #
# Creación del modelo base
# --------------------------------------------------------------------------- #
def crear_modelo_base(ruta, geometria=GEOMETRIA_BASE, vertidas=VERTIDAS_BASE, sitio=SITIO_G4):
    """Escribe un IFC4 con una fundación dividida en vertidas (unidades en metros).

    Sirve como punto de partida: después se abre y edita en Bonsai.
    """
    import ifcopenshell.util.geolocation as geoloc
    from ifcopenshell import api
    from ifcopenshell.api import aggregate, context, geometry, material, project, pset, root, spatial, unit  # noqa: F401

    f = api.project.create_file(version="IFC4")
    proyecto = api.root.create_entity(f, ifc_class="IfcProject", name="Fundacion aerogenerador - Grupo 4")
    unidades = [api.unit.add_si_unit(f, unit_type=u) for u in ("LENGTHUNIT", "AREAUNIT", "VOLUMEUNIT")]
    api.unit.assign_unit(f, units=unidades)
    modelo = api.context.add_context(f, context_type="Model")
    cuerpo = api.context.add_context(
        f, context_type="Model", context_identifier="Body", target_view="MODEL_VIEW", parent=modelo
    )

    site = api.root.create_entity(f, ifc_class="IfcSite", name=sitio["nombre"])
    site.RefLatitude = geoloc.dd2dms(sitio["latitud"], use_us=True)
    site.RefLongitude = geoloc.dd2dms(sitio["longitud"], use_us=True)
    edificio = api.root.create_entity(f, ifc_class="IfcBuilding", name="Aerogenerador")
    nivel = api.root.create_entity(f, ifc_class="IfcBuildingStorey", name="Fundacion")
    api.aggregate.assign_object(f, products=[site], relating_object=proyecto)
    api.aggregate.assign_object(f, products=[edificio], relating_object=site)
    api.aggregate.assign_object(f, products=[nivel], relating_object=edificio)

    hormigon = api.material.add_material(f, name="Hormigon C35", category="concrete")
    origen_2d = f.createIfcAxis2Placement2D(f.createIfcCartesianPoint((0.0, 0.0)), None)

    for v in vertidas:
        if v["forma"] == "losa":
            lado = geometria["lado_losa"]
            perfil = f.createIfcRectangleProfileDef("AREA", None, origen_2d, lado, lado)
        else:
            radio = geometria["diametro_pedestal"] / 2.0
            perfil = f.createIfcCircleProfileDef("AREA", None, origen_2d, radio)
        alto = v["z1"] - v["z0"]
        solido = f.createIfcExtrudedAreaSolid(
            perfil,
            f.createIfcAxis2Placement3D(f.createIfcCartesianPoint((0.0, 0.0, 0.0)), None, None),
            f.createIfcDirection((0.0, 0.0, 1.0)),
            alto,
        )
        rep = f.createIfcShapeRepresentation(cuerpo, "Body", "SweptSolid", [solido])

        elem = api.root.create_entity(f, ifc_class="IfcFooting", predefined_type="PAD_FOOTING", name=v["nombre"])
        matriz = np.eye(4)
        matriz[2, 3] = v["z0"]  # el elemento se apoya en su cota inferior
        api.geometry.edit_object_placement(f, product=elem, matrix=matriz)
        api.geometry.assign_representation(f, product=elem, representation=rep)
        api.spatial.assign_container(f, products=[elem], relating_structure=nivel)
        api.material.assign_material(f, products=[elem], type="IfcMaterial", material=hormigon)

        pset = api.pset.add_pset(f, product=elem, name=PSET_VERTIDA)
        api.pset.edit_pset(f, pset=pset, properties={**PROPS_COMUNES, **v["props"]})

    ruta = Path(ruta)
    ruta.parent.mkdir(parents=True, exist_ok=True)
    f.write(str(ruta))
    return ruta


# --------------------------------------------------------------------------- #
# Lectura
# --------------------------------------------------------------------------- #
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
    ajustes = ifcopenshell.geom.settings()
    ajustes.set("use-world-coords", True)  # coordenadas globales (aplica el placement)

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
    lineas = [f"Archivo: {resumen['archivo']}  (schema {resumen['schema']})",
              f"Elementos con geometría: {len(resumen['elementos'])}", ""]
    for gid, e in resumen["elementos"].items():
        dims = np.subtract(e["bbox_max"], e["bbox_min"])
        lineas.append(f"- {e['nombre']}  [{e['clase']}]  GlobalId={gid}")
        lineas.append(f"    bbox: {dims[0]:.2f} x {dims[1]:.2f} x {dims[2]:.2f} m, "
                      f"z = {e['bbox_min'][2]:.2f} .. {e['bbox_max'][2]:.2f} m, volumen = {e['volumen_m3']:.2f} m3")
        for pset, props in e["psets"].items():
            lineas.append(f"    {pset}: " + ", ".join(f"{k}={v}" for k, v in props.items()))
    return "\n".join(lineas)


# --------------------------------------------------------------------------- #
# Del IFC al solver 1D
# --------------------------------------------------------------------------- #
def columna_desde_ifc(ruta, x: float = 0.0, y: float = 0.0, tol: float = 1e-3) -> list[dict]:
    """Vertidas que atraviesa la vertical (x, y), ordenadas de abajo hacia arriba.

    Es la "voxelización" más simple posible: un único rayo vertical. Usa el
    bounding box de cada elemento, lo que es exacto para las losas
    rectangulares y para el eje del pedestal cilíndrico.
    Devuelve [{"nombre", "VertidaID", "z0", "z1", "T_colocacion"}, ...].
    """
    resumen = resumir_ifc(ruta)
    columna = []
    for e in resumen["elementos"].values():
        props = e["psets"].get(PSET_VERTIDA)
        if props is None:
            continue
        (x0, y0, z0), (x1, y1, z1) = e["bbox_min"], e["bbox_max"]
        if x0 - tol <= x <= x1 + tol and y0 - tol <= y <= y1 + tol:
            columna.append({
                "nombre": e["nombre"],
                "VertidaID": props.get("VertidaID"),
                "z0": z0,
                "z1": z1,
                "T_colocacion": float(props["TemperaturaColocacion_C"]),
            })
    columna.sort(key=lambda c: c["z0"])
    if not columna:
        raise ValueError(f"Ninguna vertida con {PSET_VERTIDA} pasa por (x={x}, y={y})")
    for abajo, arriba in zip(columna, columna[1:]):
        if abs(arriba["z0"] - abajo["z1"]) > tol:
            raise ValueError(f"Hueco o traslape entre {abajo['nombre']} y {arriba['nombre']} "
                             f"(z = {abajo['z1']:.3f} vs {arriba['z0']:.3f})")
    return columna
