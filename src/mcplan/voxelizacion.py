"""Voxelización del IFC: de las vertidas (mallas trianguladas) a una grilla de celdas cúbicas.

Criterio de pertenencia: un vóxel pertenece a una vertida si su CENTRO está dentro
del sólido. Para decidirlo se lanza un rayo vertical desde el centro hacia +z y se
cuentan los cruces con la superficie (triángulos): número impar = adentro (regla
de paridad, teorema de la curva de Jordan). Vale para cualquier sólido cerrado,
convexo o no.

Se calcula por columnas: para cada (x, y) de la grilla se buscan los triángulos
cuya proyección en planta contiene el punto (coordenadas baricéntricas) y la altura
z del cruce; cada centro cuenta cuántos cruces tiene por encima. El rayo se corre
una fracción ínfima de celda para que nunca pase justo por una arista o un vértice
(o por un centro que caiga exactamente en una junta entre vertidas).

La frontera queda "en escalera" (staircase): las caras inclinadas de los troncos
se representan con escalones de tamaño dx, así que el volumen y la posición de la
superficie tienen un error del orden de dx que baja al refinar la malla.
"""

from dataclasses import dataclass, field
from math import ceil, floor

import numpy as np

from .ifc_io import VertidaIFC, leer_vertidas

VACIO = 0  # vertida_id de una celda sin hormigón (aire o fuera del modelo)
_EPS = (1.23456789e-7, 2.71828183e-7, 3.14159265e-7)  # corrimiento del rayo, en fracciones de dx


def volumen_malla(vertices, triangulos) -> float:
    """Volumen exacto de una malla cerrada (teorema de la divergencia: suma de tetraedros con el origen)."""
    t = np.asarray(vertices, dtype=float)[np.asarray(triangulos)]
    return float(abs(np.einsum("ij,ij->i", t[:, 0], np.cross(t[:, 1], t[:, 2])).sum()) / 6.0)


def _cruces_por_columna(vertices, triangulos, xs, ys):
    """Cruces de los rayos verticales (una columna por cada (x, y) de la grilla) con la malla.

    Devuelve (indice_columna, z) de cada cruce; índice de columna = i * ny + j.
    """
    tri = np.asarray(vertices, dtype=float)[np.asarray(triangulos)]
    columnas, alturas = [], []
    for a, b, c in tri:
        det = (b[1] - c[1]) * (a[0] - c[0]) + (c[0] - b[0]) * (a[1] - c[1])
        if abs(det) < 1e-14:
            continue  # cara vertical: su proyección no tiene área y un rayo vertical no la cruza
        x0, x1 = min(a[0], b[0], c[0]), max(a[0], b[0], c[0])
        y0, y1 = min(a[1], b[1], c[1]), max(a[1], b[1], c[1])
        i = np.nonzero((xs >= x0) & (xs <= x1))[0]
        j = np.nonzero((ys >= y0) & (ys <= y1))[0]
        if len(i) == 0 or len(j) == 0:
            continue
        X, Y = np.meshgrid(xs[i], ys[j], indexing="ij")
        l1 = ((b[1] - c[1]) * (X - c[0]) + (c[0] - b[0]) * (Y - c[1])) / det
        l2 = ((c[1] - a[1]) * (X - c[0]) + (a[0] - c[0]) * (Y - c[1])) / det
        l3 = 1.0 - l1 - l2
        dentro = (l1 >= 0) & (l2 >= 0) & (l3 >= 0)
        if not dentro.any():
            continue
        I, J = np.meshgrid(i, j, indexing="ij")
        columnas.append((I * len(ys) + J)[dentro])
        alturas.append((l1 * a[2] + l2 * b[2] + l3 * c[2])[dentro])
    if not columnas:
        return np.zeros(0, int), np.zeros(0)
    return np.concatenate(columnas), np.concatenate(alturas)


def dentro_del_solido(vertices, triangulos, xs, ys, zs) -> np.ndarray:
    """Máscara (nx, ny, nz): True si el punto (xs[i], ys[j], zs[k]) está dentro de la malla cerrada.

    Paridad de los cruces del rayo vertical que sube desde cada punto. Para no
    recorrer cada vóxel, se usa un arreglo de diferencias: cada cruce a altura zc
    suma 1 a todos los centros de su columna que están por debajo de zc.
    """
    nx, ny, nz = len(xs), len(ys), len(zs)
    columna, zc = _cruces_por_columna(vertices, triangulos, xs, ys)
    debajo = np.searchsorted(zs, zc)  # cuántos centros quedan bajo cada cruce
    diferencias = np.zeros((nx * ny, nz + 1), dtype=np.int32)
    np.add.at(diferencias, (columna, np.zeros_like(columna)), 1)
    np.add.at(diferencias, (columna, debajo), -1)
    cruces_arriba = np.cumsum(diferencias, axis=1)[:, :nz]
    return (cruces_arriba % 2 == 1).reshape(nx, ny, nz)


@dataclass
class Voxelizacion:
    """Grilla de celdas cúbicas con la vertida de cada celda."""

    dx: float
    origen: np.ndarray  # esquina inferior de la grilla [m]
    vertida_id: np.ndarray  # (nx, ny, nz) int8: 0 = vacío, 1..n = vertida
    nombres: list  # nombres[v - 1] = nombre de la vertida v
    propiedades: list  # propiedades[v - 1] = pset UANDES_HormigonMasivo
    volumen_ifc: list  # volumen exacto de cada malla [m3]
    superposiciones: int = 0  # celdas que caían en dos vertidas (debería ser 0)
    info: dict = field(default_factory=dict)

    @property
    def forma(self):
        return self.vertida_id.shape

    @property
    def n_vertidas(self) -> int:
        return len(self.nombres)

    def centros(self, eje: int) -> np.ndarray:
        return self.origen[eje] + (np.arange(self.forma[eje]) + 0.5) * self.dx

    def celdas(self, v: int) -> int:
        return int(np.count_nonzero(self.vertida_id == v))

    def volumen_voxel(self, v: int) -> float:
        return self.celdas(v) * self.dx**3

    def error_volumen(self, v: int) -> float:
        """Error relativo del volumen voxelizado de la vertida v respecto del IFC."""
        return self.volumen_voxel(v) / self.volumen_ifc[v - 1] - 1.0

    @property
    def error_volumen_total(self) -> float:
        return (np.count_nonzero(self.vertida_id) * self.dx**3) / sum(self.volumen_ifc) - 1.0


def voxelizar(vertidas: list[VertidaIFC], dx: float, margen: int = 1) -> Voxelizacion:
    """Voxeliza las vertidas en una grilla de celdas de lado dx.

    La grilla se alinea con múltiplos de dx (la base z = 0 cae justo en una cara)
    y deja `margen` celdas vacías alrededor, para que toda superficie exterior
    tenga una vecina vacía.
    """
    if dx <= 0:
        raise ValueError("dx debe ser positivo")
    todos = np.vstack([v.vertices for v in vertidas])
    lo = np.array([floor(c / dx + 1e-9) - margen for c in todos.min(axis=0)]) * dx
    hi = np.array([ceil(c / dx - 1e-9) + margen for c in todos.max(axis=0)]) * dx
    n = np.round((hi - lo) / dx).astype(int)
    xs, ys, zs = (lo[e] + (np.arange(n[e]) + 0.5 + _EPS[e]) * dx for e in range(3))

    vertida_id = np.zeros(tuple(n), dtype=np.int8)
    superposiciones = 0
    for numero, v in enumerate(vertidas, start=1):
        dentro = dentro_del_solido(v.vertices, v.triangulos, xs, ys, zs)
        superposiciones += int(np.count_nonzero(dentro & (vertida_id != VACIO)))
        vertida_id[dentro & (vertida_id == VACIO)] = numero
    return Voxelizacion(
        dx=dx, origen=lo, vertida_id=vertida_id, nombres=[v.nombre for v in vertidas],
        propiedades=[dict(v.propiedades) for v in vertidas],
        volumen_ifc=[volumen_malla(v.vertices, v.triangulos) for v in vertidas],
        superposiciones=superposiciones,
    )


def voxelizar_ifc(ruta, dx: float) -> Voxelizacion:
    return voxelizar(leer_vertidas(ruta), dx)


# ---------------------------------------------------------------------------
# Superficies: exteriores e interfaces entre vertidas
# ---------------------------------------------------------------------------

def caras(vertida_id: np.ndarray, activas=None) -> dict:
    """Cuenta las caras de las celdas según lo que hay del otro lado.

    activas: máscara de celdas que existen (por defecto, todas las con hormigón).
    Una cara es EXTERIOR si separa una celda activa de una vacía o inactiva
    (ahí va la condición de borde: Robin con el aire, o el suelo en la base);
    es INTERFAZ si separa celdas activas de vertidas distintas (junta entre
    vertidas: el calor pasa por conducción, no hay condición de borde).
    Devuelve el número de caras por tipo y dirección.
    """
    activas = (vertida_id != VACIO) if activas is None else activas
    ids = np.where(activas, vertida_id, VACIO)
    resultado = {"exterior": {}, "interfaz": {}}
    nombres_eje = ("x", "y", "z")
    for eje in range(3):
        a = np.moveaxis(ids, eje, 0)
        relleno = np.pad(a, [(1, 1)] + [(0, 0)] * 2)
        abajo, arriba = relleno[:-1], relleno[1:]  # pares de celdas vecinas a lo largo del eje
        sale_menos = (arriba != VACIO) & (abajo == VACIO)  # cara "-eje" expuesta
        sale_mas = (abajo != VACIO) & (arriba == VACIO)  # cara "+eje" expuesta
        resultado["exterior"][f"-{nombres_eje[eje]}"] = int(np.count_nonzero(sale_menos))
        resultado["exterior"][f"+{nombres_eje[eje]}"] = int(np.count_nonzero(sale_mas))
        distintas = (abajo != VACIO) & (arriba != VACIO) & (abajo != arriba)
        for p, q in zip(abajo[distintas], arriba[distintas]):
            clave = tuple(sorted((int(p), int(q))))
            resultado["interfaz"][clave] = resultado["interfaz"].get(clave, 0) + 1
    return resultado


def area_planta_octagono(entre_caras: float) -> float:
    """Área de un octágono regular dado el ancho entre caras opuestas (para comparar con la base voxelizada)."""
    r = entre_caras / 2
    return 8 * r**2 * np.tan(np.pi / 8)


# ---------------------------------------------------------------------------
# Construcción por vertidas: activación temporal
# ---------------------------------------------------------------------------

@dataclass
class PlanVertidas:
    """Cuándo existe cada celda: la vertida v se activa en Inicio_h (pset del IFC).

    Antes de su Inicio_h, las celdas de una vertida no existen para el solver
    (no conducen ni generan calor) y las caras de las vertidas ya puestas que las
    tocan son superficies expuestas al aire. Al activarse, sus celdas parten con
    T = TemperaturaColocacion_C y edad equivalente te = 0; desde ahí generan calor
    de hidratación y conducen con sus vecinas (la junta con la vertida anterior es
    una interfaz más: mismo material, conducción con la media armónica de k).
    """

    vox: Voxelizacion

    @property
    def inicio_h(self) -> np.ndarray:
        return np.array([float(p.get("Inicio_h", 0.0)) for p in self.vox.propiedades])

    @property
    def T_colocacion(self) -> np.ndarray:
        return np.array([float(p.get("TemperaturaColocacion_C", 20.0)) for p in self.vox.propiedades])

    def vertidas_activas(self, t_h: float) -> list[int]:
        return [v for v in range(1, self.vox.n_vertidas + 1) if self.inicio_h[v - 1] <= t_h + 1e-9]

    def celdas_activas(self, t_h: float) -> np.ndarray:
        return np.isin(self.vox.vertida_id, self.vertidas_activas(t_h))

    def eventos(self) -> list[tuple]:
        """(Inicio_h, vertida, nombre, T_colocación, celdas) en orden de colocación."""
        orden = np.argsort(self.inicio_h, kind="stable")
        return [(float(self.inicio_h[i]), i + 1, self.vox.nombres[i], float(self.T_colocacion[i]),
                 self.vox.celdas(i + 1)) for i in orden]

    def activar(self, T: np.ndarray, edad: np.ndarray, t_antes_h: float, t_ahora_h: float) -> list[int]:
        """Inicializa las vertidas que se colocan en (t_antes_h, t_ahora_h]: T de colocación y te = 0.

        Modifica T y edad en el lugar y devuelve las vertidas activadas.
        """
        nuevas = [v for v in self.vertidas_activas(t_ahora_h) if self.inicio_h[v - 1] > t_antes_h + 1e-9]
        for v in nuevas:
            celdas = self.vox.vertida_id == v
            T[celdas] = self.T_colocacion[v - 1]
            edad[celdas] = 0.0
        return nuevas
