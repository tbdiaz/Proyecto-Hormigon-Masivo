"""Tubos de enfriamiento embebidos en el hormigón (post-enfriamiento con agua).

Un tubo de 25-32 mm es mucho más delgado que una celda (0,11-0,275 m), así que no se
modela como vóxeles "de tubo": cada celda por la que pasa recibe un **sumidero de calor**
hacia el agua, y el agua se va calentando a lo largo del recorrido.

**Resistencia por metro** entre la celda (a su temperatura T_c) y el agua (T_w), en serie [m·K/W]:
- R_agua = 1 / (h_w π D_i): convección forzada interior; Nu = 0,023 Re^0,8 Pr^0,4 si Re > 2300 y
  4,36 si es laminar.
- R_pared = ln(D_o/D_i) / (2π k_tubo): aquí se distinguen PVC, HDPE, acero y cobre.
- R_hormigón = ln(r_eq/r_o) / (2π k_c), con r_eq = 0,198 Δx: el radio equivalente de Peaceman, que
  relaciona la temperatura de una celda cuadrada con la de un tubo (o pozo) dentro de ella.

U_L = 1/(R_agua + R_pared + R_hormigón) [W/(m K)].

**Agua a lo largo del tubo.** En el tramo de largo L_i dentro de la celda i:
    T_w,i+1 = T_c,i + (T_w,i - T_c,i) exp(-U_L L_i / (ṁ c_w))
y el calor que esa celda entrega al agua es ṁ c_w (T_w,i+1 - T_w,i).

**Supuestos declarados:**
- El tubo existe y conduce sólo en celdas ya colocadas (fuera de ellas no hay intercambio).
- El caudal es constante mientras el circuito está activo. La temperatura de entrada es constante o sigue un
  programa por tramos (`programa`: [{hasta_h, T_agua_C}, ...], horas desde que parte el circuito), p. ej. para
  seguir enfriando el núcleo con agua cada vez más fría antes de retirar la aislación.
- La capacidad térmica del agua dentro del tubo es despreciable: régimen cuasi estacionario en cada paso.
"""

from dataclasses import dataclass, field

import numpy as np

RHO_W, C_W, K_W, MU_W = 998.0, 4186.0, 0.60, 1.0e-3  # agua a ~20 °C
PR_W = C_W * MU_W / K_W
K_TUBO = {"pvc": 0.19, "hdpe": 0.45, "acero": 45.0, "cobre": 385.0}  # W/(m K)


@dataclass
class Circuito:
    """Un tubo con su recorrido (polilínea en m) y su operación."""

    nombre: str
    puntos: np.ndarray  # (n, 3) recorrido en orden de flujo [m]
    material: str = "hdpe"  # o un número: k de la pared [W/mK]
    diametro_mm: float = 32.0  # exterior
    espesor_mm: float = 2.0
    caudal_l_min: float = 20.0
    T_agua_C: float = 15.0
    vertida: int | None = None  # el circuito parte cuando empieza esta vertida (None: desde t0)
    desde_h: float = 0.0  # después del inicio de esa vertida (o de t0)
    durante_h: float = 168.0
    programa: list = field(default_factory=list)  # [(hasta_h, T_agua_C)] desde que parte; vacío = T_agua_C constante
    celdas: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.int64))  # recorrido en celdas
    largos: np.ndarray = field(default_factory=lambda: np.zeros(0))  # m de tubo dentro de cada celda

    @property
    def k_pared(self) -> float:
        return float(K_TUBO[self.material.lower()]) if isinstance(self.material, str) else float(self.material)

    @property
    def m_c(self) -> float:
        """ṁ c_w [W/K]: capacidad del caudal de agua."""
        return RHO_W * self.caudal_l_min / 60000.0 * C_W

    def tramos(self) -> list[tuple[float, float, float]]:
        """(desde_h, hasta_h, T_entrada) de cada tramo de operación, en horas desde que parte el circuito."""
        if not self.programa:
            return [(0.0, float(self.durante_h), float(self.T_agua_C))]
        out, desde = [], 0.0
        for hasta, T in self.programa:
            out.append((desde, float(hasta), float(T)))
            desde = float(hasta)
        return out

    @property
    def largo_m(self) -> float:
        return float(self.largos.sum())

    def resistencias(self):
        """(R_agua, R_pared) [m K/W] y el número de Reynolds."""
        D_o = self.diametro_mm / 1000.0
        D_i = D_o - 2 * self.espesor_mm / 1000.0
        if D_i <= 0:
            raise ValueError(f"{self.nombre}: el espesor no deja diámetro interior")
        v = self.caudal_l_min / 60000.0 / (np.pi * D_i**2 / 4)
        Re = RHO_W * v * D_i / MU_W
        Nu = 0.023 * Re**0.8 * PR_W**0.4 if Re > 2300 else 4.36
        h_w = Nu * K_W / D_i
        return 1.0 / (h_w * np.pi * D_i), np.log(D_o / D_i) / (2 * np.pi * self.k_pared), Re

    def U_L(self, k_hormigon, dx):
        """W/(m K) por metro de tubo, para celdas de conductividad k_hormigon (escalar o arreglo)."""
        R_agua, R_pared, _ = self.resistencias()
        r_o = self.diametro_mm / 2000.0
        r_eq = 0.198 * dx  # Peaceman (1978), celda cuadrada
        R_h = np.maximum(np.log(r_eq / r_o), 0.0) / (2 * np.pi * np.asarray(k_hormigon, dtype=float))
        return 1.0 / (R_agua + R_pared + R_h)


def recorrer(dom, puntos, paso_frac=0.1):
    """Celdas por las que pasa la polilínea, en orden de flujo, y el largo de tubo dentro de cada una."""
    puntos = np.asarray(puntos, dtype=float)
    if puntos.ndim != 2 or puntos.shape[1] != 3 or len(puntos) < 2:
        raise ValueError("Un circuito necesita al menos dos puntos (x, y, z)")
    paso = dom.dx * paso_frac
    muestras, largo = [], []
    for a, b in zip(puntos[:-1], puntos[1:]):
        L = float(np.linalg.norm(b - a))
        n = max(1, int(np.ceil(L / paso)))
        s = (np.arange(n) + 0.5) / n
        muestras.append(a + s[:, None] * (b - a))
        largo.append(np.full(n, L / n))
    P, Ls = np.vstack(muestras), np.concatenate(largo)
    ijk = np.floor((P - dom.origen) / dom.dx).astype(np.int64)
    if np.any(ijk < 0) or np.any(ijk >= np.array(dom.forma)):
        raise ValueError("El recorrido del tubo sale de la grilla")
    plano = np.ravel_multi_index(ijk.T, dom.forma)
    # Agrupar muestras consecutivas en la misma celda (conservando el orden de flujo).
    corte = np.flatnonzero(np.diff(plano)) + 1
    grupos = np.split(np.arange(len(plano)), corte)
    celdas = np.array([plano[g[0]] for g in grupos], dtype=np.int64)
    largos = np.array([Ls[g].sum() for g in grupos])
    return celdas, largos


def serpentin(dom, vertida, z, separacion=1.0, margen=0.5, direccion="x"):
    """Polilínea en zigzag a la altura z, dentro de la planta de la vertida a esa altura.

    Filas paralelas a `direccion` separadas `separacion` m, recortadas `margen` m antes del borde,
    unidas en sus extremos (un solo circuito continuo).
    """
    from .simulacion3d import HORMIGON

    k = int(np.floor((z - dom.origen[2]) / dom.dx))
    capa = (dom.vertida_id[:, :, k] == vertida) & (dom.material[:, :, k] == HORMIGON)
    if not capa.any():
        raise ValueError(f"No hay hormigón de la Vertida_{vertida:02d} a z = {z} m")
    x, y = dom.centros(0), dom.centros(1)
    if direccion == "y":
        capa, x, y = capa.T, y, x
    filas = np.flatnonzero(capa.any(axis=0))
    y0, y1 = y[filas[0]] + margen, y[filas[-1]] - margen
    puntos = []
    for n, yy in enumerate(np.arange(y0, y1 + 1e-9, separacion)):
        j = int(np.argmin(np.abs(y - yy)))
        dentro = np.flatnonzero(capa[:, j])
        if len(dentro) == 0:
            continue
        xa, xb = x[dentro[0]] + margen, x[dentro[-1]] - margen
        if xb <= xa:
            continue
        fila = [(xa, yy), (xb, yy)] if n % 2 == 0 else [(xb, yy), (xa, yy)]
        puntos += fila
    if len(puntos) < 2:
        raise ValueError(f"El serpentín no cabe en la Vertida_{vertida:02d} a z = {z} m")
    P = np.array([(a, b, z) for a, b in puntos])
    return P[:, [1, 0, 2]] if direccion == "y" else P


def leer_programa(d: dict) -> list[tuple[float, float]]:
    """`programa: [{hasta_h: 168, T_agua_C: 25}, {hasta_h: 336, T_agua_C: 20}]` -> [(168, 25), (336, 20)].

    Las horas se cuentan desde que parte el circuito (igual que `durante_h`), crecen estrictamente y el último
    tramo marca el fin. Con programa, `T_agua_C` y `durante_h` no se dan aparte: salen del programa.
    """
    tramos = d.get("programa")
    if not tramos:
        return []
    if "T_agua_C" in d or "durante_h" in d:
        raise ValueError("Un circuito con `programa` no lleva T_agua_C ni durante_h: salen del programa")
    out, antes = [], 0.0
    for n, t in enumerate(tramos, start=1):
        try:
            hasta, T = float(t["hasta_h"]), float(t["T_agua_C"])
        except (KeyError, TypeError, ValueError):
            raise ValueError(f"Programa del tubo, tramo {n}: necesita hasta_h y T_agua_C numéricos") from None
        if hasta <= antes:
            raise ValueError(f"Programa del tubo, tramo {n}: hasta_h = {hasta:g} debe ser mayor que {antes:g}")
        if not -5.0 <= T <= 60.0:
            raise ValueError(f"Programa del tubo, tramo {n}: T_agua_C = {T:g} °C fuera de rango (-5 a 60)")
        out.append((hasta, T))
        antes = hasta
    return out


def circuito_desde_dict(dom, d: dict, nombre_defecto="tubo") -> Circuito:
    """Lee un circuito del YAML: `puntos` explícitos o un `serpentin` {z, separacion, margen, direccion}."""
    from .simulacion3d import HORMIGON

    vertida = d.get("vertida")
    if "puntos" in d:
        puntos = np.array(d["puntos"], dtype=float)
    elif "serpentin" in d:
        s = d["serpentin"]
        if vertida is None:
            raise ValueError("Un serpentín necesita `vertida`")
        z = s.get("z")
        if z is None:  # a media altura de la vertida
            capas = np.flatnonzero(np.any(dom.vertida_id == int(vertida), axis=(0, 1)))
            z = float(dom.centros(2)[capas].mean())
        puntos = serpentin(dom, int(vertida), float(z), float(s.get("separacion", 1.0)), float(s.get("margen", 0.5)),
                           s.get("direccion", "x"))
    else:
        raise ValueError("Un circuito necesita `puntos` o `serpentin`")
    material = d.get("material", "hdpe")
    if isinstance(material, str) and material.lower() not in K_TUBO:
        raise ValueError(f"Material de tubo desconocido: {material} (usa {', '.join(K_TUBO)} o un k en W/mK)")
    programa = leer_programa(d)
    c = Circuito(nombre=str(d.get("nombre", nombre_defecto)), puntos=puntos, material=material,
                 diametro_mm=float(d.get("diametro_mm", 32)), espesor_mm=float(d.get("espesor_mm", 2)),
                 caudal_l_min=float(d.get("caudal_l_min", 20)),
                 T_agua_C=programa[0][1] if programa else float(d.get("T_agua_C", 15)),
                 vertida=None if vertida is None else int(vertida), desde_h=float(d.get("desde_h", 0)),
                 durante_h=programa[-1][0] if programa else float(d.get("durante_h", 168)), programa=programa)
    c.celdas, c.largos = recorrer(dom, puntos)
    fuera = dom.material.reshape(-1)[c.celdas] != HORMIGON
    if fuera.any():
        raise ValueError(f"{c.nombre}: {int(fuera.sum())} tramos del tubo quedan fuera del hormigón")
    c.resistencias()
    return c


def agua_a_lo_largo(a, Tc, T_entrada, bloque=256):
    """Temperatura del agua a la salida de cada celda: T_{i+1} = a_i T_i + (1 - a_i) Tc_i (vectorizado por bloques)."""
    Tw = np.empty(len(a))
    T0 = T_entrada
    for i in range(0, len(a), bloque):
        ai, bi = a[i:i + bloque], (1 - a[i:i + bloque]) * Tc[i:i + bloque]
        P = np.cumprod(ai)
        if P[-1] > 1e-200:
            Tw[i:i + bloque] = P * (T0 + np.cumsum(bi / P))
        else:  # intercambio muy fuerte en el bloque: recurrencia directa (rara vez se usa)
            for j in range(len(ai)):
                T0 = ai[j] * T0 + bi[j]
                Tw[i + j] = T0
        T0 = Tw[i + len(ai) - 1]
    return Tw
