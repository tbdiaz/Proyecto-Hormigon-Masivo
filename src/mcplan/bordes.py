"""Condiciones de borde: cómo entra el ambiente a la física (E10).

Cada cara de celda expuesta (hormigón o suelo con aire del otro lado) intercambia
calor con el ambiente a través de una cadena de resistencias en serie:

    celda --(media celda: dx/2k)--> superficie del hormigón --(capa: e/k)--> superficie exterior --(1/h_ext)--> ambiente

    q [W/m²] = U (T_eq - T_celda),     U = 1 / (1/h_ext + R_capa + dx/(2k))

R_capa es la resistencia del moldaje o de la aislación (0 si la cara está desnuda).
Dos modelos de ambiente:

- "combinado" (E08/E09): h_ext fijo, que ya incluye convección y radiación, y
  T_eq = T_aire. No se agrega radiación aparte, porque eso la contaría dos veces.
- "fisico" (E10 en adelante): cada mecanismo entra una sola vez.
    convección:  h_c = 5,6 + 4,0 v  (v <= 5 m/s);  7,2 v^0,78  (v > 5)   [Jürges / McAdams]
    onda larga:  h_r = 4 ε σ T_m³ (linealización de ε σ (T_s⁴ - T_rad⁴)), con
                 T_rad = F_cielo T_cielo + (1 - F_cielo) T_aire
    onda corta:  q_sol = α G  sólo en caras horizontales hacia arriba (G horizontal)
    h_ext = h_c + h_r,   T_eq = (h_c T_aire + h_r T_rad + α G) / h_ext   ("temperatura sol-aire")

T_eq es la temperatura del ambiente "equivalente": con q = h_ext (T_eq - T_s) se
recupera exactamente la suma convección + radiación + sol. El viento sólo cambia
h_c; el sol sólo suma α G; la temperatura del aire entra en la convección y en
la parte del cielo/suelo que ve la cara.
"""

from dataclasses import dataclass, field

import numpy as np

SIGMA = 5.670374419e-8  # Stefan-Boltzmann [W/(m² K⁴)]
KELVIN = 273.15
SUP, LAT, INF = 0, 1, 2  # orientación de una cara expuesta: +z, ±x/±y, -z
NOMBRE_ORIENTACION = ("superior", "lateral", "inferior")
F_CIELO = np.array([1.0, 0.5, 0.0])  # fracción de cielo que ve la cara (el resto ve suelo/entorno a T_aire)
F_SOL = np.array([1.0, 0.0, 0.0])  # G es horizontal: sólo la cara superior lo recibe (simplificación declarada)


def h_conveccion(v):
    """Coeficiente de convección forzada [W/(m² K)] con el viento a 10 m v [m/s] (Jürges / McAdams)."""
    v = np.maximum(np.asarray(v, dtype=float), 0.0)
    return np.where(v <= 5.0, 5.6 + 4.0 * v, 7.2 * np.power(v, 0.78))


def T_cielo(T_aire_C, nubosidad=0.0):
    """Temperatura efectiva del cielo [°C]: Swinbank (cielo despejado) corregida por nubosidad (0-1).

    T_despejado = 0,0552 T_aire^1,5 [K];  T_cielo⁴ = (1 - c) T_despejado⁴ + c T_aire⁴
    (una nube gruesa emite casi como un cuerpo negro a la temperatura del aire).
    """
    Ta = np.asarray(T_aire_C, dtype=float) + KELVIN
    c = np.clip(np.asarray(nubosidad, dtype=float), 0.0, 1.0)
    despejado = 0.0552 * Ta**1.5
    return ((1 - c) * despejado**4 + c * Ta**4) ** 0.25 - KELVIN


@dataclass(frozen=True)
class Capa:
    """Moldaje o aislación: una resistencia térmica e/k en serie (su capacidad térmica se desprecia)."""

    nombre: str
    espesor: float  # m
    k: float  # W/(m K)
    absortividad: float = 0.6  # de su cara exterior (onda corta)
    emisividad: float = 0.9  # de su cara exterior (onda larga)

    @property
    def R(self) -> float:
        return self.espesor / self.k


@dataclass
class Clima:
    """Serie horaria del ambiente; t_h = 0 es la colocación de la Vertida_01."""

    t_h: np.ndarray
    T_aire: np.ndarray  # °C a 2 m
    viento: np.ndarray  # m/s a 10 m
    G: np.ndarray  # W/m² onda corta global horizontal
    nubosidad: np.ndarray  # 0-1
    origen: list = field(default_factory=list)  # de dónde sale cada hora: "pronóstico", "transición", "climatología"
    metadata: dict = field(default_factory=dict)

    def en(self, t_h):
        return tuple(np.interp(t_h, self.t_h, s) for s in (self.T_aire, self.viento, self.G, self.nubosidad))


@dataclass
class Ambiente:
    """Ambiente de la fundación: aire, viento, sol y cielo (o un h fijo), más el fondo del suelo."""

    T_amb: object = 20.0  # °C: número o función de t [s] (si no hay clima)
    h_superior: float = 27.8  # modelo combinado: W/(m² K) en caras ±z
    h_lateral: float = 27.8  # modelo combinado: W/(m² K) en caras ±x, ±y
    T_fondo: float = 15.0  # °C en la cara inferior del suelo (Dirichlet)
    modelo: str = "combinado"  # "combinado" o "fisico"
    clima: Clima | None = None
    absortividad: float = 0.65  # hormigón gris (onda corta)
    emisividad: float = 0.9  # hormigón (onda larga)

    def __post_init__(self):
        if self.modelo not in ("combinado", "fisico"):
            raise ValueError(f"Modelo de ambiente desconocido: {self.modelo}")

    def temperatura(self, t_s: float) -> float:
        if self.clima is not None:
            return float(np.interp(t_s / 3600.0, self.clima.t_h, self.clima.T_aire))
        return float(self.T_amb(t_s)) if callable(self.T_amb) else float(self.T_amb)

    def _meteo(self, t_s):
        if self.clima is not None:
            return self.clima.en(t_s / 3600.0)
        return self.temperatura(t_s), 0.0, 0.0, 0.0

    def exterior(self, t_s, orientacion, absortividad, emisividad):
        """(h_ext [W/m²K], T_eq [°C]) para cada categoría de cara (arreglos del mismo largo)."""
        orientacion = np.asarray(orientacion)
        if self.modelo == "combinado":
            h = np.where(orientacion == LAT, self.h_lateral, self.h_superior).astype(float)
            return h, np.full(h.shape, self.temperatura(t_s))
        T_a, v, G, c = self._meteo(t_s)
        T_rad = F_CIELO[orientacion] * T_cielo(T_a, c) + (1 - F_CIELO[orientacion]) * T_a
        h_c = np.full(orientacion.shape, float(h_conveccion(v)))
        h_r = 4 * np.asarray(emisividad) * SIGMA * ((T_a + T_rad) / 2 + KELVIN) ** 3
        q_sol = np.asarray(absortividad) * G * F_SOL[orientacion]
        h = h_c + h_r
        return h, (h_c * T_a + h_r * T_rad + q_sol) / h

    def h_maximo(self, orientacion, emisividad):
        """Cota superior de h_ext en toda la serie: fija el Δt crítico de una etapa (Euler explícito)."""
        orientacion = np.asarray(orientacion)
        if self.modelo == "combinado":
            return np.where(orientacion == LAT, self.h_lateral, self.h_superior).astype(float)
        if self.clima is None:
            v_max, T_max = 0.0, self.temperatura(0.0)
        else:
            v_max, T_max = float(np.max(self.clima.viento)), float(np.max(self.clima.T_aire))
        return float(h_conveccion(v_max)) + 4 * np.asarray(emisividad) * SIGMA * (T_max + KELVIN + 5) ** 3

    def describir(self) -> str:
        if self.modelo == "combinado":
            return f"h fijo {self.h_superior:g}/{self.h_lateral:g} W/m²K (sup./lat.), incluye radiación"
        return (f"convección h_c(v) de Jürges + onda larga al cielo (ε = {self.emisividad:g}) + sol α G "
                f"(α = {self.absortividad:g}) en caras superiores")


def coeficiente_U(h_ext, R, dx, k):
    """U [W/m²K] de la cadena ambiente -> capa -> media celda; 0 si h_ext = 0 (cara adiabática)."""
    h_ext = np.asarray(h_ext, dtype=float)
    with np.errstate(divide="ignore"):
        return np.where(h_ext > 0, 1.0 / (1.0 / np.where(h_ext > 0, h_ext, 1.0) + R + dx / (2 * np.asarray(k))), 0.0)


# ---------------------------------------------------------------------------
# Moldaje y aislación por vertida, con su retiro
# ---------------------------------------------------------------------------

@dataclass
class Cubierta:
    """Una capa (moldaje o aislación) sobre ciertas caras de ciertas vertidas, desde la colocación hasta el retiro."""

    tipo: str  # "moldaje" o "aislacion"
    capa: Capa
    orientaciones: tuple  # (LAT,) para moldaje; (SUP,) para una manta superior
    vertidas: tuple  # vertidas cubiertas (1..n)
    retiro_h: dict  # vertida -> horas desde su colocación hasta el retiro (inf = no se retira)


@dataclass
class PlanCubiertas:
    cubiertas: list = field(default_factory=list)

    def eventos(self, inicio_h, base_retiro_h=None) -> list[tuple]:
        """(t_h, descripción) de cada retiro, en orden. El retiro se cuenta desde base_retiro_h
        (fin de la colocación de cada vertida) o, si no se indica, desde inicio_h."""
        base = inicio_h if base_retiro_h is None else base_retiro_h
        ev = []
        for c in self.cubiertas:
            for v in c.vertidas:
                t = float(base[v - 1]) + float(c.retiro_h.get(v, np.inf))
                if np.isfinite(t):
                    ev.append((t, f"retiro de {c.tipo} ({c.capa.nombre}) de la Vertida_{v:02d}"))
        return sorted(ev)

    def tabla(self, t_h, inicio_h, n_vertidas, base_retiro_h=None):
        """Capas activas en t_h: (combo[v, orientación] -> índice, lista de combinaciones de capas).

        Fila 0 = sin vertida (suelo). Una combinación es una tupla de capas en serie,
        de adentro hacia afuera; la tupla vacía es una cara desnuda.
        """
        combos = [()]
        tabla = np.zeros((n_vertidas + 1, 3), dtype=np.int64)
        base = inicio_h if base_retiro_h is None else base_retiro_h
        for v in range(1, n_vertidas + 1):
            t0 = float(inicio_h[v - 1])
            if not np.isfinite(t0) or t_h < t0 - 1e-9:
                continue
            tb = float(base[v - 1])
            for o in (SUP, LAT, INF):
                capas = tuple(c.capa for c in self.cubiertas if v in c.vertidas and o in c.orientaciones
                              and t_h < tb + float(c.retiro_h.get(v, np.inf)) - 1e-9)
                if capas not in combos:
                    combos.append(capas)
                tabla[v, o] = combos.index(capas)
        return tabla, combos

    def describir(self, inicio_h) -> list[str]:
        lineas = []
        for c in self.cubiertas:
            caras = "/".join(NOMBRE_ORIENTACION[o] for o in c.orientaciones)
            retiros = ", ".join(f"V{v}: {c.retiro_h.get(v, np.inf):g} h" for v in c.vertidas)
            lineas.append(f"{c.tipo}: {c.capa.nombre} (e = {c.capa.espesor * 1000:g} mm, k = {c.capa.k:g} W/mK, "
                          f"R = {c.capa.R:.3f} m²K/W) en caras {caras} de las vertidas {list(c.vertidas)}; "
                          f"retiro desde la colocación: {retiros}")
        return lineas or ["sin moldaje ni aislación (todas las caras desnudas)"]


# ---------------------------------------------------------------------------
# Tabla de condiciones de borde (para la consola, el panel y los informes)
# ---------------------------------------------------------------------------

def _rango_h(amb: Ambiente, orientacion, absortividad, emisividad):
    """h_ext y T_eq en todas las horas del clima (o en t = 0 si no hay clima)."""
    horas = amb.clima.t_h if amb.clima is not None else np.array([0.0])
    hs, Teq = [], []
    for t in horas:
        h, T = amb.exterior(t * 3600.0, np.array([orientacion]), np.array([absortividad]), np.array([emisividad]))
        hs.append(h[0])
        Teq.append(T[0])
    return np.array(hs), np.array(Teq)


def filas_bordes(amb: Ambiente, cubiertas: PlanCubiertas, dx, k_hormigon, k_suelo, areas=None, inicio_h=None):
    """Una fila por tipo de superficie: tipo, h, R de capa, β = h/k, U con media celda, T_eq y área.

    β = h'/k [1/m], con h' = 1/(1/h_ext + R_capa) el coeficiente entre la superficie del
    hormigón y el ambiente. Neumann adiabático: β = 0; Dirichlet: β = ∞.
    """
    areas = areas or {}
    filas = []

    def robin(nombre, clave, orientacion, k, capa=None, nota=""):
        alfa = capa.absortividad if capa else amb.absortividad
        eps = capa.emisividad if capa else amb.emisividad
        h, Teq = _rango_h(amb, orientacion, alfa, eps)
        R = capa.R if capa else 0.0
        h_sup = 1.0 / (1.0 / h + R)
        U = coeficiente_U(h, R, dx, k)
        filas.append({
            "superficie": nombre, "tipo": "Robin", "capa": capa.nombre if capa else "desnuda", "R": R,
            "h": [float(h.min()), float(h.mean()), float(h.max())],
            "beta": [float(h_sup.min() / k), float(h_sup.mean() / k), float(h_sup.max() / k)],
            "h_ef": [float(U.min()), float(U.mean()), float(U.max())], "k": k,
            "T_eq": [float(Teq.min()), float(Teq.mean()), float(Teq.max())], "nota": nota,
            "area": areas.get(clave),
        })

    tipos = {c.tipo: c for c in cubiertas.cubiertas}
    sol = "aire + sol + cielo" if amb.modelo == "fisico" else "T_aire"
    robin("Hormigón: caras superiores al aire", "sup", SUP, k_hormigon, nota=sol)
    if "aislacion" in tipos and SUP in tipos["aislacion"].orientaciones:
        c = tipos["aislacion"]
        robin(f"Hormigón: caras superiores con {c.capa.nombre}", "sup", SUP, k_hormigon, c.capa,
              nota=f"vertidas {list(c.vertidas)} hasta el retiro")
    if "moldaje" in tipos:
        c = tipos["moldaje"]
        robin(f"Hormigón: caras laterales con moldaje ({c.capa.nombre})", "lat_mol", LAT, k_hormigon, c.capa,
              nota=f"vertidas {list(c.vertidas)} hasta el desmolde")
    robin("Hormigón: caras laterales desnudas", "lat", LAT, k_hormigon,
          nota="troncos (cara inclinada) y tras el desmolde")
    robin("Suelo: superficie al aire (z = 0)", "suelo_sup", SUP, k_suelo, nota=sol)
    filas.append({"superficie": "Suelo: fondo", "tipo": "Dirichlet", "capa": "—", "R": 0.0, "h": None, "beta": None,
                  "h_ef": [2 * k_suelo / dx] * 3, "k": k_suelo, "T_eq": [amb.T_fondo] * 3,
                  "nota": f"T = {amb.T_fondo:.1f} °C", "area": areas.get("fondo")})
    filas.append({"superficie": "Suelo: caras laterales", "tipo": "Neumann", "capa": "—", "R": 0.0, "h": [0.0] * 3,
                  "beta": [0.0] * 3, "h_ef": [0.0] * 3, "k": k_suelo, "T_eq": None, "nota": "q = 0 (suelo lejano)",
                  "area": areas.get("suelo_lat")})
    return filas
