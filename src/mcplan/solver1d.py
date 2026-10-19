"""Solver térmico 1D (dirección vertical) por balance de energía en celdas.

Es la versión 1D de la ecuación del curso:

    rho_i cp_i V_i dT_i/dt = sum_j G_ij (T_j - T_i) + Q_i V_i + Q_borde,i

Trabajamos por unidad de área en planta (A = 1 m2), así que V_i = dz_i y:

    C_i  = rho_i * cp_i * dz_i                       capacidad de la celda [J/(m2 K)]
    G_ij = 1 / (dz_i/(2 k_i) + dz_j/(2 k_j))         conductancia entre vecinas [W/(m2 K)]

G_ij son dos resistencias en serie (media celda i + media celda j). Así una
interfaz entre materiales distintos (hormigón-suelo, vertida-vertida) queda
bien representada sin tratamiento especial.

Integración en el tiempo: Euler explícito,

    T_i^{n+1} = T_i^n + dt / C_i * [ sum_j G_ij (T_j^n - T_i^n) + Q_i^n dz_i + q_borde,i^n ]

Estabilidad: el coeficiente que multiplica a T_i^n debe ser >= 0,

    1 - dt * sum(G_i) / C_i >= 0   =>   dt <= C_i / sum(G_i)   para toda celda i.

En una malla uniforme lejos de los bordes, sum(G_i) = 2k/dz y el criterio se
reduce al clásico  Fo = alpha*dt/dz^2 <= 1/2.

Convención: z apunta hacia arriba, con z = 0 en la base del dominio.
"""

from dataclasses import dataclass, field
from math import ceil
from typing import Callable, Optional, Union

import numpy as np

from .hidratacion import HidratacionJuLei2019

# Un valor de borde puede ser constante o una función del tiempo t [s].
Valor = Union[float, Callable[[float], float]]


def evaluar(valor: Valor, t: float) -> float:
    return float(valor(t)) if callable(valor) else float(valor)


@dataclass
class Material:
    nombre: str
    k: float  # conductividad térmica [W/(m K)]
    rho: float  # densidad [kg/m3]
    cp: float  # calor específico [J/(kg K)]
    es_hormigon: bool = False
    hidratacion: Optional[HidratacionJuLei2019] = None

    @property
    def alpha(self) -> float:
        """Difusividad térmica [m2/s]: qué tan rápido se propaga un cambio de temperatura."""
        return self.k / (self.rho * self.cp)


@dataclass
class Capa:
    material: Material
    espesor: float  # [m]
    T_inicial: float  # [°C]
    nombre: str = ""


@dataclass
class Borde:
    """Condición de borde en una cara del dominio.

    dirichlet: temperatura impuesta en la cara           T_cara = T
    neumann:   flujo de calor impuesto (+ entra al sólido) q [W/m2]
    robin:     convección con un fluido                    q = h (T_inf - T_cara)
    """

    tipo: str
    T: Valor = 0.0  # dirichlet: temperatura de la cara; robin: temperatura del fluido
    q: Valor = 0.0  # neumann: flujo entrante [W/m2]
    h: float = 0.0  # robin: coeficiente de transferencia [W/(m2 K)]

    def __post_init__(self):
        if self.tipo not in ("dirichlet", "neumann", "robin"):
            raise ValueError(f"Tipo de borde desconocido: {self.tipo!r}")
        if self.tipo == "robin" and self.h <= 0:
            raise ValueError("Un borde robin necesita h > 0")

    @classmethod
    def dirichlet(cls, T: Valor) -> "Borde":
        return cls("dirichlet", T=T)

    @classmethod
    def neumann(cls, q: Valor = 0.0) -> "Borde":
        return cls("neumann", q=q)

    @classmethod
    def robin(cls, h: float, T_inf: Valor) -> "Borde":
        return cls("robin", T=T_inf, h=h)

    def conductancia(self, k: float, dz: float) -> float:
        """Conductancia entre el centro de la celda de borde y el exterior [W/(m2 K)]."""
        if self.tipo == "dirichlet":
            return 2.0 * k / dz  # media celda de conducción
        if self.tipo == "robin":
            return 1.0 / (1.0 / self.h + dz / (2.0 * k))  # convección + media celda en serie
        return 0.0

    def flujo(self, t: float, T_celda: float, k: float, dz: float) -> float:
        """Flujo que entra a la celda de borde [W/m2]."""
        if self.tipo == "neumann":
            return evaluar(self.q, t)
        return self.conductancia(k, dz) * (evaluar(self.T, t) - T_celda)

    def temperatura_cara(self, t: float, T_celda: float, k: float, dz: float) -> float:
        """Temperatura en la superficie, reconstruida imponiendo continuidad de flujo."""
        if self.tipo == "dirichlet":
            return evaluar(self.T, t)
        if self.tipo == "neumann":
            return T_celda + evaluar(self.q, t) * dz / (2.0 * k)
        g = 2.0 * k / dz
        return (g * T_celda + self.h * evaluar(self.T, t)) / (g + self.h)

    def temperatura_exterior(self, t: float) -> float:
        return np.nan if self.tipo == "neumann" else evaluar(self.T, t)


@dataclass
class Malla1D:
    capas: list[Capa]
    dz: np.ndarray  # espesor de cada celda [m]
    z: np.ndarray  # cota del centro de cada celda [m]
    z_caras: np.ndarray  # cotas de las n+1 caras [m]
    indice_capa: np.ndarray  # a qué capa pertenece cada celda
    k: np.ndarray = field(init=False)
    rho: np.ndarray = field(init=False)
    cp: np.ndarray = field(init=False)

    def __post_init__(self):
        mats = [self.capas[i].material for i in self.indice_capa]
        self.k = np.array([m.k for m in mats])
        self.rho = np.array([m.rho for m in mats])
        self.cp = np.array([m.cp for m in mats])

    @property
    def n(self) -> int:
        return len(self.dz)

    @property
    def capacidad(self) -> np.ndarray:
        """C_i = rho cp dz [J/(m2 K)]."""
        return self.rho * self.cp * self.dz

    @property
    def conductancias(self) -> np.ndarray:
        """G entre celdas i e i+1 (n-1 valores) [W/(m2 K)]."""
        r = self.dz / (2.0 * self.k)  # resistencia de media celda
        return 1.0 / (r[:-1] + r[1:])

    @property
    def es_hormigon(self) -> np.ndarray:
        return np.array([self.capas[i].material.es_hormigon for i in self.indice_capa])


def crear_malla(capas: list[Capa], dz_objetivo: float) -> Malla1D:
    """Divide cada capa en un número entero de celdas de tamaño cercano a dz_objetivo.

    Las capas se listan de abajo hacia arriba. Cada interfaz entre capas cae
    exactamente sobre una cara de celda.
    """
    if dz_objetivo <= 0:
        raise ValueError("dz debe ser positivo")
    dz, indice = [], []
    for i, capa in enumerate(capas):
        n = max(1, round(capa.espesor / dz_objetivo))
        dz += [capa.espesor / n] * n
        indice += [i] * n
    dz = np.array(dz)
    z_caras = np.concatenate(([0.0], np.cumsum(dz)))
    z = 0.5 * (z_caras[:-1] + z_caras[1:])
    return Malla1D(capas, dz, z, z_caras, np.array(indice))


def paso_critico(malla: Malla1D, borde_inf: Borde, borde_sup: Borde) -> float:
    """Paso crítico del esquema explícito: min_i C_i / sum(G_i) [s].

    Con dt menor o igual a este valor, T_i^{n+1} es un promedio ponderado (con
    pesos >= 0) de T_i^n, sus vecinas y el exterior: no puede aparecer una
    temperatura fuera del rango de los datos (principio del máximo) y no hay
    oscilaciones. Es una condición suficiente: si se excede un poco, la solución
    puede oscilar sin explotar; si se excede claramente, diverge.
    """
    G = malla.conductancias
    suma_G = np.zeros(malla.n)
    suma_G[:-1] += G
    suma_G[1:] += G
    suma_G[0] += borde_inf.conductancia(malla.k[0], malla.dz[0])
    suma_G[-1] += borde_sup.conductancia(malla.k[-1], malla.dz[-1])
    with np.errstate(divide="ignore"):
        return float(np.min(malla.capacidad / suma_G))


@dataclass
class Resultado:
    malla: Malla1D
    t: np.ndarray  # instantes guardados [s]
    T: np.ndarray  # (nt, n) temperatura en los centros de celda [°C]
    T_cara_inf: np.ndarray  # (nt,) temperatura de la cara inferior del dominio [°C]
    T_cara_sup: np.ndarray  # (nt,) temperatura de la cara superior del dominio [°C]
    T_exterior_sup: np.ndarray  # (nt,) temperatura exterior sobre la cara superior [°C]
    te: np.ndarray  # (nt, n) edad equivalente [h] (0 donde no hay hidratación)
    dt: float
    dt_critico: float
    n_pasos: int
    diverge: bool = False

    @property
    def t_h(self) -> np.ndarray:
        return self.t / 3600.0


def simular(
    malla: Malla1D,
    borde_inf: Borde,
    borde_sup: Borde,
    t_final: float,
    dt: Optional[float] = None,
    guardar_cada: float = 3600.0,
    factor_seguridad: float = 0.9,
    dt_max: Optional[float] = None,
    forzar: bool = False,
) -> Resultado:
    """Integra el campo de temperatura desde t = 0 hasta t_final [s].

    dt = None elige automáticamente el mayor paso estable (con factor de
    seguridad), sin superar dt_max, que además divide exacto a `guardar_cada`.
    dt_max existe por PRECISIÓN, no por estabilidad: con mallas gruesas el paso
    estable puede ser > 1 h, demasiado largo para seguir el calor de hidratación.
    forzar = True permite un dt inestable (útil para mostrar qué pasa).
    """
    dt_crit = paso_critico(malla, borde_inf, borde_sup)
    if dt is None:
        dt_objetivo = factor_seguridad * dt_crit
        if dt_max is not None:
            dt_objetivo = min(dt_objetivo, dt_max)
        n_sub = ceil(guardar_cada / dt_objetivo)
        dt = guardar_cada / n_sub
    elif dt > dt_crit and not forzar:
        raise ValueError(
            f"dt = {dt:.1f} s supera el paso crítico {dt_crit:.1f} s: el esquema explícito "
            "sería inestable. Reduce dt, aumenta dz o usa forzar=True para verlo diverger."
        )

    n_pasos = int(round(t_final / dt))
    cada = max(1, int(round(guardar_cada / dt)))

    C = malla.capacidad
    G = malla.conductancias
    k, dz = malla.k, malla.dz
    T = np.array([malla.capas[i].T_inicial for i in malla.indice_capa], dtype=float)
    te = np.zeros(malla.n)

    # Celdas con calor de hidratación, agrupadas por capa (cada capa tiene un material).
    fuentes = []
    for i, capa in enumerate(malla.capas):
        modelo = capa.material.hidratacion
        if modelo is not None:
            fuentes.append((malla.indice_capa == i, modelo, capa.material.rho, capa.material.cp))

    historia = {"t": [], "T": [], "inf": [], "sup": [], "ext": [], "te": []}

    def guardar(t):
        historia["t"].append(t)
        historia["T"].append(T.copy())
        historia["te"].append(te.copy())
        historia["inf"].append(borde_inf.temperatura_cara(t, T[0], k[0], dz[0]))
        historia["sup"].append(borde_sup.temperatura_cara(t, T[-1], k[-1], dz[-1]))
        historia["ext"].append(borde_sup.temperatura_exterior(t))

    t = 0.0
    guardar(t)
    diverge = False
    for paso in range(1, n_pasos + 1):
        # 1) Conducción entre celdas vecinas: F > 0 si el calor va de i+1 hacia i.
        F = G * (T[1:] - T[:-1])
        balance = np.zeros(malla.n)
        balance[:-1] += F
        balance[1:] -= F

        # 2) Bordes del dominio.
        balance[0] += borde_inf.flujo(t, T[0], k[0], dz[0])
        balance[-1] += borde_sup.flujo(t, T[-1], k[-1], dz[-1])

        # 3) Calor de hidratación (con T y te del paso anterior: explícito).
        for m, modelo, rho, cp in fuentes:
            balance[m] += modelo.potencia(te[m], T[m], rho, cp) * dz[m]
            te[m] += dt / 3600.0 * modelo.factor_arrhenius(T[m])

        # 4) Actualización de Euler explícito.
        T = T + dt / C * balance
        t = paso * dt

        if not np.all(np.abs(T) < 1e4):  # la solución explotó (sólo posible con forzar)
            diverge = True
            guardar(t)
            break
        if paso % cada == 0 or paso == n_pasos:
            guardar(t)

    return Resultado(
        malla=malla,
        t=np.array(historia["t"]),
        T=np.array(historia["T"]),
        T_cara_inf=np.array(historia["inf"]),
        T_cara_sup=np.array(historia["sup"]),
        T_exterior_sup=np.array(historia["ext"]),
        te=np.array(historia["te"]),
        dt=dt,
        dt_critico=dt_crit,
        n_pasos=n_pasos,
        diverge=diverge,
    )
