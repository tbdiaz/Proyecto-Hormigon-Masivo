"""Simulación térmica 3D de la fundación por vertidas: IFC -> vóxeles -> solver por etapas.

Es el mismo balance de energía por celda de solver_nd (volúmenes finitos, Euler
explícito, media armónica de k en las caras, Robin con media celda en serie),
pero sobre la grilla voxelizada del IFC y con celdas que **aparecen en el tiempo**:

- cada celda tiene un material: VACIO_M (aire), SUELO o HORMIGON;
- una celda de hormigón existe sólo desde el Inicio_h de su vertida: antes no
  conduce, no genera calor y no tiene edad;
- una cara entre dos celdas que existen conduce; una cara entre una celda que
  existe y una que no (aire, o una vertida que todavía no se coloca) es una
  superficie expuesta: q = U (T_eq - T), con el ambiente, el moldaje y la aislación
  del módulo bordes.py (que también decide cuándo se retiran);
- el suelo ocupa toda la planta de la grilla bajo z = 0, con T_fondo (Dirichlet)
  en su cara inferior y caras laterales adiabáticas (dominio "lejano").

La secuencia constructiva NO está en el solver: llega como un PlanVertidas
(Inicio_h y temperatura de colocación de cada vertida), que se arma con el pset
del IFC y lo que el YAML o la línea de comandos cambien (`aplicar_cambios`).
El solver sólo integra por tramos entre eventos (colocación o retiro de moldaje o
aislación) y recalcula las conductancias cada vez que cambia la geometría.
"""

import math
from dataclasses import dataclass, field

import numpy as np

from .bordes import INF, LAT, SUP, Ambiente, Capa, Cubierta, PlanCubiertas, coeficiente_U
from .voxelizacion import VACIO, PlanVertidas, Voxelizacion

VACIO_M, SUELO, HORMIGON = 0, 1, 2  # material de cada celda


@dataclass
class Material:
    k: float  # W/(m K)
    rho: float  # kg/m3
    cp: float  # J/(kg K)

    @property
    def capacidad(self) -> float:
        return self.rho * self.cp


@dataclass
class Dominio:
    """Grilla del cálculo: vóxeles del IFC + capa de suelo + aire alrededor."""

    dx: float
    origen: np.ndarray
    material: np.ndarray  # (nx, ny, nz) int8: VACIO_M, SUELO, HORMIGON
    vertida_id: np.ndarray  # (nx, ny, nz) int8: 0 o 1..n
    vox: Voxelizacion  # la voxelización original (volúmenes, nombres, pset)

    @property
    def forma(self):
        return self.material.shape

    def centros(self, eje: int) -> np.ndarray:
        return self.origen[eje] + (np.arange(self.forma[eje]) + 0.5) * self.dx

    def indice(self, x: float, y: float, z: float) -> tuple[int, int, int]:
        """Celda que contiene el punto (x, y, z)."""
        i = [int(math.floor((c - self.origen[e]) / self.dx)) for e, c in enumerate((x, y, z))]
        if not all(0 <= i[e] < self.forma[e] for e in range(3)):
            raise ValueError(f"El punto ({x}, {y}, {z}) está fuera de la grilla")
        return tuple(i)


def construir_dominio(vox: Voxelizacion, profundidad_suelo: float = 3.0, margen_suelo: float = 2.0) -> Dominio:
    """Agrega a la voxelización una capa de suelo bajo z = 0 y un margen lateral.

    La voxelización ya deja una celda de aire alrededor. Aquí se agregan
    `margen_suelo` metros de grilla a cada lado (aire arriba de z = 0, suelo
    abajo) y celdas de suelo hasta z = -profundidad_suelo. Con profundidad 0 no
    hay suelo: la base queda expuesta (útil para comparar con solver_nd).
    """
    dx = vox.dx
    m = int(round(margen_suelo / dx)) if profundidad_suelo > 0 else 0
    k0 = int(round(-vox.origen[2] / dx))  # capas de la voxelización bajo z = 0
    p = max(int(round(profundidad_suelo / dx)) - k0, 0) if profundidad_suelo > 0 else 0
    vid = np.pad(vox.vertida_id, ((m, m), (m, m), (p, 0)))
    origen = vox.origen - np.array([m, m, p]) * dx
    material = np.where(vid != VACIO, HORMIGON, VACIO_M).astype(np.int8)
    if profundidad_suelo > 0:
        zc = origen[2] + (np.arange(vid.shape[2]) + 0.5) * dx
        material[:, :, zc < 0] = SUELO
    return Dominio(dx=dx, origen=origen, material=material, vertida_id=vid, vox=vox)


# ---------------------------------------------------------------------------
# Cambios a la secuencia constructiva (inputs, no solver)
# ---------------------------------------------------------------------------

def aplicar_cambios(vox: Voxelizacion, cambios: dict | None) -> list[str]:
    """Sobrescribe el pset de las vertidas con lo que diga el YAML o la línea de comandos.

    cambios = {"Vertida_02": {"Inicio_h": 12, "TemperaturaColocacion_C": 25, "activa": False}, ...}
    (también acepta el número de vertida como clave). Modifica vox.propiedades y
    devuelve una descripción de cada cambio, para dejar constancia en la salida.
    """
    registro = []
    for clave, valores in (cambios or {}).items():
        v = int(clave) if str(clave).isdigit() else vox.nombres.index(str(clave)) + 1
        props = vox.propiedades[v - 1]
        for nombre, valor in (valores or {}).items():
            nombre = {"inicio_h": "Inicio_h", "T_colocacion": "TemperaturaColocacion_C"}.get(nombre, nombre)
            if nombre not in ("Inicio_h", "TemperaturaColocacion_C", "activa"):
                raise ValueError(f"Cambio no reconocido para {vox.nombres[v - 1]}: {nombre}")
            antes = props.get(nombre, True if nombre == "activa" else None)
            props[nombre] = bool(valor) if nombre == "activa" else float(valor)
            registro.append(f"{vox.nombres[v - 1]}: {nombre} {antes} -> {props[nombre]}")
    return registro


# ---------------------------------------------------------------------------
# Solver por etapas
# ---------------------------------------------------------------------------

@dataclass
class Activacion:
    """Cuándo aparece cada celda de hormigón y a qué temperatura.

    Con un PlanVertidas, cada vertida del IFC aparece entera en su Inicio_h. Un plan de obra
    (planes.py) puede, en cambio, colocar cada etapa por capas a la tasa de la bomba: las
    celdas aparecen de abajo hacia arriba y cada una parte con su propia edad 0.
    """

    t_act: np.ndarray  # forma del dominio: hora de colocación de cada celda (inf = no se coloca)
    T_col: np.ndarray  # forma del dominio: temperatura de colocación [°C]
    inicio_h: np.ndarray  # por vertida del IFC: primera celda colocada (inf si no se construye)
    fin_h: np.ndarray  # por vertida del IFC: última celda colocada (desde aquí se cuenta el retiro del moldaje)
    T_col_v: np.ndarray  # por vertida del IFC: temperatura de colocación

    @classmethod
    def desde_plan(cls, dom: "Dominio", plan: PlanVertidas) -> "Activacion":
        vid = dom.vertida_id.astype(np.int64)
        ini = np.concatenate([[np.inf], plan.inicio_h])
        Tc = np.concatenate([[np.nan], plan.T_colocacion])
        return cls(t_act=ini[vid], T_col=Tc[vid], inicio_h=np.asarray(plan.inicio_h, float),
                   fin_h=np.asarray(plan.inicio_h, float), T_col_v=np.asarray(plan.T_colocacion, float))

    def vertidas_activas(self, t_h) -> list[int]:
        return [v + 1 for v, t in enumerate(self.inicio_h) if t <= t_h + 1e-9]


_NORMALES = np.array([[-1, 0, 0], [1, 0, 0], [0, -1, 0], [0, 1, 0], [0, 0, -1], [0, 0, 1]], dtype=float)


def _caras_expuestas(activas):
    """Caras entre una celda que existe y una que no: (celda plana, dirección 0..5, orientación SUP/LAT/INF).

    Las caras en el límite de la grilla no se cuentan (adiabáticas): el hormigón nunca
    llega ahí, y el suelo es "lejano" a los lados.
    """
    indice = np.arange(activas.size).reshape(activas.shape)
    celdas, direcciones = [], []
    for eje in range(3):
        a = [slice(None)] * 3
        b = [slice(None)] * 3
        a[eje], b[eje] = slice(None, -1), slice(1, None)
        a, b = tuple(a), tuple(b)
        hacia_mas = activas[a] & ~activas[b]  # la celda de abajo expone su cara +eje
        hacia_menos = activas[b] & ~activas[a]  # la celda de arriba expone su cara -eje
        celdas += [indice[a][hacia_mas], indice[b][hacia_menos]]
        direcciones += [np.full(np.count_nonzero(hacia_mas), 2 * eje + 1),
                        np.full(np.count_nonzero(hacia_menos), 2 * eje)]
    celda, direccion = np.concatenate(celdas), np.concatenate(direcciones)
    orientacion = np.where(direccion == 5, SUP, np.where(direccion == 4, INF, LAT))
    return celda, direccion, orientacion


@dataclass
class _Etapa:
    """Lo que no cambia entre dos eventos (colocación o retiro de una capa)."""

    activas: np.ndarray
    G: list  # conductancia por volumen [W/(m3 K)] entre vecinas, por eje
    superficie: np.ndarray  # celdas (índice plano) con alguna cara expuesta
    N: np.ndarray  # (n_superficie, n_categorías): cuántas caras de cada categoría tiene cada celda
    cat_orientacion: np.ndarray
    cat_R: np.ndarray  # resistencia de las capas [m²K/W]
    cat_k: np.ndarray  # k de la celda (hormigón o suelo)
    cat_alfa: np.ndarray
    cat_eps: np.ndarray
    cat_capas: list  # nombre de las capas de cada categoría
    caras_h: tuple  # (celda, categoría, dirección, plana) de cada cara expuesta de hormigón (para Tmin a 25 mm)
    fondo: np.ndarray  # celdas del fondo del suelo
    g_fondo: float  # conductancia por volumen al fondo (media celda Dirichlet)
    expuesta: np.ndarray  # máscara: hormigón con alguna cara expuesta con h > 0
    dt_critico: float
    combos: tuple

    def coeficientes(self, ambiente, t_s, dx):
        h, T_eq = ambiente.exterior(t_s, self.cat_orientacion, self.cat_alfa, self.cat_eps)
        return coeficiente_U(h, self.cat_R, dx, self.cat_k), T_eq


def _etapa(dom: Dominio, activas, k, C, amb: Ambiente, tabla, combos, k_h, k_s) -> _Etapa:
    dx = dom.dx
    G = []
    suma = np.zeros(dom.forma)
    for eje in range(3):
        a = [slice(None)] * 3
        b = [slice(None)] * 3
        a[eje], b[eje] = slice(None, -1), slice(1, None)
        a, b = tuple(a), tuple(b)
        ambas = activas[a] & activas[b]
        g = np.where(ambas, 2 * k[a] * k[b] / np.maximum(k[a] + k[b], 1e-30), 0.0) / dx**2
        G.append(g)
        suma[a] += g
        suma[b] += g

    celda, direccion, orientacion = _caras_expuestas(activas)
    material = dom.material.reshape(-1)[celda]
    combo = tabla[dom.vertida_id.reshape(-1)[celda].astype(np.int64), orientacion]
    claves, categoria = np.unique(np.stack([orientacion, combo, material], axis=1), axis=0, return_inverse=True)
    categoria = categoria.reshape(-1)
    capas = [combos[c] for c in claves[:, 1]]
    cat_k = np.where(claves[:, 2] == HORMIGON, k_h, k_s).astype(float)
    cat_R = np.array([sum(c.R for c in cs) for cs in capas], dtype=float)
    cat_alfa = np.array([cs[-1].absortividad if cs else amb.absortividad for cs in capas], dtype=float)
    cat_eps = np.array([cs[-1].emisividad if cs else amb.emisividad for cs in capas], dtype=float)
    superficie, fila = np.unique(celda, return_inverse=True)
    N = np.zeros((len(superficie), len(claves)))
    np.add.at(N, (fila.reshape(-1), categoria), 1.0)

    U_max = coeficiente_U(amb.h_maximo(claves[:, 0], cat_eps), cat_R, dx, cat_k)
    suma.reshape(-1)[superficie] += N @ (U_max / dx)
    fondo = np.flatnonzero((activas & (dom.material == SUELO))[:, :, 0])
    fondo = np.ravel_multi_index(np.unravel_index(fondo, dom.forma[:2]) + (np.zeros_like(fondo),), dom.forma)
    g_fondo = 2 * k_s / dx**2
    suma.reshape(-1)[fondo] += g_fondo

    expuesta = np.zeros(dom.forma, dtype=bool)
    con_h = (N @ (U_max > 0)) > 0
    expuesta.reshape(-1)[superficie[con_h]] = True
    expuesta &= dom.material == HORMIGON
    de_hormigon = material == HORMIGON
    # Cara "plana": su celda tiene una sola cara expuesta (no es arista ni vértice de la escalera de vóxeles).
    plana = (N.sum(axis=1) == 1)[fila.reshape(-1)]
    coef = np.where(activas, suma / C, 0.0)
    return _Etapa(
        activas=activas.copy(), G=G, superficie=superficie, N=N, cat_orientacion=claves[:, 0], cat_R=cat_R,
        cat_k=cat_k, cat_alfa=cat_alfa, cat_eps=cat_eps, cat_capas=[[c.nombre for c in cs] for cs in capas],
        caras_h=(celda[de_hormigon], categoria[de_hormigon], direccion[de_hormigon], plana[de_hormigon]), fondo=fondo,
        g_fondo=g_fondo, expuesta=expuesta, dt_critico=1.0 / coef.max() if coef.max() > 0 else math.inf,
        combos=tuple(tuple(c.nombre for c in cs) for cs in combos),
    )


@dataclass
class Resultado3D:
    dominio: Dominio
    plan: PlanVertidas
    t_h: np.ndarray  # instantes guardados [h]
    Tmax: np.ndarray  # máx. en el hormigón existente [°C]
    Tmax_pos: np.ndarray  # (n, 3) posición de Tmax [m]
    Tmin_sup: np.ndarray  # mín. en la superficie del hormigón [°C] (a `profundidad_tmin` o en el centro de la celda)
    Tmax_vertida: np.ndarray  # (n, n_vertidas) [°C], NaN si no existe
    Tmedia_vertida: np.ndarray  # (n, n_vertidas)
    puntos: dict  # nombre -> T(t) [°C] (NaN mientras la celda no existe)
    posiciones: dict  # nombre -> (x, y, z) del centro de la celda usada
    T_amb: np.ndarray  # temperatura del aire [°C]
    corte_y0: np.ndarray  # (n, nx, nz) corte vertical y = 0, NaN fuera del sólido
    campos: dict  # t_h -> campo 3D float32 (NaN donde no hay material) en los instantes elegidos y el de la Tmax
    edades: dict  # t_h -> te [h] (float32) en los mismos instantes
    dts: list = field(default_factory=list)  # (t_inicio_h, dt [s], dt_critico [s]) por tramo
    etapas: list = field(default_factory=list)  # (t_h, vertidas colocadas, dt_critico [s]) al cambiar las celdas
    eventos: list = field(default_factory=list)  # inicio de colocación (t_h, v, nombre, T_colocación)
    retiros: list = field(default_factory=list)  # (t_h, descripción) de moldajes y aislaciones retirados
    Tmin_pos: np.ndarray | None = None  # (n, 3) dónde está Tmin [m]
    profundidad_tmin: float | None = None  # None = centro de la celda de superficie
    superficie_tmin: str = "todas"  # "plana": sólo caras de celdas con una única cara expuesta
    Tmin_todas: np.ndarray | None = None  # Tmin contando aristas y vértices (dato complementario)
    Tmin_todas_pos: np.ndarray | None = None
    balance: dict = field(default_factory=dict)  # balance de energía [J] (ver simular_3d)

    @property
    def dT(self) -> np.ndarray:
        return self.Tmax - self.Tmin_sup

    @property
    def criterio_Tmin(self) -> str:
        if self.profundidad_tmin is None:
            return "Tmin en el centro de las celdas de superficie"
        donde = "caras planas (sin aristas ni vértices)" if self.superficie_tmin == "plana" else "todas las caras"
        return f"Tmin a {self.profundidad_tmin * 1000:g} mm bajo las superficies expuestas, en {donde}"


def simular_3d(
    dom: Dominio,
    plan: PlanVertidas,
    hormigon: Material,
    suelo: Material,
    hidratacion,
    ambiente: Ambiente,
    duracion_h: float,
    T_suelo_inicial=15.0,
    guardar_cada_h: float = 1.0,
    dt_max_s: float | None = 600.0,
    factor_seguridad: float = 0.9,
    puntos: dict | None = None,
    instantes_3d_h=(),
    informar=None,
    observador=None,
    cubiertas: PlanCubiertas | None = None,
    profundidad_tmin: float | None = None,
    superficie_tmin: str = "todas",
    activacion: Activacion | None = None,
) -> Resultado3D:
    """Integra la fundación por etapas: entre eventos (colocación o retiro de una capa) la geometría no cambia.

    T_suelo_inicial: número o función de z [m] (perfil inicial del suelo).
    cubiertas: moldaje y aislación con su retiro (bordes.PlanCubiertas).
    profundidad_tmin: si se indica (p. ej. 0,025 m), Tmin se evalúa a esa profundidad bajo cada
    cara expuesta, interpolando entre la superficie (por continuidad del flujo) y el centro de la celda.
    superficie_tmin: "plana" toma Tmin sólo en caras de celdas con una única cara expuesta (puntos
    representativos de una superficie); "todas" incluye aristas y vértices. Tmin_todas siempre se registra.
    observador(t_h, T, activas), si se indica, se llama en cada instante guardado (ver panel.Recolector).
    activacion: hora y temperatura de colocación de cada celda (por defecto, la de `plan`).

    Balance de energía (res.balance): E_final - E_inicial = colocación + hidratación + bordes + fondo,
    con E = suma de C T V sobre las celdas que existen. Euler explícito lo cumple a precisión de máquina;
    el residuo relativo se informa como evidencia de conservación.
    """
    forma = dom.forma
    dx = dom.dx
    cubiertas = cubiertas or PlanCubiertas()
    es_h = dom.material == HORMIGON
    es_s = dom.material == SUELO
    k = np.where(es_h, hormigon.k, np.where(es_s, suelo.k, 1.0))
    C = np.where(es_h, hormigon.capacidad, np.where(es_s, suelo.capacidad, 1.0))
    zc = dom.centros(2)
    perfil = T_suelo_inicial(zc) if callable(T_suelo_inicial) else np.full(zc.shape, float(T_suelo_inicial))
    T = np.where(es_s, np.broadcast_to(perfil, forma), np.nan)
    edad = np.zeros(forma)
    if superficie_tmin not in ("plana", "todas"):
        raise ValueError("superficie_tmin debe ser 'plana' o 'todas'")
    if profundidad_tmin is not None and not 0 <= profundidad_tmin <= dx / 2:
        raise ValueError(f"profundidad_tmin debe estar entre 0 y dx/2 = {dx / 2:g} m")

    vox_id = dom.vertida_id
    activas = es_s.copy()
    ids_h = np.flatnonzero(es_h)  # sólo estas celdas hidratan: se trabaja sobre ellas
    rho_cp = hormigon.capacidad

    puntos = puntos or {}
    idx_puntos = {n: dom.indice(*p) for n, p in puntos.items()}
    posiciones = {n: tuple(float(dom.centros(e)[i[e]]) for e in range(3)) for n, i in idx_puntos.items()}
    instantes_3d = sorted(float(t) for t in instantes_3d_h)
    j0 = int(np.argmin(np.abs(dom.centros(1))))
    nv = dom.vox.n_vertidas
    mascaras_v = [vox_id == v for v in range(1, nv + 1)]
    centros = [dom.centros(e) for e in range(3)]

    reg = {k_: [] for k_ in ("t", "Tmax", "pos", "Tmin", "pos_min", "Tmin_t", "pos_min_t", "Tmaxv", "Tmedv", "amb",
                             "corte")}
    reg_p = {n: [] for n in puntos}
    campos, edades, dts, eventos, etapas, retiros = {}, {}, [], [], [], []
    mayor = {"T": -np.inf}

    def posicion(c):
        i = np.unravel_index(c, forma)
        return np.array([centros[e][i[e]] for e in range(3)])

    def superficie_fria(t_h):
        """(Tmin, posición) del criterio y (Tmin, posición) contando aristas, en el hormigón que existe."""
        nada = (np.nan, [np.nan] * 3)
        if profundidad_tmin is None:
            sup = activas & es_h & etapa.expuesta
            if not sup.any():
                return nada, nada
            Ts = np.where(sup, T, np.inf)
            c = int(np.argmin(Ts))
            r = (float(Ts.reshape(-1)[c]), posicion(c).tolist())
            return r, r
        celda, cat, direccion, plana = etapa.caras_h
        if len(celda) == 0:
            return nada, nada
        U, T_eq = etapa.coeficientes(ambiente, t_h * 3600.0, dx)
        Tc = T.reshape(-1)[celda]
        # Perfil lineal en la media celda: flujo q = U (T_eq - Tc) entra por la cara; a la profundidad p:
        # T(p) = Tc + q (dx/2 - p) / k.  (p = 0 es la superficie del hormigón, bajo el moldaje si lo hay.)
        T25 = Tc + U[cat] * (T_eq[cat] - Tc) * (dx / 2 - profundidad_tmin) / hormigon.k

        def en(f):
            return float(T25[f]), (posicion(celda[f]) + _NORMALES[direccion[f]] * (dx / 2 - profundidad_tmin)).tolist()

        todas = en(int(np.argmin(T25)))
        if superficie_tmin == "plana" and plana.any():
            return en(int(np.flatnonzero(plana)[np.argmin(T25[plana])])), todas
        return todas, todas

    def registrar(t_h):
        h_act = activas & es_h
        Th = np.where(h_act, T, -np.inf)
        if h_act.any():
            c = int(np.argmax(Th))
            reg["Tmax"].append(float(Th.reshape(-1)[c]))
            reg["pos"].append(posicion(c).tolist())
        else:
            reg["Tmax"].append(np.nan)
            reg["pos"].append([np.nan] * 3)
        (Tmin, pos_min), (Tmin_t, pos_min_t) = superficie_fria(t_h)
        reg["Tmin"].append(Tmin)
        reg["pos_min"].append(pos_min)
        reg["Tmin_t"].append(Tmin_t)
        reg["pos_min_t"].append(pos_min_t)
        reg["Tmaxv"].append([float(T[m & activas].max()) if (m & activas).any() else np.nan for m in mascaras_v])
        reg["Tmedv"].append([float(T[m & activas].mean()) if (m & activas).any() else np.nan for m in mascaras_v])
        reg["t"].append(t_h)
        reg["amb"].append(ambiente.temperatura(t_h * 3600))
        reg["corte"].append(np.where(activas[:, j0, :], T[:, j0, :], np.nan).astype(np.float32))
        for n, i in idx_puntos.items():
            reg_p[n].append(float(T[i]) if activas[i] else np.nan)
        if observador:
            observador(t_h, T, activas)
        guardar = any(abs(ti - t_h) < 1e-6 for ti in instantes_3d)
        if reg["Tmax"][-1] > mayor["T"]:  # el campo del instante de la Tmax también se guarda (al final)
            mayor.update(T=reg["Tmax"][-1], t=t_h, campo=np.where(activas, T, np.nan).astype(np.float32),
                         edad=np.where(activas & es_h, edad, np.nan).astype(np.float32))
        if guardar:
            campos[t_h] = np.where(activas, T, np.nan).astype(np.float32)
            edades[t_h] = np.where(activas & es_h, edad, np.nan).astype(np.float32)

    # Instantes de guardado, colocaciones y retiros: la integración avanza de uno a otro.
    n_g = max(1, math.ceil(duracion_h / guardar_cada_h - 1e-9))
    guardados = [min(i * guardar_cada_h, duracion_h) for i in range(n_g + 1)]
    act = activacion if activacion is not None else Activacion.desde_plan(dom, plan)
    t_act = np.where(es_h, act.t_act, np.inf)
    colocaciones = {float(t) for t in np.unique(t_act[np.isfinite(t_act)]) if 0 <= t <= duracion_h}
    retiros_plan = [(t, d) for t, d in cubiertas.eventos(act.inicio_h, act.fin_h) if 0 <= t <= duracion_h]
    V = dx**3
    bal = {"hidratacion": 0.0, "bordes": 0.0, "fondo": 0.0, "colocacion": 0.0}
    E_inicial = float(np.sum(np.where(activas, C * np.nan_to_num(T), 0.0))) * V
    hitos = sorted(set(guardados) | colocaciones | {t for t, _ in retiros_plan})

    etapa = None
    tabla_actual = None
    t_prev = -1.0
    for n_hito, t_h in enumerate(hitos):
        if n_hito > 0:  # integrar desde el hito anterior
            t0 = hitos[n_hito - 1]
            intervalo = (t_h - t0) * 3600
            dt_obj = factor_seguridad * etapa.dt_critico
            if dt_max_s is not None:
                dt_obj = min(dt_obj, dt_max_s)
            pasos = max(1, math.ceil(intervalo / dt_obj - 1e-9))
            dt = intervalo / pasos
            if not dts or dts[-1][1] != dt:
                dts.append((t0, dt, etapa.dt_critico))
            # Coeficientes del tramo (las celdas que existen no cambian dentro de él).
            # Las celdas inexistentes se llevan a 0 °C dentro del bucle: todas sus conductancias
            # son 0 y dt/C también, así que no influyen ni cambian; al final vuelven a NaN.
            dt_C = np.where(activas, dt / C, 0.0)
            T = np.where(activas, T, 0.0)
            Tv, Ev, var_v = T.reshape(-1), edad.reshape(-1), np.empty(T.size)
            var = var_v.reshape(forma)
            sup, fondo = etapa.superficie, etapa.fondo
            ids = ids_h[activas.reshape(-1)[ids_h]]  # hormigón que existe
            difs = [np.empty([n - (e == eje) for e, n in enumerate(forma)]) for eje in range(3)]
            cortes = []
            for eje in range(3):
                a = [slice(None)] * 3
                b = [slice(None)] * 3
                a[eje], b[eje] = slice(None, -1), slice(1, None)
                cortes.append((tuple(a), tuple(b)))
            for paso in range(pasos):
                var.fill(0.0)
                # Bordes al ambiente: q = U (T_eq - T) por cada cara expuesta, agrupadas por categoría.
                U, T_eq = etapa.coeficientes(ambiente, t0 * 3600 + paso * dt, dx)
                q = etapa.N @ (U * T_eq / dx) - (etapa.N @ (U / dx)) * Tv[sup]
                var_v[sup] += q
                bal["bordes"] += float(q.sum()) * V * dt
                # Fondo del suelo: Dirichlet a media celda.
                q = etapa.g_fondo * (ambiente.T_fondo - Tv[fondo])
                var_v[fondo] += q
                bal["fondo"] += float(q.sum()) * V * dt
                # Conducción entre vecinas: flujo g (Tb - Ta), igual y opuesto.
                for (a, b), g, d in zip(cortes, etapa.G, difs):
                    np.subtract(T[b], T[a], out=d)
                    d *= g
                    var[a] += d
                    var[b] -= d
                # Hidratación (explícita, con T y te del paso anterior), sólo en el hormigón que existe.
                Th, eh = Tv[ids], Ev[ids]
                f = hidratacion.factor_arrhenius(Th)
                if hasattr(hidratacion, "dtheta_dte"):  # Ju y Lei: Q = rho cp dtheta/dte dte/dt (sin repetir Arrhenius)
                    q = rho_cp * hidratacion.dtheta_dte(eh) * f / 3600.0
                else:
                    q = hidratacion.potencia(eh, Th, hormigon.rho, hormigon.cp)
                var_v[ids] += q
                bal["hidratacion"] += float(np.sum(q)) * V * dt
                Ev[ids] = eh + dt / 3600.0 * f
                var *= dt_C
                T += var
            edad = Ev.reshape(forma)
            T = np.where(activas, T, np.nan)

        # Celdas que se colocan en (t_prev, t_h]: T de colocación y te = 0 (misma regla que PlanVertidas.activar).
        celdas = (t_act > t_prev + 1e-9) & (t_act <= t_h + 1e-9)
        nuevas = bool(celdas.any())
        if nuevas:
            T[celdas] = act.T_col[celdas]
            edad[celdas] = 0.0
            activas |= celdas
            bal["colocacion"] += float(np.sum(C[celdas] * act.T_col[celdas])) * V
        for v in range(1, nv + 1):
            if t_prev + 1e-9 < act.inicio_h[v - 1] <= t_h + 1e-9:
                eventos.append((t_h, v, dom.vox.nombres[v - 1], float(act.T_col_v[v - 1])))
        for t_r, d in retiros_plan:
            if t_prev + 1e-9 < t_r <= t_h + 1e-9:
                retiros.append((t_r, d))
                if informar:
                    informar(f"t = {t_h:6.1f} h: {d}")
        tabla, combos = cubiertas.tabla(t_h, act.inicio_h, nv, act.fin_h)
        firma = (tabla.tobytes(), tuple(tuple(c.nombre for c in cs) for cs in combos))
        if nuevas or etapa is None or firma != tabla_actual:
            etapa = _etapa(dom, activas, k, C, ambiente, tabla, combos, hormigon.k, suelo.k)
            tabla_actual = firma
            etapas.append((t_h, act.vertidas_activas(t_h), float(etapa.dt_critico)))
            if informar and (any(e[0] == t_h for e in eventos) or n_hito == 0):
                informar(f"t = {t_h:6.1f} h: activas {act.vertidas_activas(t_h)}, dt crítico = {etapa.dt_critico:.0f} s")
        t_prev = t_h
        if t_h in guardados:
            registrar(t_h)

    E_final = float(np.sum(np.where(activas, C * np.nan_to_num(T), 0.0))) * V
    entradas = sum(bal.values())
    bal = {**bal, "E_inicial": E_inicial, "E_final": E_final,
           "residuo_relativo": abs(E_final - E_inicial - entradas) / max(abs(bal["hidratacion"]), abs(entradas), 1.0)}
    if "campo" in mayor:
        campos.setdefault(mayor["t"], mayor["campo"])
        edades.setdefault(mayor["t"], mayor["edad"])
        campos, edades = dict(sorted(campos.items())), dict(sorted(edades.items()))
    return Resultado3D(
        dominio=dom, plan=plan, t_h=np.array(reg["t"]), Tmax=np.array(reg["Tmax"]), Tmax_pos=np.array(reg["pos"]),
        Tmin_sup=np.array(reg["Tmin"]), Tmax_vertida=np.array(reg["Tmaxv"]), Tmedia_vertida=np.array(reg["Tmedv"]),
        puntos={n: np.array(v) for n, v in reg_p.items()}, posiciones=posiciones, T_amb=np.array(reg["amb"]),
        corte_y0=np.stack(reg["corte"]), campos=campos, edades=edades, dts=dts, eventos=eventos, etapas=etapas,
        retiros=retiros, Tmin_pos=np.array(reg["pos_min"]), profundidad_tmin=profundidad_tmin,
        superficie_tmin=superficie_tmin if profundidad_tmin is not None else "todas",
        Tmin_todas=np.array(reg["Tmin_t"]), Tmin_todas_pos=np.array(reg["pos_min_t"]), balance=bal,
    )


# ---------------------------------------------------------------------------
# Caso YAML: todo el pipeline IFC -> grilla -> solver desde un archivo de inputs
# ---------------------------------------------------------------------------

@dataclass
class Caso3D:
    nombre: str
    ifc: object  # Path
    dx: float
    hormigon: Material
    suelo: Material
    hidratacion: object
    profundidad_suelo: float
    margen_suelo: float
    T_suelo_inicial: float
    ambiente: Ambiente
    cambios: dict  # Vertida_XX -> {Inicio_h, TemperaturaColocacion_C, activa}
    duracion_h: float
    hora_inicio: float
    factor_seguridad: float
    dt_max_s: float | None
    guardar_cada_h: float
    puntos: dict
    instantes_3d_h: list
    salida: object  # Path
    limite_Tmax: float = 70.0
    limite_dT: float = 20.0
    T_amb_texto: str = ""  # cómo se declaró T ambiente en el YAML
    profundidad_tmin: float | None = None  # criterios.profundidad_superficie_m (None = centro de la celda)
    superficie_tmin: str = "todas"  # criterios.superficie_tmin: "plana" o "todas"
    fecha_inicio: str | None = None  # hora local de la colocación de la Vertida_01 (con clima real)
    clima_cfg: dict | None = None  # sección `clima` del YAML (rutas ya resueltas)
    superficie_cfg: dict = field(default_factory=dict)
    capas: dict = field(default_factory=dict)  # catálogo de moldajes y aislaciones
    cubiertas_cfg: dict = field(default_factory=dict)  # secciones `moldaje` y `aislacion`
    suelo_cfg: dict = field(default_factory=dict)
    cubiertas: PlanCubiertas | None = None  # se arma en preparar() (necesita el pset del IFC)
    clima_info: dict = field(default_factory=dict)


def _texto_T_amb(d) -> str:
    if isinstance(d, (int, float)):
        return f"T_amb = {d:g} °C"
    return f"T_amb = {d['media']:g} ± {d['amplitud']:g} °C (sinusoide diaria, máx. a las {d['hora_max']:g} h)"


def _ruta_o_ultimo(valor, base):
    """'ultimo' (el snapshot más reciente), null, o una ruta relativa al YAML."""
    import os
    from pathlib import Path

    if valor in (None, "ninguno", False):
        return None
    if valor == "ultimo":
        return "ultimo"
    return Path(os.path.normpath(base / valor))


def cargar_caso_3d(ruta) -> Caso3D:
    """Lee un caso 3D (config/caso_3d.yaml, config/caso_14d.yaml); rutas relativas a la carpeta del YAML."""
    import os
    from pathlib import Path

    import yaml

    from .caso import temperatura_ambiente
    from .hidratacion import HidratacionJuLei2019

    ruta = Path(ruta)
    cfg = yaml.safe_load(ruta.read_text(encoding="utf-8"))
    base = ruta.parent
    mats = cfg["materiales"]
    h = dict(mats["hormigon"].get("hidratacion", {"modelo": "ju_lei_2019"}))
    if h.pop("modelo") != "ju_lei_2019":
        raise ValueError("Modelo de hidratación no implementado")
    tiempo, suelo, sal = cfg["tiempo"], cfg.get("suelo", {}), cfg.get("salida", {})
    amb = cfg.get("ambiente", {})
    fecha = cfg.get("fecha_inicio")
    hora_inicio = float(tiempo.get("hora_inicio", int(str(fecha)[11:13]) if fecha else 0.0))
    crit = cfg.get("criterios", {})
    clima_cfg = None
    if "clima" in cfg:
        if not fecha:
            raise ValueError("Con `clima` hay que indicar fecha_inicio (hora local de la primera vertida)")
        c = dict(cfg["clima"])
        clima_cfg = {"pronostico": _ruta_o_ultimo(c.get("pronostico", "ultimo"), base),
                     "climatologia": _ruta_o_ultimo(c.get("climatologia", "ultimo"), base),
                     "carpeta": Path(os.path.normpath(base / c.get("carpeta", "../data/clima/snapshots"))),
                     "transicion_h": float(c.get("transicion_h", 24)), "ventana_dias": int(c.get("ventana_dias", 3))}
    superficie = dict(cfg.get("superficie", {}))
    capas = {n: Capa(n, float(d["espesor"]), float(d["k"]), float(d.get("absortividad", 0.6)),
                     float(d.get("emisividad", 0.9))) for n, d in (cfg.get("capas") or {}).items()}
    T_suelo = suelo.get("T_inicial", 15.0)
    T_fondo = suelo.get("T_fondo", 15.0)
    if "T_amb" in amb:
        ambiente = Ambiente(temperatura_ambiente(amb["T_amb"], hora_inicio), float(amb.get("h_superior", 27.8)),
                            float(amb.get("h_lateral", 27.8)), float(T_fondo) if T_fondo != "era5" else 15.0)
        texto = _texto_T_amb(amb["T_amb"])
    else:
        ambiente, texto = Ambiente(20.0), "clima real (ver clima_info)"
    return Caso3D(
        nombre=cfg.get("nombre", ruta.stem),
        ifc=Path(os.path.normpath(base / cfg["ifc"])),
        dx=float(cfg["malla"]["dx"]),
        hormigon=Material(*(float(mats["hormigon"][c]) for c in ("k", "rho", "cp"))),
        suelo=Material(*(float(mats["suelo"][c]) for c in ("k", "rho", "cp"))),
        hidratacion=HidratacionJuLei2019(**h),
        profundidad_suelo=float(suelo.get("profundidad", 3.0)),
        margen_suelo=float(suelo.get("margen", 2.0)),
        T_suelo_inicial=T_suelo if T_suelo == "era5" else float(T_suelo),
        ambiente=ambiente,
        cambios=dict(cfg.get("vertidas") or {}),
        duracion_h=float(tiempo["duracion_h"]),
        hora_inicio=hora_inicio,
        factor_seguridad=float(tiempo.get("factor_seguridad", 0.9)),
        dt_max_s=None if tiempo.get("dt_max_s") is None else float(tiempo["dt_max_s"]),
        guardar_cada_h=float(tiempo.get("guardar_cada_h", 1.0)),
        puntos={n: tuple(map(float, p)) for n, p in (cfg.get("puntos") or {}).items()},
        instantes_3d_h=[float(t) for t in sal.get("instantes_3d_h", [])],
        salida=Path(os.path.normpath(base / sal.get("carpeta", f"../resultados/{ruta.stem}"))),
        limite_Tmax=float(crit.get("Tmax_C", 70.0)),
        limite_dT=float(crit.get("dT_C", 20.0)),
        T_amb_texto=texto,
        profundidad_tmin=None if crit.get("profundidad_superficie_m") is None
        else float(crit["profundidad_superficie_m"]),
        superficie_tmin=str(crit.get("superficie_tmin", "todas")),
        fecha_inicio=None if fecha is None else str(fecha), clima_cfg=clima_cfg, superficie_cfg=superficie,
        capas=capas, cubiertas_cfg={t: cfg[t] for t in ("moldaje", "aislacion") if cfg.get(t)},
        suelo_cfg=dict(suelo),
    )


def preparar(caso: Caso3D, vertidas_ifc=None):
    """IFC -> voxelización -> cambios a la secuencia -> dominio con suelo -> plan. Devuelve (vox, dom, plan, registro)."""
    from .ifc_io import leer_vertidas
    from .voxelizacion import voxelizar

    vertidas_ifc = leer_vertidas(caso.ifc) if vertidas_ifc is None else vertidas_ifc
    if not vertidas_ifc:
        raise ValueError(f"{caso.ifc} no tiene vertidas (Name = Vertida_01 ...)")
    vox = voxelizar(vertidas_ifc, caso.dx)
    registro = aplicar_cambios(vox, caso.cambios)
    dom = construir_dominio(vox, caso.profundidad_suelo, caso.margen_suelo)
    caso.cubiertas = armar_cubiertas(caso, vox)
    if caso.clima_cfg is not None:
        preparar_clima(caso)
    return vox, dom, PlanVertidas(vox), registro


_CARAS = {"laterales": (LAT,), "superiores": (SUP,), "todas": (SUP, LAT, INF)}


def armar_cubiertas(caso: Caso3D, vox: Voxelizacion) -> PlanCubiertas:
    """Moldaje y aislación del YAML; el retiro puede venir del pset del IFC (RetiroMoldaje_h)."""
    cubiertas = []
    for tipo, d in caso.cubiertas_cfg.items():
        nombre = d.get("capa")
        if nombre in (None, "ninguna", "ninguno"):
            continue
        if nombre not in caso.capas:
            raise ValueError(f"{tipo}: la capa {nombre!r} no está en la sección `capas` del YAML")
        vertidas = tuple(int(v) for v in d.get("vertidas", range(1, vox.n_vertidas + 1)))
        retiro = d.get("retiro_h", "pset" if tipo == "moldaje" else None)
        if retiro == "pset":
            retiro_h = {v: float(vox.propiedades[v - 1].get("RetiroMoldaje_h", np.inf)) for v in vertidas}
        elif retiro in (None, "nunca"):
            retiro_h = {v: np.inf for v in vertidas}
        else:
            retiro_h = {v: float(retiro) for v in vertidas}
        caras = d.get("caras", "laterales" if tipo == "moldaje" else "superiores")
        cubiertas.append(Cubierta(tipo, caso.capas[nombre], _CARAS[caras], vertidas, retiro_h))
    return PlanCubiertas(cubiertas)


def preparar_clima(caso: Caso3D) -> None:
    """Arma la serie horaria del clima (pronóstico + climatología) y el perfil del suelo; deja todo en caso."""
    from .clima import ensamblar, leer_snapshot, perfil_suelo, ultimo_snapshot

    cfg = caso.clima_cfg
    rutas = {}
    for clave, tipo in (("pronostico", "pronostico"), ("climatologia", "reanalisis")):
        r = cfg[clave]
        rutas[clave] = ultimo_snapshot(tipo, cfg["carpeta"]) if r == "ultimo" else r
    snaps = {k: None if r is None else leer_snapshot(r) for k, r in rutas.items()}
    clima = ensamblar(caso.fecha_inicio, caso.duracion_h, snaps["pronostico"], snaps["climatologia"],
                      cfg["transicion_h"], cfg["ventana_dias"])
    sup = caso.superficie_cfg
    T_fondo = caso.suelo_cfg.get("T_fondo", 15.0)
    info_suelo = None
    if snaps["climatologia"] is not None and "era5" in (T_fondo, caso.T_suelo_inicial):
        perfil, T_prof, info_suelo = perfil_suelo(snaps["climatologia"], caso.fecha_inicio, cfg["ventana_dias"])
        if caso.T_suelo_inicial == "era5":
            caso.T_suelo_inicial = perfil
        if T_fondo == "era5":
            T_fondo = T_prof
    if caso.T_suelo_inicial == "era5" or T_fondo == "era5":
        raise ValueError("El suelo 'era5' necesita un snapshot de climatología ERA5")
    caso.ambiente = Ambiente(20.0, float(sup.get("h_superior", 27.8)), float(sup.get("h_lateral", 27.8)),
                             float(T_fondo), modelo=sup.get("modelo", "fisico"), clima=clima,
                             absortividad=float(sup.get("absortividad", 0.65)),
                             emisividad=float(sup.get("emisividad", 0.9)))
    caso.clima_info = {**clima.metadata, "suelo": info_suelo}
    T = clima.T_aire
    caso.T_amb_texto = (f"clima real: T aire {T.min():.1f}–{T.max():.1f} °C (media {T.mean():.1f}); "
                        f"{clima.metadata['horas_pronostico']} h de pronóstico, "
                        f"{clima.metadata['horas_climatologia']} h de climatología")


def simular_caso_3d(caso: Caso3D, vertidas_ifc=None, informar=None):
    vox, dom, plan, registro = preparar(caso, vertidas_ifc)
    res = simular_3d(
        dom, plan, caso.hormigon, caso.suelo, caso.hidratacion, caso.ambiente, caso.duracion_h,
        T_suelo_inicial=caso.T_suelo_inicial, guardar_cada_h=caso.guardar_cada_h, dt_max_s=caso.dt_max_s,
        factor_seguridad=caso.factor_seguridad, puntos=caso.puntos, instantes_3d_h=caso.instantes_3d_h,
        informar=informar, cubiertas=caso.cubiertas, profundidad_tmin=caso.profundidad_tmin,
        superficie_tmin=caso.superficie_tmin,
    )
    return res, registro
