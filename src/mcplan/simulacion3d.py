"""Simulación térmica 3D de la fundación por vertidas: IFC -> vóxeles -> solver por etapas.

Es el mismo balance de energía por celda de solver_nd (volúmenes finitos, Euler
explícito, media armónica de k en las caras, Robin con media celda en serie),
pero sobre la grilla voxelizada del IFC y con celdas que **aparecen en el tiempo**:

- cada celda tiene un material: VACIO_M (aire), SUELO o HORMIGON;
- una celda de hormigón existe sólo desde el Inicio_h de su vertida: antes no
  conduce, no genera calor y no tiene edad;
- una cara entre dos celdas que existen conduce; una cara entre una celda que
  existe y una que no (aire, o una vertida que todavía no se coloca) es una
  superficie expuesta: Robin con el aire, h superior en ±z, h lateral en ±x, ±y;
- el suelo ocupa toda la planta de la grilla bajo z = 0, con T_fondo (Dirichlet)
  en su cara inferior y caras laterales adiabáticas (dominio "lejano").

La secuencia constructiva NO está en el solver: llega como un PlanVertidas
(Inicio_h y temperatura de colocación de cada vertida), que se arma con el pset
del IFC y lo que el YAML o la línea de comandos cambien (`aplicar_cambios`).
El solver sólo integra por tramos entre eventos de colocación y recalcula las
conductancias cada vez que cambia el conjunto de celdas que existen.
"""

import math
from dataclasses import dataclass, field

import numpy as np

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
class Ambiente:
    """Condiciones ambientales simplificadas (declaradas en el YAML)."""

    T_amb: object  # °C: número o función de t [s]
    h_superior: float = 27.8  # W/(m2 K) en caras ±z expuestas
    h_lateral: float = 27.8  # W/(m2 K) en caras ±x, ±y expuestas (moldaje + aire)
    T_fondo: float = 15.0  # °C en la cara inferior del suelo (Dirichlet)

    def temperatura(self, t: float) -> float:
        return float(self.T_amb(t)) if callable(self.T_amb) else float(self.T_amb)


@dataclass
class _Etapa:
    activas: np.ndarray  # celdas que existen
    G: list  # conductancia por volumen [W/(m3 K)] entre vecinas, por eje
    G_aire: np.ndarray  # suma de conductancias a T_amb [W/(m3 K)]
    G_fondo: np.ndarray  # conductancia a T_fondo [W/(m3 K)]
    dt_critico: float

    @property
    def diagonal(self) -> np.ndarray:
        return self.G_aire + self.G_fondo


def _etapa(dom: Dominio, activas, k, C, amb: Ambiente) -> _Etapa:
    dx = dom.dx
    G, G_aire = [], np.zeros(dom.forma)
    for eje in range(3):
        a = [slice(None)] * 3
        b = [slice(None)] * 3
        a[eje], b[eje] = slice(None, -1), slice(1, None)
        a, b = tuple(a), tuple(b)
        ambas = activas[a] & activas[b]
        k_cara = np.where(ambas, 2 * k[a] * k[b] / np.maximum(k[a] + k[b], 1e-30), 0.0)
        G.append(k_cara / dx**2)
        # Caras expuestas: una celda activa con la vecina inexistente.
        h = amb.h_superior if eje == 2 else amb.h_lateral
        if h > 0:
            for lado, otro in ((a, b), (b, a)):
                expuesta = activas[lado] & ~activas[otro]
                h_ef = 1.0 / (1.0 / h + dx / (2 * np.where(expuesta, k[lado], 1.0)))
                G_aire[lado] += np.where(expuesta, h_ef / dx, 0.0)
        # Las caras en el límite de la grilla son adiabáticas, salvo el fondo del suelo.
    G_fondo = np.zeros(dom.forma)
    fondo = activas[:, :, 0] & (dom.material[:, :, 0] == SUELO)
    G_fondo[:, :, 0] = np.where(fondo, 2 * k[:, :, 0] / dx**2, 0.0)

    suma = G_aire + G_fondo
    for eje, g in enumerate(G):
        a = [slice(None)] * 3
        b = [slice(None)] * 3
        a[eje], b[eje] = slice(None, -1), slice(1, None)
        suma[tuple(a)] += g
        suma[tuple(b)] += g
    coef = np.where(activas, suma / C, 0.0)
    dt_crit = 1.0 / coef.max() if coef.max() > 0 else math.inf
    return _Etapa(activas, G, G_aire, G_fondo, dt_crit)


@dataclass
class Resultado3D:
    dominio: Dominio
    plan: PlanVertidas
    t_h: np.ndarray  # instantes guardados [h]
    Tmax: np.ndarray  # máx. en el hormigón existente [°C]
    Tmax_pos: np.ndarray  # (n, 3) posición de Tmax [m]
    Tmin_sup: np.ndarray  # mín. en celdas de hormigón con alguna cara expuesta [°C]
    Tmax_vertida: np.ndarray  # (n, n_vertidas) [°C], NaN si no existe
    Tmedia_vertida: np.ndarray  # (n, n_vertidas)
    puntos: dict  # nombre -> T(t) [°C] (NaN mientras la celda no existe)
    posiciones: dict  # nombre -> (x, y, z) del centro de la celda usada
    T_amb: np.ndarray  # temperatura ambiente [°C]
    corte_y0: np.ndarray  # (n, nx, nz) corte vertical y = 0, NaN fuera del sólido
    campos: dict  # t_h -> campo 3D float32 (NaN donde no hay material) en los instantes elegidos y el de la Tmax
    edades: dict  # t_h -> te [h] (float32) en los mismos instantes
    dts: list = field(default_factory=list)  # (t_inicio_h, dt [s], dt_critico [s]) por tramo
    etapas: list = field(default_factory=list)  # (t_h, vertidas colocadas, dt_critico [s]) al cambiar las celdas
    eventos: list = field(default_factory=list)

    @property
    def dT(self) -> np.ndarray:
        return self.Tmax - self.Tmin_sup


def simular_3d(
    dom: Dominio,
    plan: PlanVertidas,
    hormigon: Material,
    suelo: Material,
    hidratacion,
    ambiente: Ambiente,
    duracion_h: float,
    T_suelo_inicial: float = 15.0,
    guardar_cada_h: float = 1.0,
    dt_max_s: float | None = 600.0,
    factor_seguridad: float = 0.9,
    puntos: dict | None = None,
    instantes_3d_h=(),
    informar=None,
    observador=None,
) -> Resultado3D:
    """Integra la fundación por etapas: entre eventos de colocación las celdas que existen no cambian.

    observador(t_h, T, activas), si se indica, se llama en cada instante guardado (ver panel.Recolector).
    """
    forma = dom.forma
    es_h = dom.material == HORMIGON
    es_s = dom.material == SUELO
    k = np.where(es_h, hormigon.k, np.where(es_s, suelo.k, 1.0))
    C = np.where(es_h, hormigon.capacidad, np.where(es_s, suelo.capacidad, 1.0))
    T = np.where(es_s, T_suelo_inicial, np.nan)
    edad = np.zeros(forma)

    # La vertida de cada celda en el dominio (con el relleno de suelo y margen).
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
    expuesta = np.zeros(forma, dtype=bool)

    reg = {k_: [] for k_ in ("t", "Tmax", "pos", "Tmin", "Tmaxv", "Tmedv", "amb", "corte")}
    reg_p = {n: [] for n in puntos}
    campos, edades, dts, eventos, etapas = {}, {}, [], [], []
    mayor = {"T": -np.inf}

    def registrar(t_h):
        h_act = activas & es_h
        Th = np.where(h_act, T, -np.inf)
        if h_act.any():
            i = np.unravel_index(int(np.argmax(Th)), forma)
            reg["Tmax"].append(float(Th[i]))
            reg["pos"].append([float(dom.centros(e)[i[e]]) for e in range(3)])
            sup = h_act & expuesta
            reg["Tmin"].append(float(T[sup].min()) if sup.any() else np.nan)
        else:
            reg["Tmax"].append(np.nan)
            reg["pos"].append([np.nan] * 3)
            reg["Tmin"].append(np.nan)
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

    # Instantes de guardado y de colocación: la integración avanza de uno a otro.
    n_g = max(1, math.ceil(duracion_h / guardar_cada_h - 1e-9))
    guardados = [min(i * guardar_cada_h, duracion_h) for i in range(n_g + 1)]
    colocaciones = sorted({t for t, *_ in plan.eventos() if 0 <= t <= duracion_h})
    hitos = sorted(set(guardados) | set(colocaciones))

    etapa = None
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
            fuente_fija = etapa.G_fondo * ambiente.T_fondo
            diagonal = etapa.diagonal
            ids = ids_h[activas.reshape(-1)[ids_h]]  # hormigón que existe
            difs = [np.empty([n - (e == eje) for e, n in enumerate(forma)]) for eje in range(3)]
            cortes = []
            for eje in range(3):
                a = [slice(None)] * 3
                b = [slice(None)] * 3
                a[eje], b[eje] = slice(None, -1), slice(1, None)
                cortes.append((tuple(a), tuple(b)))
            for paso in range(pasos):
                T_amb = ambiente.temperatura(t0 * 3600 + paso * dt)
                # Bordes: Robin con el aire y Dirichlet en el fondo del suelo.
                np.multiply(diagonal, T, out=var)
                np.subtract(fuente_fija, var, out=var)
                var += etapa.G_aire * T_amb
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
                    var_v[ids] += rho_cp * hidratacion.dtheta_dte(eh) * f / 3600.0
                else:
                    var_v[ids] += hidratacion.potencia(eh, Th, hormigon.rho, hormigon.cp)
                Ev[ids] = eh + dt / 3600.0 * f
                var *= dt_C
                T += var
            edad = Ev.reshape(forma)
            T = np.where(activas, T, np.nan)

        # Vertidas que se colocan en (t_prev, t_h]: T de colocación y te = 0 (misma regla que PlanVertidas.activar).
        nuevas = [v for v in plan.vertidas_activas(t_h) if plan.inicio_h[v - 1] > t_prev + 1e-9]
        if nuevas or etapa is None:
            for v in nuevas:
                celdas = vox_id == v
                T[celdas] = plan.T_colocacion[v - 1]
                edad[celdas] = 0.0
                activas |= celdas
                eventos.append((t_h, v, dom.vox.nombres[v - 1], float(plan.T_colocacion[v - 1])))
            etapa = _etapa(dom, activas, k, C, ambiente)
            expuesta = etapa.G_aire > 0
            etapas.append((t_h, plan.vertidas_activas(t_h), float(etapa.dt_critico)))
            if informar:
                informar(f"t = {t_h:6.1f} h: activas {plan.vertidas_activas(t_h)}, dt crítico = {etapa.dt_critico:.0f} s")
        t_prev = t_h
        if t_h in guardados:
            registrar(t_h)

    if "campo" in mayor:
        campos.setdefault(mayor["t"], mayor["campo"])
        edades.setdefault(mayor["t"], mayor["edad"])
        campos, edades = dict(sorted(campos.items())), dict(sorted(edades.items()))
    return Resultado3D(
        dominio=dom, plan=plan, t_h=np.array(reg["t"]), Tmax=np.array(reg["Tmax"]), Tmax_pos=np.array(reg["pos"]),
        Tmin_sup=np.array(reg["Tmin"]), Tmax_vertida=np.array(reg["Tmaxv"]), Tmedia_vertida=np.array(reg["Tmedv"]),
        puntos={n: np.array(v) for n, v in reg_p.items()}, posiciones=posiciones, T_amb=np.array(reg["amb"]),
        corte_y0=np.stack(reg["corte"]), campos=campos, edades=edades, dts=dts, eventos=eventos, etapas=etapas,
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


def _texto_T_amb(d) -> str:
    if isinstance(d, (int, float)):
        return f"T_amb = {d:g} °C"
    return f"T_amb = {d['media']:g} ± {d['amplitud']:g} °C (sinusoide diaria, máx. a las {d['hora_max']:g} h)"


def cargar_caso_3d(ruta) -> Caso3D:
    """Lee config/caso_3d.yaml (rutas relativas a la carpeta del YAML)."""
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
    tiempo, amb, suelo, sal = cfg["tiempo"], cfg["ambiente"], cfg.get("suelo", {}), cfg.get("salida", {})
    hora_inicio = float(tiempo.get("hora_inicio", 0.0))
    crit = cfg.get("criterios", {})
    return Caso3D(
        nombre=cfg.get("nombre", ruta.stem),
        ifc=Path(os.path.normpath(base / cfg["ifc"])),
        dx=float(cfg["malla"]["dx"]),
        hormigon=Material(*(float(mats["hormigon"][c]) for c in ("k", "rho", "cp"))),
        suelo=Material(*(float(mats["suelo"][c]) for c in ("k", "rho", "cp"))),
        hidratacion=HidratacionJuLei2019(**h),
        profundidad_suelo=float(suelo.get("profundidad", 3.0)),
        margen_suelo=float(suelo.get("margen", 2.0)),
        T_suelo_inicial=float(suelo.get("T_inicial", 15.0)),
        ambiente=Ambiente(temperatura_ambiente(amb["T_amb"], hora_inicio), float(amb.get("h_superior", 27.8)),
                          float(amb.get("h_lateral", 27.8)), float(suelo.get("T_fondo", 15.0))),
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
        T_amb_texto=_texto_T_amb(amb["T_amb"]),
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
    return vox, dom, PlanVertidas(vox), registro


def simular_caso_3d(caso: Caso3D, vertidas_ifc=None, informar=None):
    vox, dom, plan, registro = preparar(caso, vertidas_ifc)
    res = simular_3d(
        dom, plan, caso.hormigon, caso.suelo, caso.hidratacion, caso.ambiente, caso.duracion_h,
        T_suelo_inicial=caso.T_suelo_inicial, guardar_cada_h=caso.guardar_cada_h, dt_max_s=caso.dt_max_s,
        factor_seguridad=caso.factor_seguridad, puntos=caso.puntos, instantes_3d_h=caso.instantes_3d_h,
        informar=informar,
    )
    return res, registro
