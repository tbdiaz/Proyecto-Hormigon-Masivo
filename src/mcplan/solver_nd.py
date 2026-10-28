"""Solver termico 2D/3D mediante volumenes finitos y Euler explicito.

Incluye materiales homogeneos y heterogeneos, fuente volumetrica,
condiciones Dirichlet, Neumann y Robin, y superficies adiabaticas.
"""

import math

import numpy as np


def paso_critico(alpha: float, espaciamientos: tuple[float, ...]) -> float:
    """Paso critico clasico de una celda interior (material homogeneo).

    dt <= 1 / (2 alpha sum(1/d^2)). Junto a una cara Dirichlet o Robin el
    limite es menor, porque la cara esta a media celda: el limite real de
    una malla con bordes lo entrega paso_critico_malla.
    """
    if alpha <= 0 or len(espaciamientos) not in (2, 3):
        raise ValueError("Alpha y dimension deben ser validos")
    if any(d <= 0 for d in espaciamientos):
        raise ValueError("Los espaciamientos deben ser positivos")

    return 1.0 / (
        2.0 * alpha * sum(1.0 / d**2 for d in espaciamientos)
    )


def _discretizar(
    T_inicial,
    alpha,
    espaciamientos,
    fuente,
    rho,
    cp,
    bordes,
    conductividad,
):
    """Valida los datos y precalcula capacidades, conductancias y caras.

    suma_coeficientes guarda, por celda, sum(G)/C [1/s]: el coeficiente
    que multiplica a T_i en el balance. Su maximo fija el paso critico.
    """
    T = np.array(T_inicial, dtype=float, copy=True)

    if T.ndim not in (2, 3):
        raise ValueError("T debe tener 2 o 3 dimensiones")
    if len(espaciamientos) != T.ndim:
        raise ValueError("Debe existir un espaciamiento por eje")
    if any(n < 2 for n in T.shape):
        raise ValueError("Cada eje requiere al menos dos celdas")
    if not np.all(np.isfinite(T)):
        raise ValueError("Temperaturas iniciales invalidas")

    if alpha <= 0 or not np.isfinite(alpha):
        raise ValueError("Alpha debe ser positivo y finito")
    if any(d <= 0 or not np.isfinite(d) for d in espaciamientos):
        raise ValueError("Espaciamientos invalidos")

    def propiedad_malla(valor, nombre):
        """Convierte una propiedad escalar o matricial a la malla."""
        try:
            arreglo = np.broadcast_to(
                np.asarray(valor, dtype=float), T.shape
            ).copy()
        except ValueError as error:
            raise ValueError(
                f"{nombre} debe ser escalar o tener forma {T.shape}"
            ) from error

        if not np.all(np.isfinite(arreglo)):
            raise ValueError(f"{nombre} contiene valores no finitos")
        return arreglo

    rho_m = propiedad_malla(rho, "rho")
    cp_m = propiedad_malla(cp, "cp")

    if np.any(rho_m <= 0) or np.any(cp_m <= 0):
        raise ValueError("rho y cp deben ser positivos")

    capacidad = rho_m * cp_m

    if conductividad is None:
        k = alpha * capacidad
    else:
        k = propiedad_malla(conductividad, "conductividad")

    if np.any(k <= 0) or not np.all(np.isfinite(k)):
        raise ValueError("La conductividad debe ser positiva")

    fuente_m = propiedad_malla(fuente, "fuente")
    bordes = {} if bordes is None else dict(bordes)

    nombres = ("x", "y", "z")[:T.ndim]
    caras_validas = {
        f"{eje}_{lado}"
        for eje in nombres
        for lado in ("min", "max")
    }

    for cara, condicion in bordes.items():
        if cara not in caras_validas:
            raise ValueError(f"Cara desconocida: {cara}")

        if not isinstance(condicion, (tuple, list)):
            raise ValueError(f"Condicion invalida: {cara}")
        if len(condicion) != 2:
            raise ValueError(f"Condicion invalida: {cara}")

        tipo, dato = condicion

        if tipo in ("dirichlet", "neumann"):
            if not np.isscalar(dato) or not np.isfinite(dato):
                raise ValueError(f"Valor invalido en {cara}")

        elif tipo == "robin":
            if not isinstance(dato, (tuple, list)) or len(dato) != 2:
                raise ValueError(f"Robin requiere (h, T_amb) en {cara}")
            h, T_amb = dato
            if h < 0 or not np.isfinite(h) or not np.isfinite(T_amb):
                raise ValueError(f"Robin invalido en {cara}")

        else:
            raise ValueError(f"Tipo desconocido en {cara}: {tipo}")

    # Precalcular los intercambios entre celdas.
    # La misma conductancia se utiliza en ambos sentidos.
    conexiones = []
    suma_coeficientes = np.zeros_like(T)

    for eje, d in enumerate(espaciamientos):
        izq = [slice(None)] * T.ndim
        der = [slice(None)] * T.ndim

        izq[eje] = slice(None, -1)
        der[eje] = slice(1, None)

        izq = tuple(izq)
        der = tuple(der)

        k_izq = k[izq]
        k_der = k[der]

        # Conductividad equivalente por resistencias en serie.
        k_equivalente = (
            2.0 * k_izq * k_der / (k_izq + k_der)
        )

        # Conductancia por unidad de volumen, para celdas uniformes.
        coef = k_equivalente / d**2

        conexiones.append((izq, der, coef))

        suma_coeficientes[izq] += coef / capacidad[izq]
        suma_coeficientes[der] += coef / capacidad[der]

    # Precalcular las condiciones en las caras exteriores.
    caras = []

    for eje, d in enumerate(espaciamientos):
        for lado in ("min", "max"):
            cara = f"{nombres[eje]}_{lado}"
            tipo, dato = bordes.get(cara, ("neumann", 0.0))

            indice = [slice(None)] * T.ndim
            indice[eje] = 0 if lado == "min" else -1
            indice = tuple(indice)

            if tipo == "dirichlet":
                # Media celda hasta la superficie.
                coef = 2.0 * k[indice] / (
                    capacidad[indice] * d**2
                )

                suma_coeficientes[indice] += coef
                caras.append((indice, tipo, coef, float(dato)))

            elif tipo == "neumann":
                coef = 1.0 / (capacidad[indice] * d)
                caras.append((indice, tipo, coef, float(dato)))

            elif tipo == "robin":
                h, T_amb = dato

                if h == 0:
                    coef = np.zeros_like(k[indice])
                else:
                    # Resistencia de media celda y conveccion en serie.
                    h_efectivo = 1.0 / (
                        1.0 / h + d / (2.0 * k[indice])
                    )
                    coef = h_efectivo / (
                        capacidad[indice] * d
                    )

                suma_coeficientes[indice] += coef
                caras.append((indice, tipo, coef, float(T_amb)))

    return T, capacidad, fuente_m, conexiones, caras, suma_coeficientes


def _dt_critico(suma_coeficientes: np.ndarray) -> float:
    # Euler explicito: el coeficiente de T_i^n es 1 - dt * sum(G)/C >= 0.
    max_coef = float(np.max(suma_coeficientes))
    return 1.0 / max_coef if max_coef > 0 else math.inf


def paso_critico_malla(
    T_inicial: np.ndarray,
    alpha: float,
    espaciamientos: tuple[float, ...],
    rho: float = 2400.0,
    cp: float = 900.0,
    bordes: dict | None = None,
    conductividad=None,
) -> float:
    """Paso critico real de la malla: dt <= min_i C_i / sum(G_i).

    A diferencia de paso_critico, considera materiales distintos y las
    caras Dirichlet o Robin. Recibe los mismos datos que simular_difusion.
    Es una condicion suficiente (positividad), no exacta.
    """
    *_, suma_coeficientes = _discretizar(
        T_inicial, alpha, espaciamientos, 0.0,
        rho, cp, bordes, conductividad,
    )
    return _dt_critico(suma_coeficientes)


def simular_difusion(
    T_inicial: np.ndarray,
    alpha: float,
    espaciamientos: tuple[float, ...],
    dt: float,
    pasos: int,
    fuente: float = 0.0,
    rho: float = 2400.0,
    cp: float = 900.0,
    bordes: dict | None = None,
    conductividad=None,
    forzar: bool = False,
) -> np.ndarray:
    """Resuelve conduccion termica en malla uniforme 2D/3D.

    rho, cp y conductividad pueden ser escalares o matrices
    con la misma forma que T_inicial.

    Si conductividad no se especifica:
        k = alpha * rho * cp

    Convencion de bordes:
        ("dirichlet", temperatura_superficie)
        ("neumann", flujo_entrante)
        ("robin", (h, temperatura_ambiente))

    El flujo Neumann positivo entra al material.
    Las caras no especificadas son adiabaticas.

    Con forzar=True se acepta un dt inestable (para mostrar la
    divergencia); la integracion se detiene si |T| supera 1e4.
    """
    if dt <= 0 or not np.isfinite(dt):
        raise ValueError("dt debe ser positivo y finito")
    if not isinstance(pasos, (int, np.integer)) or pasos < 0:
        raise ValueError("pasos debe ser un entero no negativo")

    T, capacidad, fuente_m, conexiones, caras, suma_coeficientes = (
        _discretizar(
            T_inicial, alpha, espaciamientos, fuente,
            rho, cp, bordes, conductividad,
        )
    )

    # Limite local de estabilidad de Euler explicito.
    dt_critico = _dt_critico(suma_coeficientes)

    if dt > dt_critico * (1.0 + 1e-12) and not forzar:
        raise ValueError(
            f"dt={dt} supera el paso critico "
            f"de {dt_critico:.6g} s"
        )

    # Integracion temporal.
    for _ in range(pasos):
        variacion = fuente_m / capacidad

        # Flujos entre celdas: iguales y opuestos.
        for izq, der, coef in conexiones:
            diferencia = T[der] - T[izq]
            flujo = coef * diferencia

            variacion[izq] += flujo / capacidad[izq]
            variacion[der] -= flujo / capacidad[der]

        # Flujos impuestos en las superficies exteriores.
        for indice, tipo, coef, dato in caras:
            if tipo == "dirichlet":
                variacion[indice] += (
                    coef * (dato - T[indice])
                )

            elif tipo == "neumann":
                variacion[indice] += coef * dato

            elif tipo == "robin":
                variacion[indice] += (
                    coef * (dato - T[indice])
                )

        T = T + dt * variacion

        if forzar and not np.max(np.abs(T)) < 1e4:
            break  # divergio

    return T


def simular_hasta(
    T_inicial: np.ndarray,
    alpha: float,
    espaciamientos: tuple[float, ...],
    t_final: float,
    fuente: float = 0.0,
    rho: float = 2400.0,
    cp: float = 900.0,
    bordes: dict | None = None,
    conductividad=None,
    factor_seguridad: float = 0.9,
) -> tuple[np.ndarray, float]:
    """Simula hasta t_final con dt automatico. Devuelve (T_final, dt).

    Usa el menor numero de pasos tal que dt = t_final / pasos cumpla
    dt <= factor_seguridad * dt_critico. Asi el tiempo final es exacto y
    se pueden cambiar alpha, la malla o los bordes sin recalcular dt.
    """
    if t_final <= 0 or not np.isfinite(t_final):
        raise ValueError("t_final debe ser positivo y finito")
    if not 0 < factor_seguridad <= 1:
        raise ValueError("factor_seguridad debe estar en (0, 1]")

    dt_critico = paso_critico_malla(
        T_inicial, alpha, espaciamientos,
        rho, cp, bordes, conductividad,
    )
    pasos = max(1, math.ceil(t_final / (factor_seguridad * dt_critico)))
    dt = t_final / pasos

    T = simular_difusion(
        T_inicial, alpha, espaciamientos, dt, pasos,
        fuente, rho, cp, bordes, conductividad,
    )
    return T, dt
