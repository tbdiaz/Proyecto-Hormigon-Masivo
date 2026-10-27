
"""Solver de difusion termica 2D y 3D con Euler explicito.

Incluye:
- Material homogeneo y malla uniforme.
- Fuente volumetrica uniforme (W/m3).
- Dirichlet: temperatura fija en la superficie.
- Neumann: flujo impuesto, positivo hacia el material (W/m2).
- Robin: conveccion con el ambiente.
- Bordes adiabaticos por defecto.
"""

import numpy as np


def paso_critico(alpha: float, espaciamientos: tuple[float, ...]) -> float:
    """Limite clasico de estabilidad para difusion explicita."""
    if alpha <= 0 or len(espaciamientos) not in (2, 3):
        raise ValueError("Alpha y la dimension deben ser validos")
    if any(d <= 0 for d in espaciamientos):
        raise ValueError("Los espaciamientos deben ser positivos")

    return 1.0 / (
        2.0 * alpha * sum(1.0 / d**2 for d in espaciamientos)
    )


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
) -> np.ndarray:
    """Simula conduccion termica 2D/3D.

    bordes permite especificar cada cara:
        "x_min", "x_max", "y_min", "y_max",
        "z_min", "z_max" (solo en 3D).

    Cada cara recibe una tupla:
        ("dirichlet", temperatura)
        ("neumann", flujo_entrante)
        ("robin", (h, temperatura_ambiente))

    Una cara no especificada es adiabatica.
    """

    T = np.array(T_inicial, dtype=float, copy=True)

    if T.ndim not in (2, 3):
        raise ValueError("T debe tener 2 o 3 dimensiones")
    if len(espaciamientos) != T.ndim:
        raise ValueError("Debe haber un espaciamiento por eje")
    if any(n < 2 for n in T.shape):
        raise ValueError("Cada eje necesita al menos dos celdas")
    if alpha <= 0 or any(d <= 0 for d in espaciamientos):
        raise ValueError("Alpha y espaciamientos deben ser positivos")
    if rho <= 0 or cp <= 0:
        raise ValueError("rho y cp deben ser positivos")
    if dt <= 0 or not isinstance(pasos, (int, np.integer)) or pasos < 0:
        raise ValueError("dt y pasos deben ser validos")

    if not all(
        np.isfinite(v)
        for v in (alpha, dt, fuente, rho, cp, *espaciamientos)
    ):
        raise ValueError("Los parametros deben ser finitos")
    if not np.all(np.isfinite(T)):
        raise ValueError("Las temperaturas deben ser finitas")

    bordes = {} if bordes is None else dict(bordes)

    nombres = ("x", "y", "z")[:T.ndim]
    caras_validas = {
        f"{nombre}_{lado}"
        for nombre in nombres
        for lado in ("min", "max")
    }

    for cara, condicion in bordes.items():
        if cara not in caras_validas:
            raise ValueError(f"Cara desconocida: {cara}")
        if not isinstance(condicion, (tuple, list)) or len(condicion) != 2:
            raise ValueError(f"Condicion incorrecta para {cara}")

        tipo, dato = condicion

        if tipo == "dirichlet" or tipo == "neumann":
            if not np.isscalar(dato) or not np.isfinite(dato):
                raise ValueError(f"Valor incorrecto en {cara}")

        elif tipo == "robin":
            if not isinstance(dato, (tuple, list)) or len(dato) != 2:
                raise ValueError(f"Robin necesita (h, T_amb) en {cara}")
            h, T_amb = dato
            if not np.isfinite(h) or not np.isfinite(T_amb) or h < 0:
                raise ValueError(f"Parametros Robin invalidos en {cara}")

        else:
            raise ValueError(f"Tipo de borde desconocido: {tipo}")

    # Conductividad consistente con alpha = k / (rho * cp).
    k = alpha * rho * cp

    # Coeficientes para comprobar estabilidad local.
    suma_coeficientes = np.zeros_like(T)

    for eje, d in enumerate(espaciamientos):
        coef = alpha / d**2

        indice_izq = [slice(None)] * T.ndim
        indice_der = [slice(None)] * T.ndim
        indice_izq[eje] = slice(None, -1)
        indice_der[eje] = slice(1, None)

        suma_coeficientes[tuple(indice_izq)] += coef
        suma_coeficientes[tuple(indice_der)] += coef

        for lado in ("min", "max"):
            cara = f"{nombres[eje]}_{lado}"
            tipo, dato = bordes.get(cara, ("neumann", 0.0))

            indice_borde = [slice(None)] * T.ndim
            indice_borde[eje] = 0 if lado == "min" else -1
            indice_borde = tuple(indice_borde)

            if tipo == "dirichlet":
                suma_coeficientes[indice_borde] += 2.0 * coef

            elif tipo == "robin":
                h, _ = dato
                if h > 0:
                    h_efectivo = 1.0 / (1.0 / h + d / (2.0 * k))
                    suma_coeficientes[indice_borde] += (
                        h_efectivo / (rho * cp * d)
                    )

    max_coef = float(np.max(suma_coeficientes))
    dt_critico = 1.0 / max_coef

    if dt > dt_critico * (1.0 + 1e-12):
        raise ValueError(
            f"dt={dt} supera el paso critico de "
            f"{dt_critico:.6g} s"
        )

    for _ in range(pasos):
        variacion = np.zeros_like(T)

        for eje, d in enumerate(espaciamientos):
            coef = alpha / d**2

            indice_izq = [slice(None)] * T.ndim
            indice_der = [slice(None)] * T.ndim
            indice_izq[eje] = slice(None, -1)
            indice_der[eje] = slice(1, None)

            indice_izq = tuple(indice_izq)
            indice_der = tuple(indice_der)

            # Intercambio entre celdas vecinas:
            # lo que pierde una celda lo gana la otra.
            diferencia = T[indice_der] - T[indice_izq]

            variacion[indice_izq] += coef * diferencia
            variacion[indice_der] -= coef * diferencia

            # Intercambio con el exterior.
            for lado in ("min", "max"):
                cara = f"{nombres[eje]}_{lado}"
                tipo, dato = bordes.get(cara, ("neumann", 0.0))

                indice_borde = [slice(None)] * T.ndim
                indice_borde[eje] = 0 if lado == "min" else -1
                indice_borde = tuple(indice_borde)

                if tipo == "dirichlet":
                    T_superficie = float(dato)
                    variacion[indice_borde] += (
                        2.0 * coef *
                        (T_superficie - T[indice_borde])
                    )

                elif tipo == "neumann":
                    flujo_entrante = float(dato)
                    variacion[indice_borde] += (
                        flujo_entrante / (rho * cp * d)
                    )

                elif tipo == "robin":
                    h, T_amb = dato
                    if h > 0:
                        h_efectivo = 1.0 / (
                            1.0 / h + d / (2.0 * k)
                        )
                        variacion[indice_borde] += (
                            h_efectivo *
                            (T_amb - T[indice_borde]) /
                            (rho * cp * d)
                        )

        # Fuente volumetrica en W/m3.
        variacion += fuente / (rho * cp)

        # Euler explicito.
        T = T + dt * variacion

    return T
