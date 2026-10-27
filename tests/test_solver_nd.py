
import numpy as np

from mcplan.solver_nd import simular_difusion


def test_temperatura_uniforme_2d():
    T_inicial = np.full((10, 10), 20.0)

    T_final = simular_difusion(
        T_inicial=T_inicial,
        alpha=1e-6,
        espaciamientos=(0.1, 0.1),
        dt=100.0,
        pasos=50,
    )

    assert np.allclose(T_final, 20.0)


def test_temperatura_uniforme_3d():
    T_inicial = np.full((6, 6, 6), 20.0)

    T_final = simular_difusion(
        T_inicial=T_inicial,
        alpha=1e-6,
        espaciamientos=(0.1, 0.1, 0.1),
        dt=100.0,
        pasos=50,
    )

    assert np.allclose(T_final, 20.0)


def test_conservacion_energia_2d():
    T_inicial = np.full((10, 10), 20.0)
    T_inicial[4:6, 4:6] = 60.0

    T_final = simular_difusion(
        T_inicial=T_inicial,
        alpha=1e-6,
        espaciamientos=(0.1, 0.1),
        dt=100.0,
        pasos=100,
    )

    assert np.isclose(
        np.sum(T_inicial),
        np.sum(T_final),
        rtol=1e-10,
    )



def test_fuente_uniforme_adiabatica_3d():
    T_inicial = np.full((6, 6, 6), 20.0)

    rho = 2400.0
    cp = 900.0
    fuente = 1000.0  # W/m³
    dt = 100.0
    pasos = 50

    T_final = simular_difusion(
        T_inicial=T_inicial,
        alpha=1e-6,
        espaciamientos=(0.1, 0.1, 0.1),
        dt=dt,
        pasos=pasos,
        fuente=fuente,
        rho=rho,
        cp=cp,
    )

    T_esperada = 20.0 + fuente * (dt * pasos) / (rho * cp)

    assert np.allclose(T_final, T_esperada)


def test_dirichlet_2d():
    T_inicial = np.full((10, 10), 20.0)

    T_final = simular_difusion(
        T_inicial=T_inicial,
        alpha=1e-6,
        espaciamientos=(0.1, 0.1),
        dt=100.0,
        pasos=100,
        bordes={"x_min": ("dirichlet", 60.0)},
    )

    # La superficie caliente debe calentar las celdas cercanas.
    assert np.mean(T_final[0, :]) > 20.0
    assert np.mean(T_final[0, :]) > np.mean(T_final[-1, :])


def test_neumann_2d():
    T_inicial = np.full((10, 10), 20.0)

    T_final = simular_difusion(
        T_inicial=T_inicial,
        alpha=1e-6,
        espaciamientos=(0.1, 0.1),
        dt=100.0,
        pasos=100,
        bordes={"x_min": ("neumann", 100.0)},
    )

    # Un flujo entrante positivo debe aumentar la energia.
    energia_inicial = np.sum(T_inicial)
    energia_final = np.sum(T_final)

    aumento_esperado = (
        100.0 * 100.0 * 100 / (2400.0 * 900.0 * 0.1)
        * T_inicial.shape[1]
    )

    assert np.isclose(
        energia_final - energia_inicial,
        aumento_esperado,
    )


def test_robin_3d():
    T_inicial = np.full((6, 6, 6), 60.0)

    T_final = simular_difusion(
        T_inicial=T_inicial,
        alpha=1e-6,
        espaciamientos=(0.1, 0.1, 0.1),
        dt=100.0,
        pasos=100,
        bordes={"z_max": ("robin", (10.0, 20.0))},
    )

    # Un ambiente mas frio debe enfriar el hormigon.
    assert np.mean(T_final) < 60.0
    assert np.mean(T_final[:, :, -1]) < np.mean(T_final[:, :, 0])


def test_solucion_exacta_sinusoidal_2d():
    # Parametros del problema
    nx, ny = 20, 20
    Lx, Ly = 1.0, 1.0
    dx, dy = Lx / nx, Ly / ny

    alpha = 1e-6
    T0 = 20.0
    A = 30.0
    dt = 100.0
    pasos = 1000

    # Coordenadas de los centros de las celdas
    x = (np.arange(nx) + 0.5) * dx
    y = (np.arange(ny) + 0.5) * dy
    X, Y = np.meshgrid(x, y, indexing="ij")

    # Condicion inicial sinusoidal
    forma = np.sin(np.pi * X / Lx) * np.sin(np.pi * Y / Ly)
    T_inicial = T0 + A * forma

    # Solucion numerica
    T_numerica = simular_difusion(
        T_inicial=T_inicial,
        alpha=alpha,
        espaciamientos=(dx, dy),
        dt=dt,
        pasos=pasos,
        bordes={
            "x_min": ("dirichlet", T0),
            "x_max": ("dirichlet", T0),
            "y_min": ("dirichlet", T0),
            "y_max": ("dirichlet", T0),
        },
    )

    # Solucion analitica del enunciado
    tiempo = dt * pasos
    lambd = (np.pi / Lx)**2 + (np.pi / Ly)**2

    T_exacta = T0 + A * forma * np.exp(-alpha * lambd * tiempo)

    # Error maximo permitido: 0.2 grados Celsius
    error_maximo = np.max(np.abs(T_numerica - T_exacta))

    assert error_maximo < 0.2


def test_solucion_exacta_sinusoidal_3d():
    # Parametros del problema
    nx, ny, nz = 12, 12, 12
    Lx, Ly, Lz = 1.0, 1.0, 1.0

    dx = Lx / nx
    dy = Ly / ny
    dz = Lz / nz

    alpha = 1e-6
    T0 = 20.0
    A = 30.0
    dt = 100.0
    pasos = 1000

    # Coordenadas de los centros de las celdas
    x = (np.arange(nx) + 0.5) * dx
    y = (np.arange(ny) + 0.5) * dy
    z = (np.arange(nz) + 0.5) * dz

    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")

    # Condicion inicial sinusoidal
    forma = (
        np.sin(np.pi * X / Lx)
        * np.sin(np.pi * Y / Ly)
        * np.sin(np.pi * Z / Lz)
    )

    T_inicial = T0 + A * forma

    # Solucion numerica
    T_numerica = simular_difusion(
        T_inicial=T_inicial,
        alpha=alpha,
        espaciamientos=(dx, dy, dz),
        dt=dt,
        pasos=pasos,
        bordes={
            "x_min": ("dirichlet", T0),
            "x_max": ("dirichlet", T0),
            "y_min": ("dirichlet", T0),
            "y_max": ("dirichlet", T0),
            "z_min": ("dirichlet", T0),
            "z_max": ("dirichlet", T0),
        },
    )

    # Solucion analitica oficial del E05
    tiempo = dt * pasos

    lambd = (
        (np.pi / Lx)**2
        + (np.pi / Ly)**2
        + (np.pi / Lz)**2
    )

    T_exacta = (
        T0 + A * forma * np.exp(-alpha * lambd * tiempo)
    )

    # Error maximo permitido: 0.3 grados Celsius
    error_maximo = np.max(
        np.abs(T_numerica - T_exacta)
    )

    assert error_maximo < 0.3
