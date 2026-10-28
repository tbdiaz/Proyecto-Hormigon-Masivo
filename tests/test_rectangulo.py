"""Ejercicio de la pizarra: rectángulo 4 x 2 m a 20 °C con todo el borde a 10 °C (con y sin hidratación)."""

import numpy as np
import pytest

from mcplan.cli import main
from mcplan.hidratacion import CalorPulsoCoseno, HidratacionJuLei2019
from mcplan.rectangulo import (
    en_el_borde,
    simular_rectangulo,
    solucion_exacta,
    temperatura_en,
    tiempo_estabilizacion,
)
from mcplan.solver1d import Borde, Capa, Material, crear_malla, simular
from mcplan.solver_nd import simular_difusion, simular_serie

ALPHA_GUIA = 1e-4
HIDRATACION = HidratacionJuLei2019()
K, RHO, CP = 2.222, 2400.0, 970.0  # hormigón de Ju y Lei (2019)


def _rectangulo(**kw):
    """Caso rápido: alpha de la guía (k = alpha rho cp) y malla de 0.1 m."""
    return simular_rectangulo(k=ALPHA_GUIA * 2400 * 970, n_cuadros=60, **kw)


def test_serie_equivale_a_una_sola_simulacion():
    T0 = np.full((8, 4), 20.0)
    bordes = {c: ("dirichlet", 10.0) for c in ("x_min", "x_max", "y_min", "y_max")}
    s = simular_serie(T0, 1e-4, (0.1, 0.1), 300.0, 60.0, bordes=bordes)
    assert np.allclose(s.t, [0, 60, 120, 180, 240, 300])
    assert np.array_equal(s.T[0], T0)
    assert s.te is None
    pasos = round(300.0 / s.dt)
    assert np.array_equal(s.T[-1], simular_difusion(T0, 1e-4, (0.1, 0.1), s.dt, pasos, bordes=bordes))


def test_hidratacion_2d_igual_al_solver_1d():
    """Losa de 2 m entre dos caras a 10 °C: el 2D (caras laterales adiabáticas) debe dar lo mismo que el 1D."""
    malla = crear_malla([Capa(Material("h", K, RHO, CP, True, HIDRATACION), 2.0, 20.0)], 0.1)
    borde = Borde.dirichlet(10.0)
    res1d = simular(malla, borde, borde, t_final=3 * 86400.0, dt=150.0, guardar_cada=3 * 86400.0)

    edad = np.zeros((3, 20))
    T2d = simular_difusion(np.full((3, 20), 20.0), K / (RHO * CP), (0.1, 0.1), 150.0, round(3 * 86400 / 150),
                           rho=RHO, cp=CP, bordes={"y_min": ("dirichlet", 10.0), "y_max": ("dirichlet", 10.0)},
                           hidratacion=HIDRATACION, edad=edad)
    assert np.allclose(T2d, res1d.T[-1][None, :], atol=1e-9)
    assert np.allclose(edad, res1d.te[-1][None, :], atol=1e-9)


def test_bloque_adiabatico_2d_sigue_la_curva_del_ensayo():
    """Sin pérdidas, cada celda es el ensayo adiabático: T = T0 + theta(te), uniforme."""
    s = simular_serie(np.full((4, 3), 20.0), K / (RHO * CP), (0.1, 0.1), 7 * 86400.0, 3600.0, rho=RHO, cp=CP,
                      dt_max=60.0, hidratacion=HIDRATACION)
    np.testing.assert_allclose(s.T[:, 0, 0], 20.0 + HIDRATACION.theta(s.te[:, 0, 0]), atol=0.1)
    assert np.ptp(s.T[-1]) < 1e-9 and np.ptp(s.te[-1]) < 1e-9


def test_hidratacion_solo_en_el_hormigon():
    """Mitad hormigón, mitad sin hidratación: la edad sólo avanza donde hay hormigón."""
    es_hormigon = np.zeros((6, 3), dtype=bool)
    es_hormigon[:3] = True
    edad = np.zeros((6, 3))
    T = simular_difusion(np.full((6, 3), 20.0), K / (RHO * CP), (0.1, 0.1), 150.0, 200, rho=RHO, cp=CP,
                         hidratacion=HIDRATACION, edad=edad, es_hormigon=es_hormigon)
    assert np.all(edad[3:] == 0.0) and np.all(edad[:3] > 0.0)
    assert T[:3].mean() > T[3:].mean() > 20.0  # el calor generado se conduce hacia el resto


def test_rectangulo_con_hidratacion():
    res = simular_rectangulo(dx=0.25, t_final=5 * 86400.0, hidratacion=True, n_cuadros=20)
    assert res.hidratacion and res.dt <= 150.0 and res.T_puntos_exacta is None
    n, x, y = res.posicion_Tmax()
    assert (x, y) == pytest.approx((0.0, 0.0))  # el máximo está en el centro (simetría)
    assert 60.0 < res.T.max() < 20.0 + 60.0  # sube mucho, pero nunca más que el ensayo adiabático
    assert 24 * 3600 < res.t[n] < 4 * 86400  # el peak llega en los primeros días
    assert res.T_puntos[:, 1].max() < res.T_puntos[:, 0].max()  # cerca del borde se calienta menos
    assert np.all(res.dT_t >= 0) and res.dT_t.max() > 20.0
    T = res.T[n]
    assert np.allclose(T, T[::-1, :]) and np.allclose(T, T[:, ::-1])


def test_coincide_con_la_solucion_exacta():
    res = _rectangulo()
    assert np.abs(res.T_puntos - res.T_puntos_exacta).max() < 0.05
    assert res.error_campo[-1] < 1e-3
    assert res.dt <= res.dt_critico


def test_paso_critico_del_ejercicio():
    res = _rectangulo()
    assert res.dt_critico_interior == pytest.approx(0.1**2 / (4 * ALPHA_GUIA))  # 25 s
    assert res.dt_critico == pytest.approx(0.1**2 / (6 * ALPHA_GUIA))  # 16.7 s: la esquina manda


def test_fisica_del_enfriamiento():
    res = _rectangulo()
    centro, cerca_del_borde = res.T_puntos[:, 0], res.T_puntos[:, 1]
    assert np.all(np.diff(centro) <= 1e-12)  # sólo baja
    assert res.T.min() >= 10.0 - 1e-9 and res.T.max() <= 20.0 + 1e-9  # principio del máximo
    assert np.all(cerca_del_borde <= centro + 1e-12)  # el punto cerca del borde se enfría antes
    T = res.T[len(res.t) // 3]
    assert np.allclose(T, T[::-1, :]) and np.allclose(T, T[:, ::-1])  # simetría del rectángulo


def test_estabiliza_cuando_predice_el_primer_modo():
    res = _rectangulo()
    assert res.t_estabilizacion == pytest.approx(res.t_estabilizacion_estimado, rel=0.02)
    assert np.abs(res.T[-1] - 10.0).max() < res.tolerancia
    assert tiempo_estabilizacion(4, 2, ALPHA_GUIA, 20, 10, 0.1) == pytest.approx(
        np.log(16 / np.pi**2 * 10 / 0.1) / (ALPHA_GUIA * np.pi**2 * (1 / 16 + 1 / 4)))


def test_origen_en_el_centro_o_en_la_esquina():
    T = _rectangulo().T[10]
    assert temperatura_en(T, 4, 2, 10.0, (0, 0)) == pytest.approx(temperatura_en(T, 4, 2, 10.0, (2, 1), "esquina"))
    assert en_el_borde((1, 1), 4, 2)  # con origen al centro, (1, 1) está en el borde superior
    assert not en_el_borde((1, 1), 4, 2, "esquina")
    assert temperatura_en(T, 4, 2, 10.0, (1, 1)) == 10.0
    with pytest.raises(ValueError, match="fuera"):
        temperatura_en(T, 4, 2, 10.0, (3, 0))


def test_solucion_exacta_es_producto_de_dos_losas():
    x, y, t = 1.3, 0.4, 2000.0
    from mcplan.verificacion import solucion_analitica

    X = solucion_analitica(x, t, 4.0, ALPHA_GUIA, 1.0, 0.0)
    Y = solucion_analitica(y, t, 2.0, ALPHA_GUIA, 1.0, 0.0)
    assert solucion_exacta(x, y, t, 4.0, 2.0, ALPHA_GUIA, 20.0, 10.0)[0, 0] == pytest.approx(10.0 + 10.0 * X * Y)


def test_pulso_coseno_de_la_clase():
    p = CalorPulsoCoseno(Q_max=1800.0, T_h=42.0)
    assert p.Q(0.0) == 0.0 and p.Q(21.0) == pytest.approx(1800.0) and p.Q(42.0) == 0.0 and p.Q(60.0) == 0.0
    t = np.linspace(0.0, 42.0, 20001)
    Q = p.Q(t)
    integral = np.sum((Q[1:] + Q[:-1]) / 2 * np.diff(t)) * 3600.0  # trapecio [J/m3]
    assert p.energia(42.0) == pytest.approx(integral, rel=1e-6)
    assert p.energia(42.0) == pytest.approx(1800.0 * 42.0 * 3600.0 / 2)  # Qmax T / 2
    assert np.all(p.factor_arrhenius(np.array([5.0, 60.0])) == 1.0)  # no depende de T


def test_bloque_adiabatico_con_pulso_sube_Qmax_T_sobre_2_rho_cp():
    p = CalorPulsoCoseno()
    s = simular_serie(np.full((4, 3), 20.0), K / (RHO * CP), (0.1, 0.1), 3 * 86400.0, 3600.0, rho=RHO, cp=CP,
                      dt_max=150.0, hidratacion=p)
    assert np.ptp(s.T[-1]) < 1e-9
    assert s.T[-1, 0, 0] == pytest.approx(20.0 + p.Q_max * p.T_h * 3600.0 / (2 * RHO * CP), abs=1e-6)


def test_rectangulo_con_pulso_coincide_con_la_exacta_y_converge():
    errores = []
    for dx in (0.2, 0.1):
        res = simular_rectangulo(dx=dx, t_final=3 * 86400.0, hidratacion=CalorPulsoCoseno(), n_cuadros=20)
        assert res.T_puntos_exacta is not None
        errores.append(np.abs(res.T_puntos - res.T_puntos_exacta).max())
    assert errores[1] < 0.2  # °C, con picos de ~75 °C
    assert errores[0] / errores[1] > 3  # orden ~2: al dividir dx por 2, el error baja ~4 veces


def test_dx_debe_dividir_los_lados():
    with pytest.raises(ValueError, match="dividir"):
        simular_rectangulo(dx=0.3)


def test_comando_genera_graficos_y_animacion(tmp_path, capsys):
    codigo = main(["rectangulo", "--alpha", "1e-4", "--dx", "0.25", "--cuadros", "8",
                   "--punto", "0", "0", "--punto", "1", "1", "--salida", str(tmp_path)])
    assert codigo == 0
    assert "está en el borde" in capsys.readouterr().out
    for nombre in ("puntos.png", "campos.png", "evolucion.gif", "serie_puntos.csv", "resumen.json"):
        assert (tmp_path / nombre).stat().st_size > 0


@pytest.mark.parametrize("modelo", [[], ["pulso"]])
def test_comando_con_hidratacion(tmp_path, capsys, modelo):
    codigo = main(["rectangulo", "--hidratacion", *modelo, "--dx", "0.25", "--t", "172800", "--cuadros", "6",
                   "--sin-animacion", "--salida", str(tmp_path)])
    assert codigo == 0
    salida = capsys.readouterr().out
    assert "tope de precisión" in salida and "Tmax" in salida and "ΔTmax" in salida
    assert ("pulso de la clase" in salida) == bool(modelo)
    assert ("vs. exacta" in salida) == bool(modelo)  # el pulso tiene solución exacta; Ju y Lei no
    assert (tmp_path / "puntos.png").stat().st_size > 0
