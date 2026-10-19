import numpy as np
import pytest

from mcplan.hidratacion import HidratacionJuLei2019
from mcplan.solver1d import Borde, Capa, Material, crear_malla, simular

MODELO = HidratacionJuLei2019()


def test_curva_parte_en_cero_y_tiende_a_60():
    assert MODELO.theta(0.0) == 0.0
    assert MODELO.theta(1e4) == pytest.approx(60.0)


def test_curva_es_continua_en_21h():
    assert MODELO.theta(21.0 - 1e-9) == pytest.approx(MODELO.theta(21.0 + 1e-9), abs=0.05)


@pytest.mark.parametrize("te", [5.0, 15.0, 30.0, 100.0])
def test_derivada_coincide_con_diferencia_finita(te):
    h = 1e-4
    numerica = (MODELO.theta(te + h) - MODELO.theta(te - h)) / (2 * h)
    assert MODELO.dtheta_dte(te) == pytest.approx(numerica, rel=1e-5)


def test_arrhenius():
    assert MODELO.factor_arrhenius(20.0) == pytest.approx(1.0)
    assert MODELO.factor_arrhenius(40.0) > 1.0  # más caliente: hidrata más rápido
    assert MODELO.factor_arrhenius(5.0) < 1.0


def test_bloque_adiabatico_sigue_la_curva_del_ensayo():
    """Sin pérdidas, todo el calor queda en el hormigón: T(t) = T0 + theta(te(t)).
    Es exactamente el ensayo adiabático con el que se ajustó la curva."""
    mat = Material("hormigon", k=2.222, rho=2400.0, cp=970.0, es_hormigon=True, hidratacion=MODELO)
    malla = crear_malla([Capa(mat, 1.0, 20.0)], 0.1)
    adiabatico = Borde.neumann(0.0)
    res = simular(malla, adiabatico, adiabatico, t_final=7 * 86400.0, dt=60.0)
    esperado = 20.0 + MODELO.theta(res.te[:, 0])
    np.testing.assert_allclose(res.T[:, 0], esperado, atol=0.1)
    assert np.ptp(res.T[-1]) < 1e-9  # sin pérdidas, el bloque se calienta uniforme
