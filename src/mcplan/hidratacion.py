"""Calor de hidratación del hormigón como fuente volumétrica de calor.

Modelo de Ju y Lei (2019), ecuaciones (2)-(4):

    theta(te) = 0.02 * te^2.342                     si te <= 21 h
    theta(te) = 60 - 35 * exp(-0.03 * (te - 21))    si te >= 21 h

donde theta [°C] es el aumento adiabático de temperatura y te [h] la edad
equivalente (Arrhenius):

    dte/dt = exp( Eh/R * (1/T_ref - 1/T) )          (T en Kelvin)

En condiciones adiabáticas  rho*cp * dT/dt = Q, por lo que la fuente de calor es

    Q = rho * cp * dtheta/dt = rho * cp * dtheta/dte * dte/dt      [W/m3]

Notas:
- La curva viene de un ensayo adiabático de 7 días; después de 168 h es una
  extrapolación (tiende a 60 °C).
- Las dos ramas son continuas en te = 21 h (ambas valen ~25 °C), pero su
  derivada no lo es (2.79 vs 1.05 °C/h). Es una propiedad del ajuste del paper.
"""

from dataclasses import dataclass

import numpy as np

KELVIN = 273.15


@dataclass
class HidratacionJuLei2019:
    """Curva adiabática de Ju y Lei (2019) con edad equivalente."""

    Eh_R: float = 2400.0  # energía de activación / constante de gases [K] (paper)
    T_ref_C: float = 20.0  # temperatura de referencia [°C] (T0 = 293.15 K en el paper)
    te_quiebre_h: float = 21.0  # edad equivalente donde cambia la rama del ajuste [h]

    def theta(self, te_h):
        """Aumento adiabático de temperatura [°C] para una edad equivalente te [h]."""
        te = np.maximum(np.asarray(te_h, dtype=float), 0.0)
        temprana = 0.02 * te**2.342
        tardia = 60.0 - 35.0 * np.exp(-0.03 * (te - self.te_quiebre_h))
        return np.where(te <= self.te_quiebre_h, temprana, tardia)

    def dtheta_dte(self, te_h):
        """Derivada d(theta)/d(te) [°C/h]."""
        te = np.maximum(np.asarray(te_h, dtype=float), 0.0)
        temprana = 0.02 * 2.342 * te**1.342
        tardia = 35.0 * 0.03 * np.exp(-0.03 * (te - self.te_quiebre_h))
        return np.where(te <= self.te_quiebre_h, temprana, tardia)

    def factor_arrhenius(self, T_C):
        """dte/dt [-]: cuánto más rápido (o lento) envejece el hormigón a T que a T_ref."""
        T_K = np.asarray(T_C, dtype=float) + KELVIN
        return np.exp(self.Eh_R * (1.0 / (self.T_ref_C + KELVIN) - 1.0 / T_K))

    def potencia(self, te_h, T_C, rho, cp):
        """Calor generado por unidad de volumen Q [W/m3].

        rho*cp [J/(m3·K)] * dtheta/dte [K/h] * dte/dt [-] / 3600 [s/h]
        """
        return rho * cp * self.dtheta_dte(te_h) * self.factor_arrhenius(T_C) / 3600.0
