"""Calor de hidratación del hormigón como fuente volumétrica de calor.

Dos modelos con la misma interfaz (potencia y factor_arrhenius): la curva de
Ju y Lei (2019) y el pulso coseno de la clase (CalorPulsoCoseno, al final).

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


@dataclass
class CalorPulsoCoseno:
    """Calor de hidratación simplificado de la clase: un pulso coseno en el tiempo.

        Q(t) = Q_max/2 * (1 - cos(2 pi t / T))   si t < T
        Q(t) = 0                                 si t >= T

    Parte en 0, llega a Q_max en t = T/2 y vuelve a 0 en t = T. Sólo depende
    del tiempo (no de la temperatura): el problema sigue siendo lineal y tiene
    solución exacta. Energía total liberada: Q_max T / 2 [J/m3].

    Valores por defecto tomados de Ju y Lei (2019): el ritmo máximo de su curva
    es 2,79 °C/h a las 21 h, es decir Q_max = 2400*970*2.79/3600 ≈ 1800 W/m3 y
    T = 42 h; la elevación adiabática total queda Q_max T/(2 rho cp) ≈ 58 °C.

    Tiene la misma interfaz que HidratacionJuLei2019 (potencia y
    factor_arrhenius), con dte/dt = 1: la "edad" es el tiempo real en horas.
    """

    Q_max: float = 1800.0  # W/m3
    T_h: float = 42.0  # duración del pulso [h]

    def Q(self, t_h):
        """Calor generado [W/m3] en el instante t [h]."""
        t = np.asarray(t_h, dtype=float)
        pulso = 0.5 * self.Q_max * (1.0 - np.cos(2.0 * np.pi * t / self.T_h))
        return np.where((t >= 0.0) & (t < self.T_h), pulso, 0.0)

    def energia(self, t_h):
        """Calor acumulado desde t = 0 [J/m3]: integral de Q."""
        t = np.clip(np.asarray(t_h, dtype=float), 0.0, self.T_h)
        return 0.5 * self.Q_max * 3600.0 * (t - self.T_h / (2.0 * np.pi) * np.sin(2.0 * np.pi * t / self.T_h))

    def integral_amortiguada(self, a, t_s):
        """G(a, t) = int_0^t e^{-a (t - s)} Q(s) ds [J/m3], cerrada para el pulso coseno.

        Es la respuesta de un modo que decae con tasa a [1/s] a la fuente Q(t):
        dA/dt = -a A + Q(t)  ->  A(t) = G(a, t). Con a = 0 es la energía liberada.
        """
        a = np.asarray(a, dtype=float)
        T_s = self.T_h * 3600.0
        w = 2.0 * np.pi / T_s
        G = np.full(a.shape, float(self.energia(t_s / 3600.0)))  # a = 0: sin decaimiento
        m = a > 0
        b = a[m]
        if t_s <= T_s:
            G[m] = 0.5 * self.Q_max * ((1.0 - np.exp(-b * t_s)) / b
                                       - (b * np.cos(w * t_s) + w * np.sin(w * t_s) - b * np.exp(-b * t_s))
                                       / (b**2 + w**2))
        else:  # después del pulso, cada modo sólo decae desde su valor en t = T
            G[m] = (0.5 * self.Q_max * (1.0 - np.exp(-b * T_s)) * w**2 / (b * (b**2 + w**2))
                    * np.exp(-b * (t_s - T_s)))
        return G

    def factor_arrhenius(self, T_C):
        """dte/dt = 1: el pulso no depende de la temperatura."""
        return np.ones_like(np.asarray(T_C, dtype=float))

    def potencia(self, te_h, T_C, rho, cp):
        """Q [W/m3]; rho y cp no se usan (Q ya viene en W/m3), están por la interfaz común."""
        return self.Q(te_h) * np.ones_like(np.asarray(T_C, dtype=float))
