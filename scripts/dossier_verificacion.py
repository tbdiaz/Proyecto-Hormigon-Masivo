"""Figuras y métricas del dossier de verificación (E06): convergencia espacial y temporal.

Uso, desde la raíz del repositorio:

    python scripts/dossier_verificacion.py

Escribe docs/img/dossier_convergencia_espacial.png y docs/img/dossier_convergencia_temporal.png
e imprime las tablas que se copian en docs/dossier_verificacion.md.
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from mcplan.hidratacion import CalorPulsoCoseno  # noqa: E402
from mcplan.solver_nd import simular_difusion  # noqa: E402
from mcplan.verificacion import estudio_convergencia  # noqa: E402
from mcplan.verificacion_nd import convergencia_sinusoidal  # noqa: E402

# --- Convergencia espacial: 1D (serie de Fourier), 2D y 3D (sinusoidal de la guía) ---
filas_1d = estudio_convergencia([0.2, 0.1, 0.05, 0.025, 0.0125])
filas_2d = convergencia_sinusoidal((1.0, 1.0), (15, 30, 60), alpha=1e-4, A=50.0, T0=20.0, t_final=1200.0)
filas_3d = convergencia_sinusoidal((1.0, 1.0, 1.0), (8, 16, 32), alpha=1e-4, A=50.0, T0=20.0, t_final=600.0)

print("Convergencia espacial (dt automático, proporcional a h²)")
for nombre, h, e, o in (
    [("1D", f.dz, f.error_max, f.orden) for f in filas_1d]
    + [("2D", f.d[0], f.error_max, f.orden) for f in filas_2d]
    + [("3D", f.d[0], f.error_max, f.orden) for f in filas_3d]
):
    print(f"  {nombre}  h = {h:.4f} m  error = {e:.3e} °C  orden = {'' if o is None else f'{o:.2f}'}")

fig, ax = plt.subplots(figsize=(6.5, 4.5))
for nombre, h, e, marca in (
    ("1D (losa, serie de Fourier)", [f.dz for f in filas_1d], [f.error_max for f in filas_1d], "o"),
    ("2D sinusoidal", [f.d[0] for f in filas_2d], [f.error_max for f in filas_2d], "s"),
    ("3D sinusoidal", [f.d[0] for f in filas_3d], [f.error_max for f in filas_3d], "^"),
):
    ax.loglog(h, e, marca + "-", label=nombre)
h_ref = np.array([0.0125, 0.2])
ax.loglog(h_ref, 2.0 * h_ref**2, "k--", lw=1, label="pendiente 2")
ax.set_xlabel("tamaño de celda h [m]")
ax.set_ylabel("error máximo [°C]")
ax.set_title("Convergencia espacial: orden 2")
ax.grid(alpha=0.3, which="both")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig("docs/img/dossier_convergencia_espacial.png", dpi=150)
plt.close(fig)

# --- Convergencia temporal: bloque adiabático con el pulso de calor ---
# Sin flujo por los bordes el campo queda uniforme: no hay error de malla y sólo
# queda el de Euler explícito. Exacta: T = 20 + E(t) / (rho cp).
pulso = CalorPulsoCoseno()
rho, cp, k = 2400.0, 970.0, 2.222
t_final = 30 * 3600.0
exacta = 20.0 + pulso.energia(t_final / 3600.0) / (rho * cp)
print("\nConvergencia temporal (bloque adiabático con el pulso, t = 30 h)")
dts, errores = [600.0, 300.0, 150.0, 75.0, 37.5], []
for dt in dts:
    T = simular_difusion(np.full((3, 3, 3), 20.0), k / (rho * cp), (0.1,) * 3, dt, round(t_final / dt),
                         rho=rho, cp=cp, hidratacion=pulso, edad=np.zeros((3, 3, 3)))
    errores.append(abs(T[1, 1, 1] - exacta))
for i, (dt, e) in enumerate(zip(dts, errores)):
    orden = "" if i == 0 else f"{np.log2(errores[i - 1] / e):.2f}"
    print(f"  dt = {dt:6.1f} s  error = {e:.3e} °C  orden = {orden}")

fig, ax = plt.subplots(figsize=(6.5, 4.5))
ax.loglog(dts, errores, "o-", label="Euler explícito (pulso, t = 30 h)")
ax.loglog(dts, errores[0] * np.array(dts) / dts[0], "k--", lw=1, label="pendiente 1")
ax.set_xlabel("paso de tiempo Δt [s]")
ax.set_ylabel("error [°C]")
ax.set_title("Convergencia temporal: orden 1")
ax.grid(alpha=0.3, which="both")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig("docs/img/dossier_convergencia_temporal.png", dpi=150)
plt.close(fig)
