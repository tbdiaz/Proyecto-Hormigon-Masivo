"""Interfaz de línea de comandos.

    mcplan simular1d config/caso_base_1d.yaml [--dz 0.1] [--dt 600] [--forzar]
    mcplan malla     config/caso_base_1d.yaml --dz 0.2 0.1 0.05 0.025
    mcplan verificar [--dz 0.2 0.1 0.05 0.025] [--dt 600] [--forzar]
    mcplan sinusoidal --dim 2|3 [--L 1 2] [--alpha 1e-4] [--A 50] [--T0 20] [--t 1200] [--n 15 30] [--dt] [--forzar]
    mcplan fuente     --dim 2|3 [--q0 10000] [--rho 2400] [--cp 900] [--T0 20] [--t 28800]
    mcplan energia    --dim 2|3 [--dos-materiales] [--t 259200]
    mcplan rectangulo [--L 4 2] [--dx 0.1] [--alpha 1e-4] [--Ti 20] [--Tb 10] [--punto 0 0 --punto 1 0.5]
                      [--hidratacion [ju-lei|pulso]] [--Qmax 1800] [--T-calor 42]
    mcplan bloque     [--L 3 3 0.6] [--dx 0.05] [--casos dirichlet neumann robin] [--T-dirichlet 20 30] [--beta 12.5]
    mcplan bloque-3d  [--res 1280 720] [--muestras 32]   (animación 3D con Blender: MP4 y GIF)
    mcplan ifc-propiedades geometry/foundation.ifc
    mcplan ifc-info  geometry/foundation.ifc [--guardar antes.json] [--comparar antes.json]
    mcplan ifc-crear prueba.ifc            (misma fundación, generada sin Blender)
"""

import argparse
import json
import sys
import time

from . import __version__


def _simular_caso(caso, dz=None, dt=None, forzar=False):
    from .metricas import calcular_metricas
    from .solver1d import crear_malla, simular

    malla = crear_malla(caso.capas, dz or caso.dz)
    res = simular(
        malla, caso.borde_inf, caso.borde_sup,
        t_final=caso.duracion_h * 3600.0,
        dt=dt if dt is not None else caso.dt_s,
        guardar_cada=caso.guardar_cada_h * 3600.0,
        factor_seguridad=caso.factor_seguridad,
        dt_max=caso.dt_max_s,
        forzar=forzar,
    )
    met = None if res.diverge else calcular_metricas(
        res, caso.profundidad_superficie, caso.limite_Tmax, caso.limite_dT)
    return res, met


def cmd_simular1d(args):
    from .caso import cargar_caso

    caso = cargar_caso(args.caso)
    t0 = time.perf_counter()
    res, met = _simular_caso(caso, args.dz, args.dt, args.forzar)
    duracion = time.perf_counter() - t0

    print(f"Caso: {caso.nombre}")
    print("Capas (de abajo hacia arriba): " + ", ".join(
        f"{c.nombre or c.material.nombre} {c.espesor:.2f} m @ {c.T_inicial:.1f} °C" for c in caso.capas))
    print(f"Celdas: {res.malla.n}   dz hormigón ≈ {res.malla.dz[res.malla.es_hormigon].mean():.4f} m")
    print(f"dt = {res.dt:.1f} s   (dt crítico = {res.dt_critico:.1f} s, "
          f"Fourier_max = {0.5 * res.dt / res.dt_critico:.3f})   pasos = {res.n_pasos}   "
          f"tiempo de cómputo = {duracion:.2f} s")

    if res.diverge:
        print(f"\n¡LA SOLUCIÓN DIVERGIÓ en t = {res.t_h[-1]:.1f} h! dt > dt crítico -> esquema explícito inestable.")
        if not args.sin_graficos:
            from .graficos import graficar_perfiles

            caso.salida.mkdir(parents=True, exist_ok=True)
            graficar_perfiles(res, caso.salida / "perfiles_divergencia.png", dias=(0,))
            print(f"Perfil al momento de divergir: {caso.salida / 'perfiles_divergencia.png'}")
        return 2

    print()
    print(met.resumen())

    caso.salida.mkdir(parents=True, exist_ok=True)
    resumen = {
        "caso": caso.nombre, "archivo_caso": str(args.caso),
        "n_celdas": res.malla.n, "dt_s": res.dt, "dt_critico_s": res.dt_critico, "n_pasos": res.n_pasos,
        "Tmax_C": met.Tmax, "t_Tmax_h": met.t_Tmax_h, "altura_Tmax_sobre_base_m": met.z_Tmax,
        "dTmax_C": met.dTmax, "t_dTmax_h": met.t_dTmax_h,
        "cumple_Tmax": met.cumple_Tmax, "cumple_dT": met.cumple_dT,
    }
    (caso.salida / "resumen.json").write_text(json.dumps(resumen, indent=2), encoding="utf-8")
    if not args.sin_graficos:
        from .graficos import graficar_todo

        rutas = graficar_todo(res, met, caso.salida)
        print("\nGráficos: " + ", ".join(str(r) for r in rutas))
    print(f"Resumen:  {caso.salida / 'resumen.json'}")
    return 0


def cmd_malla(args):
    """Sensibilidad de Tmax y dTmax al tamaño de celda en el caso real."""
    from .caso import cargar_caso

    caso = cargar_caso(args.caso)
    print(f"{'dz [m]':>8} {'celdas':>7} {'dt [s]':>8} {'pasos':>7} {'Tmax [°C]':>10} {'dTmax [°C]':>11} {'cómputo [s]':>12}")
    for dz in sorted(args.dz, reverse=True):
        t0 = time.perf_counter()
        res, met = _simular_caso(caso, dz=dz)
        print(f"{dz:8.4f} {res.malla.n:7d} {res.dt:8.1f} {res.n_pasos:7d} "
              f"{met.Tmax:10.3f} {met.dTmax:11.3f} {time.perf_counter() - t0:12.2f}")
    return 0


def cmd_verificar(args):
    from .verificacion import estudio_convergencia

    print("Verificación: losa de 1 m, T0 = 50 °C, caras a 20 °C, Fo = 0.05 (≈ 14.6 h)")
    print(f"{'dz [m]':>8} {'dt [s]':>9} {'dt crít [s]':>11} {'pasos':>7} {'error máx [°C]':>15} {'orden':>6}")
    for f in estudio_convergencia(args.dz, dt=args.dt, forzar=args.forzar):
        orden = "" if f.orden is None else f"{f.orden:6.2f}"
        print(f"{f.dz:8.4f} {f.dt:9.1f} {f.dt_critico:11.1f} {f.n_pasos:7d} {f.error_max:15.3e} {orden:>6}")
    return 0


def cmd_sinusoidal(args):
    from math import exp, pi

    from .verificacion_nd import convergencia_sinusoidal

    L = args.L * args.dim if len(args.L) == 1 else args.L
    if len(L) != args.dim:
        raise ValueError(f"--L necesita 1 o {args.dim} valores")
    ns = args.n or ([15, 30, 60] if args.dim == 2 else [8, 16, 32])
    t = args.t or (1200.0 if args.dim == 2 else 600.0)  # tiempos de la guía
    lam = sum((pi / Li) ** 2 for Li in L)
    print(f"Caso sinusoidal {args.dim}D: L = {' x '.join(f'{Li:g}' for Li in L)} m, alpha = {args.alpha:g} m²/s, "
          f"A = {args.A:g} °C, T0 = {args.T0:g} °C, t = {t:g} s")
    print(f"Contorno Dirichlet a T0. Exacta: T0 + A·Πsin(πx_i/L_i)·exp(−αλt), λ = {lam:.4g} 1/m², "
          f"exp(−αλt) = {exp(-args.alpha * lam * t):.4f}")
    print(f"{'n/eje':>6} {'d [m]':>8} {'dt [s]':>10} {'dt crít [s]':>11} {'pasos':>7} {'error máx [°C]':>15} {'orden':>6}")
    filas = convergencia_sinusoidal(L, ns, alpha=args.alpha, A=args.A, T0=args.T0, t_final=t,
                                    dt=args.dt, forzar=args.forzar)
    for f in filas:
        error = "DIVERGE" if f.error_max == float("inf") else f"{f.error_max:.3e}"
        orden = "" if f.orden is None else f"{f.orden:6.2f}"
        print(f"{f.n:6d} {f.d[0]:8.4f} {f.dt:10.4g} {f.dt_critico:11.4g} {f.pasos:7d} {error:>15} {orden:>6}")
    return 2 if any(f.error_max == float("inf") for f in filas) else 0


def cmd_fuente(args):
    from .verificacion_nd import caso_fuente_uniforme

    r = caso_fuente_uniforme(args.dim, args.q0, args.rho, args.cp, args.T0, args.t, n=args.n)
    print(f"Fuente uniforme {args.dim}D en dominio adiabático: q0 = {args.q0:g} W/m³, rho = {args.rho:g} kg/m³, "
          f"cp = {args.cp:g} J/(kg·K), T0 = {args.T0:g} °C, t = {args.t:g} s ({args.t / 3600:.2f} h)")
    print(f"Exacta: T0 + q0·t/(rho·cp) = {r.T_exacta:.6f} °C   (dT/dt = {args.q0 / (args.rho * args.cp) * 3600:.4f} °C/h)")
    print(f"Numérica: T media = {r.T_media:.6f} °C   error máx = {r.error_max:.2e} °C   "
          f"Tmax − Tmin = {r.dispersion:.2e} °C")
    print(f"Energía: ganada = {r.energia_ganada:.6e} J   aportada q0·V·t = {r.energia_aportada:.6e} J   "
          f"error relativo = {r.error_balance:.2e}")
    print(f"dt = {r.dt:.4g} s   pasos = {r.pasos}")
    return 0


def cmd_energia(args):
    from .verificacion_nd import caso_energia_adiabatica

    r = caso_energia_adiabatica(args.dim, n=args.n, t_final=args.t, dos_materiales=args.dos_materiales)
    materiales = "hormigón (x < 0,5 m) y suelo (x > 0,5 m)" if args.dos_materiales else "hormigón"
    print(f"Dominio adiabático {args.dim}D sin fuente, 1 m de lado, {args.n} celdas por eje, {materiales}, "
          f"t = {args.t:g} s ({args.t / 3600:.1f} h)")
    print(f"Tmax: {r.Tmax_inicial:.2f} °C -> {r.Tmax_final:.2f} °C (el calor se redistribuye)")
    print(f"Energía E = Σ ρ·cp·V·T: inicial = {r.E_inicial:.10e} J   final = {r.E_final:.10e} J   "
          f"error relativo = {r.error_relativo:.2e}")
    print(f"dt = {r.dt:.4g} s   pasos = {r.pasos}")
    return 0


def cmd_rectangulo(args):
    from math import pi
    from pathlib import Path

    import numpy as np

    from .hidratacion import CalorPulsoCoseno, HidratacionJuLei2019
    from .rectangulo import en_el_borde, simular_rectangulo, unidad_tiempo

    Lx, Ly = args.L
    k = args.alpha * args.rho * args.cp if args.alpha else args.k
    puntos = args.punto or [(0.0, 0.0), (1.0, 0.5)]
    modelo = {None: None, "ju-lei": HidratacionJuLei2019(),
              "pulso": CalorPulsoCoseno(args.Qmax, args.T_calor)}[args.hidratacion]
    t0 = time.perf_counter()
    res = simular_rectangulo(Lx, Ly, args.dx, k, args.rho, args.cp, args.Ti, args.Tb, puntos, args.origen,
                             args.t, args.tolerancia, args.cuadros, hidratacion=modelo if modelo else False,
                             dt_max=args.dt_max)
    computo = time.perf_counter() - t0
    f, unidad = unidad_tiempo(res.t[-1])

    print(f"Rectángulo {Lx:g} x {Ly:g} m, Ti = {args.Ti:g} °C, borde Dirichlet a {args.Tb:g} °C, "
          f"origen en {'el centro' if args.origen == 'centro' else 'la esquina'}")
    print(f"Material: k = {k:.4g} W/(m·K), rho = {args.rho:g} kg/m³, cp = {args.cp:g} J/(kg·K) "
          f"-> alpha = {res.alpha:.4g} m²/s")
    if isinstance(modelo, CalorPulsoCoseno):
        E = modelo.energia(modelo.T_h)
        print(f"Calor de hidratación (pulso de la clase): Q(t) = Qmax/2·(1 − cos(2πt/T)) si t < T, 0 después; "
              f"Qmax = {modelo.Q_max:g} W/m³, T = {modelo.T_h:g} h")
        print(f"  máximo Qmax en t = T/2 = {modelo.T_h / 2:g} h; energía total Qmax·T/2 = {E / 1e6:.1f} MJ/m³ "
              f"-> elevación adiabática Qmax·T/(2·ρ·cp) = {E / (args.rho * args.cp):.1f} °C")
    elif res.hidratacion:
        print("Calor de hidratación: Ju y Lei (2019), θ(te) con edad equivalente (Eh/R = 2400 K, Tref = 20 °C), "
              "te = 0 en t = 0 en todas las celdas")
    print(f"Malla: dx = dy = {args.dx:g} m, {res.T.shape[1]} x {res.T.shape[2]} celdas")
    print("\nPaso crítico (Euler explícito, coeficiente de T_i ≥ 0):")
    print(f"  celda interior   dx²/(4α) = {res.dt_critico_interior:10.1f} s")
    print(f"  junto al borde   dx²/(6α) = {res.dt_critico:10.1f} s   <- manda (la cara Dirichlet está a media celda)")
    if res.dt_max is not None:
        print(f"  tope de precisión dt_max  = {res.dt_max:10.1f} s   <- con hidratación manda la precisión, "
              "no la estabilidad (la fuente es explícita)")
    print(f"  dt usado                  = {res.dt:10.1f} s   ({res.dt / res.dt_critico:.0%} del crítico), "
          f"{round(res.t[-1] / res.dt)} pasos")
    print(f"\nEstabilización (|T − {args.Tb:g}| < {res.tolerancia:g} °C en todo el rectángulo):")
    print(f"  modo más lento: λ1 = (π/Lx)² + (π/Ly)² = {(pi / Lx) ** 2 + (pi / Ly) ** 2:.4f} 1/m², "
          f"τ = 1/(α·λ1) = {res.tau / f:.2f} {unidad}")
    print(f"  estimado ln((16/π²)·ΔT/tol)·τ = {res.t_estabilizacion_estimado / f:.2f} {unidad}"
          + ("   (sin hidratación)" if res.hidratacion else ""))
    if res.t_estabilizacion is None:
        print(f"  simulado: no alcanza la tolerancia en t = {res.t[-1] / f:.2f} {unidad} (aumente --t)")
    else:
        print(f"  simulado: {res.t_estabilizacion / f:.2f} {unidad}   (se simuló hasta {res.t[-1] / f:.2f} {unidad})")

    if res.hidratacion:
        n, x, y = res.posicion_Tmax()
        j = int(np.argmax(res.dT_t))
        print(f"\nCriterios del curso (ΔT con Tmin a {res.profundidad * 1000:g} mm del borde):")
        print(f"  Tmax  = {res.T.max():.2f} °C a las {res.t[n] / 3600:.0f} h en ({x:g}, {y:g})   "
              f"límite 70 °C: {'CUMPLE' if res.T.max() <= 70 else 'NO CUMPLE'}")
        print(f"  ΔTmax = {res.dT_t[j]:.2f} °C a las {res.t[j] / 3600:.0f} h   "
              f"límite 20 °C: {'CUMPLE' if res.dT_t[j] <= 20 else 'NO CUMPLE'}")
        print("  (la cara a 10 °C fijos es una condición extrema: en la fundación real el borde es Robin)")

    hay_exacta = res.T_puntos_exacta is not None
    print("\nPuntos" + (" (numérica vs. exacta)" if hay_exacta else "") + ":")
    for j, p in enumerate(res.puntos):
        nota = "   <- está en el borde: T = Tb siempre" if en_el_borde(p, Lx, Ly, args.origen) else ""
        error = (f", error máx vs. exacta = {abs(res.T_puntos[:, j] - res.T_puntos_exacta[:, j]).max():.3f} °C"
                 if hay_exacta else "")
        if res.hidratacion:
            i = int(np.argmax(res.T_puntos[:, j]))
            print(f"  ({p[0]:g}, {p[1]:g}): máx {res.T_puntos[i, j]:.2f} °C a las {res.t[i] / 3600:.0f} h, "
                  f"T final = {res.T_puntos[-1, j]:.3f} °C{error}{nota}")
            continue
        llego = np.abs(res.T_puntos[:, j] - args.Tb) <= 0.5 * abs(args.Ti - args.Tb)
        mitad = f"{res.t[np.argmax(llego)] / f:.2f} {unidad}" if llego.any() else "—"
        print(f"  ({p[0]:g}, {p[1]:g}): T final = {res.T_puntos[-1, j]:.3f} °C, recorre la mitad del salto "
              f"a los {mitad}{error}{nota}")
    if hay_exacta:
        print(f"Error máx en todo el campo (t > 0): {res.error_campo[1:].max():.3f} °C "
              f"(en t = {res.t[1 + np.argmax(res.error_campo[1:])] / f:.2f} {unidad}, por el salto inicial en el borde); "
              f"al final: {res.error_campo[-1]:.2e} °C")

    sufijo = {None: "", "ju-lei": "_hidratacion", "pulso": "_pulso"}[args.hidratacion]
    salida = Path(args.salida or f"resultados/rectangulo_2d{sufijo}")
    salida.mkdir(parents=True, exist_ok=True)
    nombres, columnas = ["t_s"], [res.t]
    for j, p in enumerate(res.puntos):
        nombres.append(f"T({p[0]:g};{p[1]:g})_num")
        columnas.append(res.T_puntos[:, j])
        if res.T_puntos_exacta is not None:
            nombres.append(f"T({p[0]:g};{p[1]:g})_exacta")
            columnas.append(res.T_puntos_exacta[:, j])
    if res.hidratacion:
        nombres += ["Tmax_C", f"Tmin_{res.profundidad * 1000:g}mm_C", "dT_C"]
        columnas += [res.Tmax_t, res.Tmin_superficie, res.dT_t]
    np.savetxt(salida / "serie_puntos.csv", np.column_stack(columnas), delimiter=",", header=",".join(nombres),
               comments="", fmt="%.6f")
    resumen = {
        "Lx_m": Lx, "Ly_m": Ly, "dx_m": args.dx, "k": k, "rho": args.rho, "cp": args.cp, "alpha_m2_s": res.alpha,
        "Ti_C": args.Ti, "Tb_C": args.Tb, "origen": args.origen, "puntos": res.puntos,
        "hidratacion": args.hidratacion or "ninguna", "dt_s": res.dt, "dt_max_s": res.dt_max,
        "dt_critico_s": res.dt_critico, "dt_critico_interior_s": res.dt_critico_interior,
        "tau_s": res.tau, "t_estabilizacion_estimado_sin_hidratacion_s": res.t_estabilizacion_estimado,
        "t_estabilizacion_simulado_s": res.t_estabilizacion, "t_final_s": float(res.t[-1]),
    }
    if isinstance(modelo, CalorPulsoCoseno):
        resumen.update({"Qmax_W_m3": modelo.Q_max, "T_calor_h": modelo.T_h,
                        "elevacion_adiabatica_C": float(modelo.energia(modelo.T_h) / (args.rho * args.cp))})
    if res.hidratacion:
        n, x, y = res.posicion_Tmax()
        resumen.update({"Tmax_C": float(res.T.max()), "t_Tmax_h": float(res.t[n] / 3600), "posicion_Tmax_m": [x, y],
                        "dTmax_C": float(res.dT_t.max()), "t_dTmax_h": float(res.t[np.argmax(res.dT_t)] / 3600)})
    if hay_exacta:
        resumen["error_max_puntos_C"] = float(abs(res.T_puntos - res.T_puntos_exacta).max())
    (salida / "resumen.json").write_text(json.dumps(resumen, indent=2), encoding="utf-8")

    from .graficos import animar_rectangulo, graficar_campos_rectangulo, graficar_puntos_rectangulo

    rutas = [salida / "puntos.png", salida / "campos.png"]
    graficar_puntos_rectangulo(res, rutas[0])
    graficar_campos_rectangulo(res, rutas[1])
    if not args.sin_animacion:
        rutas.append(salida / "evolucion.gif")
        animar_rectangulo(res, rutas[-1], fps=args.fps)
    print(f"\nArchivos: {', '.join(str(r) for r in rutas)}, {salida / 'serie_puntos.csv'}, {salida / 'resumen.json'}")
    print(f"Cómputo de la simulación: {computo:.2f} s")
    return 0


def cmd_bloque(args):
    from pathlib import Path

    import numpy as np

    from .bloque import CasoBorde, H_HORMIGON_AIRE, simular_bloque, t_final_comun
    from .hidratacion import CalorPulsoCoseno

    L = tuple(args.L)
    k = args.k
    h = args.beta * k if args.beta is not None else H_HORMIGON_AIRE
    casos = []
    if "dirichlet" in args.casos:
        casos += [CasoBorde("dirichlet", T) for T in args.T_dirichlet]
    if "neumann" in args.casos:
        casos.append(CasoBorde("neumann"))
    if "robin" in args.casos:
        casos.append(CasoBorde("robin", args.T_aire, h))
    pulso = CalorPulsoCoseno(args.Qmax, args.T_calor)
    puntos = [tuple(p) for p in args.punto] if args.punto else [(0.0, 0.0, 0.275), (1.2, 1.2, 0.2)]
    alpha = k / (args.rho * args.cp)
    n = tuple(round(Li / args.dx) for Li in L)
    t_final = args.t or t_final_comun(casos, L, k, args.rho, args.cp, args.Ti, pulso, args.tolerancia)

    E = float(pulso.energia(pulso.T_h))
    print(f"Bloque {L[0]:g} x {L[1]:g} x {L[2]:g} m, Ti = {args.Ti:g} °C, origen en "
          f"{'el centro (x, y en ±' + f'{L[0] / 2:g}, z en ±{L[2] / 2:g})' if args.origen == 'centro' else 'una esquina de la base'}")
    print(f"Material: k = {k:g} W/(m·K), rho = {args.rho:g} kg/m³, cp = {args.cp:g} J/(kg·K) -> alpha = {alpha:.4g} m²/s")
    print(f"Calor (pulso de la clase): Q(t) = Qmax/2·(1 − cos(2πt/T)) si t < T; Qmax = {pulso.Q_max:g} W/m³, "
          f"T = {pulso.T_h:g} h -> energía {E / 1e6:.1f} MJ/m³, elevación adiabática {E / (args.rho * args.cp):.2f} °C")
    print(f"Malla: dx = dy = dz = {args.dx:g} m, {n[0]} x {n[1]} x {n[2]} = {np.prod(n)} celdas")
    print("Robin: −k ∂T/∂n = h (T − T_aire), o sea ∂T/∂n = −β (T − T_aire) con β = h/k "
          f"(h = {h:g} W/(m²K){' de Ju y Lei' if args.beta is None else ''}, β = {h / k:.3f} 1/m)")
    print(f"Tiempo simulado: {t_final / 3600:.1f} h "
          f"({'1,2 × la estabilización exacta más lenta' if args.t is None else 'indicado con --t'})")

    print("\nPaso crítico (Euler explícito: el coeficiente de T_i, 1 − dt·ΣG/C, debe ser ≥ 0):")
    print(f"  celda interior: 6 vecinas -> dx²/(6α) = {args.dx ** 2 / (6 * alpha):.1f} s")
    resultados = []
    for caso in casos:
        res = simular_bloque(caso, L, args.dx, k, args.rho, args.cp, args.Ti, pulso, puntos, args.origen, t_final,
                             args.tolerancia, 0.025, args.dt_max, n_modos=(100, 100, 60), con_exacta=not args.sin_exacta)
        resultados.append(res)
        # Celda de esquina: por eje, 1 vecina + la cara (Dirichlet: 2; Neumann: 0; Robin: h_ef dx / k).
        cara = {"dirichlet": 2.0, "neumann": 0.0,
                "robin": (1.0 / (1.0 / caso.h + args.dx / (2 * k))) * args.dx / k if caso.tipo == "robin" else 0.0}
        factor = 3 * (1 + cara[caso.tipo])
        manda = "la esquina" if factor > 6 else "el interior"
        print(f"  {caso.etiqueta(k):45s} esquina 3·(1 + {cara[caso.tipo]:.3f}) = {factor:.2f} "
              f"-> manda {manda}: dt crítico = {res.dt_critico:.1f} s")
    print(f"  tope de precisión dt_max = {args.dt_max:g} s (para seguir el pulso) -> dt usado = {resultados[0].dt:.1f} s, "
          f"{round(t_final / resultados[0].dt)} pasos por caso")

    print(f"\n{'Caso':38s} {'Tmax [°C]':>10} {'t [h]':>6} {'ΔTmax':>7} {'T eq':>6} {'estabiliza exacta/simulada [h]':>31} "
          f"{'error vs exacta: puntos / campo máx / final [°C]':>48}")
    for r in resultados:
        x, y, z = r.posicion_Tmax
        sim = "—" if r.t_estabilizacion is None else f"{r.t_estabilizacion / 3600:.1f}"
        if r.error_campo is not None:
            err_pts = max(np.abs(r.T_centro - r.T_centro_exacta)[1:].max(),
                          np.abs(r.T_puntos - r.T_puntos_exacta)[1:].max() if r.puntos else 0.0)
            errores = f"{err_pts:.3f} / {r.error_campo.max():.3f} / {r.error_campo[-1]:.1e}"
        else:
            errores = "(sin exacta)"
        print(f"{r.caso.etiqueta(k):38s} {r.Tmax:10.2f} {r.t_Tmax / 3600:6.0f} {r.dT_t.max():7.2f} {r.T_eq:6.2f} "
              f"{r.t_estabilizacion_exacta / 3600:>20.1f} / {sim:<8s} {errores:>48}")
        print(f"{'':38s} Tmax en ({x:g}, {y:g}, {z:g}); cómputo {r.computo_s:.1f} s; "
              f"límites 70/20 °C: {'CUMPLE' if r.Tmax <= 70 else 'NO CUMPLE'} / "
              f"{'CUMPLE' if r.dT_t.max() <= 20 else 'NO CUMPLE'}")

    salida = Path(args.salida)
    salida.mkdir(parents=True, exist_ok=True)
    from .graficos import animar_bloque, graficar_comparacion_bloque, graficar_cortes_bloque, graficar_puntos_bloque

    rutas = []
    resumen = {"L_m": L, "dx_m": args.dx, "k": k, "rho": args.rho, "cp": args.cp, "alpha_m2_s": alpha,
               "Ti_C": args.Ti, "Qmax_W_m3": pulso.Q_max, "T_calor_h": pulso.T_h, "h_robin": h, "beta_1_m": h / k,
               "t_final_s": t_final, "dt_s": resultados[0].dt, "puntos": puntos, "origen": args.origen, "casos": {}}
    for r in resultados:
        nombre = r.caso.nombre
        nombres, columnas = ["t_s", "T_centro"], [r.t, r.T_centro]
        if r.T_centro_exacta is not None:
            nombres.append("T_centro_exacta")
            columnas.append(r.T_centro_exacta)
        for j, p in enumerate(r.puntos):
            etiqueta = ";".join(f"{c:g}" for c in p)
            nombres.append(f"T({etiqueta})")
            columnas.append(r.T_puntos[:, j])
            if r.T_puntos_exacta is not None:
                nombres.append(f"T({etiqueta})_exacta")
                columnas.append(r.T_puntos_exacta[:, j])
        nombres += ["Tmax", "Tmin_25mm", "dT", "max_desvio_equilibrio"]
        columnas += [r.Tmax_t, r.Tmin_cascara_t, r.dT_t, r.desvio_t]
        if r.error_campo is not None:
            nombres.append("error_campo_vs_exacta")
            columnas.append(r.error_campo)
        np.savetxt(salida / f"serie_{nombre}.csv", np.column_stack(columnas), delimiter=",",
                   header=",".join(nombres), comments="", fmt="%.6f")
        rutas += [salida / f"puntos_{nombre}.png", salida / f"cortes_{nombre}.png"]
        graficar_puntos_bloque(r, rutas[-2])
        graficar_cortes_bloque(r, rutas[-1])
        resumen["casos"][nombre] = {
            "condicion": r.caso.etiqueta(k), "dt_critico_s": r.dt_critico, "Tmax_C": r.Tmax,
            "t_Tmax_h": r.t_Tmax / 3600, "posicion_Tmax_m": r.posicion_Tmax, "dTmax_C": float(r.dT_t.max()),
            "T_equilibrio_C": r.T_eq, "t_estabilizacion_exacta_h": r.t_estabilizacion_exacta / 3600,
            "t_estabilizacion_simulada_h": None if r.t_estabilizacion is None else r.t_estabilizacion / 3600,
            "error_max_campo_C": None if r.error_campo is None else float(r.error_campo.max()),
        }
    (salida / "resumen.json").write_text(json.dumps(resumen, indent=2), encoding="utf-8")
    rutas.append(salida / "comparacion.png")
    graficar_comparacion_bloque(resultados, rutas[-1])
    if not args.sin_animacion:
        rutas.append(salida / "evolucion.gif")
        animar_bloque(resultados, rutas[-1], fps=args.fps, cuadros=args.cuadros)
    print(f"\nArchivos en {salida}: " + ", ".join(r.name for r in rutas) + ", serie_<caso>.csv, resumen.json")
    return 0


def cmd_bloque_3d(args):
    import subprocess
    from pathlib import Path

    from .animacion3d import buscar_blender, crear_gif, generar_datos, script_render

    salida = Path(args.salida).resolve()  # Blender necesita rutas absolutas
    cuadros = salida / "cuadros"
    salida.mkdir(parents=True, exist_ok=True)
    for viejo in cuadros.glob("cuadro_*.png"):  # que no se mezclen cuadros de una corrida anterior
        viejo.unlink()
    datos = salida / "datos.npz"

    t0 = time.perf_counter()
    print(f"1/3 Simulando los cuatro casos con el solver 3D (dx = {args.dx:g} m)...", flush=True)
    _, n_cuadros, (T_min, T_max) = generar_datos(datos, dx=args.dx, beta=args.beta)
    print(f"    {n_cuadros} instantes, escala de colores {T_min:.1f} – {T_max:.1f} °C ({time.perf_counter() - t0:.0f} s)")
    if args.solo_datos:
        print(f"Datos: {datos}")
        return 0

    blender = buscar_blender(args.blender)
    mp4 = salida / "bloque_3d.mp4"
    print(f"2/3 Renderizando con Blender ({blender}), {args.res[0]}x{args.res[1]}, "
          f"{len(range(0, n_cuadros, args.cada))} cuadros...", flush=True)
    t1 = time.perf_counter()
    comando = [blender, "-b", "--factory-startup", "--python", str(script_render()), "--",
               str(datos), str(cuadros), "--res", str(args.res[0]), str(args.res[1]),
               "--muestras", str(args.muestras), "--cada", str(args.cada), "--fps", str(args.fps), "--mp4", str(mp4)]
    proceso = subprocess.run(comando, capture_output=True, text=True, encoding="utf-8", errors="replace")
    if proceso.returncode != 0 or not mp4.exists():
        print(proceso.stdout[-3000:], proceso.stderr[-3000:], sep="\n", file=sys.stderr)
        raise ValueError("Blender no terminó el render (ver el mensaje de arriba)")
    print(f"    listo ({time.perf_counter() - t1:.0f} s)")

    print("3/3 Armando el GIF...", flush=True)
    gif = crear_gif(cuadros, salida / "bloque_3d.gif", ancho=args.gif_ancho, fps=args.fps, cada=args.gif_cada)
    print(f"\nVideo: {mp4}\nGIF:   {gif}\nCuadros PNG: {cuadros}")
    return 0


def cmd_ifc_crear(args):
    from .ifc_io import crear_modelo_base

    ruta = crear_modelo_base(args.ruta, args.geometria)
    print(f"IFC creado con IfcOpenShell: {ruta}")
    return 0


def cmd_ifc_propiedades(args):
    from .ifc_io import PSET_VERTIDA, agregar_propiedades_curso

    asignadas = agregar_propiedades_curso(args.entrada, args.salida)
    for nombre in asignadas:
        print(f"OK: {nombre}")
    if not asignadas:
        print("No se encontró ninguna vertida (Name = Vertida_01 ... Vertida_05).")
        return 1
    print(f"{PSET_VERTIDA} escrito en {args.salida or args.entrada}")
    return 0


def cmd_ifc_info(args):
    from .ifc_io import cargar_resumen, comparar_resumenes, formatear_resumen, guardar_resumen, resumir_ifc

    resumen = resumir_ifc(args.ruta)
    print(formatear_resumen(resumen))
    if args.comparar:
        cambios = comparar_resumenes(cargar_resumen(args.comparar), resumen)
        print(f"\nComparación con {args.comparar}:")
        print("\n".join("  " + c for c in cambios) if cambios else "  Sin cambios.")
    if args.guardar:
        guardar_resumen(resumen, args.guardar)
        print(f"\nResumen guardado en {args.guardar}")
    return 0


def construir_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="mcplan", description="Simulación térmica de vertidas de hormigón masivo")
    p.add_argument("--version", action="version", version=f"mcplan {__version__}")
    sub = p.add_subparsers(dest="comando", required=True)

    s = sub.add_parser("simular1d", help="Simula un caso 1D definido en YAML")
    s.add_argument("caso", help="archivo YAML del caso")
    s.add_argument("--dz", type=float, help="sobrescribe el tamaño de celda [m]")
    s.add_argument("--dt", type=float, help="sobrescribe el paso de tiempo [s]")
    s.add_argument("--forzar", action="store_true", help="permite un dt inestable (para demostrar la divergencia)")
    s.add_argument("--sin-graficos", action="store_true")
    s.set_defaults(func=cmd_simular1d)

    s = sub.add_parser("malla", help="Sensibilidad de los resultados al tamaño de celda")
    s.add_argument("caso")
    s.add_argument("--dz", type=float, nargs="+", default=[0.2, 0.1, 0.05, 0.025])
    s.set_defaults(func=cmd_malla)

    s = sub.add_parser("verificar", help="Compara el solver con la solución analítica")
    s.add_argument("--dz", type=float, nargs="+", default=[0.2, 0.1, 0.05, 0.025])
    s.add_argument("--dt", type=float, help="paso fijo [s] (por defecto, automático)")
    s.add_argument("--forzar", action="store_true")
    s.set_defaults(func=cmd_verificar)

    s = sub.add_parser("sinusoidal", help="Caso exacto 2D/3D (E05) con orden de convergencia")
    s.add_argument("--dim", type=int, choices=(2, 3), default=2)
    s.add_argument("--L", type=float, nargs="+", default=[1.0], help="largo por eje [m] (uno para todos)")
    s.add_argument("--alpha", type=float, default=1e-4, help="difusividad [m²/s]")
    s.add_argument("--A", type=float, default=50.0, help="amplitud inicial [°C]")
    s.add_argument("--T0", type=float, default=20.0, help="temperatura del contorno [°C]")
    s.add_argument("--t", type=float, help="tiempo final [s] (por defecto 1200 en 2D y 600 en 3D)")
    s.add_argument("--n", type=int, nargs="+", help="celdas por eje (por defecto 15 30 60 en 2D, 8 16 32 en 3D)")
    s.add_argument("--dt", type=float, help="paso fijo [s] (por defecto, automático)")
    s.add_argument("--forzar", action="store_true", help="permite un dt inestable (para ver la divergencia)")
    s.set_defaults(func=cmd_sinusoidal)

    s = sub.add_parser("fuente", help="Fuente uniforme en dominio adiabático (E05)")
    s.add_argument("--dim", type=int, choices=(2, 3), default=3)
    s.add_argument("--q0", type=float, default=10000.0, help="fuente volumétrica [W/m³]")
    s.add_argument("--rho", type=float, default=2400.0, help="densidad [kg/m³]")
    s.add_argument("--cp", type=float, default=900.0, help="calor específico [J/(kg·K)]")
    s.add_argument("--T0", type=float, default=20.0, help="temperatura inicial [°C]")
    s.add_argument("--t", type=float, default=8 * 3600.0, help="tiempo final [s]")
    s.add_argument("--n", type=int, default=5, help="celdas por eje")
    s.set_defaults(func=cmd_fuente)

    s = sub.add_parser("energia", help="Conservación de energía en dominio adiabático (E05)")
    s.add_argument("--dim", type=int, choices=(2, 3), default=3)
    s.add_argument("--t", type=float, default=3 * 86400.0, help="tiempo final [s]")
    s.add_argument("--n", type=int, default=10, help="celdas por eje")
    s.add_argument("--dos-materiales", action="store_true", help="mitad hormigón, mitad suelo")
    s.set_defaults(func=cmd_energia)

    s = sub.add_parser("rectangulo", help="Rectángulo 2D a Ti con el borde a Tb: dt crítico, animación y T(t) en puntos")
    s.add_argument("--L", type=float, nargs=2, default=[4.0, 2.0], metavar=("LX", "LY"), help="lados [m]")
    s.add_argument("--dx", type=float, default=0.1, help="tamaño de celda (dx = dy) [m]")
    s.add_argument("--k", type=float, default=2.222, help="conductividad [W/(m·K)] (hormigón de Ju y Lei)")
    s.add_argument("--rho", type=float, default=2400.0, help="densidad [kg/m³]")
    s.add_argument("--cp", type=float, default=970.0, help="calor específico [J/(kg·K)]")
    s.add_argument("--alpha", type=float, help="difusividad [m²/s]: si se indica, k = alpha·rho·cp")
    s.add_argument("--Ti", type=float, default=20.0, help="temperatura inicial [°C]")
    s.add_argument("--Tb", type=float, default=10.0, help="temperatura del borde [°C]")
    s.add_argument("--punto", type=float, nargs=2, action="append", metavar=("X", "Y"),
                   help="punto a seguir (repetible); por defecto (0, 0) y (1, 0.5)")
    s.add_argument("--origen", choices=("centro", "esquina"), default="centro", help="dónde está (0, 0)")
    s.add_argument("--t", type=float, help="tiempo final [s] (por defecto, 1,2 × el de estabilización)")
    s.add_argument("--tolerancia", type=float, default=0.1, help="|T − Tb| para considerar estable [°C]")
    s.add_argument("--cuadros", type=int, default=150, help="instantes guardados (cuadros del GIF)")
    s.add_argument("--fps", type=int, default=15)
    s.add_argument("--hidratacion", nargs="?", const="ju-lei", choices=("ju-lei", "pulso"),
                   help="calor de hidratación: sola = Ju y Lei (2019); 'pulso' = Q(t) = Qmax/2·(1 − cos(2πt/T)), t < T")
    s.add_argument("--Qmax", type=float, default=1800.0, help="pulso: calor máximo [W/m³] (≈ 2,79 °C/h de Ju y Lei)")
    s.add_argument("--T-calor", type=float, default=42.0, help="pulso: duración T [h] (máximo en T/2 = 21 h)")
    s.add_argument("--dt-max", type=float, help="tope de dt por precisión [s] (por defecto 150 con hidratación)")
    s.add_argument("--sin-animacion", action="store_true")
    s.add_argument("--salida", help="carpeta (por defecto resultados/rectangulo_2d[_hidratacion])")
    s.set_defaults(func=cmd_rectangulo)

    s = sub.add_parser("bloque", help="Bloque 3D con el pulso de calor: Dirichlet 20/30 °C, Neumann y Robin")
    s.add_argument("--L", type=float, nargs=3, default=[3.0, 3.0, 0.6], metavar=("LX", "LY", "LZ"), help="lados [m]")
    s.add_argument("--dx", type=float, default=0.05, help="tamaño de celda (dx = dy = dz) [m]")
    s.add_argument("--casos", nargs="+", choices=("dirichlet", "neumann", "robin"),
                   default=["dirichlet", "neumann", "robin"])
    s.add_argument("--T-dirichlet", type=float, nargs="+", default=[20.0, 30.0], help="T de las caras [°C]")
    s.add_argument("--beta", type=float, help="Robin: β = h/k [1/m] (por defecto h = 27,8 W/(m²K) de Ju y Lei)")
    s.add_argument("--T-aire", type=float, default=20.0, help="Robin: temperatura del aire [°C]")
    s.add_argument("--Ti", type=float, default=20.0, help="temperatura inicial [°C]")
    s.add_argument("--Qmax", type=float, default=1800.0, help="pulso: calor máximo [W/m³]")
    s.add_argument("--T-calor", type=float, default=42.0, help="pulso: duración T [h]")
    s.add_argument("--k", type=float, default=2.222)
    s.add_argument("--rho", type=float, default=2400.0)
    s.add_argument("--cp", type=float, default=970.0)
    s.add_argument("--punto", type=float, nargs=3, action="append", metavar=("X", "Y", "Z"),
                   help="punto a seguir además del centro (repetible); por defecto (0, 0, 0.275) y (1.2, 1.2, 0.2)")
    s.add_argument("--origen", choices=("centro", "esquina"), default="centro")
    s.add_argument("--t", type=float, help="tiempo final [s] (por defecto 1,2 × la estabilización exacta más lenta)")
    s.add_argument("--tolerancia", type=float, default=0.1, help="|T − T_eq| para considerar estable [°C]")
    s.add_argument("--dt-max", type=float, default=150.0, help="tope de dt por precisión [s]")
    s.add_argument("--cuadros", type=int, default=100)
    s.add_argument("--fps", type=int, default=12)
    s.add_argument("--sin-animacion", action="store_true")
    s.add_argument("--sin-exacta", action="store_true", help="no calcula la solución exacta (más rápido)")
    s.add_argument("--salida", default="resultados/bloque_3d")
    s.set_defaults(func=cmd_bloque)

    s = sub.add_parser("bloque-3d", help="Animación 3D (Blender) de los cuatro casos del bloque: MP4 y GIF")
    s.add_argument("--dx", type=float, default=0.05, help="celda de la simulación [m]")
    s.add_argument("--beta", type=float, help="Robin: β = h/k [1/m] (por defecto h = 27,8 W/(m²K))")
    s.add_argument("--res", type=int, nargs=2, default=[1280, 720], metavar=("ANCHO", "ALTO"))
    s.add_argument("--muestras", type=int, default=32, help="muestras de EEVEE por cuadro (calidad)")
    s.add_argument("--cada", type=int, default=1, help="renderiza 1 de cada N instantes (1 instante = 1 h)")
    s.add_argument("--fps", type=int, default=15)
    s.add_argument("--gif-ancho", type=int, default=720)
    s.add_argument("--gif-cada", type=int, default=2, help="el GIF usa 1 de cada N cuadros (más liviano)")
    s.add_argument("--blender", help="ruta a blender.exe (por defecto la busca)")
    s.add_argument("--solo-datos", action="store_true", help="sólo simula y guarda datos.npz (sin Blender)")
    s.add_argument("--salida", default="resultados/bloque_3d_render")
    s.set_defaults(func=cmd_bloque_3d)

    s = sub.add_parser("ifc-crear", help="Genera la fundación oficial con IfcOpenShell (sin Blender)")
    s.add_argument("ruta", help="archivo .ifc de salida")
    s.add_argument("--geometria", default="geometry/course_geometry.json")
    s.set_defaults(func=cmd_ifc_crear)

    s = sub.add_parser("ifc-propiedades", help="Agrega/actualiza UANDES_HormigonMasivo en las vertidas")
    s.add_argument("entrada", help="archivo .ifc (se modifica si no se indica --salida)")
    s.add_argument("--salida")
    s.set_defaults(func=cmd_ifc_propiedades)

    s = sub.add_parser("ifc-info", help="Lee un IFC: elementos, propiedades y geometría")
    s.add_argument("ruta")
    s.add_argument("--guardar", help="guarda la lectura en JSON (para comparar después)")
    s.add_argument("--comparar", help="JSON de una lectura anterior: muestra qué cambió")
    s.set_defaults(func=cmd_ifc_info)
    return p


def main(argv=None) -> int:
    # La consola de Windows no siempre usa UTF-8: evita errores con tildes y "°".
    for flujo in (sys.stdout, sys.stderr):
        if hasattr(flujo, "reconfigure"):
            flujo.reconfigure(encoding="utf-8", errors="replace")
    args = construir_parser().parse_args(argv)
    try:
        return args.func(args)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
