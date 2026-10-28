"""Interfaz de línea de comandos.

    mcplan simular1d config/caso_base_1d.yaml [--dz 0.1] [--dt 600] [--forzar]
    mcplan malla     config/caso_base_1d.yaml --dz 0.2 0.1 0.05 0.025
    mcplan verificar [--dz 0.2 0.1 0.05 0.025] [--dt 600] [--forzar]
    mcplan sinusoidal --dim 2|3 [--L 1 2] [--alpha 1e-4] [--A 50] [--T0 20] [--t 1200] [--n 15 30] [--dt] [--forzar]
    mcplan fuente     --dim 2|3 [--q0 10000] [--rho 2400] [--cp 900] [--T0 20] [--t 28800]
    mcplan energia    --dim 2|3 [--dos-materiales] [--t 259200]
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
