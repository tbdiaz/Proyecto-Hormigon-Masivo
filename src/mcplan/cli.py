"""Interfaz de línea de comandos.

    mcplan simular1d config/caso_base_1d.yaml [--dz 0.1] [--dt 600] [--forzar]
    mcplan malla     config/caso_base_1d.yaml --dz 0.2 0.1 0.05 0.025
    mcplan verificar [--dz 0.2 0.1 0.05 0.025] [--dt 600] [--forzar]
    mcplan ifc-crear [models/fundacion_base.ifc]
    mcplan ifc-info  models/fundacion_base.ifc [--guardar antes.json] [--comparar antes.json]
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


def cmd_ifc_crear(args):
    from .ifc_io import crear_modelo_base

    ruta = crear_modelo_base(args.ruta)
    print(f"IFC creado: {ruta}")
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

    s = sub.add_parser("ifc-crear", help="Genera el IFC base de la fundación")
    s.add_argument("ruta", nargs="?", default="models/fundacion_base.ifc")
    s.set_defaults(func=cmd_ifc_crear)

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
