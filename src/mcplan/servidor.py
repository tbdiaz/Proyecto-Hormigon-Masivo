"""Editor local: cambiar los datos de un caso 3D y simular desde el navegador, sin usar la terminal.

Un HTML abierto como archivo no puede ejecutar Python, así que `mcplan editor` levanta un
servidor sólo para este computador (127.0.0.1). La página muestra los valores del caso
(YAML + pset del IFC); al simular, el servidor compara el formulario con esos valores y
corre `mcplan simular3d` con los cambios como argumentos, en un proceso aparte. Por eso una
corrida del editor es igual a la misma corrida en consola, y su panel muestra el comando.

Los YAML y el IFC no se modifican: cada corrida queda en resultados/editor/<fecha_hora>/.
"""

import json
import math
import os
import re
import subprocess
import sys
import threading
import time
import webbrowser
from datetime import datetime
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from types import SimpleNamespace

PUERTO = 8765
_ID_CORRIDA = re.compile(r"^\d{8}_\d{6}$")
_FECHA = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}$")
_MAX_CUERPO = 200_000


def casos_disponibles(raiz: Path) -> list[str]:
    """Casos 3D: los YAML de config/ que tienen sección `ifc`."""
    import yaml

    casos = []
    for ruta in sorted((raiz / "config").glob("*.yaml")):
        cfg = yaml.safe_load(ruta.read_text(encoding="utf-8"))
        if isinstance(cfg, dict) and "ifc" in cfg:
            casos.append(ruta.name)
    return casos


def opciones_caso(raiz: Path, nombre: str) -> dict:
    """Valores actuales de un caso: lo que el formulario muestra y contra lo que se comparan los cambios."""
    from .ifc_io import leer_vertidas
    from .simulacion3d import aplicar_cambios, cargar_caso_3d

    if nombre not in casos_disponibles(raiz):
        raise ValueError(f"Caso desconocido: {nombre}")
    caso = cargar_caso_3d(raiz / "config" / nombre)
    vertidas = leer_vertidas(caso.ifc)
    vox = SimpleNamespace(nombres=[v.nombre for v in vertidas], propiedades=[dict(v.propiedades) for v in vertidas])
    aplicar_cambios(vox, caso.cambios)

    moldaje = aislacion = None
    if caso.capas:
        m = caso.cubiertas_cfg.get("moldaje") or {}
        retiro = m.get("retiro_h", "pset")
        if retiro == "pset":
            retiro = vox.propiedades[0].get("RetiroMoldaje_h")
        moldaje = {"usar": m.get("capa") not in (None, "ninguna", "ninguno"), "capa": m.get("capa"), "retiro_h": retiro}
        a = caso.cubiertas_cfg.get("aislacion") or {}
        capa = a.get("capa") if a.get("capa") not in (None, "ninguno") else "ninguna"
        retiro_a = a.get("retiro_h")
        aislacion = {"capa": capa, "caras": a.get("caras", "superiores"),
                     "retiro_h": None if retiro_a in (None, "nunca") else float(retiro_a)}
    con_clima = caso.clima_cfg is not None
    sup = caso.superficie_cfg
    amb = caso.ambiente
    modelo = sup.get("modelo", "fisico") if con_clima else "combinado"
    return {
        "caso": nombre, "nombre": caso.nombre, "dx": caso.dx, "duracion_h": caso.duracion_h, "dt_max_s": caso.dt_max_s,
        "fecha_inicio": caso.fecha_inicio if con_clima else None, "con_clima": con_clima, "modelo": modelo,
        "h_superior": float(sup.get("h_superior", 27.8)) if con_clima else amb.h_superior,
        "h_lateral": float(sup.get("h_lateral", 27.8)) if con_clima else amb.h_lateral,
        "absortividad": float(sup.get("absortividad", 0.65)), "emisividad": float(sup.get("emisividad", 0.9)),
        "k_hormigon": caso.hormigon.k, "k_suelo": caso.suelo.k,
        "profundidad_suelo": caso.profundidad_suelo, "margen_suelo": caso.margen_suelo,
        "vertidas": [{"id": i + 1, "nombre": vox.nombres[i], "inicio_h": float(p.get("Inicio_h", 0.0)),
                      "T_colocacion": float(p.get("TemperaturaColocacion_C", 20.0)), "activa": bool(p.get("activa", True))}
                     for i, p in enumerate(vox.propiedades)],
        "moldaje": moldaje, "aislacion": aislacion, "capas": sorted(caso.capas),
        "cobertura": _cobertura(caso) if con_clima else None,
    }


def _cobertura(caso) -> str:
    """Qué fechas cubren los snapshots disponibles, para orientar al elegir la fecha de inicio."""
    from .clima import leer_snapshot, ultimo_snapshot

    partes = []
    carpeta = caso.clima_cfg["carpeta"]
    ruta = ultimo_snapshot("pronostico", carpeta)
    if ruta is not None:
        h = leer_snapshot(ruta)["metadata"]["horizonte"]["ultima_hora_con_dato"]["temperature_2m"]
        partes.append(f"pronóstico hasta {h[:10]}")
    ruta = ultimo_snapshot("reanalisis", carpeta)
    if ruta is not None:
        m = leer_snapshot(ruta)["metadata"]
        partes.append(f"climatología ERA5 de {m['horizonte']['desde'][5:10]} a {m['horizonte']['hasta'][5:10]}")
    return "; ".join(partes)


def _numero(form, clave, minimo, maximo, etiqueta):
    try:
        x = float(form[clave])
    except (KeyError, TypeError, ValueError):
        raise ValueError(f"{etiqueta}: ingresa un número") from None
    if not math.isfinite(x) or not minimo <= x <= maximo:
        raise ValueError(f"{etiqueta}: debe estar entre {minimo:g} y {maximo:g}")
    return x


def _distinto(a, b) -> bool:
    return b is None or abs(a - float(b)) > 1e-9


def construir_argumentos(form: dict, base: dict) -> list[str]:
    """Valida el formulario y devuelve los argumentos de `mcplan simular3d`: sólo lo que difiere del caso."""
    args = []

    def opcion(clave, bandera, minimo, maximo, etiqueta):
        x = _numero(form, clave, minimo, maximo, etiqueta)
        if _distinto(x, base[clave]):
            args.extend([bandera, f"{x:g}"])

    opcion("dx", "--dx", 0.05, 2.0, "Δx")
    opcion("duracion_h", "--duracion", 1, 720, "Duración")
    opcion("dt_max_s", "--dt-max", 10, 3600, "Δt máximo")
    opcion("k_hormigon", "--k-hormigon", 0.1, 10, "k del hormigón")
    opcion("k_suelo", "--k-suelo", 0.05, 10, "k del suelo")
    if base["modelo"] == "combinado":
        opcion("h_superior", "--h-superior", 0, 200, "h superior")
        opcion("h_lateral", "--h-lateral", 0, 200, "h lateral")
    else:
        opcion("absortividad", "--absortividad", 0, 1, "Absortividad")
        opcion("emisividad", "--emisividad", 0, 1, "Emisividad")

    if base["con_clima"]:
        fecha = str(form.get("fecha_inicio", ""))
        if not _FECHA.match(fecha):
            raise ValueError("Fecha de inicio: usa el formato AAAA-MM-DDTHH:MM")
        try:
            datetime.fromisoformat(fecha)
        except ValueError:
            raise ValueError(f"Fecha de inicio inválida: {fecha}") from None
        if fecha != base["fecha_inicio"]:
            args.extend(["--fecha-inicio", fecha])

    if len(form.get("vertidas", [])) != len(base["vertidas"]):
        raise ValueError("El formulario no coincide con las vertidas del caso")
    inicios, temperaturas, ocultar, mostrar = [], [], [], []
    for f, b in zip(form["vertidas"], base["vertidas"]):
        v = b["id"]
        activa = bool(f.get("activa"))
        if activa != b["activa"]:
            (mostrar if activa else ocultar).append(str(v))
        if activa:
            inicio = _numero(f, "inicio_h", 0, 720, f"{b['nombre']}: inicio")
            T = _numero(f, "T_colocacion", -5, 45, f"{b['nombre']}: T de colocación")
            if _distinto(inicio, b["inicio_h"]):
                inicios.append(f"{v}={inicio:g}")
            if _distinto(T, b["T_colocacion"]):
                temperaturas.append(f"{v}={T:g}")
    for bandera, valores in (("--inicio", inicios), ("--T-colocacion", temperaturas), ("--ocultar", ocultar),
                             ("--mostrar", mostrar)):
        if valores:
            args.extend([bandera, *valores])
    if not any(bool(f.get("activa")) for f in form["vertidas"]):
        raise ValueError("Debe quedar al menos una vertida activa")

    if base["moldaje"] is not None:
        m, bm = form.get("moldaje") or {}, base["moldaje"]
        if bm["usar"] and not m.get("usar"):
            args.append("--sin-moldaje")
        elif m.get("usar"):
            retiro = _numero(m, "retiro_h", 0, 2000, "Retiro del moldaje")
            if _distinto(retiro, bm["retiro_h"]):
                args.extend(["--retiro-moldaje", f"{retiro:g}"])
    if base["aislacion"] is not None:
        a, ba = form.get("aislacion") or {}, base["aislacion"]
        capa = a.get("capa", "ninguna")
        if capa != "ninguna" and capa not in base["capas"]:
            raise ValueError(f"Aislación desconocida: {capa}")
        if capa != ba["capa"]:
            args.extend(["--aislacion", capa])
        if capa != "ninguna":
            if a.get("caras") not in ("superiores", "laterales", "todas"):
                raise ValueError("Caras de la aislación inválidas")
            if a["caras"] != ba["caras"]:
                args.extend(["--aislacion-caras", a["caras"]])
            if a.get("retiro_h") not in (None, ""):
                retiro = _numero(a, "retiro_h", 0, 2000, "Retiro de la aislación")
                if _distinto(retiro, ba["retiro_h"]):
                    args.extend(["--retiro-aislacion", f"{retiro:g}"])
    return args


class Trabajos:
    """Una simulación a la vez, en un proceso aparte; guarda su salida para mostrar el avance."""

    def __init__(self, raiz: Path):
        self.raiz = raiz
        self.candado = threading.Lock()
        self.proceso = None
        self.estado = "libre"
        self.id = None
        self.comando = None
        self.lineas = []
        self.avance = None

    def iniciar(self, caso: str, argumentos: list[str], archivos: bool) -> dict:
        with self.candado:
            if self.proceso is not None:
                raise ValueError("Ya hay una simulación corriendo")
            self.id = datetime.now().strftime("%Y%m%d_%H%M%S")
            salida = f"resultados/editor/{self.id}"
            extra = ["--vtk"] if archivos else ["--sin-figuras"]
            pedido = ["simular3d", f"config/{caso}", *argumentos, *extra, "--salida", salida]
            self.comando = "mcplan " + " ".join(pedido)
            self.lineas, self.avance, self.estado = [], None, "corriendo"
            entorno = {**os.environ, "PYTHONIOENCODING": "utf-8", "PYTHONUNBUFFERED": "1"}
            self.proceso = subprocess.Popen(
                [sys.executable, "-m", "mcplan", *pedido, "--progreso"], cwd=self.raiz, env=entorno, text=True,
                encoding="utf-8", errors="replace", stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            threading.Thread(target=self._leer, args=(self.proceso,), daemon=True).start()
            return self.resumen()

    def _leer(self, proceso):
        for linea in proceso.stdout:
            linea = linea.rstrip()
            if linea.startswith("PROGRESO "):
                _, t, total = linea.split()
                self.avance = [float(t), float(total)]
            else:
                self.lineas.append(linea)
        codigo = proceso.wait()
        with self.candado:
            if self.estado == "corriendo":
                self.estado = "listo" if codigo == 0 else "error"
            self.proceso = None

    def cancelar(self):
        with self.candado:
            if self.proceso is not None:
                self.estado = "cancelado"
                self.proceso.terminate()

    def resumen(self) -> dict:
        return {"estado": self.estado, "id": self.id, "comando": self.comando, "avance": self.avance,
                "salida": self.lineas[-200:]}


def corridas(raiz: Path) -> list[dict]:
    """Corridas anteriores del editor, de la más nueva a la más antigua."""
    filas = []
    carpeta = raiz / "resultados" / "editor"
    if not carpeta.is_dir():
        return filas
    for d in sorted(carpeta.iterdir(), reverse=True):
        if not (_ID_CORRIDA.match(d.name) and (d / "panel.html").is_file() and (d / "resumen.json").is_file()):
            continue
        r = json.loads((d / "resumen.json").read_text(encoding="utf-8"))
        filas.append({"id": d.name, "caso": Path(r["caso"]).stem, "dx": r["dx_m"], "duracion_h": r["duracion_h"],
                      "Tmax": r["Tmax_C"], "dT": r["dTmax_C"]})
    return filas


def crear_servidor(raiz: Path, puerto: int = PUERTO) -> ThreadingHTTPServer:
    trabajos = Trabajos(raiz)
    pagina = (Path(__file__).parent / "editor.html").read_bytes()

    class Manejador(BaseHTTPRequestHandler):
        def log_message(self, *_):
            pass

        def _responder(self, codigo, cuerpo: bytes, tipo="application/json; charset=utf-8"):
            self.send_response(codigo)
            self.send_header("Content-Type", tipo)
            self.send_header("Content-Length", str(len(cuerpo)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(cuerpo)

        def _json(self, datos, codigo=200):
            self._responder(codigo, json.dumps(datos, ensure_ascii=False).encode("utf-8"))

        def _origen_propio(self) -> bool:
            """Sólo este computador: evita que otra página web le pida simulaciones al servidor."""
            puerto_real = self.server.server_address[1]
            propios = {f"127.0.0.1:{puerto_real}", f"localhost:{puerto_real}"}
            origen = self.headers.get("Origin")
            return self.headers.get("Host") in propios and (
                origen is None or origen.split("://", 1)[-1] in propios)

        def do_GET(self):
            if not self._origen_propio():
                return self._json({"error": "Origen no permitido"}, 403)
            ruta, _, consulta = self.path.partition("?")
            try:
                if ruta == "/":
                    return self._responder(200, pagina, "text/html; charset=utf-8")
                if ruta == "/api/casos":
                    return self._json({"casos": casos_disponibles(raiz)})
                if ruta == "/api/opciones":
                    return self._json(opciones_caso(raiz, consulta.removeprefix("caso=")))
                if ruta == "/api/estado":
                    return self._json(trabajos.resumen())
                if ruta == "/api/corridas":
                    return self._json({"corridas": corridas(raiz)})
                m = re.fullmatch(r"/resultados/(\d{8}_\d{6})/panel\.html", ruta)
                if m and (raiz / "resultados/editor" / m[1] / "panel.html").is_file():
                    cuerpo = (raiz / "resultados/editor" / m[1] / "panel.html").read_bytes()
                    return self._responder(200, cuerpo, "text/html; charset=utf-8")
                return self._json({"error": "No encontrado"}, 404)
            except ValueError as e:
                return self._json({"error": str(e)}, 400)

        def do_POST(self):
            if not self._origen_propio():
                return self._json({"error": "Origen no permitido"}, 403)
            largo = int(self.headers.get("Content-Length") or 0)
            if largo > _MAX_CUERPO:
                return self._json({"error": "Pedido demasiado grande"}, 413)
            try:
                if self.path == "/api/cancelar":
                    trabajos.cancelar()
                    return self._json(trabajos.resumen())
                if self.path == "/api/simular":
                    form = json.loads(self.rfile.read(largo).decode("utf-8"))
                    base = opciones_caso(raiz, str(form.get("caso")))
                    argumentos = construir_argumentos(form, base)
                    return self._json(trabajos.iniciar(base["caso"], argumentos, bool(form.get("archivos"))))
                return self._json({"error": "No encontrado"}, 404)
            except (ValueError, json.JSONDecodeError) as e:
                return self._json({"error": str(e)}, 400)

    try:
        return ThreadingHTTPServer(("127.0.0.1", puerto), Manejador)
    except OSError:  # puerto ocupado (p. ej. otro editor abierto): el sistema asigna uno libre
        return ThreadingHTTPServer(("127.0.0.1", 0), Manejador)


def raiz_proyecto() -> Path:
    """La carpeta del repositorio: la actual si tiene config/, o la que contiene src/mcplan."""
    if (Path.cwd() / "config").is_dir():
        return Path.cwd()
    return Path(__file__).resolve().parents[2]


def ejecutar(puerto: int = PUERTO, abrir: bool = True) -> int:
    servidor = crear_servidor(raiz_proyecto(), puerto)
    url = f"http://127.0.0.1:{servidor.server_address[1]}/"
    print(f"Editor en {url}\nDeja esta ventana abierta mientras lo usas; ciérrala para apagarlo.", flush=True)
    if abrir:
        threading.Timer(0.5, webbrowser.open, args=(url,)).start()
    try:
        servidor.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        servidor.server_close()
        time.sleep(0.1)
    return 0
