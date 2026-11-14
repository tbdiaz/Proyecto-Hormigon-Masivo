"""Editor local: el formulario se traduce a los mismos argumentos que la línea de comandos."""

import copy
import json
import threading
import urllib.error
import urllib.request
from pathlib import Path

import pytest

from mcplan.cli import construir_parser
from mcplan.servidor import casos_disponibles, construir_argumentos, crear_servidor, opciones_caso

RAIZ = Path(".").resolve()


@pytest.fixture(scope="module")
def base():
    return opciones_caso(RAIZ, "caso_14d.yaml")


def _form(base):
    """El formulario tal como lo manda el editor sin cambios: los valores del caso, sin tubos."""
    f = copy.deepcopy(base)
    f["tubos"] = []
    return f


def test_casos_y_valores_del_ifc(base):
    assert {"caso_3d.yaml", "caso_14d.yaml"} <= set(casos_disponibles(RAIZ))
    assert [v["inicio_h"] for v in base["vertidas"]] == [0, 6, 12, 18, 24]
    assert base["moldaje"]["retiro_h"] == 72.0 and base["aislacion"]["capa"] == "ninguna"
    with pytest.raises(ValueError):
        opciones_caso(RAIZ, "../../secreto.yaml")


def test_sin_cambios_no_hay_argumentos(base):
    assert construir_argumentos(_form(base), base) == []


def test_cambios_se_traducen_a_opciones_validas_del_cli(base):
    f = _form(base)
    f.update(dx=0.5, duracion_h=48, fecha_inicio="2026-11-08T08:00", absortividad=0.5, k_hormigon=2.5)
    f["vertidas"][1]["inicio_h"] = 12
    f["vertidas"][2]["T_colocacion"] = 10
    f["vertidas"][4]["activa"] = False
    f["moldaje"]["retiro_h"] = 24
    f["aislacion"].update(capa="manta_5cm", caras="todas", retiro_h=300)
    args = construir_argumentos(f, base)
    a = construir_parser().parse_args(["simular3d", "config/caso_14d.yaml", *args])
    assert (a.dx, a.duracion, a.fecha_inicio, a.absortividad, a.k_hormigon) == (0.5, 48, "2026-11-08T08:00", 0.5, 2.5)
    assert a.inicio == ["2=12"] and a.T_colocacion == ["3=10"] and a.ocultar == [5]
    assert a.retiro_moldaje == 24 and a.aislacion == "manta_5cm"
    assert a.aislacion_caras == "todas" and a.retiro_aislacion == 300


@pytest.mark.parametrize("cambio", [
    {"dx": 0.01}, {"dx": "abc"}, {"duracion_h": 0}, {"fecha_inicio": "ayer; rm -rf"},
    {"aislacion": {"capa": "--vtk", "caras": "todas"}},
])
def test_datos_invalidos_se_rechazan(base, cambio):
    with pytest.raises(ValueError):
        construir_argumentos({**_form(base), **cambio}, base)


def test_servidor_solo_acepta_su_propio_origen():
    servidor = crear_servidor(RAIZ, 0)
    threading.Thread(target=servidor.serve_forever, daemon=True).start()
    url = f"http://127.0.0.1:{servidor.server_address[1]}"
    try:
        with urllib.request.urlopen(f"{url}/api/casos") as r:
            assert "caso_14d.yaml" in json.load(r)["casos"]
        ajeno = urllib.request.Request(f"{url}/api/simular", data=b"{}", method="POST",
                                       headers={"Origin": "https://sitio-malicioso.example"})
        with pytest.raises(urllib.error.HTTPError) as e:
            urllib.request.urlopen(ajeno)
        assert e.value.code == 403
        with pytest.raises(urllib.error.HTTPError) as e:
            urllib.request.urlopen(f"{url}/resultados/../../README.md")
        assert e.value.code == 404
    finally:
        servidor.shutdown()
        servidor.server_close()


def test_tubos_curva_y_k_por_vertida_desde_el_editor(base):
    f = _form(base)
    f["vertidas"][2]["k"] = 4.0
    circuito = {"vertida": 3, "z": "", "material": "cobre", "direccion": "x", "separacion": 0.8, "diametro": 32,
                "caudal": 25, "T_agua": 12, "dias": 10}
    f["tubos"] = [circuito, {**circuito, "z": 1.4, "direccion": "y", "material": "pvc"}]
    f["curva"] = "M4_22A"
    args = construir_argumentos(f, base)
    a = construir_parser().parse_args(["simular3d", "config/caso_14d.yaml", *args])
    assert a.k_vertida == ["3=4"] and a.curva == "M4_22A" and len(a.tubo) == 2
    assert "z=1.4" in a.tubo[1] and "material=pvc" in a.tubo[1] and not any(p.startswith("z=") for p in a.tubo[0])
    for malo in ({"material": "oro"}, {"z": 0.2}, {"direccion": "z"}, {"vertida": 9}):
        f["tubos"] = [{**circuito, **malo}]
        with pytest.raises(ValueError):
            construir_argumentos(f, base)
    f["tubos"], f["curva"] = [], "inventada"
    with pytest.raises(ValueError):
        construir_argumentos(f, base)


def test_simular3d_con_tubo_y_curva_medida(tmp_path):
    """De punta a punta: un circuito a una altura dada, una curva del profesor y el panel con sus series."""
    import json as js

    from mcplan.cli import main

    salida = tmp_path / "corrida"
    assert main(["simular3d", "config/caso_14d.yaml", "--dx", "0.5", "--duracion", "30", "--sin-figuras",
                 "--tubo", "vertida=3", "z=1.6", "material=cobre", "T_agua=10", "dias=1",
                 "--curva", "M4_22A", "--salida", str(salida)]) == 0
    r = js.loads((salida / "resumen.json").read_text(encoding="utf-8"))
    assert list(r["tubos"]) == ["tubo1_V3"] and r["tubos"]["tubo1_V3"]["energia_MJ"] > 0
    html = (salida / "panel.html").read_text(encoding="utf-8")
    d = js.loads(html.split('id="datos">')[1].split("</script>")[0])
    assert d["tubos"][0]["puntos"][0][2] == pytest.approx(1.6) and d["tubos_series"]["circuitos"][0]["potencia_kW"]
    assert "M4_22A" in d["hidratacion"]["descripcion"]
