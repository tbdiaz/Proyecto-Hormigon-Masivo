# Proyecto de Hormigón Masivo
# Grupo 4

Planificación de vertidas de una fundación de hormigón masivo (aerogenerador,
caso Ju y Lei 2019) mediante un solver térmico propio, BIM/IFC y clima real.

## Integrantes
- Trinidad Diaz
- Pablo Arancibia
- Eduardo Vergara

## Sitio climático asignado
Litueche, Región de O'Higgins — **lat −34.105457, lon −71.649539**
(coordenada oficial del curso, en [`config/sitio.yaml`](config/sitio.yaml)).

## Registro inicial
E00 — Kickoff y registro de repositorio.
Fecha: 15 de octubre de 2026.

---

## Instalación

Requiere Python 3.11 o superior (probado con 3.12 en Windows 11).

```powershell
# Windows (PowerShell)
py -3.12 -m venv .venv
.venv\Scripts\python -m pip install -r requirements.txt
.venv\Scripts\python -m pip install -e .
```

```bash
# Linux / macOS / cluster
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -m pip install -e .
```

Luego se activa el entorno (`.venv\Scripts\Activate.ps1` o `source .venv/bin/activate`) o se
llama directo a `.venv\Scripts\mcplan`.

## Uso rápido

```powershell
pytest                                            # 30 tests automatizados
mcplan verificar                                  # solver vs. solución analítica (orden de convergencia)
mcplan simular1d config/caso_base_1d.yaml         # caso base 1D: métricas + gráficos en resultados/
mcplan simular1d config/caso_ifc_1d.yaml          # mismo caso, pero el hormigón se lee desde el IFC
mcplan malla config/caso_base_1d.yaml             # sensibilidad de Tmax y ΔT al tamaño de celda
mcplan ifc-info models/fundacion_base.ifc         # elementos, property sets y geometría del IFC
```

Para cambiar Δz o Δt sin editar el YAML:

```powershell
mcplan simular1d config/caso_base_1d.yaml --dz 0.1
mcplan simular1d config/caso_base_1d.yaml --dt 1200
mcplan simular1d config/caso_base_1d.yaml --dt 3000            # se rechaza: dt > dt crítico
mcplan simular1d config/caso_base_1d.yaml --dt 3000 --forzar   # se ve la divergencia
```

## Round-trip Bonsai → IFC → Python

```powershell
mcplan ifc-info models/fundacion_base.ifc --guardar resultados/ifc_antes.json
# ... abrir models/fundacion_base.ifc en Blender + Bonsai, modificar una propiedad
#     (MCPlan_Vertida) o la geometría, y guardar (Ctrl+S) ...
mcplan ifc-info models/fundacion_base.ifc --comparar resultados/ifc_antes.json
mcplan simular1d config/caso_ifc_1d.yaml      # la simulación toma el cambio
```

El modelo base se regenera con `mcplan ifc-crear models/fundacion_base.ifc`.

## Estructura

```
config/              casos YAML (decisiones y parámetros) y sitio del grupo
models/              modelo IFC de la fundación (editable en Bonsai)
src/mcplan/
  solver1d.py        balance de energía por celda, Euler explícito, bordes Dirichlet/Neumann/Robin
  hidratacion.py     calor de hidratación: curva adiabática Ju y Lei (2019) + edad equivalente
  metricas.py        Tmax y ΔT = Tmax − Tmin(25 mm bajo la superficie)
  caso.py            lectura de casos YAML
  ifc_io.py          creación/lectura del IFC, detección de cambios, columna IFC → solver
  verificacion.py    solución analítica y estudio de convergencia
  graficos.py        figuras (Matplotlib)
  cli.py             comando `mcplan`
tests/               tests automatizados (pytest)
docs/modelo_1d.md    ecuaciones, discretización, estabilidad, verificación y supuestos
```

## Arquitectura

```
Bonsai / Blender  →  IFC  →  IfcOpenShell  →  geometría + propiedades  →  voxelización  →  solver
```

El solver no depende de Blender: sólo lee el archivo `.ifc`. En esta versión, la
"voxelización" es una columna vertical en el eje de la fundación (`columna_desde_ifc`).

## Estado actual (semana 1)

| Componente | Estado |
|---|---|
| Solver 1D (conducción, multicapa, Dirichlet/Neumann/Robin, hidratación) | ✅ verificado (orden 2) |
| Modelo IFC con 3 vertidas y property set `MCPlan_Vertida` | ✅ probado en Bonsai (Blender 5.2) |
| Lectura IFC + detección de cambios + columna IFC → solver | ✅ |
| Tests automatizados | ✅ 30 tests |
| Activación temporal de vertidas | ⏳ pendiente |
| Clima real (Open-Meteo, snapshots) | ⏳ pendiente: hoy el ambiente es una sinusoide ilustrativa |
| Solver 3D + voxelización | ⏳ pendiente |
| Monte Carlo / HPC / PyVista | ⏳ pendiente |

## Supuestos vigentes (a validar)

| Supuesto | Valor | Dónde |
|---|---|---|
| División en 3 vertidas de 1,10 m; la tercera es el pedestal Ø7 m | geometría | `ifc_io.py` — **pendiente geometría oficial del curso** |
| Temperatura ambiente sinusoidal | media 16 °C, amplitud 7 °C | `config/*.yaml` — se reemplaza con Open-Meteo |
| Suelo modelado / temperatura | 5 m, 15 °C | `config/*.yaml` |
| Temperatura de colocación | 20 °C | decisión del plan |
| Propiedades térmicas y curva adiabática | Ju y Lei (2019) Tabla 2, ec. (2)-(3) | `config/*.yaml`, `hidratacion.py` |

## Referencias

- Ju, Y. y Lei, H. (2019). *Actual Temperature Evolution of Thick Raft Concrete Foundations and
  Cracking Risk Analysis*. Advances in Materials Science and Engineering, 7029671.
  https://doi.org/10.1155/2019/7029671
- ICH (2019). *Especificación para Hormigón Masivo Estructural*.
- IfcOpenShell: https://docs.ifcopenshell.org/ · Bonsai: https://docs.bonsaibim.org/
