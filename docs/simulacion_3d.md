# Pipeline IFC → grilla → simulación 3D por vertidas (E08)

Código: [`src/mcplan/simulacion3d.py`](../src/mcplan/simulacion3d.py) (solver por etapas),
[`voxelizacion.py`](../src/mcplan/voxelizacion.py) (IFC → vóxeles), [`vista3d.py`](../src/mcplan/vista3d.py) (vista 3D).
Caso: [`config/caso_3d.yaml`](../config/caso_3d.yaml). Tests: [`tests/test_simulacion3d.py`](../tests/test_simulacion3d.py).

```powershell
mcplan simular3d                      # caso del YAML: dx = 0,275 m, 14 días
```

## 1. Pipeline

| Paso | Qué hace | Dónde |
|---|---|---|
| 1. IFC | IfcOpenShell lee las 5 `Vertida_XX`: malla triangulada en coordenadas globales y pset `UANDES_HormigonMasivo` | `ifc_io.leer_vertidas` |
| 2. Voxelización | celda de lado dx; pertenece a la vertida que contiene su centro (paridad de un rayo vertical); `vertida_id` 0…5 | `voxelizacion.voxelizar` |
| 3. Error de volumen | (celdas · dx³) contra el volumen exacto de cada malla | `Voxelizacion.error_volumen` |
| 4. Dominio | se agrega suelo bajo z = 0 (3 m, 2 m de margen) y aire alrededor | `simulacion3d.construir_dominio` |
| 5. Secuencia | `Inicio_h` y `TemperaturaColocacion_C` del IFC, con los cambios del YAML y de la línea de comandos | `aplicar_cambios`, `PlanVertidas` |
| 6. Solver | volúmenes finitos + Euler explícito, por etapas entre colocaciones, con hidratación Ju y Lei | `simulacion3d.simular_3d` |
| 7. Salidas | Tmax(t), ΔT(t), T en puntos, cortes, plantas, vista 3D (PNG y HTML interactivo), `serie.csv`, `resumen.json` | `graficos.py`, `vista3d.py` |

## 2. La secuencia está en los inputs, no en el solver

`simular_3d` no conoce horas ni temperaturas de colocación: recibe un `PlanVertidas` y, en cada instante, le pregunta
qué vertidas existen. Las horas 0/6/12/18/24 y los 18 °C salen del pset del IFC. Hay tres niveles de input, y cada uno
manda sobre el anterior:

1. **IFC** (Bonsai): `Inicio_h` y `TemperaturaColocacion_C` de cada vertida.
2. **YAML**, sección `vertidas`: `Vertida_02: {Inicio_h: 12}`, `Vertida_05: {activa: false}`.
3. **Línea de comandos**: `--inicio 2=12`, `--T-colocacion 3=10`, `--ocultar 5`, `--mostrar 5`.

La salida imprime la secuencia usada y cada cambio, con su valor anterior y el nuevo. `test_la_secuencia_esta_en_los_inputs`
cambia el inicio y la temperatura de V2 y oculta V5 sin tocar el solver.

## 3. Modelo de cada etapa

- **Celdas que existen:** el suelo, desde t = 0; cada vertida, desde su `Inicio_h`. Al colocarse, una vertida parte
  con T = T de colocación y edad equivalente te = 0.
- **Cara entre dos celdas que existen:** conducción con la media armónica de k. Así funcionan la junta entre vertidas
  y el contacto hormigón–suelo.
- **Cara entre una celda que existe y una que no** (aire, o una vertida aún no colocada): Robin con el aire,
  $h_{ef} = 1/(1/h + \Delta x/2k)$, con $h_{sup}$ en ±z y $h_{lat}$ en ±x, ±y. Cuando llega la vertida siguiente, la
  cara superior de la anterior deja de ser Robin y pasa a ser junta.
- **Suelo:** T_fondo (Dirichlet) en la cara inferior; las caras laterales son adiabáticas (suelo lejano).
- **Hidratación:** $Q = \rho c_p\, d\theta/dt_e \cdot dt_e/dt$ en cada celda de hormigón que existe, con su propia
  $t_e$ (Arrhenius).
- **Paso de tiempo:** se recalcula en cada etapa como $\Delta t \le \min(0{,}9\,\Delta t_{crít},\ \Delta t_{máx})$, con
  $\Delta t_{crít} = \min_i C_i/\sum G_i$. Con dx = 0,11 m resulta $\Delta t_{crít}$ ≈ 2100 s; se usa
  $\Delta t_{máx}$ = 600 s por precisión (Tmax difiere unos 0,2 °C de lo que da 150 s).

**Verificación** (`tests/test_simulacion3d.py`):

| Test | Resultado |
|---|---|
| Caja con Robin en todas las caras + hidratación, contra `solver_nd` (ya verificado) | igual (< 1e-4 °C, campos en float32) |
| Dos vertidas: V2 se coloca a las 5 h a 25 °C, contra `solver_nd` en dos tramos | igual: la cara superior de V1 es Robin y luego junta |
| Bloque adiabático | T = T₀ + θ(t_e) del ensayo (< 0,1 °C) |
| Dominio | suelo bajo z = 0, mismo número de celdas de hormigón que la voxelización |
| Inputs | cambiar inicio/T y ocultar vertidas cambia el resultado sin tocar el solver |

## 4. Guion de la demo (lo que pide el ayudante)

| Pedido | Comando o acción | Qué mirar |
|---|---|---|
| Ejecutar el caso | `mcplan simular3d` | salida [1]…[5] y `resultados/caso_3d/` |
| **Cambiar Δx** | `mcplan simular3d --dx 0.5 --duracion 72` (rápido) · `--dx 0.11` (fino) | tabla [2], celdas, Δt crítico |
| **Reportar error de volumen** | sale siempre en [2]; para varios dx: `--estudio-volumen 0.5 0.275 0.2 0.11 0.1` | por vertida y total |
| **Ocultar / mostrar una vertida** | `--ocultar 5`; o `Vertida_05: {activa: false}` en el YAML y `--mostrar 5` para volver | [3] dice "OCULTA"; en `fundacion_3d.html`, clic en la leyenda |
| **Cambiar hora de inicio** | `--inicio 2=12` (o `todas=0`: vertida monolítica) | [3], líneas punteadas de `tmax.png` |
| **Cambiar T de colocación** | `--T-colocacion todas=10` o `3=25` | `puntos.png`: la curva parte en esa T |
| **Visualizar el efecto** | `--comparar-con resultados/caso_3d/serie.csv --salida resultados/demo` | `tmax.png` superpone la corrida base (línea discontinua) |
| Cambio en el BIM | editar el pset en Bonsai, guardar el IFC y volver a ejecutar | [3] muestra el valor nuevo |

Tiempos en este computador (14 días): dx = 0,5 m ≈ 5 s; dx = 0,275 m ≈ 30 s (demo); dx = 0,11 m ≈ 7 min (línea base E09).

**Ejemplo de efecto:** colocar todo a las 0 h (`--inicio todas=0`) o colocar más frío (`--T-colocacion todas=10`)
cambia Tmax y el instante en que ocurre; `tmax.png` lo muestra contra la base. Los números de la línea base están en
[baseline_3d.md](baseline_3d.md).

## 5. Limitaciones declaradas

- **Ambiente simplificado:** sinusoide 16 ± 7 °C, h constante, sin sol, sin viento variable, sin moldaje, aislación ni
  retiro de moldaje.
- **Suelo:** uniforme a 15 °C al inicio.
- **Superficies en escalera:** en los troncos el área expuesta queda sobreestimada en ~30 % (ver
  [modelo_vertidas.md](modelo_vertidas.md)); se podría ponderar cada cara por el área real.
- **ΔT:** se mide como Tmax − Tmin de las celdas de superficie, cuyo centro queda a dx/2 de la cara; no exactamente
  a 25 mm, como en el caso 1D.
