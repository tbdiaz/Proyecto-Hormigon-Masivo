# Ambiente, moldaje y clima real (E10–E11)

Código: [`bordes.py`](../src/mcplan/bordes.py) (condiciones de borde), [`clima.py`](../src/mcplan/clima.py)
(Open-Meteo y snapshots), [`simulacion3d.py`](../src/mcplan/simulacion3d.py) (solver),
[`vtk_io.py`](../src/mcplan/vtk_io.py) (VTK y PyVista). Caso: [`config/caso_14d.yaml`](../config/caso_14d.yaml).
Tests: `test_bordes.py`, `test_clima.py`, `test_vtk.py`.

| Lo que pide E10/E11 | Dónde |
|---|---|
| Módulo de condiciones de borde | `bordes.py`: `Ambiente`, `coeficiente_U`, `filas_bordes` |
| Tratamiento del suelo | `construir_dominio` (3 m de suelo bajo la base y 2 m a los lados) y perfil ERA5 (`clima.perfil_suelo`) |
| Moldaje, aislación y desmolde | `Capa`, `Cubierta`, `PlanCubiertas` en `bordes.py`; secciones `capas`, `moldaje` y `aislacion` del YAML |
| Cliente meteorológico | `clima.descargar_pronostico`, `clima.descargar_era5` (`mcplan clima pronostico`, `mcplan clima era5`) |
| Esquema de snapshots | `data/clima/snapshots/*.json` + `data/clima/indice.csv` (sección 4) |
| Snapshot real con metadata | `data/clima/snapshots/` (`mcplan clima info ARCHIVO`) |
| Serie que alimenta la BC | `clima.ensamblar` (`mcplan clima serie`) |
| Exportación VTK y PyVista | `mcplan simular3d ... --vtk`, `mcplan vtk resultados/.../vtk --recorte y` |

## 1. Cómo entra el ambiente a la física

Cada cara de celda que separa hormigón (o suelo) de aire es un borde. El flujo que entra a la celda es

$$q = U\,(T_{eq} - T_{celda}), \qquad U = \frac{1}{1/h_{ext} + R_{capa} + \Delta x/(2k)}$$

Es una cadena de resistencias en serie: ambiente → capa (moldaje o aislación) → media celda de hormigón.

| Término | Qué es | Dónde entra |
|---|---|---|
| $T_{aire}$ | temperatura a 2 m (clima) | convección y parte de la onda larga, dentro de $T_{eq}$ |
| $h_c(v)$ | convección forzada: $5{,}6 + 4{,}0\,v$ si $v \le 5$ m/s; $7{,}2\,v^{0,78}$ si $v > 5$ (Jürges/McAdams) | en $h_{ext} = h_c + h_r$ |
| $h_r$ | onda larga linealizada: $4\varepsilon\sigma T_m^3$ (de $\varepsilon\sigma(T_s^4 - T_{rad}^4)$) | en $h_{ext}$ y en $T_{eq}$ |
| $T_{rad}$ | lo que "ve" la cara: $F_{cielo} T_{cielo} + (1 - F_{cielo}) T_{aire}$, con $F_{cielo}$ = 1 arriba y 0,5 en los costados | en $T_{eq}$ |
| $T_{cielo}$ | Swinbank, $0{,}0552\,T_a^{1,5}$ [K], corregida por nubosidad $c$: $T^4 = (1-c)T_{desp}^4 + c\,T_a^4$ | en $T_{rad}$ |
| $\alpha G$ | onda corta absorbida; G = radiación global horizontal; sólo en caras hacia arriba | en $T_{eq}$ |
| $R_{capa}$ | $e/k$ del moldaje o de la aislación | en $U$ |

La **temperatura sol-aire** (equivalente) junta todo en un solo número:

$$T_{eq} = \frac{h_c T_{aire} + h_r T_{rad} + \alpha G}{h_c + h_r}$$

Con ella, $h_{ext}(T_{eq} - T_s)$ es exactamente convección + onda larga + sol. Un test lo verifica para varias
$T_s$.

**Dos modelos** (`superficie.modelo` en el YAML):
- **`combinado`** (E08/E09): $h$ fijo de 27,8 W/m²K (Ju y Lei) y $T_{eq} = T_{aire}$. Ese $h$ ya incluye la
  radiación, así que no se suma radiación aparte.
- **`fisico`** (E10 en adelante): cada mecanismo entra una vez. El viento cambia sólo $h_c$, el sol suma sólo
  $\alpha G$ y el cielo sólo cambia $T_{rad}$.

Con el clima de Litueche (27 nov–11 dic), $h_{ext}$ varía entre 18,8 y 35,3 W/m²K, con media cercana a 27.
$T_{eq}$ en la cara superior va de unos 7 °C de noche (cielo despejado) a unos 40 °C a mediodía.

**Suelo.**
- **Geometría:** 3 m de suelo bajo la base y 2 m de margen lateral.
- **Superficie al aire (z = 0, fuera de la planta):** recibe el mismo clima.
- **Fondo:** Dirichlet a la temperatura profunda de ERA5 (capa 100–255 cm, 13,9 °C).
- **Costados:** Neumann $q = 0$ (suelo lejano).
- **Perfil inicial:** la climatología ERA5 de sus 4 capas a la hora de inicio (13,5 / 16,6 / 17,2 / 13,5 °C a
  0,04 / 0,18 / 0,64 / 1,78 m).
- **Contacto hormigón–suelo:** no es un borde: es conducción con la media armónica de $k$.

**Tabla de condiciones de borde** (sale en consola, en `resumen.json` y en `panel.html`):

| Superficie | Tipo | β = h′/k [1/m] | Comentario |
|---|---|---|---|
| Hormigón, caras superiores | Robin, sol-aire | 8,5–15,9 | toda cara expuesta hacia arriba; incluye las escaleras de los troncos |
| Hormigón, laterales con moldaje | Robin + $R$ = 0,15 m²K/W | 2,2–2,5 | V1, V2 y V5 hasta el desmolde |
| Hormigón, laterales desnudas | Robin | 8,5–16,0 | troncos, y V1, V2 y V5 después del desmolde |
| Suelo al aire | Robin, sol-aire | 13,5–25,4 | |
| Fondo del suelo | Dirichlet | ∞ | 13,9 °C |
| Costados del suelo | Neumann | 0 | adiabático |

Aquí $h' = 1/(1/h_{ext} + R_{capa})$ es el coeficiente entre la superficie del hormigón y el ambiente.

## 2. Moldaje, aislación y desmolde

- **Modelo:** una capa es una resistencia $R = e/k$ en serie; su capacidad térmica se desprecia (18 mm de terciado
  frente a 0,11 m de celda).
- **Catálogo** (`capas` en el YAML):
  - terciado 18 mm, k = 0,12 W/mK → R = 0,15 m²K/W;
  - metálico 5 mm (R ≈ 0);
  - manta 5 cm, k = 0,04 → R = 1,25 m²K/W.
- **Qué caras lleva cada una:**
  - moldaje en las caras laterales de **V1, V2 y V5**, que son verticales;
  - V3 y V4 son troncos con la cara inclinada (~22°) libre;
  - la aislación va en las caras superiores.
- **Retiro:** el moldaje se retira `RetiroMoldaje_h` horas (pset del IFC, 72 h) después de colocar cada vertida.
  Cada retiro es un **evento** del solver, igual que una colocación: se recalculan las caras. Con $h_{ext}$ ≈ 27
  W/m²K, $h'$ pasa de unos 5 a unos 27 W/m²K; con la media celda de 0,11 m, $U$ pasa de 4,7 a 16 W/m²K.
- **Pset sin decidir:** el IFC trae `Moldaje` y `Aislacion` = "TBD", así que la decisión queda explícita en el YAML.

| Opción | Comando |
|---|---|
| Plan base (moldaje de terciado, desmolde a 72 h, sin aislación) | `mcplan simular3d config/caso_14d.yaml` |
| Desmoldar antes o después | `--retiro-moldaje 24` / `--retiro-moldaje 168` |
| Sin moldaje | `--sin-moldaje` |
| Manta sobre las caras superiores | `--aislacion manta_5cm [--retiro-aislacion 168]` |

**Al desmoldar**, el costado pasa de $h' \approx 5$ a $\approx 27$ W/m²K: la superficie se enfría de golpe y ΔT
sube, porque el núcleo sigue caliente. Por eso se suele desmoldar tarde o poner una manta al desmoldar.

## 3. ΔT con el criterio del curso (25 mm bajo la superficie)

`criterios.profundidad_superficie_m: 0.025`. Para cada cara expuesta de hormigón se supone un perfil lineal en la
media celda, con el flujo que entra por la cara:

$$T(p) = T_{celda} + \frac{q\,(\Delta x/2 - p)}{k}, \qquad q = U(T_{eq} - T_{celda})$$

Con $p = 0$ se obtiene la superficie del hormigón (bajo el moldaje, si lo hay) y con $p = \Delta x/2$ el centro de
la celda. Un test verifica que superficie < 25 mm < centro.

**Puntos representativos** (`criterios.superficie_tmin: plana`). Tmin es el mínimo sobre las caras de celdas con
**una sola** cara expuesta. En una grilla de vóxeles, las aristas y vértices (el borde del octágono y los cantos de
la escalera de los troncos) se enfrían por dos o tres caras a la vez y exageran el frío de una esquina real. El
enunciado pide puntos "representativos ~25 mm bajo las superficies", así que esas celdas no entran al criterio.
Igual se informan aparte: columna `Tmin_aristas_C` de `serie.csv` y `dTmax_con_aristas_C` de `resumen.json`.

## 4. Pipeline meteorológico y snapshots

**Pronóstico, reanálisis y observación no son lo mismo:**
- **Pronóstico (forecast):** lo que un modelo predice desde su emisión. Cambia con cada emisión.
- **Reanálisis (ERA5):** el pasado reconstruido por un modelo que asimila observaciones; continuo y sin huecos.
- **Observación:** la medición de una estación o un sensor en un punto. Todavía no la usamos: no hay estación en
  el sitio.

**Fuentes (Open-Meteo, punto del curso para Litueche: −34,105457; −71,649539):**
- **Pronóstico:** ECMWF IFS 0,25° (`models=ecmwf_ifs025`), 15 días con dato. El punto de grilla del modelo es
  (−34,0; −71,75), a 337 m.
- **Reanálisis:** ERA5 de la misma ventana del calendario (20 nov–18 dic) entre 2016 y 2025.

**Esquema del snapshot** (`data/clima/snapshots/<tipo>_<modelo>_<AAAAMMDDTHHMMZ>.json`): `metadata` más
`respuesta`, que es la respuesta de la API tal como llegó.

| Campo | Ejemplo |
|---|---|
| `tipo`, `fuente`, `endpoint`, `url`, `parametros` | pronostico, Open-Meteo Forecast API, URL completa de la consulta |
| `modelo` | `ecmwf_ifs025` / `era5` |
| `descargado_utc` | 2026-11-07T16:19:17Z |
| `sitio`, `punto_modelo` | coordenadas pedidas y las del punto de grilla con su elevación |
| `zona_horaria`, `utc_offset_s` | America/Santiago, −10800 |
| `unidades`, `variables` | °C, m/s, W/m², % |
| `horizonte` | desde/hasta y última hora con dato de cada variable |
| `sha256_respuesta` | huella de `respuesta`: si alguien la edita, `leer_snapshot` falla |

**Reglas:**
- **Nunca se reemplaza un snapshot.** El archivo se abre en modo "x" y falla si ya existe.
- **Índice:** cada descarga agrega una línea a `data/clima/indice.csv`.
- **Trazabilidad de la corrida:** la simulación guarda en `resumen.json` qué snapshots usó (nombre, fecha de
  descarga y SHA-256) y cuántas horas vinieron de cada fuente.

**Horizonte no cubierto.** El pronóstico descargado el 7 nov llega al 21 nov. La ventana de la obra
(27 nov–11 dic) queda **entera fuera**. `clima.ensamblar` arma la serie hora a hora en tres tramos:
1. pronóstico donde lo hay;
2. después, la anomalía de su última hora (pronóstico − climatología) decae linealmente a 0 en 24 h;
3. el resto es la **climatología ERA5**: la media, a la misma hora del día, de ±3 días alrededor de cada fecha en
   10 años, es decir 70 muestras por hora.

Así la serie conserva el ciclo diario típico de esos días. Lo que pierde son los eventos (un frente frío, una
racha de viento): la climatología los promedia. Eso se tratará con incertidumbre (ensemble o continuación
estocástica). Cada hora queda marcada con su origen: columna `origen_clima` de `serie.csv` y sombreado de
`clima.png`.

**En la fecha final** se descarga un pronóstico nuevo (`mcplan clima pronostico`); `caso_14d.yaml` usa el más
reciente (`pronostico: ultimo`). Los anteriores se conservan para ver cómo cambió el pronóstico al acercarse la
fecha y si la decisión es robusta a esos cambios.

## 5. VTK y PyVista

`mcplan simular3d config/caso_14d.yaml --vtk` escribe en `vtk/` un `.vti` (VTK ImageData, binario comprimido)
por instante, más `fundacion.pvd`, que es la serie de tiempo para ParaView.
- **Arreglos por celda:** `T`, `te`, `vertida_id`, `material` (0 aire, 1 suelo, 2 hormigón) y `superficie`
  (hormigón con una cara al aire).
- **Para mirar** se usa `mcplan vtk`: slice (`--corte`) o clipping (`--recorte`) por un plano x, y o z. Marca el
  **núcleo** (celda más caliente) y la **superficie** (celda expuesta más fría), e imprime sus coordenadas.

**Driver de video.** En el Intel Iris de este computador, VTK pinta mal los colores por celda (todo oscuro y la
barra de color repetida). `vtk_io` convierte la superficie a colores por punto, sin compartir puntos entre caras,
y dibuja su propia barra de color.

## 6. Guion de la demo E11

| Pedido | Comando |
|---|---|
| Mostrar el archivo meteorológico | `mcplan clima listar` y abrir `data/clima/snapshots/pronostico_ecmwf_ifs025_20261107T1619Z.json` |
| Mostrar la metadata | `mcplan clima info data/clima/snapshots/pronostico_ecmwf_ifs025_20261107T1619Z.json` |
| Graficar la temperatura ambiente | `mcplan clima serie --fecha-inicio 2026-11-08T08:00 --salida resultados/clima_demo` → `clima.png` |
| Alimentar una BC con esa serie | la misma serie entra a `q = U(T_eq − T)`; la tabla [4] de `simular3d` muestra h, β y U |
| Simulación corta | `mcplan simular3d config/caso_14d.yaml --fecha-inicio 2026-11-08T08:00 --dx 0.5 --duracion 72 --vtk --salida resultados/demo11` (~5 s) |
| Abrir resultados 3D | `mcplan vtk resultados/demo11/vtk --interactivo` (o abrir `fundacion.pvd` en ParaView) |
| Slice / clipping | `mcplan vtk resultados/demo11/vtk --t 48 --recorte y --interactivo` · `--corte z --pos 1.0` |
| Núcleo y superficie | los imprime `mcplan vtk` y aparecen marcados en la vista; también `panel.html` |

Si el ayudante pregunta por la ecuación exacta del viento y la radiación: `bordes.py`,
`Ambiente.exterior` ($h_c(v)$, $h_r$, $T_{rad}$ y $\alpha G$ dentro de $T_{eq}$); se usan en
`simulacion3d.simular_3d`, línea `var_v[sup] += etapa.N @ (U * T_eq / dx) - ...`.

## 7. Preguntas de la revisión E10

- **¿Qué término usa la temperatura ambiente?** El flujo de borde $q = U(T_{eq} - T)$. $T_{aire}$ entra en
  $T_{eq}$ con peso $h_c$ (convección) y, a través de $T_{rad}$, en la parte de onda larga que no ve cielo.
- **¿Cómo entra el viento?** Sólo en $h_c(v) = 5{,}6 + 4{,}0\,v$. Más viento da más $h_{ext}$ y más $U$, pero no
  cambia $T_{eq}$. Con moldaje pesa poco, porque domina $R$.
- **¿Cómo se modela la resistencia del moldaje?** $R = e/k$ en serie: $U = 1/(1/h + e/k + \Delta x/2k)$. Con
  terciado de 18 mm, $R$ = 0,15 m²K/W, y $h' = 1/(1/h + R)$ baja de ~27 a ~5 W/m²K.
- **¿Qué cambia al desmoldar?** Desaparece $R$ en esas caras desde ese instante, que es un evento del solver.
  La superficie se enfría rápido y ΔT sube; también pasa a recibir sol si mira hacia arriba.
- **¿Qué es un flujo solar?** Potencia de onda corta por unidad de área, $G$ [W/m²] (hasta ~1000 W/m² a mediodía
  en verano). La superficie absorbe $\alpha G$, con $\alpha$ ≈ 0,65 para hormigón gris.
- **¿Cómo evitar doble-contar radiación y convección?** Cada mecanismo con su término: $h_c$ sólo convección,
  $h_r$ sólo onda larga hacia $T_{rad}$ y $\alpha G$ sólo onda corta. Si se usa un $h$ "combinado" medido (como
  27,8 W/m²K), ese ya incluye la radiación y no se suma $h_r$. Por eso hay dos modelos y nunca se mezclan.
- **¿Forecast, reanálisis y observación?** Ver la sección 4: predicción emitida, pasado reconstruido por modelo
  y medición puntual.
- **¿Qué pasa si el forecast no cubre aún los 14 días posteriores al 27?** Se completa con climatología ERA5 de
  los mismos días de 2016–2025, con una transición de 24 h. Queda declarado hora a hora, y hoy son 337 de 337
  horas. La incertidumbre de ese tramo se tratará con ensemble o continuación estocástica.
- **¿Por qué conservar forecasts antiguos?** Por trazabilidad (qué dato sustentó cada decisión, verificable con
  el SHA-256), para medir el error del pronóstico según su anticipación y para ver si la decisión cambia cuando
  llega información nueva. Un archivo que se reemplaza en silencio borra esa historia.

## 8. Limitaciones declaradas

- El sol sólo llega a caras horizontales: sin sombra ni orientación de las caras verticales.
- El viento a 10 m se usa sin corregirlo a la altura de la superficie.
- Las capas no tienen capacidad térmica.
- La escalera de los troncos sobrestima su área expuesta en ~30 %.
- La climatología no trae eventos extremos.
