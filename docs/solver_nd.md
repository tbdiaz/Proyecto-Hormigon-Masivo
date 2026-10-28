# Solver térmico 2D/3D y batería de verificación

Código: [`src/mcplan/solver_nd.py`](../src/mcplan/solver_nd.py) (solver) y
[`src/mcplan/verificacion_nd.py`](../src/mcplan/verificacion_nd.py) (casos exactos).
Es la misma formulación del solver 1D ([modelo_1d.md](modelo_1d.md)), con más vecinos por celda.

## 1. Balance de un vóxel

Cada celda (vóxel) de lados $\Delta x, \Delta y, \Delta z$ tiene volumen $V$ y su temperatura en el centro:

$$C_i \frac{dT_i}{dt} = \sum_{j \in \text{vecinas}} G_{ij}(T_j - T_i) + \dot q_i V + \sum_{\text{caras de borde}} Q_{cara}$$

| Término | Significado | Unidad |
|---|---|---|
| $C_i = \rho c_p V$ | capacidad térmica de la celda | J/K |
| $G_{ij} = k_{eq} A / \Delta x$ | conductancia con la vecina ($A$ = área de la cara común) | W/K |
| $\dot q_i$ | fuente volumétrica (hidratación) | W/m³ |
| $Q_{cara}$ | calor que entra por una cara exterior (Dirichlet, Neumann o Robin) | W |

- **Signo:** $G(T_j - T_i) > 0$ cuando la vecina está más caliente, así que el calor va de lo caliente a lo frío.
- **Conservación:** lo que sale de $i$ hacia $j$ entra igual a $j$. Las dos se actualizan con el mismo flujo
  (`variacion[izq] += flujo`, `variacion[der] -= flujo`), así que sin fuente ni bordes la energía
  $E = \sum \rho c_p V T$ no cambia.
- **Malla uniforme:** al dividir por $V$, la conductancia por eje queda $k_{eq}/\Delta x^2$, que es el `coef` del código.
- **2D y 3D** usan el mismo balance: 4 vecinas en 2D y 6 en 3D.

## 2. Dos materiales

Entre dos celdas hay dos medias celdas en serie, $R = \Delta x/(2k_1) + \Delta x/(2k_2)$. Por eso:

$$k_{eq} = \frac{2 k_1 k_2}{k_1 + k_2} \quad\text{(media armónica)}$$

El material que conduce menos manda: con $k_1 = 2$ y $k_2 = 0{,}5$, $k_{eq} = 0{,}8$, y no 1,25 como daría la
media aritmética. El test `test_dos_materiales_en_serie_estacionario_exacto` reproduce exacto el perfil
lineal por tramos y falla si se usa la media aritmética. La conservación de energía, en cambio, se
cumple con cualquier media: no basta para verificar la interfaz.

## 3. Condiciones de borde

| Tipo | Dato | Aporte a $dT_i/dt$ de la celda de borde | Uso |
|---|---|---|---|
| Dirichlet | $T_s$ | $\dfrac{2k}{\rho c_p \Delta x^2}(T_s - T_i)$ | la cara está a media celda del centro; suelo profundo |
| Neumann | $q''$ [W/m²], positivo entra | $\dfrac{q''}{\rho c_p \Delta x}$ | $q'' = 0$: **adiabática** (simetría, aislación perfecta) |
| Robin | $h$, $T_\infty$ | $\dfrac{h_{ef}}{\rho c_p \Delta x}(T_\infty - T_i)$, con $h_{ef} = \dfrac{1}{1/h + \Delta x/2k}$ | superficie al aire, moldaje, manta |

- **Robin** es convección más media celda de conducción en serie. Una manta agrega $e/k_a$ a la serie.
- **Adiabática:** en el solver toda cara sin condición es adiabática. En la fundación, un plano de simetría
  permite simular medio o un cuarto del dominio.

**Fuente volumétrica y flujo superficial:**

| | Fuente $\dot q$ | Flujo $q''$ |
|---|---|---|
| Unidad | W/m³ | W/m² |
| Actúa en | todo el volumen | sólo las celdas de la cara |
| Efecto en la celda | $dT/dt = \dot q/(\rho c_p)$, no depende de $\Delta x$ | $dT/dt = q''/(\rho c_p \Delta x)$; la energía total $q'' A$ no depende de $\Delta x$ |
| En el proyecto | calor de hidratación | radiación solar, convección (Robin) |

**Por qué los bordes pueden dominar el resultado:** el calor de hidratación es el mismo en cualquier plan de
vertidas. Lo que cambia entre alternativas es cuánto calor sale por las superficies: $h$, $T_\infty$, manta y
moldaje. Además, el criterio ΔT se mide a 25 mm bajo la superficie, y la temperatura de ese punto la
fija casi por completo la condición de borde.

## 4. Estabilidad (Euler explícito)

$$T_i^{n+1} = T_i^n\left(1 - \Delta t \frac{\sum G_i}{C_i}\right) + \frac{\Delta t}{C_i}\left(\sum G_{ij} T_j^n + \dots\right)$$

Si el coeficiente de $T_i^n$ es negativo, una celda caliente pasa a quedar más fría que sus vecinas y aparecen
oscilaciones que crecen. La condición es $\Delta t \le \min_i C_i/\sum G_i$:

| Celda | $\sum G / C$ | $\Delta t$ crítico (malla cúbica) | Función |
|---|---|---|---|
| Interior 2D | $4\alpha/\Delta x^2$ | $\Delta x^2/(4\alpha)$ | `paso_critico` |
| Interior 3D | $6\alpha/\Delta x^2$ | $\Delta x^2/(6\alpha)$ | `paso_critico` |
| Esquina 2D con Dirichlet | $(3+3)\alpha/\Delta x^2$ | $\Delta x^2/(6\alpha)$ | `paso_critico_malla` |
| Esquina 3D con Dirichlet | $(3+3+3)\alpha/\Delta x^2$ | $\Delta x^2/(9\alpha)$ | `paso_critico_malla` |

- Junto a una cara Dirichlet la celda suma 1 vecina más 2 (la cara a media celda), es decir 3 en vez de 2:
  el límite baja a 2/3.
- `paso_critico_malla` calcula el límite real (materiales, bordes) y `simular_hasta` elige $\Delta t$ solo.
- El criterio es **suficiente**, no exacto. En el caso 2D de la guía con 20 celdas, el crítico es 4,17 s
  (6,25 s en el interior). Con 5 s el resultado no diverge, y con 9 s sí.
- **Costo en 3D:** si $\Delta x$ se reduce a la mitad, hay 8× más celdas y 4× más pasos, unas 32× más cómputo.

**Método de Schmidt** (lectura del Anexo 4 del ICH): en 1D, con $Fo = \alpha\Delta t/\Delta x^2 = 1/2$ exacto,
el esquema explícito queda

$$T_i^{n+1} = \frac{T_{i-1}^n + T_{i+1}^n}{2}$$

La temperatura nueva es el promedio de las vecinas, lo que permite construirla a mano o gráficamente. Es
justo el límite de estabilidad de la sección anterior. En 2D el equivalente es $Fo = 1/4$: el promedio de
las 4 vecinas.

## 5. Batería de verificación (E05)

| E05 | Tests | Comando en vivo | Resultado |
|---|---|---|---|
| Solución 1D exacta | `test_solver1d.py` (serie de Fourier) | `mcplan verificar` | orden 2 |
| Solución 2D sinusoidal | `test_sinusoidal_exacta`, `test_sinusoidal_converge_con_orden_2`, `test_solucion_exacta_sinusoidal_2d` | `mcplan sinusoidal --dim 2` | error 8,1e-3 °C con 30×30, orden 1,99 |
| Solución 3D sinusoidal | los mismos, en 3D | `mcplan sinusoidal --dim 3` | error 9,6e-3 °C con 32³, orden 1,98 |
| Energía adiabática | `test_energia_constante_en_dominio_adiabatico` (2D/3D, uno y dos materiales), `test_conservacion_energia_2d` | `mcplan energia --dim 3 --dos-materiales` | error relativo ~1e-16 |
| Fuente uniforme | `test_fuente_uniforme_adiabatica_exacta`, `test_fuente_uniforme_adiabatica_3d` | `mcplan fuente --dim 3` | 153,333333 °C exacto, balance q0·V·t ~1e-16 |
| Neumann | `test_neumann_2d` (el flujo entra completo), `test_neumann_cero_es_adiabatico` | | exacto |
| Robin | `test_robin_estacionario_exacto` (2D/3D), `test_robin_el_calor_va_de_caliente_a_frio`, `test_robin_3d` | | error < 1e-8 °C |
| Dos materiales | `test_dos_materiales_en_serie_estacionario_exacto`, `test_interfaz_dos_materiales_conserva_energia_2d` | | error < 1e-8 °C |
| Estabilidad | `test_paso_critico_con_y_sin_dirichlet`, `test_dt_inestable_se_rechaza_y_forzado_diverge` | `mcplan sinusoidal --dim 2 --n 20 --dt 9 [--forzar]` | se rechaza o diverge |

Parámetros de la guía: $L = 1$ m, $\alpha = 10^{-4}$ m²/s, $A = 50$ °C, $T_0 = 20$ °C, $t = 1200$ s en 2D y 600 s en 3D.
La fuente usa $\dot q = 10\,000$ W/m³, $\rho = 2400$, $c_p = 900$ y 8 h.

**¿La batería detecta errores?** Se introdujeron a propósito cuatro errores típicos en una copia del solver,
y cada uno hace fallar al menos un test:

- Dirichlet a una celda entera en vez de media: falla
- Robin sin la media celda de conducción: falla
- Media aritmética en la interfaz: falla
- Signo de Neumann invertido: falla

## 6. Cambiar parámetros en vivo

```powershell
pytest tests/test_verificacion_nd.py -v                      # toda la batería 2D/3D
pytest -k "sinusoidal and not converge" -v                   # sólo los casos sinusoidales
mcplan sinusoidal --dim 3 --L 1 2 1.5 --alpha 2e-4 --A 30 --T0 15 --t 900 --n 8 16
mcplan fuente --dim 2 --q0 2500 --rho 2300 --cp 1000 --T0 10 --t 259200
mcplan energia --dim 2 --t 86400 --n 20
mcplan sinusoidal --dim 2 --n 20 --dt 9            # se rechaza: dt > dt crítico (4,17 s)
mcplan sinusoidal --dim 2 --n 20 --dt 9 --forzar   # se ve la divergencia
```

En los tests, cada caso es una fila de `parametrize` (por ejemplo, `CASOS_2D` y `CASOS_3D`): se agrega o
edita una fila y $\Delta t$ se recalcula solo. Desde Python:

```python
from mcplan.solver_nd import simular_hasta
T, dt = simular_hasta(T_inicial, alpha, (dx, dy, dz), t_final=3600.0, fuente=q,
                      bordes={"z_max": ("robin", (27.8, 15.0)), "z_min": ("dirichlet", 15.0)})
```

## 7. Limitaciones actuales

- Dominio rectangular. La fundación octagonal requiere voxelizar el IFC: una máscara de celdas activas por vertida.
- Bordes y fuente constantes en el tiempo. Faltan el clima horario y la hidratación ($\dot q$ depende de $T$
  y de la edad equivalente), que ya existen en el solver 1D.
- Malla uniforme en cada eje.
