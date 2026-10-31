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

## 7. Ejercicio de la pizarra: rectángulo que se enfría por el borde

> Las figuras y animaciones de los ejercicios (secciones 7 y 8) no se guardan en el repositorio:
> se regeneran con `mcplan rectangulo` y `mcplan bloque` en `resultados/`.

**Enunciado.** Rectángulo de 4 m × 2 m que resuelve $\rho c_p\,\partial T/\partial t - k\nabla^2 T = 0$:

- inicio a $T(x,y,0) = 20$ °C;
- todo el borde $\Gamma$ a 10 °C para todo $t$ (Dirichlet, "condición de borde esencial");
- celdas cuadradas, $\Delta x = \Delta y$.

Se pide el $\Delta t$ crítico, simular hasta estabilizar, la animación de la evolución y $T(t)$ en dos puntos.

```powershell
mcplan rectangulo                                        # hormigón de Ju y Lei, puntos (0, 0) y (1, 0.5)
mcplan rectangulo --alpha 1e-4                           # difusividad de la guía (horas en vez de días)
mcplan rectangulo --punto 0 0 --punto 1 1 --dx 0.05      # otros puntos y otra malla
mcplan rectangulo --origen esquina --punto 2 1 --punto 3 1.5
```

Genera los siguientes archivos en `resultados/rectangulo_2d/`:

- `evolucion.gif`: la animación;
- `puntos.png`: T(t) numérica contra la exacta;
- `campos.png`: cuatro instantes del campo;
- `serie_puntos.csv` y `resumen.json`.

**Coordenadas.** Como en la pizarra, (0, 0) es el centro del rectángulo. Con esa convención, **(1, 1) cae en el
borde superior** y vale 10 °C siempre. Por eso el segundo punto por defecto es (1, 0.5). El comando avisa
cuando un punto está en el borde.

**Solución exacta.** Por separación de variables, es el producto de dos losas 1D (la serie de `mcplan verificar`):

$$\frac{T - T_b}{T_i - T_b} = X(x,t)\,Y(y,t), \qquad X = \sum_{m\ \text{impar}} \frac{4}{m\pi}\sin\frac{m\pi x}{L_x}\,e^{-\alpha (m\pi/L_x)^2 t}$$

A tiempos largos manda el primer modo, que vale $\tfrac{16}{\pi^2}(T_i - T_b)$ en el centro y decae con
$e^{-t/\tau}$, donde $\tau = 1/(\alpha\lambda_1)$ y $\lambda_1 = (\pi/L_x)^2 + (\pi/L_y)^2$. De ahí sale el tiempo
para quedar a menos de 0,1 °C de $T_b$ en todo el rectángulo:

$$t_{est} = \tau \ln\frac{(16/\pi^2)\,|T_i - T_b|}{0{,}1} \approx 5{,}1\,\tau$$

**Resultados** con $\Delta x = 0{,}1$ m (40 × 20 celdas):

| | Hormigón ($\alpha = 9{,}5\times10^{-7}$ m²/s) | Guía ($\alpha = 10^{-4}$ m²/s) |
|---|---|---|
| $\Delta t$ crítico, celda interior $\Delta x^2/4\alpha$ | 2619 s | 25,0 s |
| $\Delta t$ crítico, junto al borde $\Delta x^2/6\alpha$ (**manda**) | 1746 s | 16,7 s |
| $\Delta t$ usado (se guarda al menos cada 1 h) | 1198 s | 14,7 s |
| Constante de tiempo $\tau$ | 3,93 días | 0,90 h |
| Estabilización estimada / simulada | 20,0 / 20,0 días | 4,58 / 4,58 h |
| Error máx. en los puntos vs. exacta | 0,03 °C | 0,03 °C |

- El centro se queda unas horas en 20 °C, porque el frío todavía no llega desde los bordes, y después cae
  hacia 10 °C sin alcanzarlo nunca. El punto (1, 0.5) está más cerca del borde y se enfría antes: recorre
  la mitad del salto en 2,4 días, contra 4,4 días del centro.
- El error del campo es máximo en la primera hora (0,65 °C en las celdas junto al borde, por el salto de 20 a
  10 °C en $t = 0$) y al final baja a 4e-5 °C.
- En tiempo adimensional ($t/\tau$), las curvas de los dos materiales son iguales: sólo cambia la escala de tiempo.

### 7.1 El mismo rectángulo con calor de hidratación

```powershell
mcplan rectangulo --hidratacion          # resultados/rectangulo_2d_hidratacion/
```

Ahora la ecuación es $\rho c_p\,\partial T/\partial t - k\nabla^2 T = Q$, con el calor de Ju y Lei (2019). Es el
mismo modelo del solver 1D (`hidratacion.py`):

$$Q = \rho c_p \frac{d\theta}{dt_e}\frac{dt_e}{dt}, \qquad \frac{dt_e}{dt} = \exp\left[\frac{E_h}{R}\left(\frac{1}{T_{ref}} - \frac{1}{T}\right)\right]$$

- **Edad equivalente por celda.** Cada celda (vóxel) lleva su propia edad equivalente $t_e$: el centro, más
  caliente, envejece más rápido que las celdas junto al borde.
- **Esquema explícito, igual que en 1D.** Por paso: $Q$ se evalúa con $T^n$ y $t_e^n$, luego
  $t_e^{n+1} = t_e^n + \Delta t \cdot dt_e/dt$, y por último $T^{n+1}$.
- **Paso de tiempo.** La fuente no cambia el $\Delta t$ crítico, que sigue siendo 1746 s. Pero se limita a
  **150 s por precisión**, por el estudio de $\Delta t$ del modelo 1D: aquí el centro llega a subir 4 °C/h
  (a las 16 h).
- **Tiempo para estabilizar.** Ya no hay fórmula cerrada. El comando simula por bloques de 30 días hasta que
  todo el rectángulo queda a ±0,1 °C, y después simula 1,2 veces ese tiempo.

**Verificación:**

- `test_hidratacion_2d_igual_al_solver_1d`: una losa de 2 m entre dos caras a 10 °C, simulada en 2D con las
  caras laterales adiabáticas, da **lo mismo que el solver 1D** (diferencia de 1e-14 °C).
- `test_bloque_adiabatico_2d_sigue_la_curva_del_ensayo`: sin pérdidas, cada celda reproduce el ensayo
  adiabático $T = T_0 + \theta(t_e)$.
- `test_hidratacion_solo_en_el_hormigon`: la edad sólo avanza en las celdas de hormigón (máscara `es_hormigon`).

**Resultados** (hormigón, $\Delta x = 0{,}1$ m, $\Delta t = 150$ s):

| | Valor |
|---|---|
| Tmax | **71,9 °C** a las 46 h, en el centro. Supera 70 °C entre las 38 y las 57 h: **no cumple** |
| Máximo en (1, 0.5) | 58,6 °C a las 39 h: cerca del borde se calienta menos y antes |
| ΔTmax (Tmax − Tmin a 25 mm del borde) | **61,6 °C** a las 46 h. Sobre 20 °C entre las 14 y las 184 h: **no cumple** |
| A los 14 días | Tmax 14,1 °C, ΔT 4,1 °C |
| Estabilización (±0,1 °C) | 28,7 días, contra 20,0 sin hidratación: el calor generado retrasa el enfriamiento casi 9 días |

- Con 2 m de espesor el calor no alcanza a salir mientras se genera: el centro sube 52 °C en 2 días.
  En el ensayo adiabático subiría hasta 60 °C.
- El borde fijo a 10 °C es una condición extrema: la superficie no se calienta nada, por eso ΔT es tan
  grande. En la fundación real la superficie es Robin (convección), se calienta y ΔT baja. Es el mismo
  mensaje de la sección 3: **la condición de borde domina el ΔT**.

### 7.2 El mismo rectángulo con el pulso de calor de la clase

```powershell
mcplan rectangulo --hidratacion pulso                     # Qmax = 1800 W/m³, T = 42 h
mcplan rectangulo --hidratacion pulso --Qmax 2500 --T-calor 36
```

Es la fórmula de la pizarra:

$$Q(t) = \frac{Q_{max}}{2}\left(1 - \cos\frac{2\pi t}{T}\right) \quad (t < T), \qquad Q(t) = 0 \quad (t \ge T)$$

- Parte en 0, llega a $Q_{max}$ en $t = T/2$ y vuelve a 0 en $t = T$.
- Libera en total $Q_{max} T / 2$ [J/m³], es decir, una elevación adiabática de $Q_{max} T / (2\rho c_p)$.
- Sólo depende del tiempo: no hay Arrhenius ni edad equivalente. El código lo implementa como
  `CalorPulsoCoseno` en `hidratacion.py`, con la misma interfaz que Ju y Lei y $dt_e/dt = 1$.

**Valores por defecto.** La pizarra no da $Q_{max}$ ni $T$, así que se tomaron de la curva de Ju y Lei:

- su ritmo máximo es 2,79 °C/h a las 21 h, lo que da $Q_{max} = \rho c_p \cdot 2{,}79/3600 \approx 1800$ W/m³
  y $T = 2 \cdot 21 = 42$ h;
- la elevación adiabática total queda en 58,5 °C, cerca de los 60 °C del paper.

Con `--Qmax` y `--T-calor` se cambian.

**Solución exacta.** Como $Q$ no depende de $T$, el problema sigue siendo lineal: la solución es la del
rectángulo sin calor más la respuesta a la fuente. La fuente uniforme se descompone en los modos
$\sin(m\pi x/L_x)\sin(n\pi y/L_y)$, con $m, n$ impares y coeficientes $16/(mn\pi^2)$. Cada modo cumple

$$\frac{dA}{dt} = -a A + \frac{16}{mn\pi^2}\frac{Q(t)}{\rho c_p}, \qquad a = \alpha\pi^2\left(\frac{m^2}{L_x^2} + \frac{n^2}{L_y^2}\right)$$

y la integral $\int_0^t e^{-a(t-s)} Q(s)\,ds$ con el pulso coseno tiene forma cerrada. Se usan 200 × 200 modos,
con error de truncamiento menor que 1e-3 °C.

**Verificación:**

- **Bloque adiabático:** sube exactamente $Q_{max}T/(2\rho c_p)$ (error < 1e-6 °C).
- **Contra la exacta:** el error en los puntos es de 0,13 °C en el centro y 0,18 °C en (1, 0.5), con peaks de 76 °C.
- **Convergencia en malla:** a las 42 h, el error máximo del campo baja 0,37 → 0,10 → 0,026 °C con
  $\Delta x$ = 0,1 → 0,05 → 0,025 m, es decir, **orden 2**. Se concentra cerca de las esquinas; en el centro
  es de 0,03 °C. No cambia al bajar $\Delta t$ de 150 a 30 s: es error de malla.

**Resultados** (hormigón, $\Delta x = 0{,}1$ m, $\Delta t = 150$ s):

| | Pulso de la clase | Ju y Lei (sección 7.1) |
|---|---|---|
| Tmax en el centro | **76,2 °C a las 38 h** | 71,9 °C a las 46 h |
| Máximo en (1, 0.5) | 66,6 °C a las 35 h | 58,6 °C a las 39 h |
| ΔTmax (25 mm del borde) | 65,9 °C | 61,6 °C |
| Estabilización (±0,1 °C) | 28,4 días | 28,7 días |

El pulso calienta más y antes que Ju y Lei por dos razones:

- **Todo el calor sale en 42 h.** En Ju y Lei la curva tiene una cola que dura días.
- **El ritmo no depende de la temperatura.** En Ju y Lei, las celdas frías junto al borde hidratan más lento
  (Arrhenius).

Las dos fuentes liberan parecida energía total, así que el enfriamiento final y la estabilización casi coinciden.

## 8. Ejercicio de la pizarra: bloque 3D con el pulso y cuatro condiciones de borde

**Enunciado.** Bloque de **3 m × 3 m × 0,6 m** de hormigón que resuelve $\rho c_p\,\partial T/\partial t - k\nabla^2 T = Q(t)$:

- inicio a $T = 20$ °C;
- $Q(t)$ es el pulso de la sección 7.2: $Q_{max} = 1800$ W/m³ y $T = 42$ h;
- la misma condición de borde en las seis caras, en cuatro casos:

| Caso | Condición | Significado |
|---|---|---|
| Dirichlet 20 °C | $T = 20$ °C | la cara se mantiene a 20 °C |
| Dirichlet 30 °C | $T = 30$ °C | la cara se mantiene a 30 °C: el bloque además se calienta desde afuera |
| Neumann | $\partial T/\partial n = 0$ | adiabático: no sale calor (el ensayo adiabático) |
| Robin | $\partial T/\partial n = -\beta\,(T - 20)$ | convección con aire a 20 °C |

**Qué es β.** No tiene que ver con Euler explícito, que es *cómo se avanza en el tiempo*. β es el
coeficiente de la condición Robin: $-k\,\partial T/\partial n = h\,(T - T_{aire})$, de donde
$\beta = h/k$ [1/m]. Mide qué tan fuerte la superficie entrega calor al aire comparado con lo que conduce el
hormigón.

- **Valor usado.** La pizarra no lo da, así que se usa $h = 27{,}8$ W/(m²K), el de hormigón–aire de Ju y Lei:
  $\beta = 27{,}8/2{,}222 = 12{,}5$ 1/m. Se cambia con `--beta`.
- **Signo.** Con la normal hacia afuera, el signo menos hace que salga calor cuando $T > T_{aire}$. El test
  `test_exacta_robin_cumple_la_condicion_de_borde` lo verifica en la solución exacta.

```powershell
mcplan bloque                                      # los 4 casos (~2 min): resultados/bloque_3d/
mcplan bloque --beta 5 --casos robin               # otro β
mcplan bloque --T-dirichlet 15 25 --casos dirichlet
mcplan bloque --dx 0.1 --sin-animacion             # más rápido (6 celdas en el espesor)
```

Genera, por caso, `puntos_<caso>.png` (T(t) numérica y exacta, Tmax, Tmin a 25 mm, ΔT), `cortes_<caso>.png`
(planos horizontal y vertical centrales) y `serie_<caso>.csv`. Para los cuatro juntos genera `comparacion.png`,
`evolucion.gif` y `resumen.json`.

### 8.1 Solución exacta de los cuatro casos

Como $Q$ sólo depende de $t$, el problema es lineal. Con el origen en el centro, cada eje (medio lado $l$)
tiene modos simétricos $\cos(\mu s/l)$, donde $\mu \tan\mu = Bi = hl/k$:

| Borde | Bi | Autovalores μ |
|---|---|---|
| Dirichlet | ∞ | $(2k-1)\pi/2$ |
| Neumann | 0 | sólo $\mu = 0$: el campo queda uniforme |
| Robin | $hl/k$ (3,75 en z; 18,8 en x e y) | raíces de $\mu\tan\mu = Bi$, una en cada $[(k-1)\pi,\,(k-1)\pi + \pi/2]$ (bisección) |

La temperatura inicial y la fuente son uniformes, y el 1 se descompone como
$1 = \sum_k b_k \cos(\mu_k s/l)$, con $b_k = 2\sin\mu_k / (\mu_k + \sin\mu_k\cos\mu_k)$. En 3D cada modo
$(i, j, m)$ decae con $a = \alpha(\mu_i^2/l_x^2 + \mu_j^2/l_y^2 + \mu_m^2/l_z^2)$:

$$T = T_{ref} + \sum_{i,j,m} b_i b_j b_m \cos\frac{\mu_i x}{l_x}\cos\frac{\mu_j y}{l_y}\cos\frac{\mu_m z}{l_z}\left[(T_i - T_{ref})\,e^{-at} + \frac{G(a,t)}{\rho c_p}\right]$$

donde $T_{ref}$ es la temperatura de la cara o del aire, y $G(a,t) = \int_0^t e^{-a(t-s)}Q(s)\,ds$ es cerrada para
el pulso. Se usan 100 × 100 × 60 modos, con error de truncamiento menor que 2e-3 °C. Robin con $h \to \infty$
da Dirichlet, y Neumann da $T = 20 + \int Q\,dt/(\rho c_p)$: los tres bordes salen de la misma fórmula.

### 8.2 Paso crítico según el borde

En la celda de esquina, cada eje aporta 1 vecina más lo que aporte la cara (en unidades de $\alpha/\Delta x^2$):

- Dirichlet: 2 (la cara a media celda);
- Neumann: 0;
- Robin: $h_{ef}\Delta x/k$, con $h_{ef} = 1/(1/h + \Delta x/2k)$.

Con $\Delta x = 0{,}05$ m (60 × 60 × 12 = 43.200 celdas):

| Caso | Esquina | Interior | Manda | Δt crítico |
|---|---|---|---|---|
| Dirichlet | 3·(1 + 2) = 9 | 6 | la esquina | $\Delta x^2/9\alpha$ = **291 s** |
| Neumann | 3·(1 + 0) = 3 | 6 | el interior | $\Delta x^2/6\alpha$ = 436 s |
| Robin (β = 12,5) | 3·(1 + 0,48) = 4,4 | 6 | el interior | $\Delta x^2/6\alpha$ = 436 s |

En Robin, el Biot de una celda ($h\Delta x/k = 0{,}63$) es chico: la convección saca menos calor que una
vecina, así que el borde no restringe. Igual se usa $\Delta t$ = 150 s, porque **la precisión manda**
(hay que seguir el pulso). Son 3720 pasos por caso.

### 8.3 Resultados

Con $\Delta x = 0{,}05$ m, simulando 154 h (1,2 veces la estabilización más lenta, la de Robin):

| Caso | Tmax | ΔTmax (25 mm) | Equilibrio | Estabiliza (±0,1 °C), exacta / simulada | Error vs. exacta |
|---|---|---|---|---|---|
| Dirichlet 20 °C | 46,9 °C a las 29 h | 26,5 °C | 20 °C | 92,7 / 92,9 h | 0,23 °C |
| Dirichlet 30 °C | 56,0 °C a las 29 h | 25,7 °C | 30 °C | 92,5 / 92,7 h | 0,23 °C |
| Neumann | **78,5 °C** al terminar el pulso | 0 | 78,45 °C (se queda) | 39,3 / 39,4 h | 0,06 °C |
| Robin β = 12,5 | 53,9 °C a las 31 h | **31,1 °C** | 20 °C | 128,5 / 128,8 h | 0,15 °C |

En los cuatro casos Tmax está en el centro (0, 0, 0). El error es el máximo en todo el campo y en todo el
tiempo; al final baja a 1e-5 °C.

**Lectura física:**

- **Mismo bloque y mismo calor: el borde decide todo.** Tmax va de 47 a 78 °C según la condición de borde.
- **Neumann es el ensayo adiabático.** Nada sale, todo el bloque sube $Q_{max}T/(2\rho c_p) = 58{,}45$ °C de
  forma uniforme y se queda ahí. Es la cota superior y **no cumple 70 °C**, pero ΔT = 0 porque no hay gradiente.
- **Dirichlet 30 °C es Dirichlet 20 °C más 10 °C que entran desde la cara (superposición).** En el peak la
  diferencia es 9,2 °C, porque a las 29 h esos 10 °C todavía no llegan del todo al centro.
- **Robin queda entre Dirichlet y Neumann:**
  - la superficie se calienta (deja de estar fija), así que sale menos calor y el centro llega más alto que
    con Dirichlet 20 °C;
  - la superficie se calienta poco, así que **ΔT es el mayor** de los cuatro casos;
  - tarda más en estabilizarse, porque el aire saca el calor más lento que una cara fija.
- **El espesor de 0,6 m manda sobre los 3 m de planta.** La constante de tiempo de Dirichlet es
  $\tau \approx 10$ h y el gradiente es casi todo vertical, como se ve en los cortes.

**Error de Euler en el tiempo (Neumann).** El campo es uniforme, así que no hay error de malla: los
0,06 °C son el atraso de medio paso del esquema explícito, $(\Delta t/2)\,Q/(\rho c_p) = 75 \cdot 7{,}7\times10^{-4}$.
Al terminar el pulso desaparece (1e-8 °C), porque la suma de Riemann del pulso completo es exacta.

**Convergencia en malla** (error máximo desde las 6 h, hasta las 36 h):

| Δx [m] | celdas en el espesor | Dirichlet 20: centro / campo | Robin: centro / campo |
|---|---|---|---|
| 0,2 | 3 | 2,48 / 2,71 °C | 1,42 / 1,76 °C |
| 0,1 | 6 | 0,28 / 0,82 °C | 0,37 / 0,51 °C |
| 0,05 | 12 | 0,075 / 0,23 °C | 0,098 / 0,14 °C |
| orden (0,1 → 0,05) | | 1,9 / 1,9 | 1,9 / 1,8 |

Con 3 celdas en el espesor la malla todavía no está en el régimen asintótico. De 6 a 12 celdas, el error se
divide por ~4: **orden 2**.

**Tests** (`tests/test_bloque.py`, 19):

- autovalores de los tres bordes y descomposición del 1;
- Robin cumple $\partial T/\partial n = -\beta(T - 20)$;
- Dirichlet vale $T_b$ en la cara;
- Robin con $h \to \infty$ da Dirichlet;
- Neumann es la elevación adiabática;
- flujo en la cara Robin e interpolación trilineal;
- simetría;
- Δt crítico por borde;
- estabilización simulada igual a la exacta;
- convergencia de orden 2;
- el comando completo.

Los lentos van marcados: `pytest -m "not lento"` corre la batería rápida.

## 9. Limitaciones actuales

- Dominio rectangular (rectángulo, bloque). La fundación octagonal requiere voxelizar el IFC: una máscara de
  celdas activas por vertida.
- Bordes constantes en el tiempo: falta el clima horario, que ya existe en el solver 1D. La hidratación ya
  está en 2D/3D (sección 7.1).
- Malla uniforme en cada eje.
