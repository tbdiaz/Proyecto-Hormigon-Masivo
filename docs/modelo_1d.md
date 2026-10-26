# Modelo térmico 1D

## 1. Por qué una columna vertical

Ju y Lei (2019, sec. 3.1.3) midieron temperaturas casi iguales en sensores ubicados a la
misma altura: con 3,3 m de espesor y 18 m de lado, el gradiente relevante es el vertical.
Por eso el primer modelo es una columna en el eje de la fundación (suelo + hormigón).
Como no considera pérdidas laterales, es **conservador para Tmax** en el núcleo.

## 2. Ecuación

Conducción con fuente (calor de hidratación):

$$\rho c_p \frac{\partial T}{\partial t} = \frac{\partial}{\partial z}\left(k \frac{\partial T}{\partial z}\right) + Q$$

| Símbolo | Significado | Unidad |
|---|---|---|
| $k$ | conductividad: cuánto calor conduce el material por unidad de gradiente | W/(m·K) |
| $\rho$ | densidad | kg/m³ |
| $c_p$ | calor específico: energía para subir 1 K un kg | J/(kg·K) |
| $\rho c_p$ | capacidad térmica volumétrica | J/(m³·K) |
| $\alpha = k/(\rho c_p)$ | difusividad: qué tan rápido se propaga un cambio de temperatura | m²/s |
| $Q$ | calor generado por la hidratación | W/m³ |

## 3. Discretización: balance de energía por celda

La columna se divide en celdas (por unidad de área en planta). Para cada celda $i$:

$$C_i \frac{dT_i}{dt} = \sum_j G_{ij}(T_j - T_i) + Q_i\,\Delta z_i + q_{borde,i}$$

- $C_i = \rho_i c_{p,i} \Delta z_i$: capacidad térmica de la celda [J/(m²·K)]
- $G_{ij} = \dfrac{1}{\Delta z_i/(2k_i) + \Delta z_j/(2k_j)}$: conductancia entre vecinas [W/(m²·K)].
  Son dos resistencias en serie (media celda de cada lado), así que la interfaz hormigón–suelo
  (o entre vertidas) queda bien representada sin tratamiento especial.

Es la misma ecuación del enunciado (punto 5) con $V_i = \Delta z_i$, y pasar a 3D sólo
agrega vecinos en x e y.

## 4. Condiciones de borde

| Tipo | Qué se impone | Flujo hacia la celda de borde | Uso en el proyecto |
|---|---|---|---|
| Dirichlet | temperatura de la cara $T_b$ | $\dfrac{2k}{\Delta z}(T_b - T_i)$ | suelo profundo |
| Neumann | flujo $q$ [W/m²] | $q$ | simetría o adiabático ($q=0$) |
| Robin | convección $q = h(T_\infty - T_{cara})$ | $\dfrac{1}{1/h + \Delta z/(2k)}(T_\infty - T_i)$ | superficie expuesta al aire |

Una aislación de espesor $e$ y conductividad $k_a$ se suma como resistencia en serie:
$h_{eq} = 1/(1/h + e/k_a)$.

La temperatura de la superficie se reconstruye imponiendo continuidad de flujo, lo que permite
evaluar el punto a 25 mm bajo la superficie aunque quede entre la cara y el primer centro de celda.

## 5. Integración en el tiempo y estabilidad

Euler explícito:

$$T_i^{n+1} = T_i^n + \frac{\Delta t}{C_i}\left[\sum_j G_{ij}(T_j^n - T_i^n) + Q_i^n \Delta z_i + q_{borde,i}^n\right]$$

Agrupando términos, el coeficiente de $T_i^n$ es $1 - \Delta t \sum G_i / C_i$. Si es negativo,
un punto caliente queda más frío que sus vecinos en el paso siguiente: aparecen oscilaciones
que pueden crecer sin límite. Por eso:

$$\Delta t \le \min_i \frac{C_i}{\sum G_i}$$

- En una celda interior con malla uniforme: $\sum G = 2k/\Delta z$, lo que da $\Delta t \le \Delta z^2/(2\alpha)$, es decir **Fo ≤ 1/2**.
- Junto a un borde Dirichlet: $\sum G = 3k/\Delta z$, lo que da $\Delta t \le \Delta z^2/(3\alpha)$ (más restrictivo).
- Consecuencia práctica: **si Δz se reduce a la mitad, Δt debe reducirse a un cuarto** (4× más pasos y 2× más celdas, unas 8× más cómputo).

Ejemplo (hormigón, $\alpha = 9{,}5\times10^{-7}$ m²/s, Δz = 0,05 m): Δt crítico ≈ 1310 s.

Además de la estabilidad hay un límite por **precisión**: con mallas gruesas el paso estable
supera 1 h, demasiado para seguir el calor de hidratación (el núcleo llega a subir ~5 °C/h).
Euler explícito es de primer orden en el tiempo: al dividir Δt por 2, el error se divide por 2.
Estudio de Δt en el caso base (Δz = 0,05 m):

| Δt [s] | 600 | 300 | 150 | 60 | 30 |
|---|---|---|---|---|---|
| Tmax [°C] | 78,04 | 77,95 | 77,91 | 77,83 | 77,85 |
| ΔTmax [°C] | 58,74 | 58,66 | 58,61 | 58,54 | 58,55 |

Por eso los casos usan `dt_max_s: 150` (error ~0,07 °C frente a Δt = 30 s). En 1D cuesta menos de
1 s; las diferencias de ±0,02 °C entre 60 y 30 s vienen de muestrear la salida cada 1 h.

## 6. Calor de hidratación (Ju y Lei 2019, ec. 2–4)

Curva de aumento adiabático en función de la edad equivalente $t_e$ [h]:

$$\theta(t_e) = \begin{cases} 0{,}02\,t_e^{2{,}342} & t_e \le 21\text{ h} \\ 60 - 35\,e^{-0{,}03(t_e - 21)} & t_e \ge 21\text{ h}\end{cases}$$

Edad equivalente (Arrhenius; $E_h/R = 2400$ K, $T_0 = 293{,}15$ K):

$$\frac{dt_e}{dt} = \exp\left[\frac{E_h}{R}\left(\frac{1}{T_0} - \frac{1}{T}\right)\right]$$

Fuente: $Q = \rho c_p \dfrac{d\theta}{dt_e}\dfrac{dt_e}{dt}$. Cada celda de hormigón lleva su propia $t_e$:
las zonas más calientes hidratan más rápido.

Observaciones:
- $\theta \to 60$ °C: un bloque perfectamente aislado colocado a 20 °C llegaría a ~80 °C.
- La curva es continua en 21 h pero su derivada salta de 2,79 a 1,05 °C/h. Esto produce el
  quiebre que se ve en la curva de Tmax cerca de las 20 h.
- El ensayo duró 7 días; después de 168 h la curva es una extrapolación.

## 7. Conversión de unidades del paper (Tabla 2)

| Parámetro | Paper | SI |
|---|---|---|
| $k$ hormigón | 8,0 kJ/(m·h·°C) | 8000/3600 = 2,222 W/(m·K) |
| $c_p$ hormigón | 0,97 kJ/(kg·°C) | 970 J/(kg·K) |
| $\rho$ hormigón | 2400 kg/m³ | 2400 kg/m³ |
| $k$ suelo | 5,0 kJ/(m·h·°C) | 1,389 W/(m·K) |
| $c_p$ suelo | 1,273 kJ/(kg·°C) | 1273 J/(kg·K) |
| $\rho$ suelo | 1800 kg/m³ | 1800 kg/m³ |
| $h$ hormigón–aire | 100 kJ/(m²·h·°C) | 27,8 W/(m²·K) |
| $h$ moldaje metálico–aire | 95 kJ/(m²·h·°C) | 26,4 W/(m²·K) |

## 8. Verificación

**Solución analítica** (`mcplan verificar`): losa de 1 m a 50 °C con ambas caras llevadas a 20 °C.

| Δz [m] | Δt [s] | error máx [°C] | orden |
|---|---|---|---|
| 0,2 | 10477 | 2,7e-1 | |
| 0,1 | 3082 | 5,1e-2 | 2,39 |
| 0,05 | 782 | 1,3e-2 | 1,94 |
| 0,025 | 196 | 3,4e-3 | 1,99 |
| 0,0125 | 49 | 8,4e-4 | 2,00 |

Orden 2: al dividir Δz por 2 el error se divide por 4 (con Δt ∝ Δz², el error temporal
O(Δt) también es O(Δz²)).

**Batería de tests** (`pytest`), que cubre en 1D los tests que recomienda la guía G05:

| Guía G05 | Test |
|---|---|
| 1. Solución exacta 1D con convergencia | serie de Fourier, orden ≈ 2 |
| 4. Energía adiabática constante | `sum(C_i T_i)` constante a precisión de máquina |
| 5. Fuente uniforme: dT/dt = q/(ρc) | aumento lineal exacto |
| 6. Neumann | flujo impuesto entra completo; q = 0 no pierde calor |
| 7. Robin | estado estacionario exacto con resistencias en serie |

Además: estacionario exacto con dos materiales (interfaz), principio del máximo, rechazo y
divergencia con Δt inestable, bloque adiabático que reproduce la curva del ensayo, lectura del IFC
(5 vertidas, volúmenes de la guía), detección de cambios, rayo vertical sobre los troncos y
**puente IFC → solver**: la columna leída del IFC da el mismo resultado que el caso escrito a mano.
Los tests 2D y 3D exactos vendrán con el solver 3D.

## 9. Resultado de la columna del eje (sin control térmico)

`mcplan simular1d config/caso_base_1d.yaml`: geometría oficial en el eje (4,30 m de hormigón: 3,30
de fundación + 1,00 de pedestal), colocación a 18 °C a las 8:00 h, todas las vertidas desde t = 0 y
ambiente ilustrativo de 16 ± 7 °C.

- **Tmax ≈ 77,9 °C** a las 81 h, a 2,12 m sobre la base (centro del espesor). **No cumple** 70 °C.
- **ΔTmax ≈ 58,6 °C** a las 141 h. **No cumple** 20 °C.
- No depende de la malla: entre Δz = 0,2 y 0,025 m, Tmax varía 0,02 °C.
- `caso_ifc_1d.yaml` (la misma columna leída del IFC) da exactamente lo mismo.

Columna cerca del borde (x = 8,5 m: sólo 1,30 m de hormigón): Tmax ≈ 65,5 °C (cumple) pero
ΔTmax ≈ 40,9 °C (no cumple). Una sección delgada se enfría más rápido, pero igual tiene gradiente.

Es coherente con Ju y Lei: su Fundación 1, sin control, también excedió ambos límites
(73 °C y 35 °C). Las diferencias son esperables, porque no se conocen su temperatura de
colocación ni su clima, el modelo 1D no tiene pérdidas laterales y la vertida no fue
instantánea. **No es una validación**: es un orden de magnitud razonable.

## 10. Limitaciones actuales

- Todas las vertidas existen desde t = 0 (falta la activación según `Inicio_h`).
- Ambiente sinusoidal ilustrativo (sin datos reales todavía).
- Sin radiación solar ni pérdidas laterales: la columna 1D ignora el enfriamiento por las caras
  inclinadas y laterales del octágono.
- Una sola columna por simulación; el campo completo requiere la voxelización 3D.
