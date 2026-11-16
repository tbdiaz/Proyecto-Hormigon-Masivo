# Fundación en un hoyo: menos viento y menos cielo bajo el terreno

Código: `simulacion3d.Excavacion`, `construir_dominio(..., excavacion=)`; el borde usa los factores en
`bordes.Ambiente.exterior(..., f_viento, f_cielo, f_sol)`. Tests: [`tests/test_excavacion.py`](../tests/test_excavacion.py).

## 1. Qué cambia

En obra, la fundación se construye dentro de una excavación. Hasta que se rellena, las caras bajo el terreno
natural están en un hoyo: les llega menos viento y ven menos cielo. Sin esta opción, el modelo las trata como si
estuvieran al aire libre con el viento a 10 m, lo que enfría de más.

Con la opción activa:

- **Geometría:**
  - El terreno natural queda a la cota indicada sobre la base (z = 0).
  - El hoyo tiene paredes verticales a una distancia fija del borde de la base (el sobreancho).
  - Fuera del hoyo hay suelo natural desde z = 0 hasta el terreno. Existe desde el inicio, con la temperatura
    inicial de ERA5 según su profundidad bajo el terreno.
  - El margen de la grilla crece si hace falta para dejar al menos 1 m de suelo natural alrededor.
- **Viento**, sobre el de 10 m de Open-Meteo:
  - $v = v_{10}\, f_{altura}\, f_{abrigo}$.
  - $f_{abrigo}$ actúa solo bajo el terreno natural; sobre él vale 1.
  - Entra en $h_c = 5{,}6 + 4v$ (Jürges/McAdams), como antes.
- **Cielo y sol bajo el terreno:** cada cara ve solo la abertura del hoyo. Sobre el terreno no cambian.
- **Relleno:** el relleno ([docs/relleno.md](relleno.md)) llena el hoyo; las caras que tapa dejan de ser borde.
- **Modelo de superficie:** el viento, el cielo y el sol solo cambian con el modelo "fisico". Con h fijo, la
  excavación solo agrega el suelo natural alrededor.

## 2. De dónde salen los números

| Factor | Valor por defecto | Fuente | Confianza |
|---|---|---|---|
| $f_{altura}$: viento a la altura de la cara | 0,70 bajo 2 m (categoría II) | EN 1991-1-4, perfil logarítmico, ec. 4.4 y tabla 4.1 | Alta: norma |
| $f_{abrigo}$: viento en el hoyo / viento a nivel del terreno | 0,67 | Nakamura y Oke (1988), medición en una calle real de alto ≈ ancho | Media: otra geometría |
| Fracción de cielo de cada cara | geometría | cuerdas cruzadas de Hottel, hoyo largo en 2D | Alta: geometría |

**Perfil del viento (EN 1991-1-4).**

$$\frac{v(h)}{v(10)} = \frac{\ln(\max(h, z_{min})/z_0)}{\ln(10/z_0)}$$

con $h$ sobre el terreno natural. Los valores por categoría de terreno salen de la tabla 4.1:

| Categoría | Terreno | z0 [m] | z_min [m] | v/v10 bajo z_min |
|---|---|---|---|---|
| 0 | mar o costa expuesta | 0,003 | 1 | 0,72 |
| I | lagos, sin obstáculos | 0,01 | 1 | 0,67 |
| II | campo abierto | 0,05 | 2 | 0,70 |
| III | suburbano o bosque | 0,3 | 5 | 0,80 |
| IV | urbano | 1,0 | 10 | 1,00 |

El dato de Open-Meteo es el viento a 10 m sobre terreno abierto. Bajo $z_{min}$ la norma no extrapola el
logaritmo y usa el valor de $z_{min}$. Sin excavación el modelo usa $v_{10}$ en todas las caras.

**Abrigo del hoyo.**
- **El dato usado:** Nakamura y Oke (1988) midieron en una calle de Kioto con alto/ancho ≈ 1. Encontraron que
  el viento sobre el fondo de la calle es, aproximadamente, 2/3 del viento sobre el techo, con relación lineal
  hasta 5 m/s.
- **Por qué es conservador aquí:** no hay mediciones para un hoyo de fundación. Con 3,3 m de profundidad y 1 m
  de sobreancho, el hoyo es más angosto (alto/ancho ≈ 3), y en calles más angostas el viento al fondo es menor.
  En minas a rajo abierto se informa ~1/4 en el fondo. Por eso 0,67 deja más viento del que probablemente hay,
  y se muestra la sensibilidad.

**Cielo visible.** Para una franja en un hoyo largo de ancho G (distancia a lo que la tapa) y a la profundidad
d bajo el terreno:
- **cara vertical:** $F = \tfrac12\left(1 - d/\sqrt{d^2+G^2}\right)$. En el borde (d = 0) vale 0,5, igual que
  antes;
- **cara horizontal con un lado abierto**, como los escalones de los troncos:
  $F = \tfrac12\left(1 + G/\sqrt{G^2+d^2}\right)$;
- **fondo del hoyo**, tapado por la pared y por la fundación: la suma de los dos lados con la misma fórmula.

El resto de lo que ve la cara (la pared, la fundación) se toma a la temperatura del aire, como ya se hacía con
la parte de las caras laterales que no ve cielo. El sol directo sobre caras horizontales en el hoyo se reduce en
la misma fracción; es una aproximación que no sigue la posición del sol.

## 3. Cómo se usa

**En el editor**, grupo **"Excavación"**: marcar "La fundación está en un hoyo" y fijar:
- la cota del terreno natural (sugiere 3,30 m: V1 a V4 bajo tierra y el pedestal asomando 1 m);
- el sobreancho;
- el factor de abrigo;
- la categoría de terreno;
- si se corrige el viento a la altura de la obra.

**En el panel:**
- casilla "Terreno": línea verde a trazos en la vista 3D;
- en los cortes, las paredes del hoyo se ven como suelo;
- la tabla de bordes agrega "caras laterales en el hoyo".

**En la línea de comandos:**

```powershell
mcplan simular3d config/caso_14d.yaml --excavacion 3.3                       # sobreancho 1 m, abrigo 0,67, cat. II
mcplan simular3d config/caso_14d.yaml --excavacion 3.3 --sobreancho 1.5 --factor-abrigo 0.25 --categoria-terreno 0
mcplan simular3d config/caso_14d.yaml --excavacion 3.3 --relleno 3.3          # y rellenar el hoyo al desmoldar
```

**En el YAML del caso:**

```yaml
excavacion:
  profundidad: 3.3        # m: terreno natural sobre la base
  sobreancho: 1.0         # m: del borde de la base a la pared del hoyo
  factor_abrigo: 0.67     # viento en el hoyo / viento a nivel del terreno [Nakamura y Oke 1988]
  categoria_terreno: II   # EN 1991-1-4
  perfil_viento: true
```

Sin excavación el cálculo no cambia: una corrida sin ella da un `serie.csv` idéntico bit a bit al de antes.

## 4. Qué tanto ayuda

Caso de 14 días con clima real, Δx = 0,275 m, 168 h, Ju y Lei, k = 4,11 W/mK, sin tubos. Las caras laterales
llevan moldaje de terciado (retiro a 72 h) y terciado como aislación (retiro a 168 h). La excavación tiene 3,3 m
y 1 m de sobreancho.

| Caso | ΔT máx [°C] | Tmin en el ΔT máx |
|---|---|---|
| Sin excavación | 59,0 | 17,5 °C, cara lateral de V2 |
| Excavación solo con la fracción de cielo | 58,9 | 17,5 °C, ídem |
| + viento a la altura de la obra (× 0,70) | 58,1 | 18,6 °C, ídem |
| **+ abrigo 0,67 (por defecto)** | **56,9** | 19,9 °C, ídem |
| + abrigo 0,25 | 55,6 | 19,6 °C, **cara superior del pedestal** |
| abrigo 0,67 y relleno hasta el terreno a las 96 h | 55,6 | 19,6 °C, cara superior del pedestal |
| abrigo 0,67 y relleno hasta el terreno a las 30 h | 55,5 | 19,7 °C, cara superior del pedestal |

Tmax: 77,8–78,1 °C en todos.

- **Pesa más el viento que el cielo.** Con el terciado puesto, la cara lateral pierde poco por radiación, así que
  la fracción de cielo casi no cambia nada.
- **Mejora unos 2 °C con el abrigo por defecto y hasta ~3,5 °C con más abrigo.** Es poco porque las caras
  laterales ya estaban protegidas por el moldaje y el terciado: la capa (R = 0,15 m²K/W cada una) domina sobre
  el viento.
- **La cara crítica se mueve.** Con más abrigo, o con el hoyo rellenado, la superficie más fría pasa a ser la
  cara superior del pedestal, que queda sobre el terreno y al aire. Desde ahí, abrigar más abajo ya no mejora el
  ΔT: hay que proteger el pedestal (aislación superior) o bajar Tmax.

**Fuentes:**
- CEN (2005). EN 1991-1-4, *Eurocode 1: Actions on structures — Part 1-4: Wind actions*, sección 4.3.2 y
  tabla 4.1.
- Nakamura, Y. y Oke, T. R. (1988). *Wind, temperature and stability conditions in an east-west oriented urban
  canyon*. Atmospheric Environment 22(12), 2691–2700. El factor 2/3 se tomó de citas secundarias; conviene
  verificarlo en el original.
- Hottel, H. C. (método de las cuerdas cruzadas para factores de forma 2D); ver cualquier texto de
  transferencia de calor por radiación, p. ej. Modest, *Radiative Heat Transfer*, cap. 4.
