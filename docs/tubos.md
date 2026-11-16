# Tubos de enfriamiento y hormigón distinto por vertida

Código: [`tubos.py`](../src/mcplan/tubos.py) y `simulacion3d.simular_3d(..., tubos=, materiales_vertida=)`.
Tests: [`tests/test_tubos.py`](../tests/test_tubos.py).

## 1. Modelo

Un tubo de 25–32 mm es mucho más delgado que una celda (0,11–0,275 m). Convertir celdas enteras en "tubo" quitaría
hormigón y exageraría el enfriamiento. Por eso cada celda por la que pasa el tubo recibe un **sumidero de calor**
hacia el agua.

**Resistencias por metro de tubo, en serie [m·K/W]:**

| Resistencia | Fórmula | Depende de |
|---|---|---|
| Agua (convección interior) | $1/(h_w \pi D_i)$, con Nu = 0,023 Re^0,8 Pr^0,4 (turbulento) o 4,36 (laminar) | caudal, diámetro |
| Pared del tubo | $\ln(D_o/D_i)/(2\pi k_{tubo})$ | **material**: PVC 0,19; HDPE 0,45; acero 45; cobre 385 W/mK |
| Hormigón alrededor | $\ln(r_{eq}/r_o)/(2\pi k_c)$, con $r_{eq} = 0{,}198\,\Delta x$ (radio de Peaceman) | **k del hormigón** de esa vertida, malla |

$U_L = 1/(R_{agua} + R_{pared} + R_{hormigón})$ [W/(m K)]. Para Ø32 mm, 20 L/min y Δx = 0,275 m:

- **PVC:** $U_L$ ≈ 4,9 W/mK. Su pared (0,11 m·K/W) opone más resistencia que el hormigón que la rodea (0,09).
- **HDPE:** $U_L$ ≈ 7,2 W/mK.
- **Cobre:** $U_L$ ≈ 10,8 W/mK.

**El agua se calienta a lo largo del recorrido.** En el tramo de largo $L_i$ dentro de la celda $i$:

$$T_{w,i+1} = T_{c,i} + (T_{w,i} - T_{c,i})\,e^{-U_L L_i/(\dot m c_w)}$$

La celda entrega $\dot m c_w (T_{w,i+1} - T_{w,i})$. El paso de tiempo crítico incluye este intercambio, y el balance
de energía lo cuenta: lo que pierde el hormigón es exactamente lo que gana el agua (test).

**Supuestos declarados:**
- El tubo sólo intercambia calor con celdas ya colocadas.
- Caudal y temperatura de entrada son constantes mientras el circuito opera.
- No se modela la inercia térmica del agua dentro del tubo.
- No se calcula la tensión alrededor del tubo. Se informa el **ΔT hormigón–agua máximo** porque en la práctica se
  limita (del orden de 20 °C) para no fisurar alrededor del tubo; conviene verificar el valor en ACI 207.4R.

**Hormigón distinto por vertida.** Cada celda toma k, ρ y c_p de su vertida. En la junta entre vertidas se usa la
media armónica de k, y las caras expuestas usan el k de su celda. Un test lo compara con `solver_nd` usando dos k.

## 2. Cómo se edita

**En el editor** (`abrir_editor.bat`):
- **Tabla de vertidas:** columna **k [W/mK]** de cada vertida.
- **Grupo "Hidratación":** la curva de Ju y Lei o una de las 8 curvas medidas del profesor, con su aumento
  adiabático equivalente.
- **Grupo "Tubos de enfriamiento":** "+ Agregar circuito" crea un serpentín. En cada uno se eligen la vertida, la
  altura (vacía = media altura; se muestra el rango válido), el material, la dirección de las filas, la separación,
  el diámetro, el caudal, la T del agua y los días de operación. Puede haber varios por vertida.
- **En el panel:** los tubos se ven en celeste en la vista 3D (la casilla "Tubos" los muestra u oculta) y en los
  cortes, y la tarjeta "Tubos de enfriamiento"
  grafica el calor extraído, la temperatura del agua a la salida y el ΔT hormigón–agua contra 20 °C. La tabla
  marca ⚠ cuando ese ΔT se supera.

**En la línea de comandos:**

```powershell
# un circuito por opción --tubo (se puede repetir), con la curva medida M4_22A y k = 4 en la Vertida_03
mcplan simular3d config/caso_14d.yaml --tubo vertida=3 z=1.4 material=cobre direccion=x separacion=1 diametro=32 `
       caudal=20 T_agua=10 dias=7 --tubo vertida=4 material=hdpe --k-vertida 3=4.0 --curva M4_22A
# atajo: un serpentín a media altura de cada vertida indicada, todos iguales
mcplan simular3d config/caso_14d.yaml --tubos 2 3 4 --tubos-material cobre --tubos-T-agua 10
```

**En un YAML** (caso o plan del estudio), con control total del recorrido:

```yaml
materiales:
  hormigon:
    ...
    por_vertida:                  # sólo en el caso; en un plan del estudio va como `materiales: {3: {k: 4.0}}`
      Vertida_03: {k: 4.0}        # hormigón con árido recubierto de grafito [Chen et al. 2024: 4,11 W/mK, sección 3]

tubos:
  - nombre: serpentin_V3
    vertida: 3                    # parte cuando empieza esta vertida
    serpentin: {z: 1.65, separacion: 1.0, margen: 0.5, direccion: x}   # z omitido = media altura
    material: cobre               # pvc, hdpe, acero, cobre, o un número (k de la pared en W/mK)
    diametro_mm: 32
    espesor_mm: 2
    caudal_l_min: 20
    T_agua_C: 10
    desde_h: 0
    durante_h: 168
  - nombre: anillo_pedestal       # o un recorrido libre: lista de puntos (x, y, z) en metros
    vertida: 5
    puntos: [[-2.5, -2.5, 3.8], [2.5, -2.5, 3.8], [2.5, 2.5, 3.8], [-2.5, 2.5, 3.8], [-2.5, -2.4, 3.8]]
```

El recorrido se valida: si un tramo sale del hormigón, la corrida se detiene con un error.

**Dónde se ven los resultados:**
- la consola y `resumen.json`: potencia máxima, energía extraída, T de salida y ΔT hormigón–agua de cada circuito;
- `tubos.png`;
- `panel.html`: los tubos en celeste sobre la vista 3D y en los cortes;
- en el estudio, las columnas `E_tubos_MJ` y `dT_agua_max` de `escenarios.csv`.

## 3. ¿Tubos de PVC o de cobre? ¿Hormigón con grafito?

- **Tubos:** son la medida estándar de enfriamiento en hormigón masivo, y en obra se usan sobre todo HDPE y acero.
  - El cobre conduce mejor que el plástico, pero es caro.
  - Con plástico, la pared del tubo pesa en la resistencia total; el modelo lo cuantifica.
  - Lo que más importa es la **separación entre tubos** (cuánto hormigón enfría cada metro) y que el agua no esté
    demasiado fría (el ΔT hormigón–agua).
- **Grafito:** subir k ayuda a que el calor llegue a los tubos y reduce el gradiente interno.

### Respaldo: Chen et al. (2024)

Recubren árido grueso reciclado con una mezcla 1:1 de grafito y cemento (su "método B", 80 kg de grafito por m³).

| Dato | Valor | Fuente en el paper |
|---|---|---|
| k del hormigón sin grafito (RB) | 2,16 ± 0,08 W/mK | Tabla 4, 28 días |
| k con grafito suelto en el mortero (RG, 96 kg/m³) | 3,14 ± 0,04 W/mK | Tabla 4 |
| k con árido recubierto (RA11B, 80 kg/m³) | **4,11 ± 0,04 W/mK (1,9 veces)** | Tabla 4, conclusión 2 |
| Resistencia a compresión / tracción por hendimiento | −1,7 % / −2,3 % | Tabla 5 |
| Costo | ≈ +10 % por m³ (precios de China) | sección 4.3 |

- **Qué respalda:** aplicada a nuestro k base (2,222 W/mK [JL19]), la razón 1,9 da ≈ 4,2 W/mK. El k = 4,0 de los
  planes "TG" queda un poco por debajo, del lado conservador.
- **Resistencia:** el método importa. Con el árido recubierto la resistencia casi no cambia; con grafito suelto en
  la mezcla, otros autores que ellos citan informan caídas de hasta 90 %.

**Reparos que siguen en pie:**
- **k en seco:** midieron cubos de 100 mm secados a 60 °C, a 25 °C. El hormigón en obra está húmedo, así que la
  razón 1,9 se traslada mejor que el valor absoluto. Además k sube un poco con la edad: a los 3 días es algo menor
  que a los 28.
- **Mezcla de laboratorio:** razón a/c 0,6, 500 kg/m³ de cemento y árido de 5–16 mm. No es una mezcla de hormigón
  masivo. El árido se recubre y se cura 7 días antes de usarlo, lo que exige planificar la logística para todo el
  volumen de la fundación.
- **Calor de hidratación:** no lo midieron. Citan que el grafito puede retardar la hidratación, y el modelo usa la
  misma curva en todas las vertidas.
- **Durabilidad:** no la ensayaron. La esperan aceptable porque el grafito es inerte, pero advierten que el grafito
  suelto aumenta la porosidad abierta.
- **ρ y c_p:** no los informan; se mantienen los de [JL19].
- **Uso estructural:** en una fundación de aerogenerador no es habitual. Con este paper el argumento de resistencia
  queda respaldado en laboratorio; falta la durabilidad.

Chen, C., Liu, X., Liu, Z., Shao, L., Chang, H., Liu, Q., Miao, C. y Feng, P. (2024). *Developing heat-conductive
concrete with graphite-modified recycled aggregates*. Composites Part B 284, 111721.
https://doi.org/10.1016/j.compositesb.2024.111721
