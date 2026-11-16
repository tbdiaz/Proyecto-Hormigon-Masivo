# Enterrar la fundación: relleno de tierra

Código: `simulacion3d.Relleno`, `construir_dominio(..., relleno_z=)` y `simular_3d(..., relleno=)`.
Tests: [`tests/test_relleno.py`](../tests/test_relleno.py).

## 1. Modelo

- **Qué se rellena:** desde el terreno (z = 0) hasta la altura z elegida, todo el aire de la grilla fuera del
  hormigón pasa a ser tierra, en toda la planta del dominio (la fundación más el margen de suelo, 2 m por lado).
  - Con z menor que la cara superior, el relleno entierra la fundación en parte; la consola y el panel dicen qué
    vertidas toca.
  - Con z mayor, la fundación queda enterrada entera. La grilla crece hacia arriba para dejar al menos una capa
    de aire sobre el relleno.
- **Cuándo aparece:** a la hora de colocación. Hasta entonces esas celdas no existen, igual que una vertida que
  todavía no se coloca. Desde ese momento conducen como el suelo (k, ρ, c_p del caso, o los que se indiquen).
- **Bordes:**
  - La cara superior del relleno recibe el mismo clima que el terreno: viento, sol y cielo.
  - Sus caras en el límite de la grilla son adiabáticas, como el suelo lejano.
  - Las caras del hormigón que quedan en contacto con el relleno dejan de ser bordes y conducen hacia la tierra
    con la media armónica de k.
- **Criterio de ΔT:** las caras enterradas siguen siendo la superficie del hormigón. Tmin se mide a 25 mm bajo
  la junta hormigón–relleno, interpolando entre la junta y el centro de la celda. En la junta,
  $T_j = (k_h T_h + k_r T_r)/(k_h + k_r)$, porque el flujo es el mismo a ambos lados. Sin esta regla, una
  fundación enterrada entera "no tendría superficie" y el criterio quedaría vacío. La base sobre el terreno
  natural sigue fuera del criterio, como antes.
- **Valores por defecto:**
  - **Hora:** cuando las vertidas que toca el relleno ya están colocadas y se retiraron sus capas (moldaje o
    aislación) con retiro finito.
  - **Temperatura:** la media del aire en las 24 h anteriores, como tierra acopiada al aire libre.
- **Avisos:** si el relleno se coloca con una capa todavía puesta en una vertida que toca, el modelo pone el
  relleno en contacto directo con el hormigón y la capa deja de contar. La consola y el panel lo avisan.
- **Sin relleno el cálculo no cambia.** Se verificó que una corrida sin relleno da un `serie.csv` idéntico bit a
  bit al de antes, así que las corridas del estudio siguen vigentes.

**Supuestos declarados:**
- **Sin excavación**, el relleno cubre la planta completa de la grilla, como si la tierra de relleno y la de las
  paredes del hoyo fueran la misma.
- **Con excavación** ([docs/excavacion.md](excavacion.md)), el relleno llena solo el hoyo, y el suelo natural de
  alrededor tiene su propia temperatura. Si el relleno sube más que el terreno, cubre también el suelo natural.
- **En todos los casos:** superficie horizontal, colocado de una vez, sin compactación ni humedad (k constante).

## 2. Cómo se usa

**En el editor** (`abrir_editor.bat`), grupo **"Enterrar la fundación"**:
- marcar "Relleno de tierra alrededor";
- la altura en metros sobre el terreno, con la referencia de cada vertida (V1 0–0,55, V2 0,55–1,10, V3 1,10–2,20,
  V4 2,20–3,30, V5 3,30–4,30 m);
- la hora de colocación (vacía = automática), la temperatura (vacía = aire del día) y el k del relleno (vacío =
  el del suelo).

**En el panel:**
- casilla "Relleno" en la vista 3D: línea café a trazos con la altura del relleno desde que se coloca;
- en los cortes, el relleno se ve como suelo;
- la tarjeta "Ambiente, moldaje y aislación" describe el relleno y sus avisos;
- en las curvas, la colocación es una línea discontinua.

**En la línea de comandos:**

```powershell
mcplan simular3d config/caso_14d.yaml --relleno 1.1                     # entierra V1 y V2, hora automática
mcplan simular3d config/caso_14d.yaml --relleno 4.7 --relleno-desde 30  # entera, a las 30 h
mcplan simular3d config/caso_14d.yaml --relleno 4.7 --relleno-T 18 --relleno-k 1.8
```

**En el YAML del caso:**

```yaml
relleno:
  z: 1.1          # m sobre el terreno; más que la cara superior (4,3 m) = enterrada entera
  desde_h: 96     # opcional; omitido = cuando las vertidas que toca están colocadas y sin capas
  T: 15           # opcional; omitido = media del aire en las 24 h anteriores
  k: 1.8          # opcional (también rho, cp); omitido = el suelo del caso
```

## 3. Qué tanto ayuda

Caso de 14 días con clima real, Δx = 0,275 m, 168 h, Ju y Lei, k = 4,11 W/mK en las cinco vertidas, sin tubos.
Las caras laterales llevan moldaje de terciado (retiro a 72 h) y terciado como aislación (retiro a 168 h).

| Caso | ΔT máx 0–96 h [°C] | ΔT máx 97–168 h [°C] |
|---|---|---|
| Sin relleno | 58,1 | 59,0 |
| Relleno hasta 1,1 m (V1 y V2) a las 96 h | 58,1 | 58,7 |
| Enterrada entera a las 96 h | 58,1 | 54,7 |
| **Enterrada entera a las 30 h** | **42,1** | **42,8** |
| Sin moldaje ni aislación | 62,1 | 62,0 |
| Sin moldaje ni aislación, enterrada entera a las 30 h | 43,4 | 43,8 |

Tmax casi no cambia (77,8 °C; 77,9 °C enterrada a las 30 h): el relleno actúa sobre la superficie, no sobre el
núcleo.

- **Lo que más importa es cuándo se entierra.** El ΔT llega a su máximo entre las 90 y las 120 h. Si el relleno
  llega a las 96 h, el pico ya ocurrió. Enterrada a las 30 h, el ΔT baja unos 16 °C.
- **Enterrar en parte ayuda poco.** Con 1,1 m de relleno, la superficie más fría pasa a estar en las caras que
  siguen al aire (troncos y pedestal).
- **La tierra fría también enfría al principio.** Recién colocado, el relleno a ~14 °C saca calor de la
  superficie; un semiespacio de suelo equivale a un h = k/√(π α t) de ~17 W/m²K a la hora y ~3 W/m²K a las 24 h. Por eso
  enterrar se parece a poner una aislación moderada, más que a aislar del todo. Su ventaja es que corta el viento
  y el enfriamiento nocturno hacia el cielo.
- **En obra:** para enterrar a las 30 h hay que desmoldar antes, y luego rellenar y compactar en contacto con
  hormigón joven. Eso hay que justificarlo frente a la resistencia temprana. El modelo avisa cuando el relleno
  llega antes del desmolde.
- **No alcanza solo.** Incluso enterrada temprano, el ΔT sigue cerca del doble del límite de 20 °C, porque el
  núcleo llega a 78 °C. Hay que combinarlo con medidas que bajen Tmax (tubos, menos cemento, T de colocación
  más baja).
