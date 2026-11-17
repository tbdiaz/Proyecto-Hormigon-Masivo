# Artículo F01

`Articulo_F01_Grupo4.docx` (Word) y su PDF se arman **desde los resultados** con `scripts/articulo_f01.py`. Cada
número, tabla y figura sale de las corridas. Si cambia un dato, se vuelven a correr los estudios afectados y el
script, sin editar el texto a mano.

## Cómo regenerarlo

Desde la raíz del repositorio, con el entorno instalado (`pip install -r requirements.txt`, que incluye
`python-docx`):

```powershell
# 1. Clima (una vez): 30 temporadas ERA5 alrededor del 19 nov y el ensemble más reciente
mcplan clima era5 --fecha-inicio 2026-11-19T08:00 --anios 1996 2025 --dias 40 --margen 7   # sólo cuentan el día y la hora
mcplan clima ensemble

# 2. Estudios (cada uno retoma lo que falte si se corta)
mcplan estudio config/estudio_era5.yaml            # familia de planes (3 terciles climáticos)
mcplan estudio config/retiro_tubos_era5.yaml       # retiro con tubos extendidos: primera ronda
mcplan estudio config/retiro_tubos2_era5.yaml      # retiro con tubos extendidos: segunda ronda
mcplan estudio config/confiabilidad_era5.yaml      # propuesta final: 8 curvas x 30 temporadas
mcplan estudio config/estres_era5.yaml             # propuesta final con +20 % de calor
mcplan estudio config/sensibilidad_era5.yaml       # rho cp de las probetas, colocación a 8 y 12 °C, R de la manta
mcplan estudio config/verificacion_ensemble.yaml   # propuesta final con el clima del ensemble

# 3. Malla fina. Los dos escenarios que gobiernan cada fase, con tres mallas, paso de tiempo y profundidad:
P=unica-M10-T10-N25-R28
mcplan estudio config/confiabilidad_era5.yaml --procesos 1 --sensibilidad $P M4_SS ERA5_2010 --dx-variantes 0.275 0.1833 0.11
mcplan estudio config/confiabilidad_era5.yaml --procesos 1 --sensibilidad $P M4_SS ERA5_2013 --dx-variantes 0.275 0.1833 0.11
#    ... y los 10 escenarios siguientes más cercanos al límite, sólo con 0,11 m (uno por línea; se pueden correr en paralelo):
for e in "M4_SS ERA5_2022" "M4_SS ERA5_2023" "M2_SS ERA5_2010" "M2_SS ERA5_2022" "M3_SS ERA5_2010" \
         "M2_SS ERA5_2023" "M3_SS ERA5_2022" "M4_SS ERA5_2024" "M4_SS ERA5_2020" "M3_SS ERA5_2023"; do
  mcplan estudio config/confiabilidad_era5.yaml --procesos 1 --sensibilidad $P $e \
    --dx-variantes 0.275 0.11 --dt-variantes 600 --profundidades 0.025
done
#    ... y la corrida de detalle del peor escenario (cortes y VTK)
mcplan estudio config/confiabilidad_era5.yaml --detalle $P M4_SS ERA5_2010

# 4. El artículo (con --pdf, Word lo exporta a PDF, cuenta las páginas y el script informa las palabras).
#    --profesor pone su nombre en la portada (no está en el repositorio).
python scripts/articulo_f01.py --pdf --profesor "Nombre Apellido"
```

En el paso 3, M4_SS con la temporada 2010 es la fila con mayor `M_D` de `resultados/confiabilidad_era5/escenarios.csv`
(su máximo llega después del retiro, a las 718 h), y M4_SS con la temporada 2013 es el mayor `M_D` antes del retiro
(pico de las 142 h). Los otros 10 son los siguientes por `M_D` con la malla del estudio: los que quedan a menos de
0,4 °C del límite con la corrección uniforme. Si cambian los datos, se eligen de nuevo con ese mismo criterio. Cada
corrida de malla fina toma ≈35 min en un proceso; no hace falta el cluster.

## Decisiones de la pauta que el artículo respeta

- **Fecha:** t0 = jueves 19 nov 2026, 08:00 (America/Santiago), común a todos los planes. El clima de los escenarios
  son las 30 temporadas ERA5 1996–2025 alineadas a esa fecha (bloques históricos comparables); el ensemble más
  reciente es la verificación de los primeros días.
- **Retiro dentro de la evaluación:** `despues_ultimo_retiro_h: 168` en `config/estudio.yaml`. Cada escenario se
  simula hasta 7 días después del último retiro de moldaje o aislación, y M_Δ cubre todo ese horizonte.
- **Una sola probabilidad:** P_f es la enumeración ponderada de los escenarios. No hay intervalos binomiales, ni
  remuestreo de las curvas, ni distribuciones ajustadas. La energía crítica y la prueba de +20 % son sensibilidades
  deterministas.
- **Corrección de malla:** en los escenarios corridos con malla fina se usa directamente su M_Δ; a los demás se les
  suma δ, la mayor diferencia entre Δx = 0,11 y 0,275 m medida en esos escenarios. La tabla δ → P_f de la sección 5
  muestra si el error de malla puede cambiar la clasificación y da la condición de validez (δ_crít).
- **Sensibilidad "solo SS":** tiene fila propia en la tabla de sensibilidades y un párrafo en la sección 7.4.

## Qué completar a mano

- El nombre del profesor, con `--profesor` en el paso 4: no está en el repositorio y el script no lo inventa. El curso
  (Métodos Computacionales en Obras Civiles, Universidad de los Andes) ya va en la portada.
- El commit: la primera página cita el commit actual. Si dice "+ cambios sin commit", hay que hacer el commit y volver
  a correr el paso 4.
