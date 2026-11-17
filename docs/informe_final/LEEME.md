# Informe final: Plan de Control de Temperatura

`Plan_Control_Temperatura_G4.docx` (Word) y su PDF se arman **desde los resultados**, con
`scripts/informe_final.py`: cada número, tabla y figura sale de las corridas. Si cambia un dato, se vuelven a
correr los estudios afectados y el script, y el informe queda al día sin editar el texto a mano.

## Cómo regenerarlo

Desde la raíz del repositorio, con el entorno instalado (`pip install -r requirements.txt`, que incluye
`python-docx`):

```powershell
# 1. Clima real (una vez): 30 temporadas ERA5 alrededor del 19 nov, hasta 40 días después (~3 min)
mcplan clima era5 --fecha-inicio 2025-11-19T08:00 --anios 1996 2025 --dias 40 --margen 7

# 2. Estudios (~2 h en total con 7 procesos; cada uno retoma lo que falte si se corta)
mcplan estudio config/estudio_2025.yaml          # familia: 17 planes x 8 curvas x 3 historias
mcplan estudio config/confiabilidad_2025.yaml    # propuesta: 8 curvas x 30 temporadas
mcplan estudio config/estres_2025.yaml           # propuesta con +20 % de calor
mcplan estudio config/retiro_2025.yaml           # retiro de protecciones a los 14, 21 y 26 días (30 días)
mcplan estudio config/retiro_tardio_2025.yaml    # retiro el día 33 (900 h simuladas)

# 3. Peor escenario de la propuesta: malla fina y corrida de detalle (cortes, vista 3D, VTK, tubos)
mcplan estudio config/confiabilidad_2025.yaml --sensibilidad unica-M10-T10-N25 <curva> <ERA5_año>
mcplan estudio config/confiabilidad_2025.yaml --detalle unica-M10-T10-N25 <curva> <ERA5_año>

# 4. El informe (con --pdf, Word lo exporta a PDF y cuenta las páginas)
python scripts/informe_final.py --pdf
```

El peor escenario (`<curva> <ERA5_año>`) es la fila con mayor `M_D` de
`resultados/confiabilidad_2025/escenarios.csv`.

## Qué hay que completar a mano

- La portada: integrantes, curso y profesor (quedan con líneas para llenar).
- El commit: la portada muestra el commit actual; si dice "+ cambios sin commit", hacer el commit y volver a
  correr el paso 4 para que el informe cite el SHA evaluado.

## Por qué la fecha es el 19 nov 2025

El curso fija como fecha hipotética la de la defensa. El grupo la trasladó al mismo día de 2025 para que todo el
horizonte de 14 días tenga **clima real medido** (ERA5) y la propuesta se pueda verificar contra lo que ocurrió.
El flujo con pronóstico (`config/estudio.yaml`: ensemble ECMWF de 51 miembros más climatología) sigue disponible
para la fecha final: basta con cambiar `t0` y correr el estudio.

## Configuración

`config/estudio_2025.yaml` hereda de `config/estudio.yaml` (`hereda:`): toma la familia de planes, la física y
las curvas, y sólo cambia `t0`, la carpeta de salida y el clima (`fuente: era5_bloques`, 30 temporadas reales,
suelo inicial real). `confiabilidad_2025.yaml`, `estres_2025.yaml` y `retiro_2025.yaml` heredan a su vez de
`estudio_2025.yaml`.
