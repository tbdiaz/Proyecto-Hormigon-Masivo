"""mcplan: planificación de vertidas de hormigón masivo mediante simulación térmica.

Módulos principales:
    hidratacion  -- calor de hidratación (curva adiabática + edad equivalente)
    solver1d     -- solver térmico 1D por balance de energía en celdas
    metricas     -- Tmax y diferencial núcleo-superficie según criterios del curso
    caso         -- lectura de casos desde archivos YAML
    ifc_io       -- creación y lectura del modelo IFC (IfcOpenShell)
    cli          -- interfaz de línea de comandos `mcplan`
"""

__version__ = "0.1.0"
