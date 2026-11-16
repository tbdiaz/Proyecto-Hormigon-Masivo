# Datos por confirmar antes de la defensa (jueves 19 nov 2026)

Estado al 14 nov 2026. Cada punto dice qué se supone hoy, por qué importa y a quién preguntar. La prioridad
indica cuánto puede cambiar la decisión: **alta** puede cambiar qué plan cumple, **media** cambia los números sin
cambiar la conclusión, **baja** es de forma o de trazabilidad.

## 1. Curvas de hidratación (profesor)

1. **Prioridad alta. ¿Es la versión definitiva?** El archivo nuevo es distinto del anterior: otro método,
   β = 2,5 W/m²K y 61–100 MJ/m³ en vez de 193–402. Falta que el profesor confirme que es el final.
2. **Prioridad alta. ρc_p del hormigón.** El modelo usa el de Ju y Lei (2400 × 970 = 2,33 MJ/m³K); las probetas
   tienen 2476 × 1023 = 2,53. Con el del modelo el aumento adiabático sale ~9 % mayor, del lado conservador.
   ¿Cuál usar?
3. **Prioridad alta. ¿Reescalar a la dosificación de la obra?** Las probetas tienen 351 kg/m³ de cemento.
   ¿El hormigón de la fundación tiene lo mismo?
4. **Prioridad alta. Edad al inicio del archivo.** t = 0 es el "inicio del registro + 2 h". ¿Cuánto después del
   mezclado empezó el registro? Si la probeta ya tenía horas de edad, la fuente debe desplazarse (`desfase_h`).
5. **Prioridad media. Peso de la campaña 22A.** El profesor la marca de menor confianza; hoy pesa igual que SS
   (8 curvas equiprobables) y "sólo SS" se informa como sensibilidad. ¿Se acepta, o se usan sólo las SS?
6. **Prioridad media. ¿Réplicas o mezclas?** ¿M1–M4 son réplicas de una misma mezcla o mezclas distintas? Si son
   mezclas, la mezcla pasa a ser una decisión.
7. **Prioridad media. Incertidumbre de β.** El encabezado da la energía con β = 2,1 y 3,0 (×0,84 y ×1,20).
   ¿Basta informarla como sensibilidad o hay que incluirla en el modelo de incertidumbre?

## 2. Criterio del ΔT (profesor)

8. **Prioridad alta. ¿Dónde se mide Tmin?**
   - Hoy: a 25 mm bajo las caras expuestas planas, sin aristas, incluidas las caras enterradas; la base sobre el
     terreno queda fuera.
   - Con protección, el punto que manda es el borde inferior de la losa (z ≈ 0,14 m), enfriado por el suelo frío
     de abajo. ¿Ese punto cuenta?
9. **Prioridad alta. Horizonte y retiro de las mantas.** Los planes con protección dejan la manta hasta las
   600 h, después del horizonte de evaluación (14 días y al menos 7 días tras la última colocación), y el retiro
   se informa aparte. ¿Se acepta, o el retiro debe ocurrir dentro del horizonte? Retirarla de golpe a los 14 días
   haría saltar el ΔT, porque el núcleo sigue a ~50 °C.

## 3. Recursos y controles admitidos (profesor)

10. **Prioridad alta. Controles permitidos.** El plan propuesto (`unica-M10-T10-N25`) usa:
    - hormigón **preenfriado a 10 °C** (hielo o agua fría);
    - **mantas de 10 cm en todas las caras**, junto con el moldaje, por más de tres semanas;
    - **tubos HDPE sólo en V3, V4 y V5** con **agua a 25 °C** durante 14 días: hace falta agua templada, no fría.

    ¿Se admiten las tres en la defensa? Si alguna no se admite, ningún plan evaluado cumple: todos los que cumplen
    usan las tres. Con manta + 10 °C sin tubos, $P_f$ = 37,5 %; con manta + tubos a 18 °C, 62,5 %. La
    combinación 10 °C + tubos sin manta no se evaluó.
11. **Prioridad media. Tasa de colocación:** se supone una bomba a 27 m³/h (Ju y Lei vaciaron ~637 m³ en ~24 h),
    con etapas en serie y sin simultaneidad. ¿Hay una capacidad definida por el curso?
12. **Prioridad media. Moldaje:** terciado de 18 mm en V1, V2 y V5, con retiro a las 72 h (pset del IFC). Los
    planes con protección lo dejan hasta las 600 h. ¿Es admisible?

## 4. Sitio y construcción (grupo)

13. **Prioridad media. Fecha y hora de inicio.** El estudio usa t0 = jueves 19 nov 2026, 08:00 (hora de Chile),
    como fija la pauta. El caso determinista `config/caso_14d.yaml` usa 27 nov (fecha de las entregas E10–E12).
    Confirmar la hora.
14. **Prioridad media. Excavación:** ¿la fundación va enterrada? ¿A qué profundidad queda el terreno natural (se
    probó 3,3 m) y con qué sobreancho (1 m)? Cambia el abrigo del viento ([docs/excavacion.md](excavacion.md)).
15. **Prioridad media. Relleno:** ¿cuándo se rellena y con qué? Recién puesto, el relleno frío empeora el ΔT
    ([docs/relleno.md](relleno.md)).
16. **Prioridad baja. Terreno para el perfil del viento:** EN 1991-1-4, categoría II (campo abierto). ¿Litueche es
    más bien categoría 0 o I, costa expuesta?

## 5. Materiales (fichas técnicas)

17. **Prioridad media. Hormigón de la obra:** k = 2,222 W/mK de Ju y Lei; ¿hay dato propio?
18. **Prioridad media. Suelo:** k = 1,389 W/mK, ρ = 1800 kg/m³ y c_p = 1273 J/kgK de Ju y Lei; T inicial de
    ERA5; 3 m de profundidad con el fondo a temperatura fija. El suelo frío bajo la losa es lo que hoy domina el
    ΔT, así que conviene un dato del sitio o una sensibilidad.
19. **Prioridad media. Mantas:** k = 0,04 W/mK en 5 y 10 cm. Respaldar con una ficha técnica.
20. **Prioridad baja. Superficie:** absortividad 0,65 y emisividad 0,9 del hormigón (rango típico).
21. **Prioridad baja. Grafito:** k = 4,11 W/mK de Chen et al. (2024), medido en laboratorio y en seco. Faltan el
    calor de hidratación y la durabilidad ([docs/tubos.md](tubos.md)).

## 6. Clima

22. **Prioridad media. Ensemble más reciente.** El estudio usa el ensemble del 12 nov. Antes de la defensa
    conviene `mcplan clima ensemble`, para cubrir hasta ~2 dic, y volver a correr el estudio con el clima
    congelado.
23. **Prioridad baja. Factor de abrigo 0,67** en el hoyo (Nakamura y Oke 1988, tomado de citas secundarias).
    Verificar en el artículo original.

## 7. Verificación numérica (miércoles, cluster)

24. **Prioridad media. Malla fina en más escenarios.**
    - **Ya confirmado:** el peor escenario de la propuesta con Δx = 0,11 m, que subió de 18,18 a 18,33 °C (se
      mantiene bajo 20).
    - **Descartado:** el candidato anterior, `unica-M10-T10-A25`, pasó de 19,78 a 20,63 °C, porque su punto frío
      estaba bajo el escalón de V2, donde la malla importa.
    - **Pendiente para el cluster:** confirmar más escenarios de la propuesta con malla fina, y probar Δx = 0,055 m
      en el peor.
25. **Prioridad media. Datos del cluster:** si usa SLURM, cuántos núcleos asigna, cuánta memoria y qué Python
    tiene (hacen falta numpy, ifcopenshell y pyyaml). Con eso se agregan `--parte i/n` y `--fusionar` al estudio.

## 8. Entrega

26. **Prioridad baja.** Hacer el commit y declarar el SHA de la defensa en Canvas.
