# Curvas de hidratación del profesor

`source_terms.txt` trae las curvas Q(t) [W/m³] que entrega el profesor, sin modificar. Lo leen `mcplan estudio`,
`mcplan simular3d --curva <nombre>` y el editor, que muestra las curvas solo si el archivo está.

## Versión vigente: corregida (14 nov 2026)

- **Curvas:** 8, de las probetas M1–M4 en las campañas SS (13–20 abr 2025, 240,3 h) y 22A (22 abr 2025, 163,0 h).
- **Método:** inversión por volúmenes finitos con regularización de Tikhonov, con β = 2,5 W/m²K en todas las
  caras.
- **Probetas:** ρ = 2476 kg/m³, c_p = 1023 J/kgK, 350,77 kg/m³ de cemento.
- **Calor del registro:** 61–100 MJ/m³, entre 24 y 40 K adiabáticos con el ρc_p de la probeta. La lectura de
  `curvas.py` reproduce la tabla del encabezado.
- **Origen del tiempo:** inicio del registro + 2 h (`desfase_h: 2.0` en `config/estudio.yaml`); paso de 10 min.
- **Confianza:** el profesor marca 22A como de menor confianza.

## Versión retirada

Las primeras curvas (mismos nombres, recuperadas con una PINN y β de 7–15 W/m²K) estaban malas: liberaban
193–402 MJ/m³, equivalente a 83–173 °C adiabáticos. Se retiraron el 14 nov 2026 junto con los resultados del
estudio calculados con ellas. Quedan en el historial de git (commit `85516d7`).

## Formato que lee `curvas.py`

- comentarios con `#`;
- el β por curva (`#   1 NOMBRE_W_m3 [W/m^3]  beta = 9.6`) o uno común (`beta = 2.5 W/(m^2 K) on all faces`);
- una línea `# time_h NOMBRE1 NOMBRE2 ...`, con nombres PROBETA_CAMPAÑA;
- columnas: tiempo en horas y Q en W/m³ de cada curva (`NaN` donde una curva ya terminó).
