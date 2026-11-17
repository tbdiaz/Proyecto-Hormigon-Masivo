"""Elegir temporadas ERA5 por año (para filtrar variantes en las temporadas que gobiernan)."""

from types import SimpleNamespace

import numpy as np
import pytest

from mcplan.clima import elegir_temporadas


def _historias():
    return [(f"ERA5 {a}", SimpleNamespace(T_aire=np.full(24 * 8, 10.0 + i))) for i, a in enumerate(range(2018, 2024))]


def test_elige_por_anio_con_el_mismo_peso_y_la_etiqueta_del_estudio():
    el = elegir_temporadas(_historias(), [2022, 2019])
    assert [h["etiqueta"] for h in el] == ["ERA5_2022", "ERA5_2019"]
    assert [h["peso"] for h in el] == [0.5, 0.5]
    assert el[0]["indice_C"] == pytest.approx(14.0)  # 2022 es la quinta: 10 + 4


def test_temporada_sin_datos_avisa():
    with pytest.raises(ValueError, match="2030"):
        elegir_temporadas(_historias(), [2022, 2030])
