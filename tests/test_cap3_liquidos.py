"""Verificación contra los ejemplos del Cap. 3."""

import pytest

from ductos.liquidos import (Linea, PozoInyector, diametro, diametro_optimo_graper,
                             gasto)

rel = pytest.approx


def test_ejemplo_3_4_1():
    linea = Linea(d=3.937, L=3000 / 5277, rug=0.0006)
    a = linea.dp(2560, 0.9, 46)
    assert a["regimen"] == "laminar"
    assert a["NRe"] == rel(1173, abs=1)
    assert a["dp_f"] == rel(11.7, abs=0.05)
    b = linea.dp(12800, 0.9, 46)
    assert b["NRe"] == rel(5865, abs=1)
    assert b["dp_f"] == rel(192.8, abs=0.1)


def test_diametro_y_gasto_son_inversos():
    d = diametro(q=12800, gamma=0.9, mu=46, L=3000 / 5277, dp_T=192.84)
    assert d == rel(3.937, abs=1e-3)
    q = gasto(d=3.937, gamma=0.9, mu=46, L=3000 / 5277, dp_T=192.84)
    assert q == rel(12800, rel=1e-3)


def _pozo(d_linea=4.0, d_tp=2.441):
    return PozoInyector(
        linea=Linea(d=d_linea, L=2.0, dh=20.0, rug=0.0006, E=0.9),
        tp=Linea(d=d_tp, L=2.273, dh=-12000.0, rug=0.0006, E=0.9),
        gamma_linea=1.0, mu_linea=1.0, gamma_tp=1.0, mu_tp=0.8, ii=3.521)


@pytest.mark.parametrize("q, pth, pwf, pws", [
    (6000, 2997.8, 7011.3, 5307.2),
    (8000, 2932.5, 6084.4, 3812.3),
    (10000, 2851.0, 4913.1, 2073.0),
    (7250, 2959.0, 6460.7, 4401.6),
])
def test_tabla_3_1_nodos(q, pth, pwf, pws):
    r = _pozo().nodos(q, 3100)
    assert r["pth"] == rel(pth, rel=1e-3)
    assert r["pwf"] == rel(pwf, rel=5e-3)
    assert r["pws"] == rel(pws, rel=2e-2)


def test_gasto_inyeccion():
    assert _pozo().gasto_inyeccion(3100, 4402) == rel(7250, rel=3e-3)


@pytest.mark.parametrize("pp, q", [(2000, 5710), (4000, 8375), (5000, 9530)])
def test_tabla_3_2_presion_bombeo(pp, q):
    assert _pozo().gasto_inyeccion(pp, 4402) == rel(q, rel=3e-3)


@pytest.mark.parametrize("q, d, dp, bhp", [
    (15000, 5.536, 28.048, 8.957),
    (40000, 8.357, 22.978, 19.568),
    (100000, 12.267, 19.273, 41.034),
    (155000, 14.734, 17.764, 58.621),
])
def test_tabla_3_6_graper(q, d, dp, bhp):
    r = diametro_optimo_graper(q, gamma=0.7, mu=0.525)
    assert r["d"] == rel(d, abs=2e-3)
    assert r["dp_milla"] == rel(dp, rel=1e-3)
    assert r["bhp_milla"] == rel(bhp, rel=1e-3)
