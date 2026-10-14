"""Verificación contra los Ejemplos 1.1 y 1.2 (aceite y agua)."""

import pytest

from ductos.pvt_aceite import (AceiteNegro, bo_oistein, bo_standing, bo_vazquez,
                               gamma_gas_disuelto, pb_lasater, pb_oistein, pb_standing,
                               pb_vazquez, rs_lasater, rs_oistein, rs_standing, rs_vazquez)
from ductos.pvt_agua import (bw_saturada, densidad_agua_saturada, rsw,
                             tension_superficial_agua_gas, viscosidad_agua_saturada)

rel = pytest.approx
GGS = 0.85222


def _aceite():
    a = AceiteNegro(api=30, gamma_g=0.85, R=750, p_sep=120, t_sep=75)
    a.calibrar(pb_lab=3300, rsfb_lab=750, bofb_lab=1.455, t_b=220)
    return a


def test_gamma_gs_vazquez():
    assert _aceite().gamma_gs == rel(GGS, abs=1e-5)


@pytest.mark.parametrize("fn, args, esperado", [
    (pb_standing, (750, 0.85, 30, 220), 3352),
    (pb_vazquez, (750, GGS, 30, 220), 3626),
    (pb_oistein, (750, 0.85, 30, 220), 2971),
    (pb_lasater, (750, 0.85, 30, 220), 3299),
])
def test_tabla_1_pb(fn, args, esperado):
    assert fn(*args) == rel(esperado, abs=1.5)


@pytest.mark.parametrize("fn, args, esperado", [
    (rs_standing, (3300, 0.85, 30, 220), 736),
    (rs_oistein, (3300, 0.85, 30, 220), 653),
    (rs_lasater, (3300, 0.85, 30, 220), 752),
])
def test_tabla_1_rs(fn, args, esperado):
    assert fn(*args) == rel(esperado, abs=1.0)


def test_tabla_1_rs_vazquez_errata():
    # El libro reporta 610; es incoherente con su propia pb de Vázquez (3626): ver erratas
    rs = rs_vazquez(3300, GGS, 30, 220)
    assert rs == rel(676.6, abs=0.5)
    assert rs == rel(750 * (3300 / pb_vazquez(750, GGS, 30, 220)) ** 1.0937, rel=1e-6)


def test_tabla_1_bo():
    go = 141.5 / 161.5
    assert bo_standing(750, 0.85, go, 220) == rel(1.4723, abs=3e-4)
    assert bo_vazquez(750, GGS, 30, 220) == rel(1.3729, abs=2e-4)
    assert bo_oistein(750, 0.85, go, 220) == rel(1.4324, abs=2e-4)


def test_calibracion_elige_como_el_libro():
    a = _aceite()
    assert (a.corr_pb, a.corr_rs, a.corr_bo) == ("Lasater", "Lasater", "Standing")
    assert a.c_rs == rel(0.9973, abs=3e-4)
    assert a.c_pb == rel(1.0003, abs=2e-4)
    assert a.c_bo == rel(0.9882, abs=2e-4)
    assert a.c_gd == rel(1.1042, abs=1e-4)
    assert gamma_gas_disuelto(750, 30) * a.c_gd == rel(0.85, abs=1e-6)


def test_ejemplo_1_1_condiciones_de_interes():
    a = _aceite()
    rs, ggd = a.rs_y_gamma_gd(2000, 200)
    assert rs == rel(456, abs=0.5)
    assert ggd == rel(0.8847, abs=1e-4)
    assert a.gamma_gf(2000, 200) == rel(0.7962, abs=1e-4)
    assert a.viscosidad(2000, 200) == rel(0.7581, abs=1e-3)
    assert a.tension_superficial(2000, 200) == rel(6.1625, abs=1e-3)


@pytest.mark.parametrize("p, rs, ggd, ggf, bo", [
    (3000, 701.58, 0.8558, 0.7660, 1.4149),
    (2500, 580.19, 0.8701, 0.7813, 1.3504),
    (2000, 456.24, 0.8847, 0.7961, 1.2851),
    (1500, 331.02, 0.8995, 0.8109, 1.2199),
    (1000, 209.60, 0.9138, 0.8253, 1.1579),
    (500, 99.65, 0.9268, 0.8382, 1.1034),
    (100, 23.28, 0.9358, 0.8473, 1.0669),
    (15, 8.51, 0.9376, 0.8490, 1.0599),
])
def test_tabla_ejemplo_1_2(p, rs, ggd, ggf, bo):
    a = _aceite()
    r = a.tabla([p], 200)[0]
    assert r["Rs"] == rel(rs, rel=2e-3)
    assert r["gamma_gd"] == rel(ggd, abs=2e-4)
    assert r["gamma_gf"] == rel(ggf, abs=2e-3)
    assert r["Bo"] == rel(bo, abs=5e-4)


def test_agua_ejemplo_1_1():
    bw = bw_saturada(2000, 200, variante="ejemplo")
    assert bw == rel(1.0103, abs=1e-4)
    assert bw_saturada(2000, 200) == rel(1.0297, abs=1e-4)
    assert densidad_agua_saturada(bw) == rel(61.795, abs=1e-3)
    assert viscosidad_agua_saturada(200, 0.18) == rel(0.30967, abs=1e-5)
    assert tension_superficial_agua_gas(2000, 200) == rel(42.67, abs=0.01)
    assert rsw(2000, 200) == rel(12.21, abs=0.01)
    assert rsw(2000, 200, 0.18) == rel(12.12, abs=0.01)
