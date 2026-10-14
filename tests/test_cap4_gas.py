"""Verificación contra los ejemplos del Cap. 4."""

import pytest

from ductos import friccion as fr
from ductos.gas import (Gasoducto, factor_s, fraccion_loop, longitud_efectiva,
                        p2_general, presion_fondo_estatica, presion_media,
                        q_clinedinst, q_general, q_panhandle_a, q_panhandle_b,
                        q_smith, q_weymouth)
from ductos.pvt_gas import Gas

rel = pytest.approx
GAS = Gas(0.6)
Q = 1212.557807e6


def test_ejemplo_4_1_presion_estatica():
    r = presion_fondo_estatica(1920, 6320, 140, Gas(0.85), media="cubica")
    assert r["pws"] == rel(2459.3, rel=5e-4)
    r = presion_fondo_estatica(1920, 6320, 140, Gas(0.85), media="aritmetica")
    assert r["pws"] == rel(2459.3, rel=3e-3)


def test_presion_media():
    assert presion_media(1900, 1500) == rel(1707.84, abs=0.1)


def test_ejemplo_4_2_p2():
    r = p2_general(1900, Q, GAS, Gasoducto(28, 41, 150))
    assert r["p2"] == rel(1500, abs=0.5)
    assert r["Z"] == rel(0.87695, abs=1e-4)


@pytest.mark.parametrize("fn, esperado", [
    (q_general, 1212.558),
    (q_panhandle_a, 1472.509),
    (lambda *a: q_panhandle_b(*a, variante="libro"), 1721.002),
    (q_weymouth, 1126.940),
    (q_smith, 1216.592),
    (q_clinedinst, 1213.109),
])
def test_ejemplo_4_2_gastos(fn, esperado):
    q = fn(1900, 1500, GAS, Gasoducto(28, 41, 150)) / 1e6
    assert q == rel(esperado, rel=1e-3)


def test_ejemplo_4_3_inclinacion():
    s = factor_s(0.6, 50, 0.87695, 150)
    assert s == rel(0.002103, abs=1e-6)
    assert longitud_efectiva(41, s) == rel(41.043, abs=1e-3)
    d = Gasoducto(28, 41, 150, h=50)
    assert q_general(1900, 1500, GAS, d) / 1e6 == rel(1209.189, rel=1e-3)
    assert q_panhandle_a(1900, 1500, GAS, d) / 1e6 == rel(1468.912, rel=1e-3)
    assert q_weymouth(1900, 1500, GAS, d) / 1e6 == rel(1124.389, rel=1e-3)
    assert q_smith(1900, 1500, GAS, d) / 1e6 == rel(1213.879, rel=1e-3)


def test_ejemplo_4_4_factores_friccion():
    assert fr.f_rugoso(0.0006 / 14) == rel(0.010253, abs=1e-6)
    assert fr.f_rugoso(0.0006 / 16) == rel(0.010016, abs=1e-6)
    assert fr.f_rugoso(0.0006 / 18) == rel(0.009814, abs=1e-6)


def test_ejemplo_4_6_loop():
    x = fraccion_loop(Q, Q + 500e6, 28, 28, 1.0, 1.0)
    assert x == rel(0.6649, abs=1e-4)


def test_cullender_en_tablas_4_2_y_4_3():
    from ductos.gas import q_cullender_smith
    assert q_cullender_smith(1900, 1500, GAS, Gasoducto(28, 41, 150)) / 1e6 == rel(1216.335, rel=1e-3)
    assert q_cullender_smith(1900, 1500, GAS, Gasoducto(28, 41, 150, h=50)) / 1e6 == rel(1215.572, rel=1e-3)
