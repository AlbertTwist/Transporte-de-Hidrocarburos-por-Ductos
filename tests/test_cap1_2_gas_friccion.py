"""Verificación contra los valores del Ejemplo 1.1 (parte de gas) y Cap. 2."""

import pytest

from ductos import friccion as fr
from ductos.constantes import c_fric_liq, c_general_gas
from ductos.pvt_gas import Gas, pseudocriticas, wichert_aziz, z_dpr

rel = pytest.approx


def test_pseudocriticas_y_wichert_aziz():
    tpc, ppc = pseudocriticas(0.7962)
    assert tpc == rel(419.13, abs=0.01)
    assert ppc == rel(662.69, abs=0.01)
    t2, p2, eps = wichert_aziz(tpc, ppc, 0.0287, 0.2327)
    assert eps == rel(29.04, abs=0.01)
    assert t2 == rel(390.09, abs=0.01)
    assert p2 == rel(609.24, abs=0.01)


def test_z_bg_densidad_ejemplo_1_1():
    g = Gas(0.7962, y_co2=0.0287, y_h2s=0.2327)
    assert g.z(2000, 200) == rel(0.86067, abs=2e-4)
    assert g.bg(2000, 200) == rel(0.0080236, rel=1e-3)
    assert g.densidad(2000, 200) == rel(7.5813, rel=1e-3)


def test_viscosidad_gas_ejemplo_1_1():
    g = Gas(0.7962, y_co2=0.0287, y_h2s=0.2327)
    assert g.viscosidad(2000, 200, corregir=False, c_x=960) == rel(0.017145, rel=1e-3)
    assert g.viscosidad(2000, 200, gamma_correccion=0.85, c_x=960) == rel(0.018034, rel=1e-3)


def test_z_ideal_a_presion_cero():
    assert z_dpr(0.0, 1.5) == 1.0


def test_friccion_laminar_y_colebrook():
    assert fr.f_laminar(1173.0) == rel(0.05456, abs=1e-5)
    assert fr.f_colebrook(5865.0, 0.0006 / 3.937) == rel(0.03593, abs=1e-5)
    hist = fr.f_colebrook_iteraciones(5865.0, 0.0006 / 3.937)
    assert [round(x, 5) for x in hist] == [0.03700, 0.03576, 0.03596, 0.03592, 0.03593]


def test_friccion_rugoso_ejemplo_4_2():
    assert fr.f_rugoso(0.0006 / 28) == rel(0.0091085, abs=1e-6)


def test_jain_cerca_de_colebrook():
    for nre in (1e4, 1e5, 1e6):
        assert fr.f_jain(nre, 1e-4) == rel(fr.f_colebrook(nre, 1e-4), rel=0.02)


def test_constantes_derivadas():
    assert c_fric_liq() == rel(0.06056, abs=1e-5)       # con 5277 pie/milla (libro)
    assert c_general_gas() == rel(77.54, abs=0.05)       # literatura: 77.54
