"""Verificación contra los ejemplos del Cap. 4 (secciones 4.5-4.10) y el Apéndice D."""

import pytest

from ductos.costos_gas import (DatosCostosGas, barrido_presion_descarga, diseno_economico,
                               razon_compresion_optima)
from ductos.estranguladores import (area_a_64avos, area_estrangulador, cd_cook_dotterweich,
                                    diametro_cook_dotterweich, diametro_ec_4_134,
                                    dp_valvula_seguridad, funcion_h, k_gas, relacion_critica)
from ductos.gas import Gasoducto, q_cullender_smith
from ductos.pozo_gas import IPRGas, PozoGas
from ductos.pvt_gas import Gas, gamma_razasa_katz, gamma_vitter, peso_molecular_condensado
from ductos.raspatubos import velocidad_gas, velocidad_liquido

rel = pytest.approx


# ---------------------------------------------------------------- Apéndice D

def test_cullender_smith_ejemplo_d1():
    g = Gas(0.6)
    assert q_cullender_smith(1900, 1500, g, Gasoducto(28, 41, 150, h=50)) / 1e6 == rel(1215.43, rel=1e-3)
    assert q_cullender_smith(1900, 1500, g, Gasoducto(28, 41, 150)) / 1e6 == rel(1216.33, rel=1e-3)


# ---------------------------------------------------------------- 4.7 Ejemplo 4.7

def test_coeficientes_ejemplo_4_7():
    d = DatosCostosGas()
    assert d.K1 == rel(46.254, abs=1e-3)
    assert d.b == rel(0.21875)
    assert d.fi == rel(0.171, abs=1e-3)
    assert d.K2 == rel(4.823, abs=1e-3)
    assert d.K3 == rel(125.88, rel=1e-3)
    rc = razon_compresion_optima(d)
    assert rc == rel(1.1645, abs=1e-4)
    assert d.hp_por_mmpcd(rc) == rel(9.46, abs=0.01)


def test_diseno_economico_pd_1164():
    r = diseno_economico(1200e6, 1164.5)
    assert r["pe"] == rel(1000, abs=0.1)
    assert r["Z"] == rel(0.854, abs=1e-3)
    assert r["dif"] == rel(0.01404, abs=1e-5)
    assert r["d"] == rel(35.32, abs=0.01)
    assert r["f"] == rel(0.0087679, abs=1e-6)
    assert r["L"] == rel(42.83, rel=3e-3)
    assert r["tuberia"] == rel(58653.766, rel=1e-3)
    assert r["total"] == rel(87151.341, rel=2e-3)


@pytest.mark.parametrize("pd, d, L, ct", [
    (4000, 21.14, 37.60, 97272.69),
    (2750, 24.29, 38.66, 88396.06),
    (1750, 29.50, 40.89, 86269.38),
    (750, 42.61, 44.25, 88008.42),
])
def test_tabla_4_8(pd, d, L, ct):
    r = barrido_presion_descarga(1200e6, [pd])[0]
    assert r["d"] == rel(d, abs=0.02)
    assert r["L"] == rel(L, rel=3e-3)
    assert r["total"] == rel(ct, rel=1e-3)


# ---------------------------------------------------------------- 4.8 Ejemplo 4.8

def _pozo():
    return PozoGas(gas=Gas(0.7), ipr=IPRGas(16280, 0.681),
                   tp=Gasoducto(2.441, 1.1976, 140, h=1.1976 * 5280, rug=0.0007),
                   ld=Gasoducto(2.991, 1.865, 60, h=300, rug=0.0007), p_sep=50)


@pytest.mark.parametrize("q, pe", [(2e6, 153), (4e6, 292), (10e6, 698), (12e6, 828), (16e6, 1081)])
def test_ejemplo_4_8_linea_descarga(q, pe):
    assert _pozo().pe(q) == rel(pe, rel=5e-3)


@pytest.mark.parametrize("pws, q, pth", [
    (2000, 18, 472.0), (2000, 16, 886.0),
    (1500, 12, 650), (1500, 8, 1052),
    (1000, 8, 401), (1000, 6, 642), (1000, 2, 834),
    (500, 4, 185), (500, 3, 317), (500, 2, 385), (500, 1, 420),
])
def test_ejemplo_4_8_tabla_pth(pws, q, pth):
    # Se excluyen 1031 (pws=2000, q=12) y 979 (pws=1500, q=10): ver docs/erratas.md
    assert _pozo().pth(q * 1e6, pws) == rel(pth, rel=1.3e-2)


# ---------------------------------------------------------------- 4.9 estranguladores

def test_relaciones_criticas_tabla_4_9():
    for k, ref in [(1.2, 0.5644), (1.3, 0.5457), (1.4, 0.5283), (1.67, 0.4867)]:
        assert relacion_critica(k) == rel(ref, abs=2e-4)


def test_ejemplo_4_9():
    assert relacion_critica(1.25) == rel(0.555, abs=1e-3)
    assert funcion_h(0.8, 1.25) == rel(0.0305, abs=1e-4)
    a = area_estrangulador(1340927, 514.3, 411.5, Gas(0.69), 140, 1.25, cd=1.0, z=0.93)
    assert round(area_a_64avos(a)) == 23          # libro: 23.83 → 24/64
    assert cd_cook_dotterweich(32) == rel(0.828, abs=1e-3)


def test_ejemplo_4_9_cook_dotterweich_y_4_10():
    g = Gas(0.69)
    r = diametro_cook_dotterweich(1340927, 514.3, 411.5, g, 140, 1.25)
    assert r["d64_libro"] == rel(23.83, abs=0.02)        # área efectiva, como el libro
    assert r["d64"] == rel(26.17, abs=0.05)              # diámetro geométrico (A = CdA/Cd)
    assert diametro_ec_4_134(1340927, 514.3, 140, 0.69) == rel(23.1, abs=0.05)


def test_k_gas_figura_4_16():
    assert k_gas(0.69, 140) == rel(1.25, abs=0.01)       # el Ejemplo 4.9 usa k = 1.25 medido
    for t in (50, 100, 150, 200, 250):                    # k decrece con γ en el intervalo confiable
        ks = [k_gas(g, t) for g in (0.65, 0.7, 0.75, 0.8)]
        assert all(a > b for a, b in zip(ks, ks[1:]))
        assert all(1.15 < x < 1.32 for x in ks)


def test_ejemplo_4_11_valvula_seguridad():
    g = Gas(0.6)
    r = dp_valvula_seguridad(55e6, 1900, 150, g, 2.992, 1.0, z=0.87, con_psi4=False)
    assert r["psi"] == rel(0.334, abs=1e-3)
    assert r["dp"] == rel(907.74, abs=0.05)
    assert dp_valvula_seguridad(55e6, 1900, 150, g, 2.992, 1.0, z=0.87)["dp"] == rel(896.4, abs=0.1)
    it = dp_valvula_seguridad(55e6, 1900, 150, g, 2.992, 1.0, y=None, k=1.28)
    assert it["Y"] == rel(0.85, abs=0.01)                # valida la aproximación del API


# ---------------------------------------------------------------- 4.10 raspatubos

def test_ejemplos_4_12_y_4_13():
    assert velocidad_liquido(200000, 19.25)["v"] == rel(6.433, abs=1e-3)
    assert velocidad_gas(200e6, 1115, 915, 60, 19.25, Gas(0.6))["v"] == rel(13.93, abs=0.02)


# ---------------------------------------------------------------- 4.5 condensados

def test_peso_molecular_condensado_ec_4_53():
    # Las dos formas de la Ec. 4.53 son equivalentes
    api = 50
    go = 141.5 / (131.5 + api)
    assert peso_molecular_condensado(api=api) == rel(peso_molecular_condensado(gamma_o=go), rel=5e-3)


def test_vitter_variantes():
    assert gamma_vitter(0.65, 0.78, 16667) == rel(0.8103, abs=1e-4)
    assert gamma_vitter(0.65, 0.78, 16667, 0.78, variante="libro") == rel(0.8217, abs=1e-4)
    assert gamma_razasa_katz(0.65, 0.78, 16667, 138.0) > 0.65
