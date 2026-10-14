"""Cap. 4.8 — Predicción del comportamiento de un pozo productor de gas seco.

Sistema: yacimiento → tubería de producción (TP) → estrangulador → línea de descarga (LD)
→ separador. Nodo de solución: la cabeza del pozo.

    Yacimiento (Ec. 4.101):   q = C (pws² − pwf²)ⁿ
    TP y LD:                  ecuación general de flujo (Ec. 4.19)

Para cada gasto se calcula la presión en la cabeza antes del estrangulador (pth, desde el
yacimiento) y la necesaria después de él (pe, desde el separador). El gasto máximo, sin
estrangulador, es aquel en que pth = pe (Fig. 4.10).

Incluye además la velocidad de erosión (4.8.1) y el gasto mínimo para descargar líquidos
según Turner (4.8.2).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from scipy.optimize import brentq

from .constantes import f_a_r
from .gas import Gasoducto, p1_general, p2_general
from .pvt_gas import Gas


@dataclass
class IPRGas:
    """Ec. 4.101 — q [pie³/d] = C (pws² − pwf²)ⁿ."""
    C: float
    n: float

    def q(self, pws: float, pwf: float) -> float:
        return self.C * (pws**2 - pwf**2) ** self.n

    def pwf(self, pws: float, q: float) -> float:
        """Ec. 4.102."""
        arg = pws**2 - (q / self.C) ** (1.0 / self.n)
        if arg <= 0:
            raise ValueError("el gasto excede el potencial del yacimiento")
        return math.sqrt(arg)

    def q_max(self, pws: float) -> float:
        return self.q(pws, 0.0)


@dataclass
class PozoGas:
    gas: Gas
    ipr: IPRGas
    tp: Gasoducto            # tubería vertical: L = profundidad [mi], h = profundidad [pies]
    ld: Gasoducto            # línea de descarga
    p_sep: float             # presión de separación [psia]

    def pth(self, q: float, pws: float) -> float:
        """Presión en la cabeza antes del estrangulador (desde el yacimiento)."""
        pwf = self.ipr.pwf(pws, q)
        return p2_general(pwf, q, self.gas, self.tp)["p2"]

    def pe(self, q: float) -> float:
        """Presión corriente abajo del estrangulador (desde el separador)."""
        return p1_general(self.p_sep, q, self.gas, self.ld)["p1"]

    def nodos(self, q: float, pws: float) -> dict:
        pwf = self.ipr.pwf(pws, q)
        try:
            pth = p2_general(pwf, q, self.gas, self.tp)["p2"]
        except ValueError:
            pth = float("nan")
        return {"q": q, "pws": pws, "pwf": pwf, "pth": pth, "pe": self.pe(q), "psep": self.p_sep}

    def gasto_maximo(self, pws: float) -> float:
        """Gasto sin estrangulador: pth(q) = pe(q)."""
        def g(q):
            try:
                return self.pth(q, pws) - self.pe(q)
            except ValueError:
                return -pws
        q_hi = 0.999 * self.ipr.q_max(pws)
        return brentq(g, 1e3, q_hi, xtol=1.0)


# ------------------------------------------------------------- 4.8.1 erosión

def velocidad_erosion(p: float, t_f: float, gas: Gas, c4: float = 100.0) -> float:
    """Ecs. 4.103-4.105 — ve = c4 / √ρg  [pie/s]; c4 entre 75 y 150."""
    return c4 / math.sqrt(gas.densidad(p, t_f))


def gasto_erosion(d: float, p: float, t_f: float, gas: Gas, c4: float = 100.0,
                  p_b: float = 14.7, t_b: float = 60.0) -> float:
    """Ec. 4.106 — gasto a c.s. [pie³/d] a la velocidad de erosión en una tubería de d [pg].

    Deducido aquí a partir de ve; la constante de la fuente (1291.539×10⁶) parece tener
    el exponente equivocado (ver docs/erratas.md): para c4 = 100 el coeficiente es 1.2906×10⁶.
    """
    area_ft2 = math.pi * d**2 / 4.0 / 144.0
    z = gas.z(p, t_f)
    ve = velocidad_erosion(p, t_f, gas, c4)
    return ve * area_ft2 * 86400.0 * (p / p_b) * (f_a_r(t_b) / (z * f_a_r(t_f)))


# ------------------------------------------------------------- 4.8.2 Turner

def velocidad_turner(p: float, liquido: str = "agua") -> float:
    """Ecs. 4.108-4.109 — velocidad mínima del gas para arrastrar gotas [pie/s].

    Válidas para γg = 0.6 y T = 140 °F (ρg ≈ 0.0031 p). Si fluyen agua y condensado,
    Turner recomienda la de agua.
    """
    if liquido == "agua":
        return 5.62 * (67.0 - 0.0031 * p) ** 0.25 / (0.0031 * p) ** 0.5
    if liquido == "condensado":
        return 4.02 * (45.0 - 0.0031 * p) ** 0.25 / (0.0031 * p) ** 0.5
    raise ValueError("liquido debe ser 'agua' o 'condensado'")


def velocidad_terminal(sigma: float, rho_l: float, rho_g: float) -> float:
    """Ec. 4.107 (forma general) — σ [dina/cm], densidades [lbm/pie³]."""
    return 2.04 * (sigma * (rho_l - rho_g)) ** 0.25 / rho_g**0.5


def gasto_minimo_turner(d: float, p: float, t_f: float, gas: Gas, liquido: str = "agua",
                        p_b: float = 14.7, t_b: float = 60.0) -> float:
    """Gasto a c.s. [pie³/d] que mantiene la velocidad de Turner en la tubería de d [pg]."""
    area_ft2 = math.pi * d**2 / 4.0 / 144.0
    z = gas.z(p, t_f)
    v = velocidad_turner(p, liquido)
    return v * area_ft2 * 86400.0 * (p / p_b) * (f_a_r(t_b) / (z * f_a_r(t_f)))
