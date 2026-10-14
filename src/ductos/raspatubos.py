"""Cap. 4.10 — Raspatubos ("diablos"): velocidad de corrida.

Rangos recomendados por el libro:
    gas natural  5 a 15 mph  (7.3 a 22 pie/s)
    líquidos     2 a 10 mph  (2.9 a 14.7 pie/s)
Con la conexión lateral de la trampa abierta la velocidad baja de 3 a 5 %.
"""

from __future__ import annotations

from .constantes import f_a_r
from .pvt_gas import Gas

MPH_A_PIES_S = 5280.0 / 3600.0
RANGO_GAS_MPH = (5.0, 15.0)
RANGO_LIQUIDO_MPH = (2.0, 10.0)


def velocidad_gas(q: float, p1: float, p2: float, t_f: float, d: float, gas: Gas) -> dict:
    """Ec. 4.138 — v = 5.995×10⁻⁵ q Z T / (p̄ d²)  [pie/s].

    q [pie³/d @ c.s.], p [psia], T [°F], d [pg]; p̄ = (p1 + p2)/2 como en el Ejemplo 4.13.
    """
    pm = 0.5 * (p1 + p2)
    z = gas.z(pm, t_f)
    v = 5.995e-5 * q * z * f_a_r(t_f) / (pm * d**2)
    return {"v": v, "mph": v / MPH_A_PIES_S, "Z": z, "p_media": pm,
            "en_rango": RANGO_GAS_MPH[0] <= v / MPH_A_PIES_S <= RANGO_GAS_MPH[1]}


def velocidad_liquido(q: float, d: float) -> dict:
    """Ec. 4.139 — v = 0.01192 q / d²  [pie/s], q [bl/d], d [pg]."""
    v = 0.01192 * q / d**2
    return {"v": v, "mph": v / MPH_A_PIES_S,
            "en_rango": RANGO_LIQUIDO_MPH[0] <= v / MPH_A_PIES_S <= RANGO_LIQUIDO_MPH[1]}
