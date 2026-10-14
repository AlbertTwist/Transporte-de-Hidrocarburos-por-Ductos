"""Cap. 4.9 — Flujo de gas a través de estranguladores y válvulas de seguridad.

Flujo isoentrópico de un gas a través de una restricción:

    p*/p0 = (2/(k+1))^(k/(k−1))                       relación crítica (Ec. 4.113)
    H     = r^(2/k) − r^((k+1)/k),  r = pe/pth          (Ec. 4.130)
    q     = C · Cd · A · pth · √[ k/(k−1) · H / (γ Z T) ]   (Ecs. 4.129 y 4.131)

Si pe/pth ≤ p*/p0 el flujo es crítico y H se evalúa con r = p*/p0 (la presión corriente
abajo ya no afecta al gasto). C se deduce de primeros principios (≈ 1.2434×10⁶ con A en
pg², q en pie³/d a 14.7 psia y 60 °F); la ecuación de Cook y Dotterweich (4.131) usa
155 500·√64.34 ≈ 1.2473×10⁶, el mismo valor.
"""

from __future__ import annotations

import math

from .constantes import GC, M_AIRE, R_GAS, f_a_r
from .pvt_gas import Gas


def relacion_critica(k: float) -> float:
    """Ec. 4.113 — p*/p0."""
    return (2.0 / (k + 1.0)) ** (k / (k - 1.0))


def relacion_temperatura_critica(k: float) -> float:
    """Ec. 4.112 — T*/T0."""
    return 2.0 / (k + 1.0)


def relacion_densidad_critica(k: float) -> float:
    """Ec. 4.114 — ρ*/ρ0."""
    return (2.0 / (k + 1.0)) ** (1.0 / (k - 1.0))


def funcion_h(r: float, k: float) -> float:
    """Ec. 4.130."""
    return r ** (2.0 / k) - r ** ((k + 1.0) / k)


def es_critico(pe: float, pth: float, k: float) -> bool:
    return pe / pth <= relacion_critica(k)


def constante_flujo(p_b: float = 14.7, t_b: float = 60.0) -> float:
    """C de la ecuación de gasto, derivada (A en pg², q en pie³/d @ c.s.)."""
    c_rho = 144.0 * M_AIRE / (R_GAS * 144.0)       # ρ = c_rho p γ /(Z T) ≈ 2.70 p γ/(Z T)
    rho_b = c_rho * p_b / f_a_r(t_b)               # densidad a c.s. por unidad de γ
    # ṁ/A = p √(2 gc ·144 · c_rho · k/(k−1) · H · γ/(Z T))   [lbm/(pie² s)], A en pie²
    return math.sqrt(2.0 * GC * 144.0 * c_rho) / 144.0 * 86400.0 / rho_b


def gasto_estrangulador(area: float, pth: float, pe: float, gas: Gas, t_f: float, k: float,
                        cd: float = 1.0, z: float | None = None) -> float:
    """Gasto de gas [pie³/d] por un estrangulador de área [pg²]."""
    z = z if z is not None else gas.z(pth, t_f)
    r = max(pe / pth, relacion_critica(k))
    h = funcion_h(r, k)
    return (constante_flujo() * cd * area * pth
            * math.sqrt(k / (k - 1.0) * h / (gas.gamma * z * f_a_r(t_f))))


def area_estrangulador(q: float, pth: float, pe: float, gas: Gas, t_f: float, k: float,
                       cd: float = 1.0, z: float | None = None) -> float:
    """Área [pg²] necesaria para un gasto dado (Ec. 4.129 despejada)."""
    return q / gasto_estrangulador(1.0, pth, pe, gas, t_f, k, cd, z)


def area_a_64avos(area: float) -> float:
    """Ec. 4.133 — diámetro en 64avos de pg."""
    return 64.0 * math.sqrt(4.0 * area / math.pi)


def cd_cook_dotterweich(d64: float) -> float:
    """Ec. 4.132 — coeficiente de descarga (2/64 ≤ d ≤ 32/64); 0.828 arriba de 32/64."""
    if d64 > 32.0:
        return 0.828
    return (0.5489239 + 3.720401e-2 * d64 - 1.603191e-3 * d64**2
            + 2.387117e-5 * d64**3 - 4.942371e-8 * d64**4)


def diametro_cook_dotterweich(q: float, pth: float, pe: float, gas: Gas, t_f: float, k: float,
                              tol: float = 1e-6) -> dict:
    """Ecs. 4.131-4.133 — diámetro en 64avos, iterando Cd(d). No usa Z (gas ideal)."""
    r = max(pe / pth, relacion_critica(k))
    h = funcion_h(r, k)
    d64 = 32.0
    for i in range(200):
        cd = cd_cook_dotterweich(d64)
        area = q / (155500.0 * cd * pth * math.sqrt(64.34 * k * h / (gas.gamma * f_a_r(t_f) * (k - 1.0))))
        d_n = area_a_64avos(area)
        if abs(d_n - d64) < tol:
            return {"d64": d_n, "area": area, "Cd": cd, "iteraciones": i + 1}
        d64 = d_n
    raise RuntimeError("no convergió")


# ------------------------------------------------- 4.9.4 válvulas de seguridad

def dp_valvula_seguridad(q: float, p1: float, t1_f: float, gas: Gas, d_tp: float, d_v: float,
                         cd: float = 0.9, y: float = 0.85, z: float | None = None) -> dict:
    """Caída de presión en una válvula de seguridad subsuperficial en flujo subcrítico.

    Ecuación de orificio con factor de expansión Y (el API 14B propone Y = 0.85):
        Δp = (1 − β⁴) ṁ² / (2 gc ρ1 Cd² Y² A²),  β = dv/dtp
    Equivale a las Ecs. 4.135-4.137; con Cd ≈ 0.895 reproduce el Ejemplo 4.11.
    """
    z = z if z is not None else gas.z(p1, t1_f)
    beta = d_v / d_tp
    rho1 = gas.densidad(p1, t1_f)
    rho_b = 0.0764 * gas.gamma
    m = q * rho_b / 86400.0                                  # lbm/s
    area = math.pi * d_v**2 / 4.0 / 144.0                    # pie²
    dp = (1.0 - beta**4) * m**2 / (2.0 * GC * rho1 * cd**2 * y**2 * area**2) / 144.0
    return {"dp": dp, "beta": beta, "Z": z, "p2": p1 - dp}
