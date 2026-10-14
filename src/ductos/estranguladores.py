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
    """Ecs. 4.131-4.133 — diámetro en 64avos, iterando Cd(d). No usa Z (gas ideal).

    d64: diámetro geométrico (A = CdA/Cd).
    d64_libro: el Ejemplo 4.9 obtiene CdA = 0.1089 pg² y lo convierte directamente con la
    Ec. 4.133, es decir, reporta el diámetro de la *área efectiva* (23.83/64 → 24/64).
    """
    r = max(pe / pth, relacion_critica(k))
    h = funcion_h(r, k)
    d64 = 32.0
    for i in range(200):
        cd = cd_cook_dotterweich(d64)
        area = q / (155500.0 * cd * pth * math.sqrt(64.34 * k * h / (gas.gamma * f_a_r(t_f) * (k - 1.0))))
        d_n = area_a_64avos(area)
        if abs(d_n - d64) < tol:
            return {"d64": d_n, "area": area, "Cd": cd, "iteraciones": i + 1,
                    "d64_libro": area_a_64avos(cd * area)}
        d64 = d_n
    raise RuntimeError("no convergió")


def diametro_ec_4_134(q: float, p1: float, t_f: float, gamma: float) -> float:
    """Ec. 4.134 — estranguladores de 6 pg en flujo crítico (ajuste de campo).

    q = 99.188 dφ² p1 / ((T+460) γ)^0.5,   dφ en 64avos de pg.
    """
    return math.sqrt(q * math.sqrt(f_a_r(t_f) * gamma) / (99.188 * p1))


def diametro_ec_4_129(q: float, pth: float, pe: float, gas: Gas, t_f: float, k: float,
                      z: float | None = None, cd: float = 1.0) -> float:
    """Ec. 4.129 resuelta para dφ [64avos de pg] con la constante deducida.

    La Ec. 4.129 impresa trae Cd = 1.9083×10⁻⁴, y el Ejemplo 4.9 lo multiplica en lugar de
    dividir: obtiene dφ² = 25 y lo reporta como dφ = 25/64. Con la constante física,
    q = 238.4 dφ² pth √[k H/((k−1) γ Z T)], el resultado es 23.4/64 (ver docs/erratas.md).
    """
    return area_a_64avos(area_estrangulador(q, pth, pe, gas, t_f, k, cd, z))


# ------------------------------------------------- k = Cp/Cv (Ecs. 4.116-4.128)

def k_gas(gamma: float, t_f: float, advertir: bool = True) -> float:
    """Relación de calores específicos del gas natural (ajuste polinomial de la Fig. 4.16).

    x = 10 γ − 7.1, Y = 0.01 T − 1.5. El polinomio de grado 9 se comporta bien para
    0.65 ≤ γ ≤ 0.80 y 50 ≤ T ≤ 250 °F; fuera de ese intervalo pierde monotonía y diverge (da
    k > 1.4 para γ = 0.9 a 50 °F), así que conviene usar un valor medido.
    Coeficientes tomados del libro impreso; b5 usa 0.006673 Y⁴ (la impresión dice 0.06673,
    que hace divergir la correlación; la tesis de 1990 trae 0.006673).
    """
    if advertir and not (0.65 <= gamma <= 0.80 and 50 <= t_f <= 250):
        import warnings
        warnings.warn("k_gas fuera del intervalo confiable (0.65 ≤ γ ≤ 0.80, 50-250 °F)")
    x = 10.0 * gamma - 7.1
    Y = 0.01 * t_f - 1.5
    b = [
        1.245874 - 0.027331 * Y - 0.017771 * Y**2 + 0.00305 * Y**3 + 0.013167 * Y**4,
        -0.027336 + 0.001484 * Y - 0.015829 * Y**2 - 0.0033 * Y**3 + 0.015654 * Y**4,
        0.002485 - 0.008877 * Y + 0.020643 * Y**2 + 0.009488 * Y**3 - 0.021162 * Y**4,
        0.002334 - 0.007175 * Y + 0.025258 * Y**2 + 0.009191 * Y**3 - 0.024669 * Y**4,
        0.000717 + 0.006973 * Y - 0.017431 * Y**2 - 0.007491 * Y**3 + 0.017444 * Y**4,
        0.002823 - 0.004985 * Y**2 - 0.004656 * Y**3 + 0.006673 * Y**4,
        -0.001784 * Y + 0.003282 * Y**2 + 0.002643 * Y**3 - 0.004141 * Y**4,
        0.000029 - 0.000056 * Y - 0.000191 * Y**2 + 0.000223 * Y**3,
        -0.000041 + 0.00008 * Y + 0.000275 * Y**2 - 0.000321 * Y**3,
        0.000028 + 0.000056 * Y,
    ]
    return float(sum(bi * x**i for i, bi in enumerate(b)))


# ------------------------------------------------- 4.9.4 válvulas de seguridad

K_SSSV = 69.22e-8   # constante de la Ec. 4.135 que reproduce el Ejemplo 4.11 (impresa: 69.22×10⁻⁶)


def factor_expansion(p1: float, p2: float, psi: float, k: float) -> float:
    """Ec. 4.137 — Y = 1 − (0.41 + 0.35 ψ⁴)(p1 − p2)/(k p1)."""
    return 1.0 - (0.41 + 0.35 * psi**4) * (p1 - p2) / (k * p1)


def dp_valvula_seguridad(q: float, p1: float, t1_f: float, gas: Gas, d_tp: float, d_v: float,
                         y: float | None = 0.85, k: float | None = None, z: float | None = None,
                         con_psi4: bool = True, tol: float = 0.01) -> dict:
    """Ec. 4.135 (API 14B) — caída de presión en una válvula de seguridad subsuperficial.

        p1 − p2 = 2.7 γ p1 / (Z1 T1) · (1 − ψ⁴) · [K Z1 T1 q / (p1 dv² Y)]²,  ψ = dv/d

    y: factor de expansión fijo (el API propone 0.85). Con y=None y un valor de k, se
       itera la Ec. 4.137 como describe el libro (tolerancia 0.01 lb/pg²).
    con_psi4: el Ejemplo 4.11 omite el factor (1 − ψ⁴) al sustituir (907.74 lb/pg²);
       con el factor el resultado es 896.4 lb/pg².
    K = 69.22×10⁻⁸: la impresión dice 10⁻⁶, pero el ejemplo usa en el corchete 12.51, que
    corresponde a 10⁻⁸; equivale a un orificio con Cd ≈ 0.90.
    """
    z = z if z is not None else gas.z(p1, t1_f)
    psi = d_v / d_tp
    t1 = f_a_r(t1_f)
    pre = 2.7 * gas.gamma * p1 / (z * t1) * ((1.0 - psi**4) if con_psi4 else 1.0)

    def dp_con(yy):
        return pre * (K_SSSV * z * t1 * q / (p1 * d_v**2 * yy)) ** 2

    if y is not None:
        dp = dp_con(y)
        return {"dp": dp, "psi": psi, "Y": y, "Z": z, "p2": p1 - dp, "iteraciones": 1}
    if k is None:
        raise ValueError("para iterar Y indique k")
    p2 = 0.9 * p1
    for i in range(200):
        yy = factor_expansion(p1, p2, psi, k)
        p2_n = p1 - dp_con(yy)
        if abs(p2_n - p2) < tol:
            return {"dp": p1 - p2_n, "psi": psi, "Y": yy, "Z": z, "p2": p2_n, "iteraciones": i + 1}
        p2 = 0.5 * (p2 + p2_n)
    raise RuntimeError("no convergió")
