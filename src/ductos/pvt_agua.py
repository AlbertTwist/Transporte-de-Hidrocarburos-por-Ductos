"""Cap. 1.5-1.6 — Propiedades del agua de formación, saturada y bajosaturada.

Unidades: p [lb/pg² abs], T [°F], salinidad en % de NaCl equivalente, Rsw [pie³/bl],
ρ [lbm/pie³], μ [cp], σ [dina/cm], cw [1/(lb/pg²)].
"""

from __future__ import annotations

import math


def bw_saturada(p: float, t_f: float, variante: str = "correlacion") -> float:
    """Ec. 1.69 — factor de volumen del agua saturada.

    variante="correlacion": Bw = 1 + 1.2×10⁻⁴ (T−60) + 1.0×10⁻⁶ (T−60)² − 3.33×10⁻⁶ p
    variante="ejemplo":     sin el cuadrado, como se sustituye en el Ejemplo 1.1 (Bw = 1.0103)
    Con el cuadrado, el agua a 200 °F se expande ~3 %, que es lo físicamente esperado
    (ver docs/erratas.md).
    """
    dt = t_f - 60.0
    termino = dt**2 if variante == "correlacion" else dt
    if variante not in ("correlacion", "ejemplo"):
        raise ValueError("variante debe ser 'correlacion' o 'ejemplo'")
    return 1.0 + 1.2e-4 * dt + 1.0e-6 * termino - 3.33e-6 * p


def densidad_agua_saturada(bw: float, gamma_w: float = 1.0) -> float:
    """Ec. 1.70 [lbm/pie³]."""
    return 62.43 * gamma_w / bw


def viscosidad_agua_saturada(t_f: float, nacl: float = 0.0, p: float | None = None) -> float:
    """Ecs. 1.71-1.74 — μw = A + B/T con corrección por presión si se da p.

    Intervalo: 60 < T < 400 °F, p ≤ 10 000 lb/pg², NaCl < 26 %.
    """
    A = -0.04518 + 0.009313 * nacl - 0.000393 * nacl**2
    B = 70.634 + 0.09576 * nacl**2
    mu = A + B / t_f
    if p is not None:
        mu *= 1.0 + 3.5e-12 * p**2 * (t_f - 40.0)
    return mu


def tension_superficial_agua_gas(p: float, t_f: float) -> float:
    """Ecs. 1.75-1.77 — interpolación lineal entre 74 °F y 280 °F."""
    s280 = 52.5 - 0.006 * p
    s74 = 76.0 * math.exp(-0.00025 * p)
    return (280.0 - t_f) / 206.0 * (s74 - s280) + s280


def rsw(p: float, t_f: float, nacl: float = 0.0) -> float:
    """Ecs. 1.78-1.83 — solubilidad del gas en el agua [pie³/bl], corregida por salinidad."""
    pp = 1.0 - math.exp(-p / 2276.0)
    tc = 5.0 / 9.0 * (t_f - 32.0)
    tt = (tc - 90.0) / 10.0
    S = pp * (3.69051 + 0.08746 * tt + 0.01129 * tt**2 - 0.00647 * tt**3)
    rsw_pura = 5.6146 * S
    cs = 1.0 + (0.0001736 * t_f - 0.07703) * nacl
    return rsw_pura * cs


def compresibilidad_agua(p: float, t_f: float, rsw_: float = 0.0) -> float:
    """Ecs. 1.84-1.88 — con corrección por gas en solución.

    Intervalo: 1000 < p < 6000 lb/pg² abs, 80 < T < 250 °F, 0 < Rsw < 25 pie³/bl.
    """
    A = 3.8546 - 0.000134 * p
    B = -0.01052 + 4.77e-7 * p
    C = 3.9267e-5 - 8.8e-10 * p
    return (A + B * t_f + C * t_f**2) * 1e-6 * (1.0 + 8.9e-3 * rsw_)


def bw_bajosaturada(bwb: float, cw: float, p: float, pb: float) -> float:
    """Ec. 1.90 — Bw = Bwb exp[cw (pb − p)]."""
    return bwb * math.exp(cw * (pb - p))


def densidad_agua_bajosaturada(rho_wb: float, cw: float, p: float, pb: float) -> float:
    """Ec. 1.89 — ρw = ρwb exp[cw (p − pb)]."""
    return rho_wb * math.exp(cw * (p - pb))
