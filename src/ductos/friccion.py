"""Cap. 2 — Fundamentos de flujo a través de tuberías: factor de fricción.

Todas las funciones devuelven el factor de fricción de **Darcy** (f_Darcy = 4 f_Fanning),
que es el que usa el libro en las ecuaciones de Darcy (Ec. 2.10).

Regímenes según el libro:
    NRe < 2300          laminar           f = 64/NRe            (Ec. 2.16)
    2300 ≤ NRe ≤ 3100   crítico           interpolación         (Ec. 2.19)
    NRe > 3100          turbulento        Colebrook-White       (Ec. 2.17)
Para flujo francamente turbulento (rugosidad dominante): Ec. 2.18.
Aproximación explícita para turbulento: Ec. 2.20 (tipo Jain).
"""

from __future__ import annotations

import math

RE_LAMINAR = 2300.0
RE_TURBULENTO = 3100.0


def f_laminar(nre: float) -> float:
    """Ec. 2.16."""
    return 64.0 / nre


def f_colebrook(nre: float, rug_rel: float, f0: float = 0.03,
                tol: float = 1e-10, max_iter: int = 200) -> float:
    """Ec. 2.17 (Colebrook-White, constantes 3.715 y 2.514 del libro).

    Sustituciones sucesivas, igual que el procedimiento del Ejemplo 3.4.1(b).

    rug_rel : ε/d (adimensional)
    """
    f = f0
    for _ in range(max_iter):
        f_nuevo = (-2.0 * math.log10(rug_rel / 3.715 + 2.514 / (nre * math.sqrt(f)))) ** -2
        if abs(f_nuevo - f) < tol:
            return f_nuevo
        f = f_nuevo
    raise RuntimeError("Colebrook no convergió")


def f_colebrook_iteraciones(nre: float, rug_rel: float, f0: float = 0.03, n: int = 5) -> list[float]:
    """Devuelve la historia de iteraciones de la Ec. 2.17 (útil para clase)."""
    hist, f = [], f0
    for _ in range(n):
        f = (-2.0 * math.log10(rug_rel / 3.715 + 2.514 / (nre * math.sqrt(f)))) ** -2
        hist.append(f)
    return hist


def f_rugoso(rug_rel: float) -> float:
    """Ec. 2.18: zona francamente turbulenta, f independiente de NRe."""
    return (-2.0 * math.log10(rug_rel / 3.715)) ** -2


def f_jain(nre: float, rug_rel: float) -> float:
    """Ec. 2.20: aproximación explícita para NRe > 3100."""
    return (1.14 - 2.0 * math.log10(rug_rel + 21.25 / nre**0.9)) ** -2


def f_critico(nre: float, rug_rel: float) -> float:
    """Ec. 2.19: zona crítica 2300 < NRe < 3100 (tal como la presenta el libro).

    f = (NRe - 2300)/2300 · 1.3521 / [ln(ε/(3.715 d) + 2.514/(3100 √f_s))]² + 0.032

    donde f_s es el factor de Colebrook evaluado en NRe = 3100.
    Nota: la expresión no es continua con 64/NRe en NRe = 2300 (da 0.032 vs 0.0278);
    se reproduce como aparece en la fuente.
    """
    fs = f_colebrook(RE_TURBULENTO, rug_rel)
    x = math.log(rug_rel / 3.715 + 2.514 / (RE_TURBULENTO * math.sqrt(fs)))
    return (nre - RE_LAMINAR) / RE_LAMINAR * 1.3521 / x**2 + 0.032


def factor_friccion(nre: float, rug_rel: float, metodo: str = "colebrook") -> float:
    """Selecciona el factor de fricción según el régimen (criterio del libro).

    metodo (solo para turbulento): "colebrook" | "jain" | "rugoso"
    """
    if nre <= 0:
        raise ValueError("NRe debe ser positivo")
    if nre < RE_LAMINAR:
        return f_laminar(nre)
    if nre <= RE_TURBULENTO:
        return f_critico(nre, rug_rel)
    if metodo == "colebrook":
        return f_colebrook(nre, rug_rel)
    if metodo == "jain":
        return f_jain(nre, rug_rel)
    if metodo == "rugoso":
        return f_rugoso(rug_rel)
    raise ValueError(f"método desconocido: {metodo}")


def regimen(nre: float) -> str:
    if nre < RE_LAMINAR:
        return "laminar"
    if nre <= RE_TURBULENTO:
        return "crítico"
    return "turbulento"
