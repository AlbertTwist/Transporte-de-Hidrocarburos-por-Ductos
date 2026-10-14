"""Cap. 4 — Flujo de gas por tuberías (secciones 4.1 a 4.4).

Ecuación general (Ec. 4.19/4.20, flujo permanente, una fase, energía cinética
despreciable, f de Darcy):

    p1² − p2² · e^s = f q² (p_b/T_b)² γ Z̄ T̄ L_e / (C² d⁵)

    C = 77.556 (derivada de Darcy-Weisbach, ver constantes.c_general_gas)
    s = 0.0375 γ h / (Z̄ T̄)              (Ec. 4.23)
    L_e = L (e^s − 1)/s                  (Ec. 4.24)

Unidades: p [psia], T [°F], q [pie³/d @ c.s.], L [millas], d [pg], h [pies].
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from . import friccion as fr
from .constantes import c_general_gas, f_a_r
from .pvt_gas import Gas, integral_ppr_sobre_z

C_GEN = c_general_gas()      # ≈ 77.556
P_BASE = 14.7                # lb/pg² abs (el libro usa 14.7 en el Cap. 4)
T_BASE = 60.0                # °F


def presion_media(p1: float, p2: float) -> float:
    """Ec. 4.46 — presión media en un gasoducto."""
    if abs(p1 - p2) < 1e-9:
        return p1
    return 2.0 / 3.0 * (p1**3 - p2**3) / (p1**2 - p2**2)


def reynolds(q: float, gamma: float, d: float, mu: float) -> float:
    """Ec. 4.33 — NRe = 0.0201056 q γ / (d μ)."""
    return 0.0201056 * q * gamma / (d * mu)


def reynolds_anular(q: float, gamma: float, d_ci: float, d_te: float, mu: float) -> float:
    """Ec. 4.44 — NRe en espacio anular."""
    return 0.0201056 * q * gamma / (mu * (d_ci + d_te))


def diametro_equivalente_anular(d_ci: float, d_te: float) -> float:
    """Ec. 4.42: d⁵ equivalente = (dci − dte)³ (dci + dte)²; devuelve d."""
    return ((d_ci - d_te) ** 3 * (d_ci + d_te) ** 2) ** 0.2


def factor_s(gamma: float, h: float, z: float, t_f: float) -> float:
    """Ec. 4.23."""
    return 0.0375 * gamma * h / (z * f_a_r(t_f))


def longitud_efectiva(L: float, s: float) -> float:
    """Ec. 4.24 (pendiente constante)."""
    if abs(s) < 1e-12:
        return L
    return L * (math.exp(s) - 1.0) / s


# ------------------------------------------------------------- columna estática

def presion_fondo_estatica(p_th: float, h: float, t_f: float, gas: Gas,
                           media: str = "aritmetica", tol: float = 1e-6) -> dict:
    """Ec. 4.18 — presión de fondo de un pozo de gas cerrado (Ejemplo 4.1).

    media: "aritmetica" (algoritmo del texto, p̄ = (p1+p2)/2) o
           "cubica" (Ec. 4.46). Ambas se muestran en el notebook.
    """
    p = p_th * 1.1
    for i in range(200):
        pm = (p_th + p) / 2.0 if media == "aritmetica" else presion_media(p, p_th)
        z = gas.z(pm, t_f)
        p_n = math.sqrt(p_th**2 + 0.03756 * gas.gamma * pm**2 * h / (z * f_a_r(t_f)))
        if abs(p_n - p) < tol:
            return {"pws": p_n, "Z": z, "p_media": pm, "iteraciones": i + 1}
        p = p_n
    raise RuntimeError("no convergió")


# ------------------------------------------------------------- ecuaciones de flujo

@dataclass
class Gasoducto:
    d: float                 # diámetro interno [pg]
    L: float                 # longitud [millas]
    t_f: float               # temperatura media de flujo [°F]
    h: float = 0.0           # desnivel salida − entrada [pies]
    rug: float = 0.0006      # rugosidad [pg]
    E: float = 1.0           # eficiencia
    p_b: float = P_BASE
    t_b: float = T_BASE

    def f_rugoso(self) -> float:
        """Ec. 4.45 — el libro usa flujo francamente turbulento por defecto."""
        return fr.f_rugoso(self.rug / self.d)


def _z_medio(gas: Gas, p1: float, p2: float, t_f: float) -> float:
    return gas.z(presion_media(p1, p2), t_f)


def _terminos_inclinacion(gas: Gas, ducto: Gasoducto, z: float) -> tuple[float, float]:
    """Devuelve (e^s, L_e)."""
    s = factor_s(gas.gamma, ducto.h, z, ducto.t_f)
    return math.exp(s), longitud_efectiva(ducto.L, s)


def q_general(p1: float, p2: float, gas: Gas, ducto: Gasoducto,
              f: float | None = None, c: float = C_GEN) -> float:
    """Ec. 4.20 (horizontal o inclinada) — gasto [pie³/d @ c.s.]."""
    f = f if f is not None else ducto.f_rugoso()
    z = _z_medio(gas, p1, p2, ducto.t_f)
    es, le = _terminos_inclinacion(gas, ducto, z)
    t = f_a_r(ducto.t_f)
    tb = f_a_r(ducto.t_b)
    return (ducto.E * c * (tb / ducto.p_b)
            * math.sqrt((p1**2 - es * p2**2) * ducto.d**5 / (gas.gamma * z * t * le * f)))


def p2_general(p1: float, q: float, gas: Gas, ducto: Gasoducto,
               f: float | None = None, c: float = C_GEN, tol: float = 1e-6) -> dict:
    """Ec. 4.19 resuelta para p2 por aproximaciones sucesivas (Ejemplo 4.2)."""
    f = f if f is not None else ducto.f_rugoso()
    t = f_a_r(ducto.t_f)
    tb = f_a_r(ducto.t_b)
    p2 = 0.8 * p1
    for i in range(200):
        z = _z_medio(gas, p1, p2, ducto.t_f)
        es, le = _terminos_inclinacion(gas, ducto, z)
        dp2 = f * (q / (ducto.E * c) * ducto.p_b / tb) ** 2 * gas.gamma * z * t * le / ducto.d**5
        arg = (p1**2 - dp2) / es
        if arg <= 0:
            raise ValueError("el gasto excede la capacidad de la línea")
        p2_n = math.sqrt(arg)
        if abs(p2_n - p2) < tol:
            return {"p2": p2_n, "Z": z, "p_media": presion_media(p1, p2_n), "f": f,
                    "iteraciones": i + 1}
        p2 = p2_n
    raise RuntimeError("no convergió")


def q_weymouth(p1, p2, gas: Gas, ducto: Gasoducto) -> float:
    """Tabla 4.1/4.2 — Weymouth: f = 0.032/d^(1/3) (Darcy); constante 433.5."""
    f = 0.032 / ducto.d ** (1.0 / 3.0)
    return q_general(p1, p2, gas, ducto, f=f)


def q_panhandle_a(p1, p2, gas: Gas, ducto: Gasoducto) -> float:
    """Panhandle A — líneas troncales de gran longitud."""
    z = _z_medio(gas, p1, p2, ducto.t_f)
    es, le = _terminos_inclinacion(gas, ducto, z)
    t = f_a_r(ducto.t_f)
    tb = f_a_r(ducto.t_b)
    return (ducto.E * 435.87 * (tb / ducto.p_b) ** 1.0788
            * ((p1**2 - es * p2**2) / le) ** 0.5394
            * (1.0 / gas.gamma) ** 0.4606 * (1.0 / (z * t)) ** 0.5394 * ducto.d ** 2.6182)


def q_panhandle_b(p1, p2, gas: Gas, ducto: Gasoducto, variante: str = "estandar") -> float:
    """Panhandle B — líneas troncales.

    variante="estandar": forma de la literatura (GPSA), γ^0.961 dentro del exponente 0.51.
    variante="libro":    constante K4 de la Tabla 4.3, con (1/γ)^0.961 fuera del
                         corchete. Reproduce el valor del Ejemplo 4.2 (1721 MMpcd),
                         pero sobrestima el gasto ~27 % respecto a la forma estándar
                         (ver docs/erratas.md).
    """
    z = _z_medio(gas, p1, p2, ducto.t_f)
    es, le = _terminos_inclinacion(gas, ducto, z)
    t = f_a_r(ducto.t_f)
    tb = f_a_r(ducto.t_b)
    base = ducto.E * 737.0 * (tb / ducto.p_b) ** 1.02 * ducto.d ** 2.53
    if variante == "estandar":
        return base * ((p1**2 - es * p2**2) / (gas.gamma**0.961 * z * t * le)) ** 0.51
    if variante == "libro":
        return base * ((p1**2 - es * p2**2) / (z * t * le)) ** 0.51 * (1.0 / gas.gamma) ** 0.961
    raise ValueError("variante debe ser 'estandar' o 'libro'")


def q_smith(p1, p2, gas: Gas, ducto: Gasoducto, f: float | None = None) -> float:
    """Smith — ecuación general con la constante 77.831 de la Tabla 4.3."""
    return q_general(p1, p2, gas, ducto, f=f, c=77.831)


def q_clinedinst(p1, p2, gas: Gas, ducto: Gasoducto, f: float | None = None) -> float:
    """Clinedinst — usa ∫ ppr/Z dppr en lugar de un Z medio (horizontal)."""
    f = f if f is not None else ducto.f_rugoso()
    tpc, ppc = gas.criticas()
    tpr = f_a_r(ducto.t_f) / tpc
    i1 = integral_ppr_sobre_z(p1 / ppc, tpr)
    i2 = integral_ppr_sobre_z(p2 / ppc, tpr)
    t = f_a_r(ducto.t_f)
    tb = f_a_r(ducto.t_b)
    k = math.sqrt(2.0) * C_GEN * tb * ppc / ducto.p_b * math.sqrt(1.0 / (gas.gamma * ducto.L * t * f))
    return ducto.E * k * math.sqrt(i1 - i2) * ducto.d ** 2.5


ECUACIONES = {
    "General (4.20)": q_general,
    "Panhandle A": q_panhandle_a,
    "Panhandle B": q_panhandle_b,
    "Weymouth": q_weymouth,
    "Smith": q_smith,
    "Clinedinst": q_clinedinst,
}
# Pendiente: Cullender-Smith (Apéndice D del libro).


# ------------------------------------------------- 4.4 sistemas complejos

def longitud_equivalente_serie(tramos: list[tuple[float, float, float]],
                               d_ref: float, f_ref: float) -> float:
    """Tuberías en serie: L_E con diámetro de referencia.

    tramos: [(L_i [mi], d_i [pg], f_i), ...]
    Igualando las pérdidas p1²−p2² ∝ f L / d⁵:
        L_E = (d_ref⁵ / f_ref) Σ L_i f_i / d_i⁵
    (ver docs/erratas.md: la fuente usa exponente 2.5 en el diámetro y 1 en f).
    """
    return d_ref**5 / f_ref * sum(L * f / d**5 for L, d, f in tramos)


def diametro_equivalente_paralelo(ramas: list[tuple[float, float]], rug: float = 0.0006,
                                  tol: float = 1e-8) -> tuple[float, float]:
    """Ec. 4.49 — diámetro equivalente de n líneas paralelas de igual longitud.

    ramas: [(d_i [pg], f_i), ...]. Como q ∝ d^2.5/√f,
        d_E^2.5/√f_E = Σ d_i^2.5/√f_i,  con f_E = f_rugoso(ε/d_E).
    Devuelve (d_E, f_E).
    """
    suma = sum(d**2.5 / math.sqrt(f) for d, f in ramas)
    d = max(d for d, _ in ramas)
    for _ in range(200):
        f_e = fr.f_rugoso(rug / d)
        d_n = (suma * math.sqrt(f_e)) ** 0.4
        if abs(d_n - d) < tol:
            return d_n, fr.f_rugoso(rug / d_n)
        d = d_n
    raise RuntimeError("no convergió")


def fraccion_loop(q1: float, q2: float, d_orig: float, d_paralela: float,
                  f_orig: float, f_paralela: float) -> float:
    """Ecs. 4.50-4.51 — fracción x de la línea que debe duplicarse (loop)
    para pasar de q1 a q2 con la misma caída de presión."""
    r = math.sqrt(f_orig / f_paralela) * (d_paralela / d_orig) ** 2.5
    return (1.0 - (q1 / q2) ** 2) / (1.0 - 1.0 / (1.0 + r) ** 2)
