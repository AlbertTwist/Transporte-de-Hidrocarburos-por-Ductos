"""Cap. 4.6–4.7 — Ecuación general de costos y diseño económico de gasoductos.

Costo total anual por milla (Ec. 4.80):

    CT = K2 · P(pd) · (A d² + B d³) + Ft  +  K3 · q a Z / (10⁶ L)
         └──────── tubería ────────┘        └─── compresión ───┘

    P(p) = p (2S' − p) / [4 (S' − p)²],  S' = S·F·E·T           (de Ecs. 4.58-4.63)
    K2   = 28.2 · f_i                                           (Ec. 4.81)
    K3   = g_r X + 87.6 C2 + Fe X                               (Ec. 4.82)
    a    = K1/(b e) (Rc^b − 1) + he   [BHP por MMpcd]           (Ec. 4.69)

El diseño económico (método de multiplicadores de Lagrange, sección 4.7.3) busca
d, L, pd y pe que minimizan CT sujetos a la ecuación de flujo (Ec. 4.20).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from scipy.optimize import brentq, minimize_scalar

from . import friccion as fr
from .constantes import c_general_gas, f_a_r
from .pvt_gas import Gas


def factor_anualidad(i: float, n: float) -> float:
    """Ecs. 4.66 y 4.75 — factor de recuperación de capital."""
    return i * (1.0 + i) ** n / ((1.0 + i) ** n - 1.0)


def espesor(p: float, d: float, s_prima: float) -> float:
    """Ec. 4.62 (ASME B31) — espesor en función del diámetro interno [pg]."""
    return p * d / (2.0 * (s_prima - p))


def toneladas_por_milla(d: float, t: float) -> float:
    """Ec. 4.58 — peso de tubería de acero [ton/milla], d interno y t en pg."""
    return 28.2 * (d + t) * t


@dataclass
class DatosCostosGas:
    """Datos del Ejemplo 4.7 por defecto."""
    gas: Gas = field(default_factory=lambda: Gas(0.6))
    k: float = 1.28           # Cp/Cv
    he: float = 1.5           # pérdidas en válvulas [BHP/MMpcd]
    e: float = 0.9            # eficiencia mecánica adiabática
    ts: float = 80.0          # temperatura de succión [°F]
    t0: float = 60.0
    p0: float = 14.7
    t_flujo: float = 70.0     # temperatura de flujo en la línea [°F]
    y0: float = 220.0         # costo base de tubería no instalada [dll/ton]
    g0: float = 180.0         # costo base de colocación [dll/ton]
    im: float = 4.583         # incremento de Y por pg [dll/ton/pg]
    in_: float = 3.75         # incremento de G por pg [dll/ton/pg]
    s: float = 60000.0        # esfuerzo máximo permisible [lb/pg²]
    f_diseno: float = 0.72    # F: tipo de construcción
    e_union: float = 1.0      # E: unión longitudinal
    t_termico: float = 1.0    # T: deformación térmica
    i: float = 0.15           # tasa de interés (tubería)
    n: float = 15.0           # vida útil de la tubería [años]
    r: float = 0.15           # tasa de interés (compresión)
    m: float = 15.0           # vida útil de compresoras [años]
    x: float = 430.0          # costo de compresoras instaladas [dll/BHP]
    c2: float = 0.5           # valor del gas combustible [dll/1000 pie³]
    fe: float = 0.02          # costo fijo anual de mantenimiento de compresión (fracción de X)
    ft: float = 0.0           # costo fijo anual de la tubería [dll/milla/año]
    rug: float = 0.0006

    # --- coeficientes derivados
    @property
    def s_prima(self) -> float:
        return self.s * self.f_diseno * self.e_union * self.t_termico

    @property
    def A(self) -> float:
        return self.y0 + self.g0

    @property
    def B(self) -> float:
        return self.im + self.in_

    @property
    def fi(self) -> float:
        return factor_anualidad(self.i, self.n)

    @property
    def gr(self) -> float:
        return factor_anualidad(self.r, self.m)

    @property
    def K1(self) -> float:
        """Ec. 4.71."""
        return 3.03 * self.p0 * f_a_r(self.ts) / f_a_r(self.t0)

    @property
    def b(self) -> float:
        """Ec. 4.72."""
        return (self.k - 1.0) / self.k

    @property
    def K2(self) -> float:
        return 28.2 * self.fi

    @property
    def K3(self) -> float:
        return self.gr * self.x + 87.6 * self.c2 + self.fe * self.x

    def hp_por_mmpcd(self, rc: float) -> float:
        """Ec. 4.69 — potencia por MMpcd para una razón de compresión."""
        return self.K1 / (self.b * self.e) * (rc**self.b - 1.0) + self.he

    def factor_presion(self, p: float) -> float:
        """P(p) = p(2S'−p)/(4(S'−p)²) — "dif" en el Ejemplo 4.7."""
        sp = self.s_prima
        return p * (2.0 * sp - p) / (4.0 * (sp - p) ** 2)


def razon_compresion_optima(datos: DatosCostosGas) -> float:
    """Ec. 4.95 — razón de compresión que minimiza el costo de compresión por milla.

    Con d fijo, L ∝ pd² (1 − Rc⁻²), así que la compresión cuesta ∝ a(Rc)/(1 − Rc⁻²);
    se minimiza numéricamente (el libro lo hace por aproximaciones sucesivas).
    """
    g = lambda rc: datos.hp_por_mmpcd(rc) / (1.0 - rc**-2)
    return float(minimize_scalar(g, bounds=(1.001, 5.0), method="bounded", options={"xatol": 1e-9}).x)


def k4(d: float, z: float, datos: DatosCostosGas, c: float | None = None) -> float:
    """Constante de la Ec. 4.20 (tabla 4.3) con f de la Ec. 2.18."""
    c = c if c is not None else c_general_gas()
    f = fr.f_rugoso(datos.rug / d)
    t = f_a_r(datos.t_flujo)
    return c * f_a_r(datos.t0) / datos.p0 * math.sqrt(1.0 / (datos.gas.gamma * z * t * f))


def longitud_entre_estaciones(q: float, d: float, pd: float, pe: float, z: float,
                              datos: DatosCostosGas) -> float:
    """Ec. 4.100 — L = (K4/q)² (pd² − pe²) d⁵  [millas]."""
    return (k4(d, z, datos) / q) ** 2 * (pd**2 - pe**2) * d**5


def costo_total(q: float, d: float, L: float, pd: float, rc: float, z: float,
                datos: DatosCostosGas) -> dict:
    """Ec. 4.80 — costo total anual por milla y sus dos componentes."""
    tub = datos.K2 * datos.factor_presion(pd) * (datos.A * d**2 + datos.B * d**3) + datos.ft
    comp = datos.K3 * q * datos.hp_por_mmpcd(rc) * z / (1e6 * L)
    return {"tuberia": tub, "compresion": comp, "total": tub + comp}


def diseno_economico(q: float, pd: float, datos: DatosCostosGas | None = None,
                     rc: float | None = None, d_min: float = 1.0, d_max: float = 80.0) -> dict:
    """Pasos 4 a 7 de la sección 4.7.4 para una presión de descarga dada.

    El diámetro óptimo cumple ∂CT/∂d = 0 con L ligado a d por la ecuación de flujo
    (equivale a las Ecs. 4.96-4.98, d⁸ + α d⁷ = β con α = 2A/(3B)); como en el libro,
    f se evalúa en cada diámetro pero no se deriva.
    """
    datos = datos or DatosCostosGas()
    rc = rc if rc is not None else razon_compresion_optima(datos)
    pe = pd / rc
    z = datos.gas.z(pe, datos.t_flujo)
    a = datos.hp_por_mmpcd(rc)
    pf = datos.factor_presion(pd)

    def residuo(d):
        L = longitud_entre_estaciones(q, d, pd, pe, z, datos)
        lado_tub = datos.K2 * pf * (2.0 * datos.A * d + 3.0 * datos.B * d**2)
        lado_comp = 5.0 * datos.K3 * q * a * z / (1e6 * L * d)
        return lado_tub - lado_comp

    d = brentq(residuo, d_min, d_max, xtol=1e-8)
    L = longitud_entre_estaciones(q, d, pd, pe, z, datos)
    costos = costo_total(q, d, L, pd, rc, z, datos)
    return {"pd": pd, "pe": pe, "Rc": rc, "Z": z, "a": a, "dif": pf, "d": d,
            "f": fr.f_rugoso(datos.rug / d), "K4": k4(d, z, datos), "L": L,
            "t": espesor(pd, d, datos.s_prima), **costos}


def barrido_presion_descarga(q: float, pds, datos: DatosCostosGas | None = None) -> list[dict]:
    """Tabla 4.8 — diseño económico para varias presiones de descarga."""
    datos = datos or DatosCostosGas()
    return [diseno_economico(q, pd, datos) for pd in pds]


def presion_descarga_optima(q: float, datos: DatosCostosGas | None = None,
                            pd_min: float = 500.0, pd_max: float = 5000.0) -> dict:
    """Mínimo de CT(pd) por búsqueda acotada (alternativa a leerlo de la Fig. 4.8)."""
    datos = datos or DatosCostosGas()
    res = minimize_scalar(lambda pd: diseno_economico(q, pd, datos)["total"],
                          bounds=(pd_min, pd_max), method="bounded", options={"xatol": 0.1})
    return diseno_economico(q, float(res.x), datos)
