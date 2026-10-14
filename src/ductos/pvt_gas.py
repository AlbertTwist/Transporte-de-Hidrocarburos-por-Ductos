"""Cap. 1.2 — Propiedades del gas natural.

Incluye: propiedades pseudocríticas (Standing), corrección por gases amargos
(Wichert-Aziz), factor Z (Dranchuk-Purvis-Robinson, 8 constantes, Ecs. 1.16-1.19),
factor de volumen, densidad, viscosidad (Lee-González-Eakin) y su corrección
por N2, CO2 y H2S (Standing).

Unidades: p [lb/pg² abs], T [°F], densidad [lbm/pie³], viscosidad [cp].
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from scipy.integrate import quad
from scipy.optimize import brentq

from .constantes import f_a_r

# Constantes de Dranchuk-Purvis-Robinson (Ec. 1.19)
A1, A2, A3, A4 = 0.31506237, -1.04670990, -0.57832729, 0.53530771
A5, A6, A7, A8 = -0.61232032, -0.10488813, 0.68157001, 0.68446549


def pseudocriticas(gamma_g: float, tipo: str = "superficial") -> tuple[float, float]:
    """Propiedades pseudocríticas de Standing.

    tipo = "superficial"  Ecs. 1.12-1.13 (gas de separador / tanque)
    tipo = "humedo"       Ecs. 1.14-1.15 (gas en equilibrio con aceite, "condensados")

    Devuelve (Tpc [°R], ppc [psia]).
    """
    if tipo == "superficial":
        return 167.0 + 316.67 * gamma_g, 702.5 - 50.0 * gamma_g
    if tipo == "humedo":
        return 238.0 + 210.0 * gamma_g, 740.0 - 100.0 * gamma_g
    raise ValueError("tipo debe ser 'superficial' o 'humedo'")


def wichert_aziz(tpc: float, ppc: float, y_co2: float, y_h2s: float) -> tuple[float, float, float]:
    """Ecs. 1.20-1.22: corrección por CO2 y H2S.

    Devuelve (T'pc [°R], p'pc [psia], ε [°R]).
    """
    a = y_co2 + y_h2s
    b = y_h2s
    eps = 120.0 * (a**0.9 - a**1.6) + 15.0 * (b**0.5 - b**4)
    tpc_c = tpc - eps
    ppc_c = ppc * tpc_c / (tpc + b * (1.0 - b) * eps)
    return tpc_c, ppc_c, eps


def z_dpr(ppr: float, tpr: float, z0: float = 0.9, tol: float = 1e-8, max_iter: int = 500) -> float:
    """Factor de compresibilidad por Dranchuk-Purvis-Robinson (Ecs. 1.16-1.19).

    Se resuelve F(ρr) = Z(ρr)·ρr·Tpr − 0.27·ppr = 0 sobre la densidad reducida:
    primero Newton desde el gas ideal y, si no converge (cerca de Tpr ≈ 1),
    se acota la primera raíz física y se usa Brent. El resultado es el mismo
    que el de la sustitución directa del libro (ver z_dpr_sustitucion).
    """
    if ppr == 0:
        return 1.0

    def z_de(rho):
        return (1.0 + (A1 + A2 / tpr + A3 / tpr**3) * rho + (A4 + A5 / tpr) * rho**2
                + A5 * A6 * rho**5 / tpr
                + A7 * rho**2 / tpr**3 * (1.0 + A8 * rho**2) * math.exp(-A8 * rho**2))

    def F(rho):
        return z_de(rho) * rho * tpr - 0.27 * ppr

    rho = 0.27 * ppr / (z0 * tpr)
    for _ in range(50):
        h = 1e-7 * max(rho, 1e-6)
        dF = (F(rho + h) - F(rho)) / h
        if dF <= 0:
            break
        rho_n = rho - F(rho) / dF
        if rho_n <= 0:
            break
        if abs(rho_n - rho) < tol:
            return 0.27 * ppr / (rho_n * tpr)
        rho = rho_n

    malla = [1e-6 + i * 4.0 / 400 for i in range(401)]
    for a, b in zip(malla, malla[1:]):
        if F(a) * F(b) <= 0:
            rho = brentq(F, a, b, xtol=1e-12)
            return 0.27 * ppr / (rho * tpr)
    raise RuntimeError("DPR no convergió")


def z_dpr_sustitucion(ppr: float, tpr: float, z0: float = 0.8, tol: float = 1e-3,
                      max_iter: int = 200) -> tuple[float, int]:
    """Procedimiento literal del libro: suponer Z, calcular ρr, recalcular Z.

    Devuelve (Z, iteraciones). Tolerancia del libro: 0.001.
    """
    z = z0
    for i in range(1, max_iter + 1):
        rho = 0.27 * ppr / (z * tpr)
        z_c = (1.0 + (A1 + A2 / tpr + A3 / tpr**3) * rho + (A4 + A5 / tpr) * rho**2
               + A5 * A6 * rho**5 / tpr
               + A7 * rho**2 / tpr**3 * (1.0 + A8 * rho**2) * math.exp(-A8 * rho**2))
        if abs(z_c - z) <= tol:
            return z_c, i
        z = z_c
    raise RuntimeError("sustitución directa no convergió")


@dataclass
class Gas:
    """Gas natural caracterizado por su densidad relativa y contaminantes."""
    gamma: float
    y_co2: float = 0.0
    y_h2s: float = 0.0
    y_n2: float = 0.0
    tipo: str = "superficial"

    def criticas(self) -> tuple[float, float]:
        tpc, ppc = pseudocriticas(self.gamma, self.tipo)
        if self.y_co2 > 0 or self.y_h2s > 0:
            tpc, ppc, _ = wichert_aziz(tpc, ppc, self.y_co2, self.y_h2s)
        return tpc, ppc

    def reducidas(self, p: float, t_f: float) -> tuple[float, float]:
        tpc, ppc = self.criticas()
        return p / ppc, f_a_r(t_f) / tpc

    def z(self, p: float, t_f: float) -> float:
        ppr, tpr = self.reducidas(p, t_f)
        return z_dpr(ppr, tpr)

    def bg(self, p: float, t_f: float) -> float:
        """Ec. 1.9 — factor de volumen [pie³ @ c.esc. / pie³ @ c.s.]."""
        return 0.02825 * self.z(p, t_f) * f_a_r(t_f) / p

    def densidad(self, p: float, t_f: float) -> float:
        """Ec. 1.11 — [lbm/pie³]."""
        return 2.7044 * p * self.gamma / (self.z(p, t_f) * f_a_r(t_f))

    def viscosidad(self, p: float, t_f: float, corregir: bool = True,
                   gamma_correccion: float | None = None, c_x: float = 986.0) -> float:
        """Lee-González-Eakin (Ecs. 1.23-1.26) + corrección de Standing (1.27-1.30).

        c_x: coeficiente del término c_x/T en X. La correlación original usa 986;
        el Ejemplo 1.1 de la fuente usa 960 (ver docs/erratas.md). Con c_x=960 se
        reproduce μg = 0.017145 cp del libro.

        gamma_correccion: densidad relativa usada en las correcciones por
        contaminantes. En el Ejemplo 1.1 el libro usa la del gas producido
        (0.85) aunque la viscosidad se calcula con la del gas libre.
        """
        t = f_a_r(t_f)
        m = 28.97 * self.gamma
        k = (9.4 + 0.02 * m) * t**1.5 / (209.0 + 19.0 * m + t)
        x = 3.5 + c_x / t + 0.01 * m
        y = 2.4 - 0.2 * x
        rho_gcc = self.densidad(p, t_f) / 62.428
        mu = 1e-4 * k * math.exp(x * rho_gcc**y)
        if corregir:
            g = gamma_correccion if gamma_correccion is not None else self.gamma
            lg = math.log10(g)
            mu += self.y_n2 * (8.48e-3 * lg + 9.59e-3)
            mu += self.y_co2 * (9.08e-3 * lg + 6.24e-3)
            mu += self.y_h2s * (8.49e-3 * lg + 3.73e-3)
        return mu


def integral_ppr_sobre_z(ppr: float, tpr: float) -> float:
    """∫₀^ppr (ppr/Z) dppr, usada por la ecuación de Clinedinst (Tabla 4.1).

    El libro la toma de tablas impresas (ref. 6, Tabla A-6); aquí se integra
    numéricamente con Z de DPR.
    """
    val, _ = quad(lambda x: x / z_dpr(x, tpr), 0.0, ppr, limit=200)
    return val


# ------------------------------------------------- 4.5 presencia de condensados

def peso_molecular_condensado(api: float | None = None, gamma_o: float | None = None) -> float:
    """Ec. 4.53 — peso molecular del condensado en el tanque.

    Con °API:   Mo = 6084 / (°API − 5.9)
    Con γo:     Mo = 44.29 γo / (1.03 − γo)
    """
    if api is not None:
        return 6084.0 / (api - 5.9)
    if gamma_o is not None:
        return 44.29 * gamma_o / (1.03 - gamma_o)
    raise ValueError("indique api o gamma_o")


def gamma_razasa_katz(gamma_g: float, gamma_o: float, rgc: float, mo: float) -> float:
    """Ec. 4.52 — densidad relativa de la corriente gas + condensado.

    rgc: relación gas-condensado [pie³/bl]. Usar el resultado en lugar de γg para Z.
    """
    return (gamma_g + 4584.0 * gamma_o / rgc) / (1.0 + 132800.0 * gamma_o / (mo * rgc))


def gamma_vitter(gamma_g: float, gamma_l: float, rgl: float) -> float:
    """Ec. 4.54 — densidad relativa ponderada por agua y condensado (Vitter).

    rgl: relación gas-líquido [pie³/bl]; gamma_l: densidad relativa del líquido.
    El numerador no es legible en la tesis escaneada; se usa la forma de la
    literatura (Peffer, Miller y Hill, SPE 15655).
    """
    return (gamma_g + 4591.0 * gamma_l / rgl) / (1.0 + 1123.0 / rgl)


def equivalente_gas_condensado(gamma_o: float, mo: float) -> float:
    """Ec. 4.55 (Ikoku) — equivalente en gas del condensado [pie³ @ c.s./bl]."""
    return 133037.0 * gamma_o / mo


def gasto_total_gas(qg: float, qo: float, ge: float) -> float:
    """Ec. 4.56 — gasto de gas total incluyendo el equivalente del condensado."""
    return qg + ge * qo
