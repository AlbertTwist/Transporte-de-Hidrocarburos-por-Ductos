"""Cap. 3 — Flujo de líquidos por tuberías.

Ecuación general en unidades prácticas (Ec. 3.23):

    Δp_T = 0.433 γL Δh + C · f γL q² L / (E² d⁵)

    Δp [psi], γL [agua = 1], Δh [pies, + hacia arriba], q [bl/d], L [millas],
    d [pg], E = eficiencia de flujo (0-1), C = 0.06056 (libro).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from scipy.optimize import brentq

from . import friccion as fr
from .constantes import C_BHP, C_ELEV_LIQ, C_FRIC_LIQ_LIBRO, C_RE_LIQ


def reynolds(q: float, gamma: float, d: float, mu: float) -> float:
    """Ec. 3.22 — NRe = 92.2 q γ / (d μ)."""
    return C_RE_LIQ * q * gamma / (d * mu)


def factor_friccion(q: float, gamma: float, d: float, mu: float, rug: float,
                    metodo: str = "colebrook") -> float:
    """Factor de Darcy para las condiciones dadas (rug en pg)."""
    return fr.factor_friccion(reynolds(q, gamma, d, mu), rug / d, metodo)


def dp_elevacion(gamma: float, dh: float) -> float:
    """Ec. 3.2 — Δp por elevación [psi]."""
    return C_ELEV_LIQ * gamma * dh


def dp_friccion(q: float, gamma: float, d: float, L: float, f: float,
                E: float = 1.0, c: float = C_FRIC_LIQ_LIBRO) -> float:
    """Ec. 3.13 / 3.23 — Δp por fricción [psi]."""
    return c * f * gamma * q**2 * L / (E**2 * d**5)


@dataclass
class Linea:
    """Tramo de tubería que conduce líquido."""
    d: float            # diámetro interno [pg]
    L: float            # longitud [millas]
    dh: float = 0.0     # desnivel salida - entrada [pies]
    rug: float = 0.0006 # rugosidad absoluta [pg]
    E: float = 1.0      # eficiencia de flujo

    def dp(self, q: float, gamma: float, mu: float, metodo: str = "colebrook") -> dict:
        """Caída de presión total entrada → salida (positiva = pérdida)."""
        nre = reynolds(q, gamma, self.d, mu)
        f = fr.factor_friccion(nre, self.rug / self.d, metodo)
        dpf = dp_friccion(q, gamma, self.d, self.L, f, self.E)
        dpe = dp_elevacion(gamma, self.dh)
        return {"NRe": nre, "regimen": fr.regimen(nre), "f": f,
                "dp_f": dpf, "dp_e": dpe, "dp_T": dpf + dpe}


# ---------------------------------------------------------------- 3.4.2 Diseño

def diametro(q: float, gamma: float, mu: float, L: float, dp_T: float, dh: float = 0.0,
             rug: float = 0.0006, E: float = 1.0, d0: float = 4.0,
             tol: float = 1e-6, max_iter: int = 200) -> float:
    """Ec. 3.25 — diámetro para un gasto y caída de presión dados (sustituciones sucesivas)."""
    disp = dp_T - dp_elevacion(gamma, dh)
    if disp <= 0:
        raise ValueError("La caída disponible no alcanza a vencer el desnivel")
    d = d0
    for _ in range(max_iter):
        f = factor_friccion(q, gamma, d, mu, rug)
        d_n = (C_FRIC_LIQ_LIBRO * f * gamma * q**2 * L / (E**2 * disp)) ** 0.2
        if abs(d_n - d) < tol:
            return d_n
        d = d_n
    raise RuntimeError("diámetro no convergió")


def gasto(d: float, gamma: float, mu: float, L: float, dp_T: float, dh: float = 0.0,
          rug: float = 0.0006, E: float = 1.0, q0: float = 1000.0,
          tol: float = 1e-6, max_iter: int = 500) -> float:
    """Ec. 3.26 — gasto para un diámetro y caída de presión dados."""
    disp = dp_T - dp_elevacion(gamma, dh)
    if disp <= 0:
        return 0.0
    q = q0
    for _ in range(max_iter):
        f = factor_friccion(q, gamma, d, mu, rug)
        q_n = E * math.sqrt(d**5 * disp / (C_FRIC_LIQ_LIBRO * f * gamma * L))
        if abs(q_n - q) < tol * max(1.0, q):
            return q_n
        q = q_n
    raise RuntimeError("gasto no convergió")


# ------------------------------------------- 3.4.3-3.4.4 Pozo inyector (nodal)

@dataclass
class PozoInyector:
    """Sistema: bomba → línea superficial → tubería vertical → yacimiento."""
    linea: Linea
    tp: Linea                 # tubería vertical (dh negativo: el agua baja)
    gamma_linea: float
    mu_linea: float
    gamma_tp: float
    mu_tp: float
    ii: float                 # índice de inyectividad [bl/d/psi]

    def nodos(self, q: float, p_bomba: float) -> dict:
        """Presiones en cada nodo para un gasto supuesto (pasos b-d del procedimiento)."""
        pth = p_bomba - self.linea.dp(q, self.gamma_linea, self.mu_linea)["dp_T"]
        pwf = pth - self.tp.dp(q, self.gamma_tp, self.mu_tp)["dp_T"]
        pws = pwf - q / self.ii                      # Ec. 3.27
        return {"q": q, "pth": pth, "pwf": pwf, "pws": pws}

    def gasto_inyeccion(self, p_bomba: float, pws_real: float,
                        q_min: float = 10.0, q_max: float = 1e5) -> float:
        """Gasto que hace coincidir la Pws calculada con la real (paso e)."""
        return brentq(lambda q: self.nodos(q, p_bomba)["pws"] - pws_real, q_min, q_max, xtol=1e-6)


# ------------------------------------------------ 3.5 Optimización (Graper)

@dataclass
class DatosGraper:
    """Datos económicos del método de Graper (Tabla 3.5)."""
    ip1: float = 998.8      # pendiente del incremento de costo por pg de d.i. [dll/pg²/milla]
    ip2: float = 1421.6     # ordenada al origen [dll/pg/milla]
    xs: float = 400.0       # inversión en bombeo [dll/BHP]
    zv: float = 0.9091      # factor a valor presente de la inversión
    yp: float = 15.0        # operación y mantenimiento de tubería [dll/pg/milla/año]
    ys: float = 130.0       # costos variables de la estación [dll/BHP/año]
    yf: float = 26.0        # costos fijos de la estación [dll/BHP/año]
    w: float = 7.036        # factor a valor presente de gastos de operación
    fc: float = 0.95        # factor de carga
    e: float = 0.8          # eficiencia de bombeo

    @property
    def yx(self) -> float:
        """Costo anual de operación por BHP: Yf + Ys·Fc."""
        return self.yf + self.ys * self.fc


def potencia_por_milla(q: float, dp_milla: float, e: float = 0.8) -> float:
    """Ec. 3.36 — BHP por milla."""
    return C_BHP * q * dp_milla / e


def diametro_optimo_graper(q: float, gamma: float, mu: float, eco: DatosGraper = DatosGraper(),
                           rug: float = 0.0006, d0: float | None = None,
                           tol: float = 1e-3, max_iter: int = 200) -> dict:
    """Diámetro interno óptimo por milla (Ecs. 3.33-3.40, diagrama de la Fig. 3.9).

    Minimiza el costo total a valor presente por milla:
        Ct(d) = Zv[(Ip1 d + Ip2) d + Xs·BHP(d)] + W[Yp d + Yx·BHP(d)],
        BHP(d) = B/d⁵
    La condición dCt/dd = 0 da
        d⁶ = 5 B (Zv Xs + W Yx) / (Zv (2 Ip1 d + Ip2) + W Yp)
    que se resuelve por sustituciones sucesivas (f se actualiza con d).
    """
    d = d0 if d0 is not None else 0.04564 * math.sqrt(q)  # estimación inicial (Ec. 3.35)
    hist = []
    for _ in range(max_iter):
        f = factor_friccion(q, gamma, d, mu, rug)
        dp_m = dp_friccion(q, gamma, d, 1.0, f)
        bhp = potencia_por_milla(q, dp_m, eco.e)
        b = bhp * d**5
        num = 5.0 * b * (eco.zv * eco.xs + eco.w * eco.yx)
        den = eco.zv * (2.0 * eco.ip1 * d + eco.ip2) + eco.w * eco.yp
        d_c = (num / den) ** (1.0 / 6.0)
        hist.append(d_c)
        if abs(d_c - d) < tol:
            d = d_c
            break
        d = d_c
    f = factor_friccion(q, gamma, d, mu, rug)
    dp_m = dp_friccion(q, gamma, d, 1.0, f)
    return {"d": d, "f": f, "dp_milla": dp_m, "bhp_milla": potencia_por_milla(q, dp_m, eco.e),
            "iteraciones": len(hist)}


def costo_total_milla(d: float, q: float, gamma: float, mu: float,
                      eco: DatosGraper = DatosGraper(), rug: float = 0.0006) -> float:
    """Ec. 3.33 (forma explícita) — costo total a valor presente por milla [dll/milla]."""
    f = factor_friccion(q, gamma, d, mu, rug)
    bhp = potencia_por_milla(q, dp_friccion(q, gamma, d, 1.0, f), eco.e)
    inversion = (eco.ip1 * d + eco.ip2) * d + eco.xs * bhp
    operacion = eco.yp * d + eco.yx * bhp
    return eco.zv * inversion + eco.w * operacion


# ------------------------------------------------ 3.5.5 Crudos viscosos

def dp_arranque_gel(tau: float, L: float, d: float, factor_seguridad: float = 1.0) -> float:
    """Ec. 3.42 — presión para romper el gel de un crudo parafínico [psi].

    Δp = 4 τ L / D  en unidades prácticas:  Δp = 17.6 τ L / d
        τ [lbf/100 pie²], L [millas], d [pg]

    El libro recomienda un factor de seguridad de 1.5 a 2.0 sobre el valor predicho.
    El coeficiente se derivó aquí; el de la fuente es ilegible en la tesis
    escaneada (ver docs/erratas.md).
    """
    return factor_seguridad * 4.0 * (tau / 100.0) / 144.0 * (L * 5280.0) / (d / 12.0)
