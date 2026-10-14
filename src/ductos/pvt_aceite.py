"""Cap. 1.1, 1.3, 1.4 y 1.7 — Propiedades del aceite (aceite negro y volátil).

Correlaciones para aceite saturado (Standing, Vázquez, Oistein (Glaso), Lasater), densidad,
viscosidad (Beggs-Robinson), tensión superficial, propiedades del aceite bajosaturado
(Vázquez-Beggs) y los ajustes pVT a las condiciones de separación.

La clase `AceiteNegro` sigue el procedimiento del Ejemplo 1.1: se comparan las correlaciones
contra el laboratorio en el punto de burbujeo, se eligen las que mejor ajustan y se aplican
sus factores de corrección C = valor de laboratorio / valor calculado.

Unidades: p [lb/pg² abs], T [°F], Rs [pie³/bl], Bo [bl @ c.esc./bl @ c.s.], ρ [lbm/pie³],
μ [cp], σ [dina/cm], co [1/(lb/pg²)].
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from scipy.optimize import brentq


# --------------------------------------------------------------- 1.1-1.2 gas

def gamma_promedio_separadores(gammas: list[float], gastos: list[float]) -> float:
    """Ec. 1.1 — densidad relativa del gas producido en n etapas de separación."""
    return sum(g * q for g, q in zip(gammas, gastos)) / sum(gastos)


def gamma_gas_disuelto(rs: float, api: float) -> float:
    """Ec. 1.2 (Katz) — densidad relativa del gas disuelto."""
    return 0.25 + 0.02 * api + rs * 1e-6 * (0.6874 - 3.5864 * api)


def gamma_gas_libre(R: float, gamma_g: float, rs: float, gamma_gd: float) -> float:
    """Ec. 1.8 — densidad relativa del gas libre por balance de masa (R > Rs)."""
    return (R * gamma_g - rs * gamma_gd) / (R - rs)


def gamma_aceite(api: float) -> float:
    return 141.5 / (131.5 + api)


def gamma_gas_vazquez(gamma_gp: float, api: float, t_sep: float, p_sep: float) -> float:
    """Ec. 1.35 — γg normalizada a una separación de 100 lb/pg² man (114.7 abs)."""
    return gamma_gp * (1.0 + 5.912e-5 * api * t_sep * math.log10(p_sep / 114.7))


# --------------------------------------------------------------- 1.3.1 Standing

def pb_standing(rs: float, gamma_gd: float, api: float, t_f: float) -> float:
    """Ecs. 1.31 / 1.64."""
    return 18.0 * (rs / gamma_gd) ** 0.83 * 10.0 ** (0.00091 * t_f - 0.0125 * api)


def rs_standing(p: float, gamma_gd: float, api: float, t_f: float) -> float:
    """Ec. 1.32."""
    return gamma_gd * (p / 18.0 * 10.0 ** (0.0125 * api - 0.00091 * t_f)) ** (1.0 / 0.83)


def bo_standing(rs: float, gamma_gd: float, gamma_o: float, t_f: float) -> float:
    """Ecs. 1.33-1.34."""
    F = rs * (gamma_gd / gamma_o) ** 0.5 + 1.25 * t_f
    return 0.972 + 0.000147 * F**1.175


# --------------------------------------------------------------- 1.3.2 Vázquez

def _coef_vazquez(api: float) -> dict:
    if api <= 30.0:
        return {"C1": 0.0362, "C2": 1.0937, "C3": 25.724,
                "B1": 4.677e-4, "B2": 1.751e-5, "B3": -1.811e-8}
    return {"C1": 0.0178, "C2": 1.1870, "C3": 23.931,
            "B1": 4.670e-4, "B2": 1.100e-5, "B3": 1.337e-9}


def rs_vazquez(p: float, gamma_gs: float, api: float, t_f: float) -> float:
    """Ec. 1.36."""
    c = _coef_vazquez(api)
    return c["C1"] * gamma_gs * p ** c["C2"] * math.exp(c["C3"] * api / (t_f + 460.0))


def pb_vazquez(rs: float, gamma_gs: float, api: float, t_f: float) -> float:
    """Ec. 1.65."""
    c = _coef_vazquez(api)
    return (rs / (c["C1"] * gamma_gs * math.exp(c["C3"] * api / (t_f + 460.0)))) ** (1.0 / c["C2"])


def bo_vazquez(rs: float, gamma_gs: float, api: float, t_f: float) -> float:
    """Ec. 1.37."""
    c = _coef_vazquez(api)
    return 1.0 + c["B1"] * rs + (t_f - 60.0) * (api / gamma_gs) * (c["B2"] + c["B3"] * rs)


# --------------------------------------------------------------- 1.3.3 Oistein (Glaso)

def _p_estrella(p: float) -> float:
    """Ec. 1.38."""
    lp = math.log10(p)
    return 10.0 ** (-2.57364 + 2.35772 * lp - 0.703988 * lp**2 + 0.098479 * lp**3)


def rs_oistein(p: float, gamma_gd: float, api: float, t_f: float, a: float = 0.172) -> float:
    """Ec. 1.39 — a = 0.172 aceite negro, 0.130 aceite volátil (puede ajustarse)."""
    return gamma_gd * (_p_estrella(p) * api**0.989 / t_f**a) ** (1.0 / 0.816)


def pb_oistein(rs: float, gamma_gd: float, api: float, t_f: float, a: float = 0.130) -> float:
    """Ecs. 1.66-1.67 — pb* = (Rs/γgd)^0.816 T^a / °API^0.989.

    El libro escribe el exponente 0.130 en la Ec. 1.67 (y así resuelve el Ejemplo 1.1),
    aunque para el Rs del mismo aceite usa 0.172.
    """
    pbe = (rs / gamma_gd) ** 0.816 * t_f**a / api**0.989
    lpb = math.log10(pbe)
    return 10.0 ** (1.7669 + 1.7447 * lpb - 0.30218 * lpb**2)


def bo_oistein(rs: float, gamma_gd: float, gamma_o: float, t_f: float) -> float:
    """Ecs. 1.40-1.42."""
    be = rs * (gamma_gd / gamma_o) ** 0.526 + 0.968 * t_f
    lb = math.log10(be)
    A = -6.58511 + 2.91329 * lb - 0.27683 * lb**2
    return 1.0 + 10.0**A


# --------------------------------------------------------------- 1.3.4 Lasater

def peso_molecular_aceite_lasater(api: float) -> float:
    """Ecs. 1.46-1.47 — peso molecular del aceite en el tanque."""
    if api <= 40.0:
        return (63.506 - api) / 0.0996
    return (1048.33 / api) ** 1.6736


def _pf_de_yg(yg: float) -> float:
    """Ec. 1.44."""
    return 5.043 * yg**3 + 3.10526 * yg**2 + 1.36226 * yg + 0.119118


def _yg_de_pf(pf: float) -> float:
    """Ec. 1.49 (ajuste inverso de la Ec. 1.44).

    y_g = 419.545×10⁻⁵ pf³ − 591.428×10⁻⁴ pf² + 334.519×10⁻³ pf + 169.879×10⁻⁴
    (confirmado con el libro impreso, p. 19; reproduce el y_g = 0.685 del Ejemplo 1.1).
    """
    return 4.19545e-3 * pf**3 - 5.91428e-2 * pf**2 + 0.334519 * pf + 0.0169879


def fraccion_molar_gas(rs: float, gamma_o: float, mo: float) -> float:
    """Ec. 1.45."""
    n_g = rs / 379.3
    return n_g / (n_g + 350.0 * gamma_o / mo)


def pb_lasater(rs: float, gamma_gd: float, api: float, t_f: float) -> float:
    """Ecs. 1.43-1.45 / 1.68."""
    go = gamma_aceite(api)
    yg = fraccion_molar_gas(rs, go, peso_molecular_aceite_lasater(api))
    return _pf_de_yg(yg) * (t_f + 460.0) / gamma_gd


def rs_lasater(p: float, gamma_gd: float, api: float, t_f: float) -> float:
    """Ecs. 1.48-1.49."""
    go = gamma_aceite(api)
    mo = peso_molecular_aceite_lasater(api)
    yg = _yg_de_pf(p * gamma_gd / (t_f + 460.0))
    return 132755.0 * go * yg / ((1.0 - yg) * mo)


# --------------------------------------------------------------- 1.3.5-1.3.7

def densidad_aceite_saturado(rs: float, gamma_gd: float, gamma_o: float, bo: float) -> float:
    """Ec. 1.50 [lbm/pie³]."""
    return (62.4 * gamma_o + 0.01362 * rs * gamma_gd) / bo


def viscosidad_aceite_muerto(api: float, t_f: float) -> float:
    """Ecs. 1.54-1.57 (Beggs y Robinson).

    X = Y T^−1.163. El libro imprime T^−1.63 en la Ec. 1.55, pero el Ejemplo 1.1 usa
    −1.163, que es el exponente de la correlación original (ver docs/erratas.md).
    """
    z = 3.0324 - 0.02023 * api
    x = 10.0**z * t_f**-1.163
    return 10.0**x - 1.0


def viscosidad_aceite_saturado(rs: float, api: float, t_f: float) -> float:
    """Ecs. 1.51-1.53 (Beggs y Robinson)."""
    a = 10.715 * (rs + 100.0) ** -0.515
    b = 5.44 * (rs + 150.0) ** -0.338
    return a * viscosidad_aceite_muerto(api, t_f) ** b


def tension_superficial_aceite(p: float, t_f: float, api: float) -> float:
    """Ec. 1.58 [dina/cm]."""
    return (42.4 - 0.047 * t_f - 0.267 * api) * math.exp(-0.0007 * p)


# --------------------------------------------------------------- 1.4 bajosaturado

def compresibilidad_aceite(p: float, rs: float, t_f: float, gamma_gs: float, api: float,
                           c4: float = -1180.0) -> float:
    """Ec. 1.59 (Vázquez y Beggs).

    co = (−1433 + 5 Rs + 17.2 T + c4 γgs + 12.61 °API) / (10⁵ p),  c4 = −1180
    (confirmado con el libro impreso, p. 19).
    """
    return (-1433.0 + 5.0 * rs + 17.2 * t_f + c4 * gamma_gs + 12.61 * api) / (1e5 * p)


def viscosidad_aceite_bajosaturado(mu_ob: float, p: float, pb: float) -> float:
    """Ecs. 1.61-1.62 (Vázquez y Beggs)."""
    m = 2.6 * p**1.187 * math.exp(-11.513 - 8.98e-5 * p)
    return mu_ob * (p / pb) ** m


# --------------------------------------------------------------- 1.7 ajustes pVT

def bo_flash(bofb: float, bod: float, bodb: float) -> float:
    """Ec. 1.91 — Bof = Bofb · (Bod / Bodb) (curva de encogimiento)."""
    return bofb * bod / bodb


def rs_flash(rsfb: float, rsdb: float, rsd: float, bofb: float, bodb: float) -> float:
    """Ec. 1.92 — Rs a condiciones del separador a partir de la liberación diferencial."""
    return rsfb - (rsdb - rsd) * bofb / bodb


# --------------------------------------------------------------- procedimiento del Ej. 1.1

CORRELACIONES_PB = {
    "Standing": lambda rs, ggd, ggs, api, t: pb_standing(rs, ggd, api, t),
    "Vázquez": lambda rs, ggd, ggs, api, t: pb_vazquez(rs, ggs, api, t),
    "Oistein": lambda rs, ggd, ggs, api, t: pb_oistein(rs, ggd, api, t),
    "Lasater": lambda rs, ggd, ggs, api, t: pb_lasater(rs, ggd, api, t),
}
CORRELACIONES_RS = {
    "Standing": lambda p, ggd, ggs, api, t: rs_standing(p, ggd, api, t),
    "Vázquez": lambda p, ggd, ggs, api, t: rs_vazquez(p, ggs, api, t),
    "Oistein": lambda p, ggd, ggs, api, t: rs_oistein(p, ggd, api, t),
    "Lasater": lambda p, ggd, ggs, api, t: rs_lasater(p, ggd, api, t),
}
CORRELACIONES_BO = {
    "Standing": lambda rs, ggd, ggs, api, t: bo_standing(rs, ggd, gamma_aceite(api), t),
    "Vázquez": lambda rs, ggd, ggs, api, t: bo_vazquez(rs, ggs, api, t),
    "Oistein": lambda rs, ggd, ggs, api, t: bo_oistein(rs, ggd, gamma_aceite(api), t),
}


@dataclass
class AceiteNegro:
    """Aceite caracterizado como en el Ejemplo 1.1, con correlaciones calibradas.

    api, gamma_g (gas producido), R (relación gas-aceite producida), p_sep y t_sep de la
    primera etapa de separación. Las correlaciones elegidas y sus factores de corrección
    se fijan con `calibrar` o se indican a mano.
    """
    api: float
    gamma_g: float
    R: float
    p_sep: float = 114.7
    t_sep: float = 60.0
    corr_pb: str = "Lasater"
    corr_rs: str = "Lasater"
    corr_bo: str = "Standing"
    c_pb: float = 1.0
    c_rs: float = 1.0
    c_bo: float = 1.0
    c_gd: float = 1.0
    calibracion: list = field(default_factory=list)

    @property
    def gamma_o(self) -> float:
        return gamma_aceite(self.api)

    @property
    def gamma_gs(self) -> float:
        return gamma_gas_vazquez(self.gamma_g, self.api, self.t_sep, self.p_sep)

    # ---------------------------------------------------- calibración
    def calibrar(self, pb_lab: float, rsfb_lab: float, bofb_lab: float, t_b: float,
                 elegir: bool = True) -> list[dict]:
        """Compara las correlaciones contra el laboratorio en el punto de burbujeo (Tabla 1).

        En pb todo el gas está disuelto, así que γgd = γg. Con elegir=True selecciona,
        para cada propiedad, la correlación con el factor de corrección más cercano a 1.
        """
        ggd, ggs = self.gamma_g, self.gamma_gs
        filas = []
        for nombre in CORRELACIONES_PB:
            rs = CORRELACIONES_RS[nombre](pb_lab, ggd, ggs, self.api, t_b)
            pb = CORRELACIONES_PB[nombre](rsfb_lab, ggd, ggs, self.api, t_b)
            fila = {"correlación": nombre, "Rs": rs, "C_Rs": rsfb_lab / rs,
                    "pb": pb, "C_pb": pb_lab / pb, "Bo": None, "C_Bo": None}
            if nombre in CORRELACIONES_BO:
                bo = CORRELACIONES_BO[nombre](rsfb_lab, ggd, ggs, self.api, t_b)
                fila.update({"Bo": bo, "C_Bo": bofb_lab / bo})
            filas.append(fila)
        self.calibracion = filas
        if elegir:
            mejor = lambda clave: min((f for f in filas if f[clave] is not None),
                                      key=lambda f: abs(f[clave] - 1.0))
            f_pb, f_rs, f_bo = mejor("C_pb"), mejor("C_Rs"), mejor("C_Bo")
            self.corr_pb, self.c_pb = f_pb["correlación"], f_pb["C_pb"]
            self.corr_rs, self.c_rs = f_rs["correlación"], f_rs["C_Rs"]
            self.corr_bo, self.c_bo = f_bo["correlación"], f_bo["C_Bo"]
        # Ajuste de la Ec. 1.2 para que en pb dé γgd = γg
        self.c_gd = self.gamma_g / gamma_gas_disuelto(rsfb_lab, self.api)
        return filas

    # ---------------------------------------------------- propiedades
    def pb(self, t_f: float) -> float:
        return self.c_pb * CORRELACIONES_PB[self.corr_pb](self.R, self.gamma_g, self.gamma_gs,
                                                          self.api, t_f)

    def gamma_gd(self, rs: float) -> float:
        return min(self.c_gd * gamma_gas_disuelto(rs, self.api), 10.0)

    def rs_y_gamma_gd(self, p: float, t_f: float, tol: float = 1e-6) -> tuple[float, float]:
        """Rs y γgd a p y T por ensaye y error (pasos 1-6 del Ejemplo 1.1)."""
        if p >= self.pb(t_f):
            return self.R, self.gamma_g
        ggd = self.gamma_g
        for _ in range(200):
            rs = self.c_rs * CORRELACIONES_RS[self.corr_rs](p, ggd, self.gamma_gs, self.api, t_f)
            rs = min(rs, self.R)
            ggd_n = self.gamma_gd(rs)
            if abs(ggd_n - ggd) < tol:
                return rs, ggd_n
            ggd = ggd_n
        raise RuntimeError("Rs-γgd no convergió")

    def rs(self, p: float, t_f: float) -> float:
        return self.rs_y_gamma_gd(p, t_f)[0]

    def gamma_gf(self, p: float, t_f: float) -> float | None:
        """Ec. 1.8 — densidad relativa del gas libre (None si está bajosaturado)."""
        rs, ggd = self.rs_y_gamma_gd(p, t_f)
        if rs >= self.R:
            return None
        return gamma_gas_libre(self.R, self.gamma_g, rs, ggd)

    def bo_saturado(self, p: float, t_f: float) -> float:
        rs, ggd = self.rs_y_gamma_gd(p, t_f)
        return self.c_bo * CORRELACIONES_BO[self.corr_bo](rs, ggd, self.gamma_gs, self.api, t_f)

    def co(self, p: float, t_f: float) -> float:
        return compresibilidad_aceite(p, self.R, t_f, self.gamma_gs, self.api)

    def bo(self, p: float, t_f: float) -> float:
        """Bo saturado o bajosaturado (Ec. 1.63: Bo = Bob exp[co (pb − p)])."""
        pb = self.pb(t_f)
        if p < pb:
            return self.bo_saturado(p, t_f)
        return self.bo_saturado(pb, t_f) * math.exp(self.co(p, t_f) * (pb - p))

    def densidad(self, p: float, t_f: float) -> float:
        """Ecs. 1.50 y 1.60 [lbm/pie³]."""
        pb = self.pb(t_f)
        pp = min(p, pb)
        rs, ggd = self.rs_y_gamma_gd(pp, t_f)
        rho = densidad_aceite_saturado(rs, ggd, self.gamma_o, self.bo_saturado(pp, t_f))
        if p > pb:
            rho *= math.exp(self.co(p, t_f) * (p - pb))
        return rho

    def viscosidad(self, p: float, t_f: float) -> float:
        """Ecs. 1.51-1.57 y 1.61-1.62 [cp]."""
        pb = self.pb(t_f)
        if p <= pb:
            return viscosidad_aceite_saturado(self.rs(p, t_f), self.api, t_f)
        mu_ob = viscosidad_aceite_saturado(self.R, self.api, t_f)
        return viscosidad_aceite_bajosaturado(mu_ob, p, pb)

    def tension_superficial(self, p: float, t_f: float) -> float:
        return tension_superficial_aceite(p, t_f, self.api)

    def tabla(self, presiones, t_f: float) -> list[dict]:
        """Tabla del Ejemplo 1.2: Rs, γgd, γgf y Bo contra presión."""
        filas = []
        for p in presiones:
            rs, ggd = self.rs_y_gamma_gd(p, t_f)
            filas.append({"p": p, "Rs": rs, "gamma_gd": ggd, "gamma_gf": self.gamma_gf(p, t_f),
                          "Bo": self.bo(p, t_f)})
        return filas
