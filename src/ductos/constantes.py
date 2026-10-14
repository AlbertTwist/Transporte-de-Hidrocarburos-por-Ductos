"""Constantes y conversiones en unidades prácticas de campo.

Convención de unidades (la del libro):
    presión          lb/pg² abs (psia)
    temperatura      °F  (se convierte a °R sumando 460)
    longitud         millas para ductos, pies para desniveles
    diámetro         pg (pulgadas), diámetro interno
    gasto de líquido bl/día
    gasto de gas     pies³/día a condiciones estándar
    viscosidad       cp
    densidad rel.    agua = 1.0 (líquidos), aire = 1.0 (gases)

Notas sobre las constantes del libro
------------------------------------
El libro deriva sus coeficientes con 1 milla = 5277 pies (posible errata de 5280).
Para reproducir los ejemplos se usa el valor del libro por defecto; el valor exacto
queda disponible para comparación. La diferencia es de ~0.06 %.
"""

import math

# Condiciones estándar en México (Cap. 1)
P_STD_MEXICO = 14.69   # lb/pg² abs
T_STD = 60.0           # °F

# Constantes físicas
R_GAS = 10.7316        # lb/pg² · pie³ / (lbmol · °R)
M_AIRE = 28.97         # lbm/lbmol
GC = 32.174            # lbm·pie/(lbf·s²)
RHO_AGUA = 62.428      # lbm/pie³ a 60 °F
RANKINE = 460.0        # °F -> °R (el libro usa 460)

# Longitud de la milla
MILLA_PIES_LIBRO = 5277.0
MILLA_PIES = 5280.0

# Coeficientes del Cap. 3 (líquidos)
C_FRIC_LIQ_LIBRO = 0.06056   # Ec. 3.13  Δp = C f γ q² L / d⁵
C_ELEV_LIQ = 0.433           # Ec. 3.2   Δp = 0.433 γ Δh
C_RE_LIQ = 92.2              # Ec. 3.22  NRe = 92.2 q γ / (d μ)
C_BHP = 1.70322e-5           # Ec. 3.36  BHP = C q Δp / e


def c_fric_liq(pies_por_milla: float = MILLA_PIES_LIBRO) -> float:
    """Coeficiente de la Ec. 3.13 derivado desde Darcy-Weisbach.

    Δp[psi] = c · f · γ · q[bl/d]² · L[mi] / d[pg]⁵
    """
    v_por_q = 5.6146 / 86400.0                      # pie³/s por bl/d
    coef_por_pie = (RHO_AGUA * v_por_q**2 * 16.0 * 12.0**5
                    / (math.pi**2 * 2.0 * GC * 144.0))
    return coef_por_pie * pies_por_milla


def c_general_gas() -> float:
    """Constante de la ecuación general de flujo de gas (Ec. 4.20, Darcy).

    q[pie³/d @ c.s.] = C · (T_b/p_b) · [ (p1² - p2²) · d⁵ / (γ · Z · T · L · f) ]^0.5

    con T en °R, p en psia, L en millas, d en pg y f de Darcy.
    Derivada de primeros principios; vale ≈ 77.56 (literatura: 77.54).
    """
    c1 = M_AIRE / (R_GAS * 144.0)                    # lbm/pie³ por (psia/°R), γ=1
    num = 16.0 * c1 * MILLA_PIES * 12.0**5
    den = 86400.0**2 * GC * math.pi**2
    return math.sqrt(den / num)


def f_a_r(t_f: float) -> float:
    """°F -> °R con la convención del libro (+460)."""
    return t_f + RANKINE
