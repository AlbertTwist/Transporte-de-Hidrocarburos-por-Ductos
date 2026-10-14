# Discrepancias encontradas al reproducir los ejemplos

Se verificó cada resultado contra la tesis UNAM de 1990 (Bernal Huicochea y López Ortiz, asesor
F. Garaicochea P.), antecedente directo del libro. Algunas diferencias pueden estar ya corregidas en la
edición impresa del CIPM; **conviene cotejar cada punto con el libro** y actualizar esta tabla.

| # | Sección | Qué dice la fuente | Qué hace el código | Efecto |
|---|---|---|---|---|
| 1 | 3.1, Ec. 3.13 | Coeficiente 0.06056, derivado con 1 milla = 5277 pies | Usa 0.06056 por defecto; `constantes.c_fric_liq(5280)` da 0.06060 | 0.06 % |
| 2 | 1.2.6, Ej. 1.1 | Lee-González-Eakin con 960/T en el término X | Por defecto 986/T (correlación original); `c_x=960` reproduce 0.017145 cp | 0.7 % en μg |
| 3 | 4.1, Ej. 4.1 | El algoritmo usa p̄ aritmética, pero el resultado (2459.3) coincide con la p̄ de la Ec. 4.46 | `presion_fondo_estatica(..., media=...)` ofrece ambas | 0.25 % |
| 4 | Tabla 4.3, Ej. 4.2 | Panhandle B con (1/γ)^0.961 fuera del exponente 0.51 → 1721 MMpcd (+42 %) | `q_panhandle_b(variante="estandar")` usa la forma GPSA (1353 MMpcd); `variante="libro"` reproduce 1721 | +27 % sobre la forma estándar |
| 5 | 4.4.1, Ej. 4.4 | Longitud equivalente en serie con d^2.5 y f^1 → 27.92 millas | `longitud_equivalente_serie` usa f·L/d⁵, que conserva p1²−p2² (38.11 millas) | se comprueba en el notebook con el gasto |
| 6 | 4.4.1, Fig. 4.1 | Longitudes de tramo poco legibles en la tesis escaneada | Se usan 10, 10 y 15 millas | **verificar con el libro** |
| 7 | 4.4.2, Fig. 4.2 | Configuración del sistema en paralelo ilegible | Ejemplo con configuración ilustrativa (14 y 16 pg) | **verificar con el libro** |
| 8 | 3.4.4 | Diámetro de la línea superficial ilegible | 4 pg, que reproduce la Tabla 3.1 | coincide al 0.1 % |
| 9 | 3.5.5, Ec. 3.42 | Coeficiente ilegible | Derivado: Δp = 17.6 τ L / d | **verificar con el libro** |
| 10 | 3.5.3 | En el texto aparece 2554 BHP para la alternativa 3; la Tabla 3.9 dice 2455 | El cálculo da 2454 | tipográfica |
| 11 | 2.2.3, Ec. 2.19 | La zona crítica no empata con 64/NRe en NRe = 2300 | Se implementa tal cual | discusión en el notebook 02 |
| 12 | 3.5.1, Ec. 3.35 | Estimación inicial del diámetro ilegible | 0.04564·√q, ajustada a los valores del ejemplo; no afecta el resultado convergido | ninguno |
