# Discrepancias encontradas al reproducir los ejemplos

Se verificó cada resultado contra la tesis UNAM de 1990 (Bernal Huicochea y López Ortiz, asesor
F. Garaicochea P.), antecedente directo del libro. Algunas diferencias pueden estar ya corregidas en la
edición impresa del CIPM; **conviene cotejar cada punto con el libro** y actualizar esta tabla.

Tipos: **E** error de la fuente · **I** dato ilegible en el escaneo, reconstruido · **N** nota de método.

| # | Tipo | Sección | Qué dice la fuente | Qué hace el código | Efecto |
|---|---|---|---|---|---|
| 1 | N | 3.1, Ec. 3.13 | Coeficiente 0.06056, derivado con 1 milla = 5277 pies | Usa 0.06056; `constantes.c_fric_liq(5280)` da 0.06060 | 0.06 % |
| 2 | E | 1.2.6, Ej. 1.1 | Lee-González-Eakin con 960/T en el término X | Por defecto 986/T (original); `c_x=960` reproduce 0.017145 cp | 0.7 % en μg |
| 3 | N | 4.1, Ej. 4.1 | El algoritmo usa p̄ aritmética, pero el resultado coincide con la p̄ de la Ec. 4.46 | `presion_fondo_estatica(..., media=...)` ofrece ambas | 0.25 % |
| 4 | E | Tabla 4.3, Ej. 4.2 | Panhandle B con (1/γ)^0.961 fuera del exponente 0.51 → 1721 MMpcd (+42 %) | `variante="estandar"` (GPSA, 1353 MMpcd) y `variante="libro"` (1721) | +27 % |
| 5 | E | 4.4.1, Ej. 4.4 | L_E con d^2.5 y f^1, y error aritmético en el tramo de 18 pg (anota 35.4×10⁻⁶ en lugar de 107.1×10⁻⁶) → 27.92 millas | f·L/d⁵ → 38.11 millas; el notebook comprueba que el gasto coincide | −27 % en L_E |
| 6 | E | 4.4.2, Ej. 4.5 | Pondera ambas ramas (14 y 16 pg) con √(f14/f16) y vuelve a mezclar exponentes → d_E = 19.95 pg, L_E = 86.64 millas | √(f_E/f_i) con f_E(d_E) → d_E = 19.68 pg, L_E = 80.3 millas | 8 % en L_E |
| 7 | I | 4.4.2, Fig. 4.2 | Figura ilegible | Configuración (14 y 16 pg en paralelo, 30 + 70 millas) reconstruida desde los números del ejemplo | **verificar con el libro** |
| 8 | I | 3.4.4 | Diámetro de la línea superficial ilegible | 4 pg, que reproduce la Tabla 3.1 | 0.1 % |
| 9 | I | 3.5.5, Ec. 3.42 | Coeficiente ilegible | Derivado: Δp = 17.6 τ L / d | **verificar con el libro** |
| 10 | E | 3.5.3 | El texto dice 2554 BHP para la alternativa 3; la Tabla 3.9, 2455 | El cálculo da 2454 | tipográfica |
| 11 | N | 2.2.3, Ec. 2.19 | La zona crítica no empata con 64/NRe en NRe = 2300 | Se implementa tal cual | discusión en el notebook 02 |
| 12 | I | 3.5.1, Ec. 3.35 | Estimación inicial del diámetro ilegible | 0.04564·√q; no afecta el resultado convergido | ninguno |
| 13 | I | 4.5, Ec. 4.54 | Numerador de la ecuación de Vitter ilegible | Forma de Peffer, Miller y Hill (SPE 15655) | **verificar con el libro** |
| 14 | I | 4.7, Ej. 4.7 | Costo de compresoras ilegible ("4ll0 dll/BHP") | X = 430 dll/BHP → K3 = 125.94 (libro 125.88) | 0.05 % |
| 15 | N | 4.7, Tabla 4.8 | El mínimo se lee en una malla de 250 lb/pg² (1750) | `presion_descarga_optima` da ≈1840 lb/pg² | 0.03 % en costo |
| 16 | I | 4.8, Ej. 4.8 | No es claro a qué elemento corresponden "300 pies" | Desnivel de la línea de descarga: reproduce 292, 698, 828 y 1081 lb/pg² | **verificar con el libro** |
| 17 | E | 4.8.1, Ec. 4.106 | Coeficiente 1291.539×10⁶ | Deducido: 1.29×10⁶ (el exponente debería ser 10³) | factor 1000 |
| 18 | I | 4.9.3, Ec. 4.129 | Constante incoherente con el resultado del Ejemplo 4.9 | Deducida de primeros principios: 23.4/64 pg (libro 23.83/64; ambos → 24/64) | 2 % en d |
| 19 | I | 4.9.3, Ec. 4.134 y Ecs. 4.116-4.128 | Constante y coeficientes ilegibles | No implementadas; se usa k medido | pendiente |
| 20 | I | 4.9.4, Ecs. 4.135-4.137 | Válvula de seguridad parcialmente ilegible | Ecuación de orificio con Y = 0.85; Cd ≈ 0.895 reproduce 907.7 lb/pg² | **verificar con el libro** |
| 21 | N | Apéndice D | Cullender-Smith con 5277 pies/milla | 5280 pies/milla | 0.03 % |
