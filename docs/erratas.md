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
| 6 | E | 4.4.2, Ej. 4.5 | Pondera ambas ramas (14 y 16 pg) con √(f14/f16) y vuelve a mezclar exponentes → d_E = 19.95 pg, L_E = 86.64 millas | √(f_E/f_i) con f_E(d_E) → d_E = 19.68 pg, L_E = 80.3 millas | 8 % en L_E; confirmado con el libro impreso (p. 74) |
| 7 | N | 4.4.2, Fig. 4.2 | Ilegible en la tesis | 14 y 16 pg en paralelo a lo largo de 30 millas + 70 millas de 16 pg | confirmado con el libro impreso (p. 74) |
| 8 | N | 3.4.4 | Ilegible en la tesis | Línea superficial de 4 pg | confirmado con el libro impreso (p. 41) |
| 9 | E | 3.5.5, Ec. 3.42 | Δp = 175 900 Lτ/d, que es 4(5277)(12)/144 multiplicado por 100 en lugar de dividido (τ en lbf/100 pie²) | Δp = 17.6 τ L / d; `variante="libro"` conserva el impreso | factor 10 000; confirmado con el libro impreso (p. 57) |
| 10 | E | 3.5.3 | El texto dice 2554 BHP para la alternativa 3; la Tabla 3.9, 2455 | El cálculo da 2454 | tipográfica |
| 11 | N | 2.2.3, Ec. 2.19 | La zona crítica no empata con 64/NRe en NRe = 2300 | Se implementa tal cual | discusión en el notebook 02 |
| 12 | I | 3.5.1, Ec. 3.35 | Estimación inicial del diámetro ilegible | 0.04564·√q; no afecta el resultado convergido | ninguno |
| 13 | E | 4.5, Ec. 4.54 | Vitter con denominador 1 + 1123 γo/R_L | Por defecto 1 + 1123/R_L (Peffer, Miller y Hill, SPE 15655); `variante="libro"` conserva el impreso | ~1.4 % en γ; confirmado con el libro impreso (p. 75) |
| 14 | I | 4.7, Ej. 4.7 | Costo de compresoras ilegible ("4ll0 dll/BHP") | X = 430 dll/BHP → K3 = 125.94 (libro 125.88) | 0.05 % |
| 15 | N | 4.7, Tabla 4.8 | El mínimo se lee en una malla de 250 lb/pg² (1750) | `presion_descarga_optima` da ≈1840 lb/pg² | 0.03 % en costo |
| 16 | N | 4.8, Ej. 4.8 | Los 300 pies (dato h) corresponden al desnivel de la línea de descarga | h = 300 pies en la LD: reproduce 153, 292, 698, 828 y 1081 lb/pg² | confirmado con el libro impreso (p. 85) |
| 16b | E | 4.8, Ej. 4.8 | IPR impresa como 16 280 (pws − pwf)^0.681, sin cuadrados (la tesis: 16 290 con cuadrados) | 16 280 (pws² − pwf²)^0.681, como la Ec. 4.101; reproduce la tabla | confirmado con el libro impreso (p. 86) |
| 16c | E | 4.8, tabla del Ej. 4.8 | 1031 lb/pg² (pws = 2000, 12 MMpcd) y 979 lb/pg² (pws = 1500, 10 MMpcd) rompen la forma de sus curvas | El cálculo da 1305 y 888; los otros 16 valores de la tabla coinciden dentro del 1.2 % | probable errata de transcripción |
| 17 | E | 4.8.1, Ec. 4.106 | Coeficiente 1291.539×10⁶ | Deducido: 1.29×10⁶ (el exponente debería ser 10³) | factor 1000 |
| 18 | E | 4.9.3, Ec. 4.129 y Ej. 4.9 | Cd = 1.9083×10⁻⁴; el ejemplo lo multiplica, obtiene dφ² = 25 y lo reporta como dφ = 25/64 | Constante deducida: q = 238.4 dφ² pth √[kH/((k−1)γZT)] → 23.4/64 (`diametro_ec_4_129`) | ambos → 24/64; confirmado con el libro impreso |
| 19 | E | 4.9.1, Ecs. 4.116-4.128 | Correlación de k: b5 impreso con 0.06673 Y⁴, que la hace divergir; la tesis trae 0.006673 | `k_gas` con 0.006673; fiable solo para 0.65 ≤ γ ≤ 0.80 (el polinomio de grado 9 diverge fuera); k(0.69, 140 °F) = 1.254 vs 1.25 medido | confirmado con el libro impreso |
| 19b | N | 4.9.3, Ej. 4.9 (Cook-Dotterweich) | Convierte CdA con la Ec. 4.133 → 23.83/64 es el diámetro del área efectiva | `d64_libro` = 23.83/64 y `d64` geométrico = 26.2/64 | 10 % en d |
| 20 | E | 4.9.4, Ec. 4.135 y Ej. 4.11 | Constante impresa 69.22×10⁻⁶; el ejemplo usa en el corchete 12.51, que corresponde a 10⁻⁸, y omite (1 − ψ⁴) al sustituir | K = 69.22×10⁻⁸; `con_psi4=False` reproduce 907.74; con el factor, 896.4; Y iterado (Ec. 4.137) = 0.845 | 1.2 %; confirmado con el libro impreso |
| 21 | N | Apéndice D | Cullender-Smith con 5277 pies/milla | 5280 pies/milla | 0.03 % |
| 22 | E | 1.3.2, Ej. 1.1 (Tabla 1) | Rs de Vázquez = 610 pie³/bl, incoherente con su propia pb de Vázquez (3626) | 676.6 pie³/bl (C_Rs = 1.1085 en lugar de 1.2295) | no cambia la selección de correlaciones |
| 23 | N | 1.3.4, Ec. 1.49 | Exponentes ilegibles en la tesis | y_g = 419.545×10⁻⁵ pf³ − 591.428×10⁻⁴ pf² + 334.519×10⁻³ pf + 169.879×10⁻⁴ | confirmado con el libro impreso (p. 19) |
| 24 | N | 1.3.3, Ec. 1.67 | El exponente de T para pb es 0.130 (aceite volátil), aunque el Rs del mismo aceite negro usa 0.172 | Se reproduce el libro; ambos exponentes son parámetros | pb Oistein = 2971 |
| 25 | E | 1.5.1, Ec. 1.69 y Ej. 1.1 | El libro imprime y usa 1.0×10⁻⁶ (T − 60); la correlación original (Gould) lleva (T − 60)² → Bw = 1.0103 | Por defecto con el cuadrado (1.0297); `variante="ejemplo"` reproduce el libro | 2 %; forma impresa confirmada (p. 20) |
| 26 | N | 1.4.1, Ec. 1.59 | C4 = −1180 | −1180 | confirmado con el libro impreso (p. 19) |
| 27 | N | 1.8, Ej. 1.1 | Elige Standing para Bo, pero calcula Bo a 2000 lb/pg² con Oistein corregido (1.28363); la tabla del Ej. 1.2 sí usa Standing (1.2851) | Standing corregido; el notebook muestra ambos | 0.1 % |
| 28 | E | 1.8, Ej. 1.2 | "se supuso una temperatura media igual a 200 °C" | 200 °F | tipográfica |
| 29 | N | 1.4.2, 1.4.4, 1.6.2 y 1.6.3 (Ecs. 1.60, 1.63, 1.89 y 1.90) | ρ = ρb exp[c (p − pb)], B = Bb exp[−c (p − pb)] | Las mismas | confirmado con el libro impreso (pp. 19 y 21) |

| 30 | E | 1.3.6, Ec. 1.55 | X = Y T^−1.63 | X = Y T^−1.163 (Beggs y Robinson), el que usa el Ejemplo 1.1 | con −1.63, μo cambia por completo; tipográfica |
| 31 | E | 1.5.5, Ec. 1.80 | T' = (T* − 32)/10 | T' = (T* − 90)/10, el que usa el Ejemplo 1.1 (T' = 0.333) | con 32, Rsw = 10.4 en lugar de 12.21 pie³/bl |
