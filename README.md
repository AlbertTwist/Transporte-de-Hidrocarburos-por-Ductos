# Transporte de Hidrocarburos por Ductos

*Principios fundamentales y el estado del arte del transporte de hidrocarburos por tubería.*

Material de clase que implementa en Python los métodos del libro **Transporte de Hidrocarburos por
Ductos** del Ing. Francisco Garaicochea Petrirena (†), Ing. César Bernal Huicochea e Ing. Oscar López
Ortiz (Colegio de Ingenieros Petroleros de México, 1991; reedición 2022).

Cada ejemplo resuelto del libro se reproduce en un notebook y se **verifica numéricamente** con `pytest`
contra los valores publicados. Sobre esa base se proponen soluciones nuevas con inteligencia artificial.

> Los enunciados de los ejemplos están parafraseados y las ecuaciones se citan por su número. Este
> repositorio no sustituye al libro: adquiéralo en el [CIPM](https://cipm.org.mx/libros/transporte-de-hidrocarburos-por-ductos/).

## Estructura

```
src/ductos/          paquete con toda la lógica de cálculo
  constantes.py      unidades prácticas y constantes derivadas de primeros principios
  friccion.py        Cap. 2  factor de fricción (laminar, crítico, Colebrook, rugoso, Jain)
  pvt_gas.py         Cap. 1.2 propiedades del gas natural (Standing, Wichert-Aziz, DPR, Lee)
  liquidos.py        Cap. 3  flujo de líquidos, pozos inyectores, optimización de Graper
  gas.py             Cap. 4  flujo de gas: general, Weymouth, Panhandle A/B, Smith, Clinedinst, Cullender-Smith
  costos_gas.py      Cap. 4.6-4.7 costos y diseño económico de gasoductos (Lagrange)
  pozo_gas.py        Cap. 4.8 pozo productor de gas: nodal, erosión, Turner
  estranguladores.py Cap. 4.9 estranguladores y válvulas de seguridad
  raspatubos.py      Cap. 4.10 velocidad de corrida de diablos
notebooks/           un notebook por capítulo, ejecutado y con salidas
  ia/                extensiones con IA
tests/               verificación contra los ejemplos del libro
docs/erratas.md      discrepancias encontradas con la fuente
scripts/             utilidades (ejecutar notebooks)
```

## Instalación

```bash
git clone https://github.com/AlbertTwist/Transporte-de-Hidrocarburos-por-Ductos.git
cd Transporte-de-Hidrocarburos-por-Ductos
python -m venv .venv && source .venv/bin/activate   # en Windows: .venv\Scripts\activate
pip install -e ".[dev,ia]"
pytest                        # 53 pruebas contra los ejemplos del libro
jupyter lab notebooks/
```

## Uso rápido

```python
from ductos.pvt_gas import Gas
from ductos.gas import Gasoducto, p2_general, q_weymouth

gas = Gas(gamma=0.6)
ducto = Gasoducto(d=28, L=41, t_f=150)               # pg, millas, °F
print(p2_general(1900, 1212.56e6, gas, ducto)["p2"])  # ≈ 1500 lb/pg² abs (Ejemplo 4.2)
print(q_weymouth(1900, 1500, gas, ducto) / 1e6)       # ≈ 1127 MMpcd
```

## Avance por capítulo

| Cap. | Tema | Estado |
|---|---|---|
| 1 | Propiedades de los fluidos | Gas natural ✅ · aceite y agua pendientes |
| 2 | Fundamentos de flujo, factor de fricción | ✅ |
| 3 | Flujo de líquidos por tuberías | ✅ 3.1–3.5 (Ejemplos 3.4.1, 3.4.4 y 3.2 reproducidos) |
| 4 | Flujo de gas por tuberías y estranguladores | ✅ 4.1–4.10 y Apéndice D (Ejemplos 4.1–4.13, D-1) |
| 5 | Flujo multifásico: conceptos fundamentales | pendiente |
| 6 | Multifásico horizontal: Bertuzzi-Tek-Poettmann, Eaton, Beggs y Brill, Dukler | pendiente |
| 7 | Multifásico vertical: Poettmann-Carpenter, Orkiszewski, Beggs y Brill, Gilbert, anulares, baches | pendiente |
| 8 | Multifásico en estranguladores: Gilbert, Ros, Achong, Poettmann-Beck, Ashford, Omaña | pendiente |
| 9 | Pozos fluyentes: IPR (Vogel, Fetkovich, Standing), análisis nodal, BNC | pendiente |
| 10 | Distribución de temperatura en tuberías | pendiente |

### Resultados verificados

| Ejemplo | Magnitud | Libro | Código |
|---|---|---|---|
| 1.1 | Z del gas libre amargo | 0.86067 | 0.86066 |
| 3.4.1 b | Δp por fricción | 192.8 lb/pg² | 192.8 |
| 3.4.4 | Gasto de inyección | 7250 bl/d | 7259 (gráfica) |
| 3.2 (Graper) | d óptimo, 155 000 bl/d | 14.734 pg | 14.734 |
| 4.2 | Presión de llegada | 1500 lb/pg² | 1499.8 |
| 4.2 | Panhandle A | 1472.509 MMpcd | 1472.513 |
| 4.3 | Weymouth inclinada | 1124.389 MMpcd | 1124.572 |
| 4.6 | Fracción de loop | 0.6649 | 0.6649 |
| 4.7 | d óptimo para pd = 1750 lb/pg² | 29.50 pg | 29.50 |
| 4.7 | Razón de compresión óptima | 1.1645 | 1.1645 |
| 4.8 | Presión en la línea de descarga, 12 MMpcd | 828 lb/pg² | 827 |
| 4.12 | Velocidad del diablo en líquido | 6.433 pie/s | 6.433 |
| D-1 | Cullender-Smith inclinada | 1215.43 MMpcd | 1215.23 |

## Extensiones con IA

`notebooks/ia/IA_01_sustituto_factor_z.ipynb` entrena una red neuronal que sustituye a la correlación
de DPR y discute sus límites (extrapolación, restricciones físicas). Siguientes propuestas:

1. Sustitutos con restricciones físicas (Z → 1 cuando p → 0, monotonía).
2. PINN para flujo transitorio en gasoductos (*line pack*), validada con la Ec. 4.19 en estado estacionario.
3. Agente LLM que recomiende la ecuación de flujo (Tabla 4.4) a partir de datos de campo y cite la sección del libro.
4. Calibración bayesiana de la eficiencia de flujo E con mediciones.
5. Optimización de redes de distribución con diámetros comerciales como restricción.

## Cómo contribuir

1. Implemente la función en `src/ductos/` con el número de ecuación en el docstring.
2. Agregue una prueba en `tests/` con el valor que publica el libro.
3. Use la función en el notebook del capítulo y ejecútelo:
   `python scripts/ejecutar_notebooks.py notebooks/*.ipynb`.
4. Si el resultado no coincide, documente la diferencia en `docs/erratas.md`.

## Referencias

- Garaicochea P., F., Bernal H., C. y López O., O. *Transporte de Hidrocarburos por Ductos*. CIPM, 1991; 2022.
- Bernal Huicochea, C. A. y López Ortiz, O. *Transporte de hidrocarburos*. Tesis de licenciatura, Facultad de
  Ingeniería, UNAM, 1990. Asesor: F. Garaicochea P. <https://hdl.handle.net/20.500.14330/TES01000117387>

## Licencia

Código bajo licencia MIT (ver `LICENSE`). Los métodos y datos numéricos pertenecen a sus autores originales
y se citan con fines educativos.
