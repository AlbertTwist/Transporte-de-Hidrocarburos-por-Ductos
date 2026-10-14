"""ductos — Transporte de hidrocarburos por ductos.

Implementación en Python de los métodos del libro *Transporte de Hidrocarburos
por Ductos* (F. Garaicochea P., C. Bernal H., O. López O., CIPM), para uso
como material de clase.

Módulos
-------
constantes : unidades prácticas y constantes derivadas
friccion   : Cap. 2, factor de fricción (laminar, crítico, Colebrook, rugoso, Jain)
pvt_gas    : Cap. 1.2, propiedades del gas natural (Z de DPR, μg de Lee, etc.)
liquidos   : Cap. 3, flujo de líquidos, pozos inyectores, optimización de Graper
gas        : Cap. 4.1-4.4, flujo de gas (general, Weymouth, Panhandle, Clinedinst, Cullender-Smith)
costos_gas : Cap. 4.6-4.7, costos y diseño económico de gasoductos
pozo_gas   : Cap. 4.8, pozo productor de gas, erosión y descarga de líquidos
estranguladores : Cap. 4.9, estranguladores y válvulas de seguridad
raspatubos : Cap. 4.10, velocidad de corrida de diablos
"""

from . import (constantes, costos_gas, estranguladores, friccion, gas, liquidos, pozo_gas,
               pvt_gas, raspatubos)

__all__ = ["constantes", "costos_gas", "estranguladores", "friccion", "gas", "liquidos",
           "pozo_gas", "pvt_gas", "raspatubos"]
__version__ = "0.2.0"
