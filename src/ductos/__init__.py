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
gas        : Cap. 4, flujo de gas (general, Weymouth, Panhandle, Clinedinst, ...)
"""

from . import constantes, friccion, gas, liquidos, pvt_gas

__all__ = ["constantes", "friccion", "gas", "liquidos", "pvt_gas"]
__version__ = "0.1.0"
