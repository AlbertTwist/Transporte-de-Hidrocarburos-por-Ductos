"""Ejecuta los notebooks y guarda sus salidas (se usa en CI y antes de cada commit).

Uso:
    python scripts/ejecutar_notebooks.py notebooks/*.ipynb notebooks/ia/*.ipynb
"""
import sys
from pathlib import Path

import nbformat
from nbclient import NotebookClient

ok = True
for p in map(Path, sys.argv[1:]):
    nb = nbformat.read(p, as_version=4)
    try:
        NotebookClient(nb, timeout=600, kernel_name="python3",
                       resources={"metadata": {"path": str(p.parent)}}).execute()
        nbformat.write(nb, p)
        print("OK ", p)
    except Exception as e:  # noqa: BLE001
        ok = False
        print("ERROR", p, str(e)[-1500:])
sys.exit(0 if ok else 1)
