"""Permite ejecutar `python -m mcplan ...` además de `mcplan ...`."""

import sys

from .cli import main

if __name__ == "__main__":  # necesario para los procesos paralelos (multiprocessing en Windows)
    sys.exit(main())
