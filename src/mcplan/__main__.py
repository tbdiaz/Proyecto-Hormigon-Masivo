"""Permite ejecutar `python -m mcplan ...` además de `mcplan ...`."""

import sys

from .cli import main

sys.exit(main())
