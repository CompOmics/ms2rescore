"""Entrypoint for MS²Rescore GUI."""

import multiprocessing
import os
import sys
from pathlib import Path

from ms2rescore.gui.app import app


def main():
    """Entrypoint for MS²Rescore GUI."""
    multiprocessing.freeze_support()

    # Fix for PyInstaller windowed mode: sys.stdout/stderr can be None
    # This causes issues with libraries that try to write to stdout (e.g., Keras progress bars)
    if sys.stdout is None:
        sys.stdout = Path(os.devnull).open("w")  # noqa: SIM115 - must stay open for process lifetime
    if sys.stderr is None:
        sys.stderr = Path(os.devnull).open("w")  # noqa: SIM115 - must stay open for process lifetime

    app()


if __name__ == "__main__":
    main()
