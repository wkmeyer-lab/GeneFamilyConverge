"""
convgeno.utils.logging
~~~~~~~~~~~~~~~~~~~~~~

Centralised logging configuration.  Every CLI entry point calls
:func:`setup_logging` once at startup; library code uses the standard
``logging.getLogger(__name__)`` pattern and never configures handlers.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Optional, Union


def setup_logging(
    level: Union[int, str] = logging.INFO,
    log_file: Optional[Path] = None,
) -> None:
    """Configure the root logger for the pipeline.

    Parameters
    ----------
    level : int or str
        Logging level (e.g. ``logging.DEBUG``, ``"WARNING"``).
    log_file : Path, optional
        If provided, log messages are also written to this file.
    """
    fmt = "%(asctime)s | %(name)s | %(levelname)s | %(message)s"
    datefmt = "%Y-%m-%d %H:%M:%S"

    handlers: list[logging.Handler] = [logging.StreamHandler(sys.stderr)]
    if log_file is not None:
        log_file = Path(log_file)
        log_file.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(
        level=level,
        format=fmt,
        datefmt=datefmt,
        handlers=handlers,
        force=True,  # override any prior basicConfig
    )