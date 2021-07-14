import typer

try:
    from cupcake.__about__ import __author__, __email__, __version__
except ImportError:
    __author__ = "Elizabeth Tseng"
    __email__ = "etseng@pacb.com"
    __version__ = "unknown"

import logging

from cupcake.logger import setup_logging

cupcake_logger = setup_logging("cupcake")


def version_callback(value: bool) -> None:
    """Prints the version of the package."""
    if value:
        print(f"{__name__} version: {__version__}")
        raise typer.Exit()


def set_verbosity(v: int = 3) -> None:
    if v == 0:
        cupcake_logger.setLevel(logging.CRITICAL)
    elif v == 1:
        cupcake_logger.setLevel(logging.ERROR)
    elif v == 2:
        cupcake_logger.setLevel(logging.WARNING)
    elif v == 3:
        cupcake_logger.setLevel(logging.INFO)
    elif v == 4:
        cupcake_logger.setLevel(logging.DEBUG)
