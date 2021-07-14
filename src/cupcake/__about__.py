r"""*Single-source version number for* ``cdna_cupcake``.
"""
try:
    from importlib.metadata import metadata
except ImportError:
    # Running on pre-3.8 Python; use importlib-metadata package
    from importlib_metadata import metadata

__author__ = metadata("cupcake")["Author"]
__email__ = metadata("cupcake")["Author-email"]
__version__ = metadata("cupcake")["Version"]
