r"""*Single-source version number for* ``cdna_cupcake``.
"""

try:
    from importlib.metadata import metadata
except ImportError:
    # Running on pre-3.8 Python; use importlib-metadata package
    from importlib_metadata import metadata

__author__ = metadata["Author"]
__email__ = metadata["Author-email"]
__version__ = metadata["Version"]
