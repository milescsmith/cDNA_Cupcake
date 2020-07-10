__version__ = "8.6"

# from . import (
#     ice,
#     io,
#     tofu,
# )

try:
    with open(osp.join("cupcake", "__about__.py")) as f:
        exec(f.read())
except:
    __author__ = "Elizabeth Tseng"
    __email__ = "etseng@pacb.com"
    __version__ = "unknown"