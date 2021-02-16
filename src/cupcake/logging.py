import logging
from typing import Optional


def setup_logging(name: Optional[str] = None):
    if name:
        logger = logging.getLogger(name)
    else:
        logger = logging.getLogger(__name__)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    logger.setLevel(logging.DEBUG)
    logger.propagate = False

    if name:
        fh = logging.FileHandler(filename=name)
    else:
        fh = logging.FileHandler(filename=f"{__name__}.log")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    st = logging.StreamHandler()
    st.setLevel(logging.INFO)
    st.setFormatter(formatter)
    logger.addHandler(st)
