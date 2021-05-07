__author__ = "miles-smith@omrf.org"
import gzip
from mimetypes import guess_type
from functools import partial
from pathlib import Path
from typing import Union, Optional
from io import StringIO
from Bio import bgzf


class OpenFile(StringIO):
    """
    Convenience class to transparently handle potentially (b)gzipped files

    Partially based on code from https://stackoverflow.com/a/52839332/15270148
    and https://stackoverflow.com/a/52839332/15270148
    """
    def __init__(self, input_file: Union[str, Path], mode: Optional[str] = None):
        self.filename = Path(input_file)
        if mode is None:
            if Path(self.filename).suffix == ".gz":
                self.mode = "rt"
            else:
                self.mode = "r"
        self.file = None
        super().__init__()

    def __enter__(self):
        encoding = guess_type(self.filename)[1]
        if self.filename.suffix == ".bam":
            _open = partial(bgzf.open, mode=self.mode)
        elif encoding == "gzip":
            _open = partial(gzip.open, mode=self.mode)
        else:
            _open = partial(open, mode=self.mode)
        self.file = _open(self.filename)

    def __exit__(self, filetype, value, traceback):
        self.file.close()
        return True

    def __iter__(self):
        for line in self.file:
            yield line

    def __next__(self):
        line = self.file.readline().strip()
        if len(line) == 0:
            raise StopIteration
        return line

    def close(self):
        self.file.close()

    def read(self, n: int):
        if self.file is None:
            self.__enter__()
            return self.file.read(n)
        else:
            return self.file.read(n)
    
    def readlines(self):
        if self.file is None:
            self.__enter__()
            return [line for line in self.file]
        else:
            return [line for line in self.file]
    
    def readline(self):
        if self.file is None:
            self.__enter__()
            for line in self.file:
                return line
        else:
            for line in self.file:
                return line
