__author__ = "etseng@pacb.com"

"""
Misc APIs for reading BED format

UCSC BED format:
https://genome.ucsc.edu/FAQ/FAQformat#format1

SimpleBED --- chr, 0-based start, 1-based end
"""
from sys import stdout


class SimpleBED:
    def __init__(self, chrom, start, end, name=None):
        self.chr = chrom
        self.start = start
        self.end = end
        self.name = name

    def __str__(self):
        return f"{self.chr}:{self.start}-{self.end} (name:{self.name})"


class SimpleBEDReader:
    def __init__(self, filename, start_base=0, end_base=1):
        self.filename = filename
        self.f = open(filename)

        if start_base != 0 and start_base != 1:
            raise Exception("start_base can only be 0 or 1!")
        if end_base != 0 and end_base != 1:
            raise Exception("end_base can only be 0 or 1!")

        self.start_base = start_base
        self.end_base = end_base

    def __iter__(self):
        return self

    def __next__(self):
        return self.read()

    def read(self):
        cur = self.f.tell()
        line = self.f.readline()
        if self.f.tell() == cur:
            raise StopIteration

        raw = line.strip().split("\t")
        if len(raw) >= 4:
            name = raw[3]
        else:
            name = None
        return SimpleBED(
            raw[0],
            int(raw[1]) - self.start_base,
            int(raw[2]) + (1 - self.end_base),
            name,
        )


class SimpleBEDWriter:
    def __init__(self, handle):
        self.handle = handle

    def writerow(self, r):
        if r.name is not None:
            self.handle.write(f"{r.chrom}\t{r.start}\t{r.end}\t{r.name}\n")
        else:
            self.handle.write(f"{r.chrom}\t{r.start}\t{r.end}\n")

    def writerows(self, rows):
        for r in rows:
            self.writerow(r)


class LazyBEDPointReader(SimpleBEDReader):
    def __init__(
        self, filename, start_base=0, end_base=1, windowsize=100000, debug=False
    ):
        self.filename = filename
        self.f = open(filename)

        if start_base != 0 and start_base != 1:
            raise Exception("start_base can only be 0 or 1!")
        if end_base != 0 and end_base != 1:
            raise Exception("end_base can only be 0 or 1!")

        self.start_base = start_base
        self.end_base = end_base

        self.windowsize = windowsize
        self.pos_d = {}  # chrom --> (pos/windowsize) --> file location

        while True:
            cur = self.f.tell()
            line = self.f.readline()
            if self.f.tell() == cur:
                break

            raw = line.strip().split("\t")
            chrom, start = (raw[0], int(raw[1]) - self.start_base)

            if chrom not in self.pos_d:
                self.pos_d[chrom] = {}

            i = start / self.windowsize
            if i not in self.pos_d[chrom]:
                self.pos_d[chrom][i] = cur
                if debug:
                    print(f"**** Hashing {chrom}:{start}....", file=stdout)

    def get_pos(self, chrom, pos):
        if chrom not in self.pos_d:
            return "NA"

        i = pos / self.windowsize
        if i not in self.pos_d[chrom]:
            return "NA"

        self.f.seek(self.pos_d[chrom][i])
        while True:
            cur = self.f.tell()
            line = self.f.readline()
            if self.f.tell() == cur:
                return "NA"  # raise Exception, "EOF reached and {0}:{1} not seen!".format(chrom, pos)

            raw = line.strip().split("\t")
            _chrom, start = (raw[0], int(raw[1]) - self.start_base)

            if _chrom != chrom or start > pos:
                return "NA"  # raise Exception, "End of chrom reached and {0}:{1} not seen!".format(chrom, pos)

            if start == pos:
                return raw[3]
