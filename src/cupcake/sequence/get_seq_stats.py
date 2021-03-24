#!/usr/bin/env python

__version__ = "1.0"

import numpy as np
from Bio import SeqIO


def type_fa_or_fq(file):
    file = file.upper()
    if file.endswith(".FA") or file.endswith(".FASTA") or file.endswith("CLIPS"):
        return "fasta"
    else:
        return "fastq"


def get_seq_stats(file, binwidth):
    print("file type is:", type_fa_or_fq(file))

    f = open(f"{file}.seqlengths.txt", "w")
    lens = []
    for r in SeqIO.parse(open(file), type_fa_or_fq(file)):
        f.write(f"{r.id}	{str(len(r.seq))}\n")
        lens.append(len(r.seq))
    f.close()

    print(f"{len(lens)} sequences")
    print("min:", min(lens))
    print("max:", max(lens))
    print("avg:", sum(lens) * 1.0 / len(lens))

    # print by 1 kb bins
    print("Length Breakdown by kb range:")

    _max = (max(lens) // binwidth) + 1
    bin = [0] * _max
    for x in lens:
        bin[x // binwidth] += 1

    for i in range(0, _max):
        if binwidth == 1000:
            print(f"{i}-{i + 1} kb: {bin[i]}")
        else:
            print(f"{i * binwidth}-{(i + 1) * binwidth}: {bin[i]}")

    print("5-95% percentile:", np.percentile(lens, 5), np.percentile(lens, 95))


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser("Summarize sequence lengths in fasta/fastq")
    parser.add_argument("filename", help="Input fasta/fastq filename")
    parser.add_argument(
        "-b",
        "--binwidth",
        default=1000,
        type=int,
        help="Bin width, in bp (default: 1000 bp)",
    )
    args = parser.parse_args()

    file = args.filename

    get_seq_stats(file, args.binwidth)


if __name__ == "__main__":
    main()
