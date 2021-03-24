#!/usr/bin/env python

__version__ = "1.0"


import random
import sys

from Bio import SeqIO

from cupcake.logging import cupcake_logger as logger


def type_fa_or_fq(file):
    file = file.upper()
    if file.endswith(".FA") or file.endswith(".FASTA"):
        return "fasta"
    else:
        return "fastq"


def sep_by_primer(filename, output_prefix, sample_size):
    filetype = type_fa_or_fq(filename)

    ids = [r.id for r in SeqIO.parse(open(filename), filetype)]

    n = len(ids)
    if sample_size > n:
        logger.warning(f"WARNING: {filename} contains only {n} sequences but subsample size at {sample_size}! Simply output whole file.")

    chosen_ids = random.sample(ids, min(n, sample_size))

    with open(f"{output_prefix}.random{str(sample_size)}.{filetype}", "w") as f:
        for r in SeqIO.parse(open(filename), filetype):
            if r.id in chosen_ids:
                SeqIO.write(r, f, filetype)

        logger.info(f"Randomly selected sequences written to {f.name}.")


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser("Randomly select N sequences from fasta/fastq files")
    parser.add_argument("filename", help="Input fasta/fastq filename")
    parser.add_argument("output_prefix", help="Output file prefix")
    parser.add_argument("sample_size", type=int, help="Subsample size")
    args = parser.parse_args()

    sep_by_primer(args.filename, args.output_prefix, args.sample_size)


if __name__ == "__main__":
    main()
