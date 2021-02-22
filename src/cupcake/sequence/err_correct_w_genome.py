#!/usr/bin/env python

__version__ = "1.0"

import sys
from pathlib import Path
from typing import Dict

from Bio import SeqIO

from cupcake.sequence import BioReaders
from cupcake.sequence.coordinate_mapper import consistute_genome_seq_from_exons

###### MODIFY FILENAME BELOW #######
# genome_file = 'hg38.fa'
# sam_file = 'alz.rep.fq.sam'
# sam_file = 'isoseq_flnc.fasta.sam'
# output_err_corrected_fasta = 'alz.rep.genome_corrected.fasta'
# output_err_corrected_fasta = 'isoseq_flnc.genome_corrected.fasta'
####################################


def err_correct(
    genome_file: Path,
    sam_file: Path,
    output_err_corrected_fasta: Path,
    genome_dict: Dict[str, SeqIO.SeqRecord] = None,
) -> None:
    if genome_dict is None:
        genome_dict = {}
        for r in SeqIO.parse(genome_file.open(), "fasta"):
            genome_dict[r.name] = r
        logger.info(f"done reading {genome_file}")

    with output_err_corrected_fasta.open("w") as f:
        reader = BioReaders.GMAPSAMReader(str(sam_file), True)
        for r in reader:
            if r.sID == "*":
                continue
            seq = consistute_genome_seq_from_exons(
                genome_dict, r.sID, r.segments, r.flag.strand
            )
            f.write(f">{r.qID}\n{seq}\n")

    logger.info(f"output written to {output_err_corrected_fasta}")


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser(
        "Generate sequences using genome bases and SAM alignment file"
    )
    parser.add_argument("genome_file", help="Genome Fasta File")
    parser.add_argument("sam_file", help="GMAP SAM File")
    parser.add_argument("output_file", help="Output Fasta File")
    args = parser.parse_args()

    err_correct(Path(args.genome_file), Path(args.sam_file), Path(args.output_file))


if __name__ == "__main__":
    main()
