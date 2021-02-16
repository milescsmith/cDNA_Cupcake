#!/usr/bin/env python
__author__ = "etseng@pacb.com"

import logging
import sys
from pathlib import Path

from cupcake.logging import setup_logging
from cupcake.phasing.io.MummerSNPReader import write_snp_to_vcf
from cupcake.sequence.SeqReaders import LazyFastaReader


def main() -> None:
    setup_logging("cDNA_Cupcake.phasing.snps_to_vcf")

    from argparse import ArgumentParser

    parser = ArgumentParser(
        "Process one or more .snps_files files from dnadiff to VCF format."
    )
    parser.add_argument(
        "snps_filename", help="Filename containing list of .snps_files to process."
    )
    parser.add_argument(
        "genome_filename",
        help="Genome fasta. Chromosome IDs must agree with .snps_files files!",
    )

    args = parser.parse_args()

    snps_filename = Path(args.snps_filename)
    genome_filename = Path(args.genome_filename)

    snps_files = []
    # sanity checking of input files
    for line in snps_filename.open():
        filename = Path(line.strip())
        if not filename.suffix(".snps"):
            logging.critical(
                f"Input files listed in {snps_filename} must end with .snps_files!"
            )
            sys.exit(-1)
        if not filename.exists():
            logging.critical(f"{filename} does not exist! Abort.")
            sys.exit(-1)
        snps_files.append(filename)

    if not genome_filename.exists():
        logging.critical(f"Genome file {genome_filename} does not exist!")

    logging.info(f"Reading genome file {genome_filename}...")
    genome_d = LazyFastaReader(genome_filename)

    # quick checking if the genome chromosomes have the |arrow|arrow style suffix, if they do, process it
    keys = list(genome_d.keys())
    for k in keys:
        k2 = k.split("|")[0]
        if k2 != k and k2 not in keys:
            genome_d.d[k2] = genome_d.d[k]
            logging.info(
                f"Detected | string in chromosome ID, stripping {k} to {k2}..."
            )
    logging.info("Finished reading genome.")

    for snp_file in snps_files:
        assert snp_file.suffix(".snps")
        if snp_file.st_size == 0:
            logging.info(f"Skipping {snp_file} because empty file.")
            continue
        vcf_file = snp_file.with_suffix(".vcf")
        logging.info(f"Processing {snp_file} --> {vcf_file}")
        write_snp_to_vcf(snp_file, vcf_file, genome_filename, genome_d)


if __name__ == "__main__":
    main()
