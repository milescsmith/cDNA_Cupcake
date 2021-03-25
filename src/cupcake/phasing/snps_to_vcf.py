#!/usr/bin/env python
__author__ = "etseng@pacb.com"

from pathlib import Path

import typer
from cupcake.logging import cupcake_logger as logger
from cupcake.phasing.io.MummerSNPReader import write_snp_to_vcf
from cupcake.sequence.SeqReaders import LazyFastaReader

app = typer.Typer(
    name="cupcake.phasing.snps_to_vcf",
    help="Process one or more .snps_files files from dnadiff to VCF format.",
)


@app.command(name="")
def main(
    snps_filename: str = typer.Argument(
        ..., help="Filename containing list of .snps_files to process."
    ),
    genome_filename: str = typer.Argument(
        ..., help="Genome fasta. Chromosome IDs must agree with .snps_files files!"
    ),
) -> None:
    snps_filename = Path(snps_filename)
    genome_filename = Path(genome_filename)

    snps_files = []
    # sanity checking of input files
    for line in snps_filename.open():
        filename = Path(line.strip())
        if filename.suffix != "snps":
            raise FileNotFoundError(
                f"Input files listed in {snps_filename} must end with .snps_files!"
            )
        if not filename.exists():
            raise FileNotFoundError(f"{filename} does not exist! Abort.")
        snps_files.append(filename)

    if not genome_filename.exists():
        raise FileNotFoundError(f"Genome file {genome_filename} does not exist!")

    logger.info(f"Reading genome file {genome_filename}...")
    genome_d = LazyFastaReader(genome_filename)

    # quick checking if the genome chromosomes have the |arrow|arrow style suffix, if they do, process it
    keys = list(genome_d.keys())
    for k in keys:
        k2 = k.split("|")[0]
        if k2 != k and k2 not in keys:
            genome_d.d[k2] = genome_d.d[k]
            logger.info(f"Detected | string in chromosome ID, stripping {k} to {k2}...")
    logger.info("Finished reading genome.")

    for snp_file in snps_files:
        assert snp_file.suffix(".snps")
        if snp_file.st_size == 0:
            logger.info(f"Skipping {snp_file} because empty file.")
            continue
        vcf_file = snp_file.with_suffix(".vcf")
        logger.info(f"Processing {snp_file} --> {vcf_file}")
        write_snp_to_vcf(snp_file, vcf_file, genome_filename, genome_d)


if __name__ == "__main__":
    typer.run(main)
