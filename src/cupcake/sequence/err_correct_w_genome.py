#!/usr/bin/env python

from pathlib import Path
from typing import Dict

import typer
from Bio import SeqIO

from cupcake.__about__ import __version__
from cupcake.logging import cupcake_logger as logger
from cupcake.sequence import BioReaders
from cupcake.sequence.coordinate_mapper import consistute_genome_seq_from_exons

###### MODIFY FILENAME BELOW #######
# genome_file = 'hg38.fa'
# sam_file = 'alz.rep.fq.sam'
# sam_file = 'isoseq_flnc.fasta.sam'
# output_err_corrected_fasta = 'alz.rep.genome_corrected.fasta'
# output_err_corrected_fasta = 'isoseq_flnc.genome_corrected.fasta'
####################################


app = typer.Typer(
    name="cupcake.sequence.err_correct_w_genome",
    help="Generate sequences using genome bases and SAM alignment file",
)


def err_correct(
    genome_file: Path,
    sam_file: Path,
    output_err_corrected_fasta: Path,
    genome_dict: Dict[str, SeqIO.SeqRecord] = None,
) -> None:
    if genome_dict is None:
        genome_dict = {}
        for r in SeqIO.parse(open(genome_file, "r"), "fasta"):
            genome_dict[r.name] = r
        logger.info(f"done reading {genome_file}")

    with open(output_err_corrected_fasta, "w") as f:
        reader = BioReaders.GMAPSAMReader(str(sam_file), True)
        for r in reader:
            if r.sID == "*":
                continue
            seq = consistute_genome_seq_from_exons(
                genome_dict, r.sID, r.segments, r.flag.strand
            )
            f.write(f">{r.qID}\n{seq}\n")

    logger.info(f"output written to {output_err_corrected_fasta}")


@app.command(name="")
def main(
    genome_file: str = typer.Argument(..., help="Genome Fasta File"),
    sam_file: str = typer.Argument(..., help="GMAP SAM File"),
    output_file: str = typer.Argument(..., help="Output Fasta File"),
) -> None:
    err_correct(Path(genome_file), Path(sam_file), Path(output_file))


if __name__ == "__main__":
    typer.run(main)
