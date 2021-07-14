#!/usr/bin/env python

from collections import defaultdict
from pathlib import Path

import typer
from Bio import SeqIO

from cupcake import version_callback
from cupcake.logger import cupcake_logger as logger
from cupcake.sequence import BioReaders

app = typer.Typer(
    name="cupcake.sequence.summarize_gmap_sam",
    help="Summarize GMAP SAM file in tab-delimited file.",
)


def type_fa_or_fq(filename):
    filename = Path(filename)
    if filename.stem.upper() not in (".FA", ".FASTA"):
        return "fasta"
    else:
        return "fastq"


def summarize_GMAP_sam(input_fa_or_fq, input_sam):
    d = {
        r.id: len(r.seq)
        for r in SeqIO.parse(open(input_fa_or_fq), type_fa_or_fq(input_fa_or_fq))
    }

    map_count = defaultdict(lambda: 0)
    for r in BioReaders.GMAPSAMReader(input_sam, True):
        map_count[r.qID] += 1
    multi = [x for x in map_count if map_count[x] > 1]

    with open(f"{input_sam}.summary.txt", "w") as f:
        f.write(
            "id\tqLength\tqCoverage\tidentity\tnum_nonmatch\tnum_ins\tnum_del\tunique\n"
        )
        for r in BioReaders.GMAPSAMReader(input_sam, True, query_len_dict=d):
            if r.sID == "*":
                continue
            if r.qID in multi:
                uni = "N"
            else:
                uni = "Y"
            f.write(
                f"{r.qID}\t{d[r.qID]}\t{r.qCoverage:.4f}\t"
                f"{r.identity:.4f}\t\t{r.num_nonmatches}\t"
                f"{r.num_ins}\t{r.num_del}\t{uni}\n"
            )

        logger.info(f"Output written to: {f.name}")


@app.command(name="")
def main(
    input_fa_or_fq: str = typer.Argument(..., help="Input fasta/fastq filename"),
    sam_file: str = typer.Argument(..., help="(GMAP) SAM filename"),
    version: bool = typer.Option(
        None,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Prints the version of the SQANTI3 package.",
    ),
) -> None:
    summarize_GMAP_sam(input_fa_or_fq, sam_file)


if __name__ == "__main__":
    typer.run(main)
