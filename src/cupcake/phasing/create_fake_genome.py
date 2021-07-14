#!/usr/bin/env python
__author__ = "etseng@pacb.com"

import re
from enum import Enum
from pathlib import Path
from typing import Optional

import typer
from Bio import SeqIO
from bx.intervals.cluster import ClusterTree

import cupcake.sequence.GFF as GFF
from cupcake import version_callback
from cupcake.logger import cupcake_logger as logger

"""
Given a GFF file and selected loci (ex: PB.45),
create a fake genome that is the concatenation of all seen exons.
"""

extra_bp_around_junctions = 50  # get this much around junctions to be safe AND to not screw up GMAP who doesn't like microintrons....
__padding_before_after__ = 10  # get this much before and after the start


class strand_direction(str, Enum):
    plus = "+"
    minus = "-"


app = typer.Typer(name="cupcake.phasing.create_fake_genome")


def make_fake_genome(
    genome_filename,
    gff_filename,
    ref_chr,
    ref_start,
    ref_end,
    ref_strand,
    output_prefix,
    output_name=None,
    genome_d=None,
):
    if genome_d is None:
        logger.info(f"Reading genome file {genome_filename}...")
        d = SeqIO.to_dict(SeqIO.parse(open(genome_filename), "fasta"))
    else:
        d = genome_d

    if output_name is None:
        output_name = f"fake_{genome_filename}"

    logger.info(f"Reading GFF file {gff_filename}...")
    good = []
    reader = GFF.collapseGFFReader(gff_filename)
    for r in reader:
        if (
            r.chr == ref_chr
            and r.strand == ref_strand
            and (ref_start <= r.start < r.end <= ref_end)
            and len(r.ref_exons) > 1
        ):
            logger.info(f"Adding {r.seqid} to fake genome.")
            good.append(r)

    if len(good) == 0:
        raise RuntimeError(
            f"Did not find any transcripts strictly within {ref_chr}:{ref_start}-{ref_end} on strand {ref_strand}. Abort!"
        )

    c = ClusterTree(0, 0)
    for r in good:
        for e in r.ref_exons:
            c.insert(
                e.start - extra_bp_around_junctions,
                e.end + extra_bp_around_junctions,
                1,
            )

    regions = [(a, b) for (a, b, junk) in c.getregions()]
    regions[0] = (regions[0][0] - __padding_before_after__, regions[0][1])
    regions[-1] = (regions[-1][0], regions[-1][1] + __padding_before_after__)

    with open(output_prefix + ".fasta", "w") as f:
        f.write(">" + output_name + "\n")
        for a, b in regions:
            f.write(str(d[r.chr][a:b].seq))
        f.write("\n")

    # for mapping, write <0-based index on fake genome>, <ref chrom>, <0-based index on ref genome>
    with open(output_prefix + ".mapping.txt", "w") as f:
        i = 0
        for a, b in regions:
            for j in range(a, b):
                f.write(f"{i},{ref_chr},{j}\n")
                i += 1

        with open(output_prefix + ".pbids.txt", "w") as f:
            f.write("\n".join(r.seqid for r in good) + "\n")

    logger.info(
        f"Output written to {output_prefix}.fasta, {output_prefix}.mapping.txt, {output_prefix}.pbids.txt."
    )


@app.command(name="")
def main(
    genome_filename: str = typer.Argument(...),
    gff_filename: str = typer.Argument(...),
    locus: str = typer.Option(..., help="locus in format <chr>:<start>-<end>"),
    strand: strand_direction = typer.Option(..., help="strand of locus"),
    output_prefix: str = typer.Option(
        ..., "--output_prefix", "-o", help="Output prefix"
    ),
    output_name: Optional[str] = typer.Option(None, help="Output genome name"),
    version: bool = typer.Option(
        None,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Prints the version of the SQANTI3 package.",
    ),
) -> None:

    rex = re.compile(r"(\S+):(\d+)-(\d+)")
    m = rex.match(locus)
    if m is None:
        raise RuntimeError(f"{locus} is not a defined chr location! Abort.")

    ref_chr = m.group(1)
    ref_start = int(m.group(2)) - 1  # make it 0-based
    ref_end = int(m.group(3))  # keep it 1-based

    genome_filename = Path(genome_filename)
    gff_filename = Path(gff_filename)

    if not genome_filename.exists():
        raise FileNotFoundError(f"Genome file {genome_filename} does not exist! Abort.")
    if not gff_filename.exists():
        raise FileNotFoundError(f"GFF {gff_filename} does not exist! Abort.")

    make_fake_genome(
        genome_filename,
        gff_filename,
        ref_chr,
        ref_start,
        ref_end,
        strand,
        output_prefix,
        output_name,
    )


if __name__ == "__main__":
    typer.run(main)
