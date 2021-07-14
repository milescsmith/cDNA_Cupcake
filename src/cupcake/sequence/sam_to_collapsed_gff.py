#!/usr/bin/env python
from pathlib import Path

import typer

from cupcake import version_callback
from cupcake import cupcake_logger as logger
from cupcake.sequence.BioReaders import GMAPSAMReader
from cupcake.sequence.GFF import write_collapseGFF_format

app = typer.Typer(
    name="cupcake.sequence.sam_to_collapsed_gff",
    help="Convert SAM to collapsed GFF format",
)


@app.command(name="")
def main(
    sam_filename: str = typer.Argument(...),
    version: bool = typer.Option(
        None,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Prints the version of the SQANTI3 package.",
    ),
) -> None:
    sam_filename = Path(sam_filename)
    if sam_filename.suffix != ".sam":
        raise RuntimeError("Only accepts files ending in .sam. Abort!")

    prefix = sam_filename.stem
    output_gff = f"{prefix}.collapsed.gff"

    with open(output_gff, "w") as f:
        reader = GMAPSAMReader(sam_filename, True)
        for r in reader:
            if r.sID == "*":
                continue
            r.strand = r.flag.strand
            r.geneid = r.qID
            r.seqid = r.qID
            r.chr = r.sID
            r.ref_exons = r.segments
            r.start = r.sStart
            r.end = r.sEnd
            r.cds_exons = None
            write_collapseGFF_format(f, r)

    logger.info(f"Output written to {output_gff}.")


if __name__ == "__main__":
    typer.run(main)
