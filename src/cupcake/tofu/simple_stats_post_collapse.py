#!/usr/bin/env python

import sys
from pathlib import Path

import typer
from cupcake.logging import cupcake_logger as logger
from cupcake.sequence.GFF import collapseGFFReader

app = typer.Typer(name="cupcake.tofu.simple_stats_post_collapse")


def simple_stats_post_collapse(input_prefix: str) -> None:
    input_gff = f"{input_prefix}.gff"

    if not Path(input_gff).exists():
        logger.error(f"Looking for input GFF {input_gff} but not found! Abort!")
        sys.exit(-1)

    simple_stats = Path(f"{input_prefix}.simple_stats.txt")
    exon_stats = Path(f"{input_prefix}.exon_stats.txt")

    with simple_stats.open("w") as f1, exon_stats.open("w") as f2:
        f1.write("pbid\tlocus\tlength\tnum_exon\n")
        f2.write("pbid\texon_index\texon_size\tintron_size\n")

        for r in collapseGFFReader(input_gff):
            f1.write(f"{r.seqid}	")
            f1.write(f"{r.seqid.split('.')[1]}	")
            sum_len = 0
            for i, e in enumerate(r.ref_exons):
                exon_len = e.end - e.start
                sum_len += exon_len
                f2.write(f"{r.seqid}\t{i+1}\t{exon_len}\t")
                if i == 0:
                    f2.write("NA\n")
                else:
                    f2.write(f"{str(e.start - r.ref_exons[i - 1].end)}\n")

            f1.write(f"{str(sum_len)}\t")
            f1.write(f"{str(len(r.ref_exons))}\n")
    logger.info(f"Output written to: {simple_stats.name},{exon_stats.name}\n")


@app.command(name="")
def main(
    input_prefix: str = typer.Argument(
        ..., help="Input prefix, ex: hq.5merge.collapsed"
    )
) -> None:
    simple_stats_post_collapse(input_prefix)


if __name__ == "__main__":
    typer.run(main)
