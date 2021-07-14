#!/usr/bin/env python
__author__ = "etseng@pacb.com"
"""
Wrapper for running STARlong.
Parameters are pre-set according to:


"""
import shutil
import subprocess
import tempfile
from pathlib import Path

import typer

from cupcake import version_callback

app = typer.Typer(name="cupcake.sequence.STARwrapper", help="Wrapper for running STAR")


CMD_STARlong = "/home/UNIXHOME/etseng/software_downloads/STAR-2.5.3a/bin/Linux_x86_64/STAR --runMode alignReads --outSAMattributes NH HI NM MD --readNameSeparator space --outFilterMultimapScoreRange 1 --outFilterMismatchNmax 2000 --scoreGapNoncan -1  --scoreGapGCAG -4 --scoreGapATAC -8 --scoreDelOpen -1 --scoreDelBase -1 --scoreInsOpen -1 --scoreInsBase -1 --alignEndsType Local --seedSearchStartLmax 50 --seedPerReadNmax 100000 --seedPerWindowNmax 1000 --alignTranscriptsPerReadNmax 100000 --alignTranscriptsPerWindowNmax 10000"
CMD_STAR2_format = (
    CMD_STARlong
    + " --twopassMode None --runThreadN {c} --genomeDir {d} --readFilesIn {i}"
)


def run_STAR(in_fasta, out_sam, genome_dir, cpus):
    with tempfile.mkdtemp(prefix="STARtmp") as tmp_dir:
        in_fasta = Path(in_fasta)
        out_sam = Path(out_sam)
        cmd = CMD_STAR2_format.format(c=cpus, d=genome_dir, i=in_fasta)
        if subprocess.check_call(cmd, shell=True, cwd=tmp_dir) != 0:
            raise subprocess.CalledProcessError(f"ERROR RUNNING CMD: {cmd}")

        shutil.move(Path(tmp_dir, "Aligned.out.sam"), out_sam)


@app.command(name="")
def main(
    genome_dir: str = typer.Argument(...),
    in_fasta: str = typer.Argument(...),
    out_sam: str = typer.Argument(...),
    cpus: int = typer.Option(10, help="Number of threads (default: 10)"),
    version: bool = typer.Option(
        None,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Prints the version of the SQANTI3 package.",
    ),
) -> None:

    run_STAR(in_fasta, out_sam, genome_dir, cpus)


if __name__ == "__main__":
    typer.run(main)
