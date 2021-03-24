#!/usr/bin/env python
__author__ = "etseng@pacb.com"

"""
Demultiplex IsoSeq (SMRT Link 8.0) job output (with genome mapping)
"""

import re
import sys
from collections import Counter, defaultdict
from csv import DictReader
from pathlib import Path
from Bio import SeqIO
from typing import Tuple, Dict, Optional, Set

import typer

from cupcake.logging import cupcake_logger as logger

mapped_id_rex = re.compile(r"(PB.\d+.\d+)")


app = typer.Typer(name="cupcake.post_isoseq_cluster.demux_isoseq_with_genome")


def type_fafq(fafq: str) -> str:
    x = fafq.upper()
    if x.endswith(".FA") or x.endswith(".FASTA"):
        return "fasta"
    elif x.endswith(".FQ") or x.endswith(".FASTQ"):
        return "fastq"
    else:
        raise Exception(
            f"Mapped fasta/fastq filename must end with .fasta or .fastq! Saw {fafq} instead, abort!"
        )


def link_files(src_dir: str, out_dir=Path.cwd()) -> Tuple[Path, Path, Path, Path]:
    """
    :param src_dir: job directory
    Locate mapped.fastq, read-stat, classify report link to current directory
    """

    src_dir = Path(src_dir)
    # location for mapped fastq in IsoSeq3
    mapped_fastq = src_dir.joinpath("outputs", "collapse_isoforms.fastq")  # for <SL8
    mapped_fasta = src_dir.joinpath("outputs", "collapse_isoforms.fasta")  # SL8+ only fasta
    # mapped_gff = os.path.join(
    #     os.path.abspath(src_dir), "outputs", "collapse_isoforms.gff"
    # )
    read_stat = src_dir.joinpath("outputs", "collapse_isoforms.read_stat.txt")
    primer_csv = src_dir.joinpath("outputs", "flnc.report.csv")

    if mapped_fastq.exists():
        logger.info("Detecting IsoSeq task directories...")
        return out_dir, mapped_fastq, read_stat, primer_csv
    elif mapped_fasta.exists():
        logger.info("Detecting IsoSeq task directories...")
        return out_dir, mapped_fasta, read_stat, primer_csv
    else:
        raise FileNotFoundError("Cannot find expected files (ex: collapse_isoforms.fastq) in job directory! Does not look like a Iso-Seq job!")


def read_read_stat(read_stat: Path, classify_info: Dict[str, int]):
    """
    :return: dict of pbid --> (int) primer --> FL count
    """
    info = defaultdict(lambda: Counter())
    for r in DictReader(open(read_stat), delimiter="\t"):
        p = classify_info[r["id"]]
        info[r["pbid"]][p] += 1
    return dict(info)


def read_classify_csv(classify_csv: Path) -> Tuple[Set[str], Dict[str, str]]:
    """
    :param classify_csv: classify report csv
    :return: primer list, dict of FL id --> primer
    """
    info = {}
    primer_list = set()
    for r in DictReader(open(classify_csv), delimiter=","):
        p = r["primer"]
        primer_list.add(p)
        if r["id"] in info:
            raise Exception(f"{r['id']} listed more than once in {classify_csv}!")
        info[r["id"]] = p
    return primer_list, info


def demux_isoseq_with_genome(
    job_dir: Optional[str] = None,
    mapped_fafq: Optional[str] = None,
    read_stat: Optional[str] = None,
    classify_csv: Optional[str] = None,
    output_filename=sys.stdout,
    primer_names: Optional[str] = None,
) -> None:
    mapped_fafq = Path(mapped_fafq)
    read_stat = Path(read_stat)
    classify_csv = Path(classify_csv)

    if job_dir is not None:
        _, mapped_fafq, read_stat, classify_csv = link_files(job_dir)
    else:
        for _ in (mapped_fafq, read_stat, classify_csv):
            if not _.exists():
                raise FileNotFoundError(f"Cannot find {_.name}")

    # info: dict of hq_isoform --> primer --> FL count
    logger.info(f"Reading {classify_csv}...")
    primer_list, classify_info = read_classify_csv(classify_csv)
    logger.info(f"Reading {read_stat}...")
    info = read_read_stat(read_stat, classify_info)

    primer_list = list(primer_list)
    primer_list.sort()
    # if primer names are not given, just use as is...
    tmp_primer_names = {x: x for x in primer_list}
    if primer_names is None:
        primer_names = tmp_primer_names
    else:
        for k, v in tmp_primer_names.items():
            if k not in primer_names:
                primer_names[k] = v

    with open(output_filename, "w") as f:
        f.write(f"id,{','.join(list(primer_names.values()))}\n")
        logger.info(f"Reading {mapped_fafq}....")
        for r in SeqIO.parse(open(mapped_fafq), type_fafq(mapped_fafq)):
            m = mapped_id_rex.match(r.id)  # expected ID: PB.X.Y|xxxx.....
            if m is None:
                raise Exception(f"Expected ID format PB.X.Y but found {r.id}!")
            pbid = m.group(1)
            f.write(pbid)
            for p in primer_names:
                f.write(f",{info[pbid][p]}")
            f.write("\n")
        logger.info(f"Count file written to {f.name}.")


@app.command(name="")
def main(
    job_dir: str = typer.Argument(
        ...,
        "--job_dir",
        "-j",
        help="Job directory (if given, automatically finds required files)",
    ),
    mapped_fafq: str = typer.Argument(..., help="mapped fasta/fastq (overridden by --job_dir if given)"),
    read_stat: str = typer.Argument(..., help="read_stat txt (overridden by --job_dir if given)"),
    classify_csv: str = typer.Argument(..., help="Classify report CSV (overriden by --job_dir if given)"),
    primer_names: Optional[str] = typer.Option(None, help="Text file showing primer sample names (default: None)",),
    output: str = typer.Argument(..., "--output", "-o", help="Output count filename")
):
    if primer_names is not None:
        primer_names = {}
        for line in open(primer_names):
            index, name = line.strip().split()
            primer_names[index] = name
    else:
        primer_names = None

    demux_isoseq_with_genome(
        job_dir,
        mapped_fafq,
        read_stat,
        classify_csv,
        output,
        primer_names,
    )


if __name__ == "__main__":
    typer.run(main)
