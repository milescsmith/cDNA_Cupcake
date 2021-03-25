#!/usr/bin/env python
__author__ = "etseng@pacb.com"

import re
import sys
from collections import Counter, defaultdict
from csv import DictReader
from pathlib import Path
from typing import Any
from typing import Counter as counter
from typing import Dict, Optional, Set, Tuple, Union

import typer
from Bio import SeqIO
from cupcake.logging import cupcake_logger as logger

"""
Demultiplex IsoSeq (SMRT Link 8.0) job output (without genome mapping)
"""

# hq1_id_rex = re.compile('(i\d+_HQ_\S+\|\S+)\/f\d+p\d+\/\d+')
# hq2_id_rex = re.compile('HQ_\S+\|(\S+)\/f\d+p\d+\/\d+')
hq3_id_rex = re.compile(r"[\S]+(transcript/\d+)")  # ex: UnnamedSample_HQ_transcript/0


app = typer.Typer(name="cupcake.post_isoseq_cluster.demux_isoseq_no_genome")


def type_fafq(fafq: str) -> str:
    x = fafq.upper()
    if x.endswith(".FA") or x.endswith(".FASTA"):
        return "fasta"
    elif x.endswith(".FQ") or x.endswith(".FASTQ"):
        return "fastq"
    else:
        raise Exception(
            f"HQ fasta/fastq filename must end with .fasta or .fastq! Saw {fafq} instead, abort!"
        )


def link_files(src_dir, out_dir="./") -> Tuple[Path, Path, Path, Path]:
    """
    :param src_dir: job directory
    Locate HQ isoform, (cluster) report.csv, (classify) file.csv link to current directory
    """
    # location for mapped fastq in IsoSeq3
    src_dir = Path(src_dir)
    hq_fasta = src_dir.joinpath("output", "hq_transcript.fasta")
    hq_fastq = src_dir.joinpath("outputs", "hq_transcripts.fastq")
    primer_csv = src_dir.joinpath("outputs", "flnc.report.csv")
    cluster_csv = src_dir.joinpath("outputs", "polished.cluster_report.csv")

    if hq_fasta.exists() or hq_fastq.exists():
        logger.info("Detecting IsoSeq directories...")
    else:
        raise FileNotFoundError(
            "Cannot find hq_transcripts.fasta/fastq in job directory! Does not look like SMRTLink 8 Iso-Seq job!",
        )
    return (
        out_dir,
        hq_fastq if hq_fastq.exists() else hq_fasta,
        cluster_csv,
        primer_csv,
    )


def read_cluster_csv(cluster_csv, classify_info):
    """
    :param report_csv: cluster_report.csv
    :return: dict of cluster_id --> integer primer --> FL count
    """
    info = defaultdict(lambda: Counter())
    for r in DictReader(open(cluster_csv), delimiter=","):
        assert r["read_type"] == "FL"  # always FL for isoseq3
        p = classify_info[r["read_id"]]
        cid = r["cluster_id"]
        info[cid][p] += 1
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


def demux_isoseq_no_genome(
    job_dir: Path = None,
    hq_fafq: Path = None,
    cluster_csv: Path = None,
    classify_csv: Path = None,
    output_filename=sys.stdout,
    primer_names=None,
) -> None:
    if job_dir is not None:
        out_dir_ignore, hq_fafq, cluster_csv, classify_csv = link_files(job_dir)
    else:
        for _ in (hq_fafq, cluster_csv, classify_csv):
            if not _.exists():
                raise FileNotFoundError(f"{_.name} cannot be found!")

    # info: dict of hq_isoform --> primer --> FL count
    logger.info(f"Reading {classify_csv}....")
    primer_list, classify_csv = read_classify_csv(classify_csv)

    logger.info(f"Reading {cluster_csv}....")
    info: Union[Dict[Any, counter[Any]], Dict[Any, Any]] = read_cluster_csv(
        cluster_csv, classify_csv
    )

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
        f.write(f"id,{','.join(list(primer_names.keys()))}\n")
        logger.info(f"Reading {hq_fafq}....")
        for r in SeqIO.parse(open(hq_fafq), type_fafq(hq_fafq)):
            f.write(r.id)
            m = hq3_id_rex.match(r.id)
            cid = m.group(1)
            for p in primer_names:
                f.write(f",{info[cid][p]}")
            f.write("\n")

        logger.info(f"Count file written to {f.name}.")


@app.command(name="")
def main(
    job_dir: Optional[str] = typer.Option(
        None,
        "--job_dir",
        "-j",
        help="Job directory (if given, automatically finds required files)",
    ),
    hq_fafq: Optional[str] = typer.Option(
        None,
        help="HQ isoform fasta/fastq (overridden by --job_dir if given)",
    ),
    cluster_csv: Optional[str] = typer.Option(
        None,
        help="Cluster report CSV (overridden by --job_dir if given)",
    ),
    classify_csv: Optional[str] = typer.Option(
        None,
        help="Classify report CSV (overriden by --job_dir if given)",
    ),
    primer_names: Optional[str] = typer.Option(
        None,
        help="Text file showing primer sample names (default: None)",
    ),
    output: str = typer.Option(
        sys.stdout, "--output", "-o", help="Output count filename"
    ),
) -> None:

    job_dir = Path(job_dir) if job_dir is not None else job_dir
    hq_fafq = Path(hq_fafq) if hq_fafq is not None else hq_fafq
    cluster_csv = Path(cluster_csv) if cluster_csv is not None else cluster_csv
    classify_csv = Path(classify_csv) if classify_csv is not None else classify_csv

    if primer_names is not None:
        primer_names = {}
        for line in open(primer_names):
            index, name = line.strip().split()
            primer_names[index] = name
    else:
        primer_names = None

    demux_isoseq_no_genome(
        job_dir,
        hq_fafq,
        cluster_csv,
        classify_csv,
        output,
        primer_names,
    )


if __name__ == "__main__":
    typer.run(main)
