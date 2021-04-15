#!/usr/bin/env python
__author__ = "etseng@pacb.com"

"""
Demultiplex IsoSeq1/IsoSeq2 job output (without genome mapping)

INPUT: HQ isoform and report.csv and file.csv (alternatively, job directory)
OUTPUT: CSV file containing associated FL count for each isoform

HQ isoform ID format:
  (isoseq1) i0_HQ_sample3f1db2|c44/f3p0/2324
  (isoseq2) HQ_sampleAZxhBguy|cb1063_c22/f3p0/6697
Cluster report format:
   (isoseq1)
    cluster_id,read_id,read_type
    i0_ICE_sample3f1db2|c23,m54033_171031_152256/26476965/30_5265_CCS,FL
    i0_ICE_sample3f1db2|c43,m54033_171031_152256/32441242/25_2283_CCS,FL
    i0_ICE_sample3f1db2|c43,m54033_171031_152256/50201527/30_2323_CCS,FL
    i0_ICE_sample3f1db2|c44,m54033_171031_152256/49545897/2374_68_CCS,FL
   (isoseq2)
    cluster_id,read_id,read_type
    cb10063_c6,m54006_170729_232022/56361426/29_9138_CCS,FL
    cb10407_c1,m54006_170729_232022/16712197/2080_72_CCS,FL
    cb10467_c49,m54006_170729_232022/48104064/0_2421_CCS,NonFL
Classify report format:
   (isoseq1 and 2)
    id,strand,fiveseen,polyAseen,threeseen,fiveend,polyAend,threeend,primer,chimera
    m54033_171031_152256/27919078/31_1250_CCS,+,1,1,1,31,1250,1280,3,0
    m54033_171031_152256/27919079/31_3840_CCS,+,1,1,1,31,3840,3869,2,0
    m54033_171031_152256/27919086/29_3644_CCS,+,1,1,1,29,3644,3674,2,0
"""

import re
import sys
from collections import Counter, defaultdict
from csv import DictReader
from pathlib import Path
from typing import Optional

import typer
from Bio import SeqIO

from cupcake.logging import cupcake_logger as logger

hq1_id_rex = re.compile(r"i\d+_HQ_\S+\|(\S+)\/f\d+p\d+\/\d+")
hq2_id_rex = re.compile(r"HQ_\S+\|(\S+)\/f\d+p\d+\/\d+")


app = typer.Typer(name="cupcake.post_isoseq_cluster.demux_isoseq2_no_genome")


def link_files(src_dir, out_dir=Path.cwd()):
    """
    :param src_dir: job directory
    Locate HQ isoform, (cluster) report.csv, (classify) file.csv link to current directory
    """
    src_dir = Path(src_dir)
    # location for HQ fastq in IsoSeq1
    hq_fastq = src_dir.joinpath(
        "tasks",
        "pbtranscript.tasks.combine_cluster_bins-0",
        "hq_isoforms.fastq",
    )
    # location for HQ fastq in IsoSeq2
    hq_fastq2 = src_dir.joinpath(
        "tasks",
        "pbtranscript2tools.tasks.collect_polish-0",
        "all_arrowed_hq.fastq",
    )
    # location for cluster report in IsoSeq1
    cluster_csv = src_dir.joinpath(
        "tasks",
        "pbtranscript.tasks.combine_cluster_bins-0",
        "cluster_report.csv",
    )
    cluster_csv2 = src_dir.joinpath(
        "tasks",
        "pbtranscript2tools.tasks.collect_polish-0",
        "report.csv",
    )
    # location for classify report in IsoSeq1 and 2
    primer_csv = src_dir.joinpath("tasks", "pbcoretools.tasks.gather_csv-1", "file.csv")

    if hq_fastq.exists():
        logger.info("Detecting IsoSeq1 task directories...")
        hq_fastq.symlink_to(out_dir.joinpath("hq_isoforms.fastq"))
        cluster_csv.symlink_to(out_dir.joinpath("cluster_report.csv"))
        primer_csv.symlink_to(out_dir.joinpath("classify_report.csv"))
        isoseq_version = "1"
    else:
        logger.info("Detecting IsoSeq2 task directories...")
        hq_fastq2.symlink_to(out_dir.joinpath("hq_isoforms.fastq"))
        cluster_csv2.symlink_to(out_dir.joinpath("cluster_report.csv"))
        primer_csv.symlink_to(out_dir.joinpath("classify_report.csv"))
        isoseq_version = "2"
    return (
        out_dir,
        "hq_isoforms.fastq",
        "cluster_report.csv",
        "classify_report.csv",
        isoseq_version,
    )


def read_cluster_csv(cluster_csv, classify_info, isoseq_version):
    """
    :param report_csv: cluster_report.csv
    :return: dict of cluster_id --> integer primer --> FL count
    """
    info = defaultdict(lambda: Counter())
    for r in DictReader(open(cluster_csv), delimiter=","):
        if r["read_type"] == "FL":
            p = classify_info[r["read_id"]]
            if isoseq_version == "1":
                cid = r["cluster_id"].split("|")[1]
            else:
                cid = r["cluster_id"]
            info[cid][p] += 1
    return dict(info)


def read_classify_csv(classify_csv):
    """
    :param classify_csv: classify report csv
    :return: primer range, dict of FL/nFL id --> primer (in integer)
    """
    info = {}
    max_p = 0
    for r in DictReader(open(classify_csv), delimiter=","):
        if r["primer"] == "NA":
            continue  # skip nFL
        p = int(r["primer"])
        max_p = max(max_p, p)
        info[r["id"]] = p
    return max_p, info


def demux_isoseq2_no_genome(
    job_dir: Optional[Path] = None,
    hq_fastq: Optional[Path] = None,
    cluster_csv: Optional[Path] = None,
    classify_csv: Optional[Path] = None,
    output_filename=sys.stdout,
):

    if job_dir is not None:
        (
            out_dir_ignore,
            hq_fastq,
            cluster_csv,
            classify_csv,
            isoseq_version,
        ) = link_files(job_dir)
        assert isoseq_version in ("1", "2")
    else:
        for _ in (
            hq_fastq,
            cluster_csv,
            classify_csv,
        ):
            if not _.exists():
                raise FileNotFoundError(f"{_.name} was not found!")

    # info: dict of hq_isoform --> primer --> FL count
    logger.info(f"Reading {classify_csv}...")
    max_primer, classify_csv = read_classify_csv(classify_csv)
    logger.info(f"Reading {cluster_csv}...")
    info = read_cluster_csv(cluster_csv, classify_csv, isoseq_version)

    with open(output_filename, "w") as f:
        f.write(f"id,{','.join('primer' + str(i) for i in range(max_primer + 1))}\n")
        logger.info(f"Reading {hq_fastq}...")
        for r in SeqIO.parse(open(hq_fastq), "fastq"):
            if isoseq_version == "1":
                m = hq1_id_rex.match(r.id)
            else:
                m = hq2_id_rex.match(r.id)

            if m is None:
                raise RuntimeError(f"Unexpected HQ isoform ID format: {r.id}! Abort.")
            cid = m.group(1)
            f.write(r.id)
            for p in range(max_primer + 1):
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
    hq_fastq: Optional[str] = typer.Option(
        None, help="HQ isoform fastq (overridden by --job_dir if given)"
    ),
    cluster_csv: Optional[str] = typer.Option(
        None, help="Cluster report CSV (overridden by --job_dir if given)"
    ),
    classify_csv: Optional[str] = typer.Option(
        None, help="Classify report CSV (overriden by --job_dir if given)"
    ),
    output: Optional[str] = typer.Option(
        None, "--output", "-o", help="Output count filename"
    ),
):

    job_dir = Path(job_dir) if job_dir is not None else job_dir
    hq_fastq = Path(hq_fastq) if hq_fastq is not None else hq_fastq
    cluster_csv = Path(cluster_csv) if cluster_csv is not None else cluster_csv
    classify_csv = Path(classify_csv) if classify_csv is not None else classify_csv

    demux_isoseq2_no_genome(job_dir, hq_fastq, cluster_csv, classify_csv, output)


if __name__ == "__main__":
    typer.run(main)
