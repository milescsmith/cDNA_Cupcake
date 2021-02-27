#!/usr/bin/env python
__author__ = "etseng@pacb.com"

"""
Used for getting per-barcode Fl count information after IsoSeq cluster is run and HQ is mapped and collapsed.

Required Input:

isoseq_classify.primer_info.csv -- CSV output from classify
*.collapsed.group.txt, *.collapsed.rep.fq --- to get the PB IDs and group info

Optional input: primer names in a text file where each line is:

<primer_index>, <primer_name>

"""

from collections import Counter, defaultdict
from csv import DictReader
from pathlib import Path
from typing import List, Optional, Tuple

import typer

from cupcake.logging import cupcake_logger as logger

app = typer.Typer(name="cupcake.tofu.get_counts_by_barcode")


def read_classify_csv(csv_filename: str) -> Tuple[list, str]:
    """
    :param csv_filename: CSV file from IsoSeq classify, must have "id" and "primer" field
    :return: primer_ranges (list of indices), primer_info which is dict of read_id --> primer (all strings)
    """
    primer_info = {}
    for r in DictReader(open(csv_filename), delimiter=","):
        if r["primer"] != "NA":
            primer_info[r["id"]] = r["primer"]
    primer_ranges = list({map(int, list(primer_info.values()))})
    primer_ranges.sort()
    primer_ranges = list(map(str, primer_ranges))
    return primer_ranges, primer_info


# primer_names = ['B6-GV', 'B6-MII', 'B6-1C', 'B6-2C', 'B6-8C', 'B6-BI', 'DBA-GV', 'DBA-MII', 'DBA-1C', 'DBA-2C', 'DBA-8C', 'DBA-BI']


def get_fl_count_by_barcode(
    collapse_prefix: str,
    classify_csv: str,
    cluster_csv: str,
    primer_names: Optional[List[str]] = None,
) -> None:

    primer_ranges, primer_info = read_classify_csv(classify_csv)
    if primer_names is None:
        primer_names = {int(x): x for x in primer_ranges}

    cluster_info = defaultdict(lambda: [])
    for r in DictReader(open(cluster_csv), delimiter=","):
        cluster_info[r["cluster_id"]].append(r)

    group_filename = f"{collapse_prefix}.group.txt"
    logger.info(f"Reading {group_filename}...")

    fl = Path(f"{collapse_prefix}.fl_count_by_barcode.txt")
    with fl.open(mode="w") as f:
        f.write("pbid")
        for p in primer_ranges:
            f.write(f"\t{str(primer_names[int(p)])}")
        f.write("\n")
        for line in open(group_filename):
            pbid, members = line.strip().split("\t")
            tally = Counter()
            for m in members.split(","):
                cid = m.split("/")[0]
                for r in cluster_info[cid]:
                    if r["read_type"] == "FL":
                        tally[primer_info[r["read_id"]]] += 1
            f.write(pbid)
            for p in primer_ranges:
                f.write(f"\t{str(tally[p])}")
            f.write("\n")

    logger.info(f"Output written to: {fl.name}.")


def main(
    collapse_prefix: str = typer.Argument(
        ..., help="Collapse prefix (ex: hq_isoforms.fastq.no5merge.collapsed)"
    ),
    classify_csv: str = typer.Argument(
        ..., help="Classify output CSV (ex: classify.primer_info.csv)"
    ),
    cluster_csv: str = typer.Argument(
        ..., help="Cluster output CSV (ex: cluster_report.csv)"
    ),
) -> None:
    get_fl_count_by_barcode(collapse_prefix, classify_csv, cluster_csv)


if __name__ == "__main__":
    typer.run(main)
