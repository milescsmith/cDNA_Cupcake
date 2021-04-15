#!/usr/bin/env python

__author__ = "etseng@pacb.com"

"""
Filter away 5' degraded products (isoforms that are shorter isoforms of longer ones)
3' differences are always honored.

If input is:

Isoform A: exon 1, 2, 3, 4
Isoform B: exon 3, 4
Isoform C: exon 3, 4, 5

Then Isoform A and C are preserved. Isoform B (degraded form of A) is filtered out.

Required input: <input_prefix>.gff, .rep.fq, .abundance.txt

Example:
    filter_away_subset.py test.collapsed test.collapsed.filtered

"""

import sys
from collections import defaultdict
from csv import DictReader, DictWriter
from pathlib import Path
from typing import Dict, Tuple

import typer
from Bio import SeqIO

from cupcake.logging import cupcake_logger as logger
from cupcake.sequence import GFF
from cupcake.tofu import compare_junctions

app = typer.Typer(name="filer_away_subset")


def sanity_check_collapse_input(input_prefix: str) -> Tuple[Path, Path, Path, str]:
    """
    Check that
    1. the count, gff, rep files exist
    2. the number of records agree among the three
    """
    # group_filename = f"{input_prefix}.group.txt"
    count_filename = Path(f"{input_prefix}.abundance.txt")
    gff_filename = Path(f"{input_prefix}.gff")
    rep_filenames = [
        (f"{input_prefix}.rep.fq", "fastq"),
        (f"{input_prefix}.rep.fastq", "fastq"),
        (f"{input_prefix}.rep.fa", "fasta"),
        (f"{input_prefix}.rep.fasta", "fasta"),
    ]

    rep_filename = None
    rep_type = None
    for x, filetype in rep_filenames:
        if Path(x).exists():
            rep_filename = Path(x)
            rep_type = filetype

    if rep_filename is None:
        raise RuntimeError(
            f"Expected to find input fasta or fastq files {input_prefix}.rep.fa or {input_prefix}.rep.fq. Not found. Abort!"
        )

    if not count_filename.exists():
        raise RuntimeError(f"File {count_filename} does not exist. Abort!")
    if not gff_filename.exists():
        raise RuntimeError(f"File {gff_filename} does not exist. Abort!")

    pbids1 = set([r.id for r in SeqIO.parse(open(rep_filename), rep_type)])
    pbids2 = set([r.seqid for r in GFF.collapseGFFReader(gff_filename)])
    pbids3 = set(read_count_file(count_filename)[0].keys())

    if (
        len(pbids1) != len(pbids2)
        or len(pbids2) != len(pbids3)
        or len(pbids1) != len(pbids3)
    ):
        logger.error(
            "The number of PBID records in the files disagree! Sanity check failed."
        )
        logger.error(f"# of PBIDs in {rep_filename}: {len(pbids1)}")
        logger.error(f"# of PBIDs in {gff_filename}: {len(pbids2)}")
        logger.error(f"# of PBIDs in {count_filename}: {len(pbids3)}")
        sys.exit(-1)

    return count_filename, gff_filename, rep_filename, rep_type


def read_count_file(count_filename: Path) -> Tuple[Dict[str, str], str]:
    f = count_filename.open()
    count_header = ""
    while True:
        cur_pos = f.tell()
        line = f.readline()
        if not line.startswith("#"):
            f.seek(cur_pos)
            break
        else:
            count_header += line
    d = {r["pbid"]: r for r in DictReader(f, delimiter="\t")}
    f.close()
    return d, count_header


def can_merge(
    m: str, r1: GFF.gmapRecord, r2: GFF.gmapRecord, internal_fuzzy_max_dist: int
) -> bool:
    if m == "subset":
        r1, r2 = r2, r1  # rotate so r1 is always the longer one
    if m == "super" or m == "subset":
        n2 = len(r2.ref_exons)
        if r1.strand == "+":
            # if r2 is monoexonic, it can start after the last exon of r1's last exon
            # if r2 is multiexonic, the last start must be pretty close (fuzzy allowed)
            if n2 == 1:  # r2 is mono-exonic
                return (
                    r1.ref_exons[-1].start - r2.ref_exons[-1].start
                    <= internal_fuzzy_max_dist
                )
            else:
                return (
                    abs(r1.ref_exons[-1].start - r2.ref_exons[-1].start)
                    <= internal_fuzzy_max_dist
                    and r1.ref_exons[-n2].start
                    <= r2.ref_exons[0].start
                    < r1.ref_exons[-n2].end
                )
        else:
            if n2 == 1:
                return (
                    r1.ref_exons[0].end - r2.ref_exons[0].end
                    >= -internal_fuzzy_max_dist
                )
            else:
                return (
                    abs(r1.ref_exons[0].end - r2.ref_exons[0].end)
                    <= internal_fuzzy_max_dist
                    and r1.ref_exons[n2 - 1].start
                    <= r2.ref_exons[-1].end
                    < r1.ref_exons[n2].end
                )


def filter_out_subsets(
    recs: Dict[int, GFF.gmapRecord], internal_fuzzy_max_dist: int
) -> None:
    # recs must be sorted by start becuz that's the order they are written
    i = 0
    while i < len(recs) - 1:
        no_change = True
        j = i + 1
        while j < len(recs):
            if recs[j].start > recs[i].end:
                break
            recs[i].segments = recs[i].ref_exons
            recs[j].segments = recs[j].ref_exons
            m = compare_junctions.compare_junctions(
                recs[i], recs[j], internal_fuzzy_max_dist
            )
            if can_merge(m, recs[i], recs[j], internal_fuzzy_max_dist):
                if m == "super":  # pop recs[j]
                    recs.pop(j)
                else:
                    recs.pop(i)
                    no_change = False
            else:
                j += 1
        if no_change:
            i += 1


@app.command(name="")
def main(
    input_prefix: str = typer.Argument(
        ..., help="Input prefix (ex: test.collapsed.min_fl_2)"
    ),
    fuzzy_junction: int = typer.Option(
        5, help="Fuzzy junction max dist (default: 5bp)"
    ),
) -> None:

    output_prefix = f"{input_prefix}.filtered"

    count_filename, gff_filename, rep_filename, rep_type = sanity_check_collapse_input(
        input_prefix
    )

    recs = defaultdict(lambda: [])
    reader = GFF.collapseGFFReader(gff_filename)
    for r in reader:
        assert r.seqid.startswith("PB.")
        recs[int(r.seqid.split(".")[1])].append(r)

    good = []
    with open(f"{output_prefix}.gff", "w") as f:
        keys = list(recs.keys())
        keys.sort()
        for k in recs:
            xxx = recs[k]
            filter_out_subsets(xxx, fuzzy_junction)
            for r in xxx:
                GFF.write_collapseGFF_format(f, r)
                good.append(r.seqid)

    # read abundance first
    d, count_header = read_count_file(count_filename)

    # write output rep.fq
    with open(
        f"{output_prefix}.rep.{('fq' if rep_type == 'fastq' else 'a')}", "w"
    ) as f:
        for r in SeqIO.parse(open(rep_filename), rep_type):
            if r.name.split("|")[0] in good:
                SeqIO.write(r, f, rep_type)

    # write output to .abundance.txt
    with open(f"{output_prefix}.abundance.txt", "w") as f:
        f.write(count_header)
        writer = DictWriter(
            f,
            fieldnames=[
                "pbid",
                "count_fl",
                "count_nfl",
                "count_nfl_amb",
                "norm_fl",
                "norm_nfl",
                "norm_nfl_amb",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for k in good:
            r = d[k]
            writer.writerow(r)

    logger.info(
        f"Output written to: {output_prefix}.gff\n"
        f"Output written to: {rep_filename}\n"
        f"Output written to: {output_prefix}.abundance.txt"
    )


if __name__ == "__main__":
    typer.run(main)
