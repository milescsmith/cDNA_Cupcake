#!/usr/bin/env python
__author__ = "etseng@pacb.com"

"""
After running UMI_BC_error_correct.py with --bc_rank_file and --only_top_ranked parameters,
 take the (1) .annotated.correct.csv file, and
          (2) short read cluster (cell type) csv file
          (3) fasta, gff, [faa]

 to de-duplicate the FLNC reads by cluster (cell type).

Outputs:
(1) a de-dup table of

(UMI-ed, BC-ed, gene) — pbid — associate gene — associated transcript — category — length — cluster #

(2)  A “master” file, one sequence for each [pbid] that appeared at least once in (1)
- fasta
- gff

(3) A “per cluster” file, one sequence for each [pbid] that appeared once in each cluster
- fasta
- gff
"""

from collections import Counter, defaultdict
from csv import DictReader, DictWriter
from pathlib import Path

import typer
from Bio import SeqIO

from cupcake import version_callback
from cupcake.sequence.GFF import collapseGFFReader, write_collapseGFF_format

app = typer.Typer(name="cupcake.singlecell.dedup_FLNC_per_cluster")


CORRECTED_CSV_FILELDS = [
    "id",
    "pbid",
    "length",
    "transcript",
    "gene",
    "category",
    "ORFgroup",
    "UMI_ed",
    "BC_ed",
]
# PBID_FORMAT = re.compile(r"PB.(\d+).(\d+)")


def dedup_FLNC_per_cluster(
    corrected_csv,
    cluster_info,
    output_prefix,
    fasta_file=None,
    gff_file=None,
    faa_file=None,
):

    # read corrected CSV
    reader = DictReader(open(corrected_csv), delimiter="\t")
    for k in CORRECTED_CSV_FILELDS:
        if k not in reader.fieldnames:
            raise RuntimeError(
                "The following fields must exist in {}!\n{}".format(
                    corrected_csv, "\n".join(CORRECTED_CSV_FILELDS)
                )
            )

    per_unique = {}  # tag -> record
    per_unique_count = Counter()  # tag -> number of duplicates
    per_pbid = defaultdict(
        lambda: {"gene": None, "transcript": None, "clusters": []}
    )  # pbid --> list of clusters it is in
    for r in reader:
        tag = f"{r['BC_ed']}-{r['UMI_ed']}-{r['gene']}"
        per_unique[tag] = r
        per_unique_count[tag] += 1

    # now link barcode to cell type, also PCR dup counts
    for tag in per_unique:
        c = cluster_info[per_unique[tag]["BC_ed"]]
        rec = per_unique[tag]
        rec["cluster"] = c
        rec["num_dups"] = per_unique_count[tag]
        pbid = rec["pbid"]
        if pbid in per_pbid:
            per_pbid[pbid]["clusters"].add(c)
        else:
            per_pbid[pbid] = {
                "gene": rec["gene"],
                "transcript": rec["transcript"],
                "clusters": {c},
            }

    # write out de-dup CSV file
    with open(f"{output_prefix}.csv", "w") as f:
        writer = DictWriter(
            f,
            CORRECTED_CSV_FILELDS + ["cluster", "num_dups"],
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        keys = per_unique.keys()
        for k in sorted(keys):
            writer.writerow(per_unique[k])

    if fasta_file is not None:
        f_d = {}  # cluster --> file handle
        # writer pbid master file
        with open(f"{output_prefix}.fasta", "w") as f:
            for r in SeqIO.parse(open(fasta_file), "fasta"):
                if r.id in per_pbid:
                    newid = f"{r.id}|{per_pbid[r.id]['gene']}|{per_pbid[r.id]['transcript']}|{';'.join(per_pbid[r.id]['clusters'])}"
                    f.write(f">{newid}\n{r.seq}\n")
                    for c in per_pbid[r.id]["clusters"]:
                        if c not in f_d:
                            f_d[c] = open(f"{output_prefix}.{c}.fasta", "w")
                        f_d[c].write(f">{newid}\n{r.seq}\n")

    if faa_file is not None:
        f_d = {}  # cluster --> file handle
        # writer pbid master file
        with open(f"{output_prefix}.faa", "w") as f:
            for r in SeqIO.parse(open(faa_file), "fasta"):
                if r.id in per_pbid:
                    newid = f'{r.id}|{per_pbid[r.id]["gene"]}|{per_pbid[r.id]["transcript"]}|{";".join(per_pbid[r.id]["clusters"])}'
                    f.write(f">{newid}\n{r.seq}\n")
                    for c in per_pbid[r.id]["clusters"]:
                        if c not in f_d:
                            f_d[c] = open(f"{output_prefix}.{c}.faa", "w")
                        f_d[c].write(f">{newid}\n{r.seq}\n")
        for handle in f_d.values():
            handle.close()

    if gff_file is not None:
        f_d = {}  # cluster --> file handle
        # writer pbid master file
        with open(f"{output_prefix}.gff", "w") as f:
            for r in collapseGFFReader(gff_file):
                if r.seqid in per_pbid:
                    newid = f'{r.seqid}|{per_pbid[r.seqid]["gene"]}|{per_pbid[r.seqid]["transcript"]}|{";".join(per_pbid[r.seqid]["clusters"])}'
                    write_collapseGFF_format(f, r)
                    for c in per_pbid[r.seqid]["clusters"]:
                        if c not in f_d:
                            f_d[c] = open(f"{output_prefix}.{c}.gff", "w")
                        write_collapseGFF_format(f_d[c], r)
        for handle in f_d.values():
            handle.close()


@app.command(name="")
def main(
    corrected_csv: str = typer.Argument(
        ..., help="Annotated, error-corrected FLNC CSV file"
    ),
    cluster_file: str = typer.Argument(
        ..., help="Short read barcode to cluster CSV file"
    ),
    fasta: str = typer.Argument(
        ..., help="(Optional) Fasta file (IDs should be PB.X.Y)"
    ),
    gff: str = typer.Argument(..., help="(Optional) GFF file (IDs should be PB.X.Y)"),
    faa: str = typer.Argument(..., help="(Optional) Faa file (IDs should be PB.X.Y)"),
    version: bool = typer.Option(
        None,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Prints the version of the SQANTI3 package.",
    ),
):

    if not Path(corrected_csv).exists():
        raise FileNotFoundError(f"Input file {corrected_csv} does not exist! Abort!")

    if not Path(cluster_file).exists():
        raise FileNotFoundError(f"Input file {cluster_file} does not exist! Abort!")

    cluster_info = {}
    reader = DictReader(open(cluster_file), delimiter="\t")
    if ("cell_barcode" not in reader.fieldnames) or (
        "cluster" not in reader.fieldnames
    ):
        raise RuntimeError(
            f"Cluster file {cluster_file} must contain 'cell_barcode' and 'cluster' fields!"
        )
    for r in reader:
        if r["cluster"] != "NA":
            cluster_info[r["cell_barcode_rev"]] = r["cluster"]

    output_prefix = f"{corrected_csv[:corrected_csv.rfind('.')]}.dedup"

    dedup_FLNC_per_cluster(corrected_csv, cluster_info, output_prefix, fasta, gff, faa)


if __name__ == "__main__":
    typer.run(main)
