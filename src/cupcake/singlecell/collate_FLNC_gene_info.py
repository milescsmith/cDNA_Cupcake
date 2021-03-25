#!/usr/bin/env python
"""
Given:

   1. single cell UMI, BC file (ccs -> umi, bc)
   2. collapse group.txt  (ccs -> pbid)
   3. SQANTI classification (pbid -> transcript, isoform, category)
   4. optional ontarget file (pbid -> ontarget or not)

Output a collated infor file that is:

   <ccs>, <pbid>, <transcript>, <gene>, <category>, <ontarget Y|N|NA>, <UMI>, <BC>, <UMIrev>, <BCrev>
"""

from csv import DictReader, DictWriter
from pathlib import Path
from typing import Optional

import typer
from Bio.Seq import Seq
from cupcake.logging import cupcake_logger as logger

app = typer.Typer(name="cupcake.singlecell.collate_FLNC_gene_info")


def read_group_info(group_filename):
    """
    :return: dict of ccs -> pbid
    """
    d = {}
    for line in open(group_filename):
        pbid, members = line.strip().split("\t")
        for m in members.split(","):
            d[m] = pbid
    return d


def collate_gene_info(
    group_filename,
    csv_filename,
    class_filename,
    output_filename,
    ontarget_filename=None,
    dedup_ORF_prefix=None,
    no_extra_base=False,
    is_clustered=False,
):
    """
    <id>, <pbid>, <length>, <transcript>, <gene>, <category>, <ontarget Y|N|NA>, <ORFgroup NA|NoORF|groupID>, <UMI>, <BC>
    """
    FIELDS = [
        "id",
        "pbid",
        "length",
        "transcript",
        "gene",
        "category",
        "ontarget",
        "ORFgroup",
        "UMI",
        "UMIrev",
        "BC",
        "BCrev",
    ]

    group_info = read_group_info(group_filename)
    umi_bc_info = {r["id"]: r for r in DictReader(open(csv_filename), delimiter="\t")}
    sqanti_info = {
        r["isoform"]: r for r in DictReader(open(class_filename), delimiter="\t")
    }
    if ontarget_filename is not None:
        ontarget_info = {
            r["read_id"]: r for r in DictReader(open(ontarget_filename), delimiter="\t")
        }

    if dedup_ORF_prefix is not None:
        dedup_ORF_info = (
            {}
        )  # seqid --> which group they belong to (ex: PB.1.2 --> ORFgroup_PB.1_1)
        for line in open(f"{dedup_ORF_prefix}.group.txt"):
            group_id, members = line.strip().split("\t")
            for pbid in members.split(","):
                dedup_ORF_info[pbid] = group_id

    f = open(output_filename, "w")
    writer = DictWriter(f, FIELDS, delimiter="\t")
    writer.writeheader()

    for ccs_id, pbid in group_info.items():
        if pbid not in sqanti_info:
            logger.error(f"ignoring ID {pbid} cuz not in classification file.")
            continue

        if is_clustered:
            # id: 1-ATCGAATGT-GCTTCTTTCACCTATCGATGATGGCTCAT-m64015_200531_015713/110297924/ccs
            _index, _umi, _bc, _ccs_id = ccs_id.split("-")
            ccs_id = _ccs_id

        if no_extra_base and (
            not is_clustered and umi_bc_info[ccs_id]["extra"] != "NA"
        ):
            logger.info(f"ignoring ID {pbid} cuz extra bases.")
            continue
        rec = {"id": ccs_id, "pbid": pbid}
        rec["length"] = sqanti_info[pbid]["length"]
        rec["category"] = sqanti_info[pbid]["structural_category"]
        rec["transcript"] = sqanti_info[pbid]["associated_transcript"]
        rec["gene"] = sqanti_info[pbid]["associated_gene"]

        if is_clustered:
            rec["UMI"] = _umi
            rec["BC"] = _bc
        else:
            rec["UMI"] = umi_bc_info[ccs_id]["UMI"]
            rec["BC"] = umi_bc_info[ccs_id]["BC"]
        rec["UMIrev"] = Seq(rec["UMI"]).reverse_complement()
        rec["BCrev"] = Seq(rec["BC"]).reverse_complement()
        if ontarget_filename is None:
            rec["ontarget"] = "NA"
        else:
            rec["ontarget"] = "Y" if ontarget_info[pbid]["genes"] != "" else "N"
        if dedup_ORF_prefix is None:
            rec["ORFgroup"] = "NA"
        else:
            if pbid not in dedup_ORF_info:
                rec["ORFgroup"] = "NoORF"
            else:
                rec["ORFgroup"] = dedup_ORF_info[pbid]

        writer.writerow(rec)

    f.close()


@app.command(name="")
def main(
    group_filename: str = typer.Argument(..., help="Collapse .group.txt"),
    csv_filename: str = typer.Argument(..., help="Trimmed UMI/BC CSV info"),
    class_filename: str = typer.Argument(..., help="SQANTI classification.txt"),
    output_filename: str = typer.Argument(..., help="Output filename"),
    ontarget_filename: str = typer.Option(
        ..., "-i", help="(Optional) on target information text"
    ),
    dedup_ORF_prefix: Optional[str] = typer.Option(
        None,
        "-p",
        help="Dedup-ed ORF group prefix, must have <pre>.faa and <pre>.group.txt",
    ),
    no_extra_base: bool = typer.Option(
        False, help="Drop all reads where there are extra bases"
    ),
    is_clustered: bool = typer.Option(
        False, help="group.txt contains post-UMI clustering result"
    ),
):
    if Path(output_filename).exists():
        raise FileExistsError(f"Output file {output_filename} already exists. Abort!")

    if not Path(group_filename).exists():
        raise FileNotFoundError(f"Group file {group_filename} not found. Abort!")

    if not Path(csv_filename).exists():
        raise FileNotFoundError(f"CSV file {csv_filename} not found. Abort!")

    if not Path(class_filename).exists():
        raise FileNotFoundError(f"Class file {class_filename} not found. Abort!")

    if ontarget_filename is not None and not Path(ontarget_filename).exists():
        raise FileNotFoundError(
            f"Ontarget file {ontarget_filename} given but not found. Abort!"
        )

    if dedup_ORF_prefix is not None:
        if not Path(f"{dedup_ORF_prefix}.group.txt").exists():
            raise FileNotFoundError(
                f"Dedup {dedup_ORF_prefix}.group.txt not found. Abort!"
            )
        if not Path(f"{dedup_ORF_prefix}.faa").exists():
            raise FileNotFoundError(f"Dedup {dedup_ORF_prefix}.faa not found. Abort!")

    collate_gene_info(
        group_filename,
        csv_filename,
        class_filename,
        output_filename,
        ontarget_filename,
        dedup_ORF_prefix,
        no_extra_base,
        is_clustered,
    )


if __name__ == "__main__":
    typer.run(main)
