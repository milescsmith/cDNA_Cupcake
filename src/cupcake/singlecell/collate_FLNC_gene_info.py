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

import os
import sys
from csv import DictReader, DictWriter

from Bio.Seq import Seq


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
            print(
                f"ignoring ID {pbid} cuz not in classification file.", file=sys.stderr
            )
            continue

        if is_clustered:
            # id: 1-ATCGAATGT-GCTTCTTTCACCTATCGATGATGGCTCAT-m64015_200531_015713/110297924/ccs
            _index, _umi, _bc, _ccs_id = ccs_id.split("-")
            ccs_id = _ccs_id

        if no_extra_base and (
            not is_clustered and umi_bc_info[ccs_id]["extra"] != "NA"
        ):
            print(f"ignoring ID {pbid} cuz extra bases.", file=sys.stderr)
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


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("group_filename", help="Collapse .group.txt")
    parser.add_argument("csv_filename", help="Trimmed UMI/BC CSV info")
    parser.add_argument("class_filename", help="SQANTI classification.txt")
    parser.add_argument("output_filename", help="Output filename")
    parser.add_argument(
        "-i", "--ontarget_filename", help="(Optional) on target information text"
    )
    parser.add_argument(
        "-p",
        "--dedup_ORF_prefix",
        help="(Optional) dedup-ed ORF group prefix, must have <pre>.faa and <pre>.group.txt",
    )
    parser.add_argument(
        "--no-extra-base",
        dest="no_extra_base",
        action="store_true",
        default=False,
        help="Drop all reads where there are extra bases",
    )
    parser.add_argument(
        "--is_clustered",
        action="store_true",
        default=False,
        help="group.txt contains post-UMI clustering result",
    )

    args = parser.parse_args()

    if os.path.exists(args.output_filename):
        print(
            f"Output file {args.output_filename} already exists. Abort!",
            file=sys.stderr,
        )
        sys.exit(-1)

    if not os.path.exists(args.group_filename):
        print(f"Group file {args.group_filename} not found. Abort!", file=sys.stderr)
        sys.exit(-1)

    if not os.path.exists(args.csv_filename):
        print(f"CSV file {args.csv_filename} not found. Abort!", file=sys.stderr)
        sys.exit(-1)

    if not os.path.exists(args.class_filename):
        print(f"Class file {args.class_filename} not found. Abort!", file=sys.stderr)
        sys.exit(-1)

    if args.ontarget_filename is not None and not os.path.exists(
        args.ontarget_filename
    ):
        print(
            f"Ontarget file {args.ontarget_filename} given but not found. Abort!",
            file=sys.stderr,
        )
        sys.exit(-1)

    if args.dedup_ORF_prefix is not None:
        if not os.path.exists(args.dedup_ORF_prefix + ".group.txt"):
            print(
                f"Dedup {args.dedup_ORF_prefix}.group.txt not found. Abort!",
                file=sys.stderr,
            )
            sys.exit(-1)
        if not os.path.exists(args.dedup_ORF_prefix + ".faa"):
            print(
                f"Dedup {args.dedup_ORF_prefix}.faa not found. Abort!", file=sys.stderr
            )
            sys.exit(-1)

    collate_gene_info(
        args.group_filename,
        args.csv_filename,
        args.class_filename,
        args.output_filename,
        args.ontarget_filename,
        args.dedup_ORF_prefix,
        args.no_extra_base,
        args.is_clustered,
    )


if __name__ == "__main__":
    main()