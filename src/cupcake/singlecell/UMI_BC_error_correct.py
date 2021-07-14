#!/usr/bin/env python
import sys
from collections import Counter, defaultdict
from csv import DictReader, DictWriter
from pathlib import Path
from typing import Dict, Optional, Union

import typer
from Bio.Seq import Seq

from cupcake import version_callback
from cupcake.logger import cupcake_logger as logger

app = typer.Typer(name="cupcake.singlecell.UMI_BC_error_correct")


def read_dropseq_clean_report(report_filename: str) -> Dict[Seq, Seq]:
    """
    EXAMPLE
    # FILTER_AMBIGUOUS=true
    # MIN_UMIS_PER_CELL=20
    # UMI_BIAS_THRESHOLD=20
    # EDIT_DISTANCE=20
    #
    # TOTAL_BARCODES_TESTED=996206
    # BARCODES_COLLAPSED=12546
    # ESTIMATED_UMIS_COLLAPSED=1107912
    # AMBIGUOUS_BARCODES=8206
    # ESTIMATED_AMBIGUOUS_UMIS=427356
    # POLY_T_BIASED_BARCODES=3374
    # POLY_T_BIASED_BARRCODE_UMIS=568661
    # POLY_T_POSITION=8
    intended_barcode        neighbor_barcode        intended_size   neighbor_size   position        intended_base   neighbor_base   repaired
    GGTTGGTCGCGG    GGTTGGTCTCGG    5458    20      9       G       T       FALSE
    TGGAGCTGGCGC    TGGCGCTGGCGC    3826    190     4       A       C       TRUE
    GGCGGCACGGCC    GACGGCACGGCC    16997   69      2       G       A       FALSE
    GAAGCCAGAGGT    GAAGCCAGAAGT    7037    47      10      G       A       FALSE
    GTGACGGACGGT    GTGACCGACGGT    8959    69      6       G       C       FALSE
    """
    with open(report_filename) as f:
        while True:
            cur_pos = f.tell()
            if not f.readline().startswith("#"):
                break
        f.seek(cur_pos)

        bc_repair_dict = {}
        reader = DictReader(f, delimiter="\t")
        for r in reader:
            if r["repaired"] == "TRUE":
                seq_from = Seq(r["neighbor_barcode"]).reverse_complement()
                seq_to = Seq(r["intended_barcode"]).reverse_complement()
                bc_repair_dict[seq_from] = seq_to
    return bc_repair_dict


def read_dropseq_synthesis_report(
    report_filename: Union[str, Path], bc_repair_dict: Optional[Dict[Seq, Seq]] = None
) -> Dict[Seq, Seq]:
    """
        EXAMPLE
        intended_sequence       related_sequences       num_related     deleted_base    deleted_base_pos        non_incorporated_rate   intended
    _UMIs   related_median_UMIs     intended_TBias  related_median_TBias
        NA      CGCAGCTCTGAG    1       NA      NA      NA      NA      20      NA      1
        NA      GCTGGACTTACA    1       NA      NA      NA      NA      22      NA      0.95
        CTCCACTGGAAA    CTCCACTGAAAA:CTCCACTGAAAC:CTCCACTGAAAG:CTCCACTGAAAT     4       G       9       0.47    4851    1039    0.3     0.98
        TTTAATATGGAT    TTAATATGGATG:TTAATATGGATT
    """
    with open(report_filename) as f:
        while True:
            cur_pos = f.tell()
            if not f.readline().startswith("#"):
                break
        f.seek(cur_pos)

        if bc_repair_dict is None:
            bc_repair_dict = {}
        reader = DictReader(f, delimiter="\t")
        for r in reader:
            if r["intended_sequence"] != "NA":
                seq_to = Seq(r["intended_sequence"]).reverse_complement()
                for s in r["related_sequences"].split(":"):
                    seq_from = Seq(s).reverse_complement()
                    bc_repair_dict[seq_from] = seq_to
    return bc_repair_dict


def edit_distance(seq1, seq2) -> int:
    assert len(seq1) == len(seq2)
    diff = 0
    for i in enumerate(seq1):
        diff += seq1[i] != seq2[i]
    return diff


def error_correct_BC_or_UMI(
    records: Dict[str, str], key: str, threshold: int = 1
) -> Dict[str, str]:
    """
    :param records: should be list of records all from the same gene!
    """
    assert key in ("BC", "UMI")
    merge_map = {}
    bc_count = Counter()
    for r in records:
        bc_count[r[key]] += 1

    # most common BC, in decreasing order
    bc_sorted = [bc for bc, _ in bc_count.most_common()]

    i = 0
    while i < len(bc_sorted) - 1:
        j = i + 1
        while j < len(bc_sorted):
            if edit_distance(bc_sorted[i], bc_sorted[j]) <= threshold:
                # merge j into i
                merge_map[bc_sorted[j]] = bc_sorted[i]
                bc_sorted.pop(j)
            else:
                j += 1
        i += 1

    # if len(merge_map) > 0:
    #    print merge_map
    #    raw_input()
    return merge_map


def umi_bc_error_correct(
    csv_filename: Union[str, Path],
    output_filename: Union[str, Path],
    shortread_bc: Optional[Dict[str, str]] = None,
    only_top_ranked: bool = False,
    bc_repair_dict: Optional[Dict[Seq, Seq]] = None,
) -> None:

    shortread_bc = {} if shortread_bc is None else shortread_bc
    reader = DictReader(open(csv_filename), delimiter="\t")

    FIELDS = reader.fieldnames + ["BC_ed", "UMI_ed", "BC_match", "BC_top_rank"]
    f = open(output_filename, "w")
    writer = DictWriter(f, FIELDS, delimiter="\t")
    writer.writeheader()

    recs_by_gene = defaultdict(lambda: [])
    for r in reader:
        recs_by_gene[r["gene"]].append(r)

    # error correct BCs by gene group
    for gene in recs_by_gene:
        recs_by_bc = defaultdict(lambda: [])

        if bc_repair_dict is not None:  # has DropSeq BC cleaning report! Just use it!
            for r in recs_by_gene[gene]:
                if r["BC"] in bc_repair_dict:
                    r["BC_ed"] = bc_repair_dict[r["BC"]]
                else:
                    r["BC_ed"] = r["BC"]
                recs_by_bc[r["BC"]].append(r)
        else:
            bc_merge_map = error_correct_BC_or_UMI(recs_by_gene[gene], "BC")
            for r in recs_by_gene[gene]:
                if r["BC"] in bc_merge_map:
                    r["BC_ed"] = bc_merge_map[r["BC"]]
                else:
                    r["BC_ed"] = r["BC"]
                recs_by_bc[r["BC"]].append(r)
        # now error correct by UMI
        for bc in recs_by_bc:
            umi_merge_map = error_correct_BC_or_UMI(recs_by_bc[bc], "UMI")
            for r in recs_by_bc[bc]:
                if r["UMI"] in umi_merge_map:
                    r["UMI_ed"] = umi_merge_map[r["UMI"]]
                else:
                    r["UMI_ed"] = r["UMI"]

                BC_ed_rev = str(Seq(r["BC_ed"]).reverse_complement())

                r["BC_match"] = "Y" if BC_ed_rev in shortread_bc else "N"
                r["BC_top_rank"] = (
                    "Y"
                    if (r["BC_match"] == "Y" and shortread_bc[BC_ed_rev] == "Y")
                    else "N"
                )

                if not only_top_ranked or r["BC_top_rank"] == "Y":
                    writer.writerow(r)


@app.command(name="")
def main(
    input_csv: str = typer.Argument(..., help="Input CSV"),
    output_csv: str = typer.Argument(..., help="Output CSV"),
    bc_rank_file: str = typer.Argument(
        ..., help="Cell barcode rank file from short read data"
    ),
    only_top_ranked: bool = typer.Option(
        False,
        help="Only output those that are top-ranked. Must have --bc_rank_file.",
    ),
    dropseq_clean_report: str = typer.Option(
        ...,
        help="Output from running DetectBeadSubstitutionErrors in DropSeq cookbook (ex: star_gene_exon_tagged_clean_substitution.bam_report.txt)",
    ),
    dropseq_synthesis_report: str = typer.Option(
        ...,
        help="Output from running DetectBeadSynthesisErrors in DropSeq cookbook (ex: star_gene_exon_tagged_clean_substitution_clean2.bam_report.txt)",
    ),
    version: bool = typer.Option(
        None,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Prints the version of the SQANTI3 package.",
    ),
) -> None:
    shortread_bc = {}  # dict of cell barcode -> "Y" for top ranked
    if bc_rank_file is not None:
        reader = DictReader(open(bc_rank_file), delimiter="\t")
        for r in reader:
            shortread_bc[r["cell_barcode"]] = r["top_ranked"]
    else:
        if only_top_ranked:
            logger.error("--bc_rank_file must be given if using --only_top_ranked!")
            sys.exit(-1)

    bc_repair_dict = None
    if dropseq_clean_report is not None:
        bc_repair_dict = read_dropseq_clean_report(dropseq_clean_report)
    if dropseq_synthesis_report is not None:
        bc_repair_dict = read_dropseq_synthesis_report(
            dropseq_synthesis_report, bc_repair_dict
        )

    umi_bc_error_correct(
        input_csv,
        output_csv,
        shortread_bc,
        only_top_ranked,
        bc_repair_dict,
    )


if __name__ == "__main__":
    typer.run(main)
