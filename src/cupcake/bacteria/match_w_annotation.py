#!/usr/bin/env python

"""match_w_annotation.py: functions for categorizing read alignments
to annotations. mainly for bacteria.

Criteria for matching genes:

Single --- query must cover the ref gene by 50% and does not extend upstream
more than 300 bp or downstream more than 500 bp

Poly --- query must cover each of the ref gene by 50%

Novel ---
 (a) novel-antisense: has overlapping genes on opposite strand but not the same
    strand
 (b) novel-unannotated: no overlapping genes on either strand
 (c) novel-partial: has overlapping genes on same strand (but less than the
    "single" criterion)
"""

__copyright__ = "Copyright 2016, cDNA_Cupcake"

from collections import defaultdict, namedtuple
from csv import DictReader
from typing import List, Tuple

import typer
from Bio import SeqIO
from bx.intervals import IntervalTree
from bx.intervals.cluster import ClusterTree

from cupcake.__about__ import __author__, __email__, __version__
from cupcake.logging import cupcake_logger as logger
from cupcake.sequence import BioReaders

app = typer.Typer(
    name="cupcake.bacteria.match_w_annotation",
    help="Match alignment with annotation. Categorize and Report.",
)


AMatch = namedtuple("AMatch", "name strand start end record")


def calc_overlap(s1: int, e1: int, s2: int, e2: int) -> int:
    return min(e1, e2) - max(s1, s2)


def calc_overlap_ratio1(s1: int, e1: int, s2: int, e2: int) -> float:
    return (min(e1, e2) - max(s1, s2)) * 1.0 / (e1 - s1)


def check_multigene(
    overlaps: List[Tuple[str, int, float, float]],
    min_overlap_bp: int = 0,
    min_query_overlap: int = 0,
    min_gene_overlap: float = 0.5,
) -> str:
    """
    overlaps is a list of: (gene, overlap_bp, overlap_gene_ratio, overlap_query_ratio)
    """
    if all(
        x[1] >= min_overlap_bp
        and x[2] >= min_gene_overlap
        and x[3] >= min_query_overlap
        for x in overlaps
    ):
        new_name = f"poly-{'-'.join(x[0] for x in overlaps)}"
        return new_name
    #    elif overlaps[0][2] >= min_gene_overlap: # first gene covers 50% of query
    #        return overlaps[0][0]
    #    elif overlaps[-1][2] >= min_gene_overlap:
    #        return overlaps[-1][0]
    else:
        return "novel"


def categorize_novel(t, r, info, same_strand_overlap_gene=None) -> namedtuple:
    """
    Further categorize novel as:

    (a) anti-sense: opposite strand of any gene
    (b) unannotated: does not overlap with any gene on any strand
    (c) other: none of the above
    """
    s, e = r.sStart, r.sEnd
    matches = t[r.sID]["+" if r.flag.strand == "-" else "-"].find(s, e)
    if len(matches) > 0 and same_strand_overlap_gene is None:
        *_, gene2 = info[matches[0]]
        return AMatch(f"novel-antisense-{gene2}", r.flag.strand, s, e, r)
    else:
        if same_strand_overlap_gene is None:
            return AMatch("novel-unannotated", r.flag.strand, s, e, r)
        else:
            return AMatch(
                f"novel-partial-{same_strand_overlap_gene}", r.flag.strand, s, e, r
            )


def match_w_annotation(
    t,
    r,
    info,
    min_overlap_bp=0,
    min_query_overlap=0,
    min_gene_overlap=0.5,
    max_upstream_bp=300,
    max_downstream_bp=500,
):
    """
    Input:
       t -- dict of chr -> strand -> bx.intervals.IntervalTree
       r -- GMAP SAM record
       info --- dict of (start, end, gene_name)

    Return: AMatch(matching_name, strand, start, end, record)
    matching_name --- "novel" if it matches 0 genes;
        gene_name if it matches that 1 gene;
        gene1_gene2_... if it matches multiple genes
    """
    s, e = r.sStart, r.sEnd
    matches = t[r.sID][r.flag.strand].find(s, e)
    if len(matches) == 0:
        return categorize_novel(t, r, info, same_strand_overlap_gene=None)
    elif len(matches) == 1:
        s2, e2, gene2 = info[matches[0]]
        # check that the query overlaps >= 50% (or min_gene_overlap) of the annotation
        if calc_overlap_ratio1(s2, e2, s, e) >= min_gene_overlap and (
            (
                r.flag.strand == "+"
                and s2 - s <= max_upstream_bp
                and e - e2 <= max_downstream_bp
            )
            or (
                r.flag.strand == "-"
                and s2 - s <= max_downstream_bp
                and e - e2 <= max_upstream_bp
            )
        ):
            return AMatch(gene2, r.flag.strand, s, e, r)
        else:
            return categorize_novel(t, r, info, same_strand_overlap_gene=gene2)
    else:  # matches 2+ genes
        # pdb.set_trace()
        overlaps = (
            []
        )  # list of (gene, overlap_bp, overlap_gene_ratio, overlap_query_ratio)
        for gene2 in matches:
            s2, e2, gene2 = info[gene2]
            o = calc_overlap(s2, e2, s, e)
            overlaps.append((gene2, o, o * 1.0 / (e2 - s2), o * 1.0 / (e - s)))

        result = check_multigene_helper(
            r, s, e, overlaps, min_overlap_bp, min_query_overlap, min_gene_overlap
        )
        if result.name.startswith("poly"):
            return result
        else:
            # check if any of the single gene criterion matches
            for gene2 in matches:
                s2, e2, gene2 = info[gene2]
                if calc_overlap_ratio1(s2, e2, s, e) >= min_gene_overlap and (
                    (
                        r.flag.strand == "+"
                        and s2 - s <= max_upstream_bp
                        and e - e2 <= max_downstream_bp
                    )
                    or (
                        r.flag.strand == "-"
                        and s2 - s <= max_downstream_bp
                        and e - e2 <= max_upstream_bp
                    )
                ):
                    return AMatch(gene2, r.flag.strand, s, e, r)
            # if we reach here, we are NOVEL
            return categorize_novel(
                t, r, info, same_strand_overlap_gene=info[matches[0]][2]
            )


def check_multigene_helper(
    r, s, e, overlaps, min_overlap_bp, min_query_overlap, min_gene_overlap
):
    flag = check_multigene(
        overlaps, min_overlap_bp, min_query_overlap, min_gene_overlap
    )
    if flag.startswith("poly"):
        return AMatch(flag, r.flag.strand, s, e, r)
    else:
        # try shaving off genes at ends
        for i in range(1, len(overlaps) - 1):
            result = check_multigene_helper(
                r,
                s,
                e,
                overlaps[i:],
                min_overlap_bp,
                min_query_overlap,
                min_gene_overlap,
            )
            if result.name.startswith("poly"):
                return result
        for i in range(len(overlaps) - 1, 1, -1):
            result = check_multigene_helper(
                r,
                s,
                e,
                overlaps[:i],
                min_overlap_bp,
                min_query_overlap,
                min_gene_overlap,
            )
            if result.name.startswith("poly"):
                return result
        return AMatch("novel", r.flag.strand, s, e, r)


def categorize_aln_by_annotation(
    gene_annotation_file: str,
    input_fasta: str,
    input_sam: str,
    output_prefix: str,
    min_overlap_bp: int = 200,
    min_query_overlap: float = 0.5,
    min_gene_overlap: float = 0.8,
) -> None:

    t = defaultdict(
        lambda: {"+": IntervalTree(), "-": IntervalTree()}
    )  # chr -> strand -> IntervalTree
    info = {}

    # reader = DictReader(open('ProteinTable149_154224.txt'),delimiter='\t')
    for r in DictReader(open(gene_annotation_file), delimiter="\t"):
        if r["#Replicon Name"] != "chr":
            logger.info(f"Ignore {r}")
            continue
        info[r["Locus tag"]] = (int(r["Start"]), int(r["Stop"]), r["Locus tag"])
        t[r["Replicon Accession"]][r["Strand"]].add(
            int(r["Start"]), int(r["Stop"]), r["Locus tag"]
        )

    # pdb.set_trace()

    result = defaultdict(lambda: [])  # gene -> list of rec
    d = {r.id: len(r.seq) for r in SeqIO.parse(open(input_fasta), "fasta")}

    reader = BioReaders.GMAPSAMReader(input_sam, True, query_len_dict=d)
    for r in reader:
        # if r.qID == 'm151125_055539_42275_c100921822550000001823204305121656_s1_p0/121461/30_2108_CCS':
        #    pdb.set_trace()
        ans = match_w_annotation(
            t, r, info, min_overlap_bp, min_query_overlap, min_gene_overlap
        )
        # ans is AMatch(name, strand, start, end, record)
        result[ans.name].append(ans)

    novel_ct = defaultdict(lambda: {"+": ClusterTree(0, 0), "-": ClusterTree(0, 0)})
    novel_list = []
    novel_index = 0

    with open(f"{output_prefix}.sam", "w") as f, open(
        f"{output_prefix}.report.txt", "w"
    ) as f1:
        f.write(reader.header)
        f1.write("id\tread_group\tgene_name\tserial_number\tstrand\tstart\tend\n")
        for k, v in result.items():
            # v is: list of AMatch(name, strand, start, end, record)
            if k.startswith("novel-unannotated"):
                # write novel later, we are grouping them by loci first
                # tagRG='novel'
                for x in v:
                    novel_ct[x.record.sID][x.strand].insert(x.start, x.end, novel_index)
                    novel_index += 1
                    novel_list.append(x)
                continue
            elif k.startswith("novel-antisense"):
                tagRG = "novel-antisense"
            elif k.startswith("novel-partial"):
                tagRG = "novel-partial"
            elif k.startswith("poly-"):
                tagRG = "poly"
            else:
                tagRG = "single"
            v.sort(
                key=lambda x: (x.start, x.end),
                reverse=bool(v[0].strand == "-"),
            )  # sort by start, then end
            for i, x in enumerate(v):
                f.write(
                    f"{x.record.record_line}\tSN:Z:{i + 1:06d}\tRG:Z:{tagRG}\tgn:Z:{k}\n"
                )
                if x.strand == "+":
                    f1.write(
                        f"{x.record.qID}\t{tagRG}\t{k}\t{i + 1:06d}\t{x.strand}\t{x.start + 1}\t{x.end}\n"
                    )
                else:  # - strand, start is end, end is start
                    f1.write(
                        f"{x.record.qID}\t{tagRG}\t{k}\t{i + 1:06d}\t{x.strand}\t{x.end}\t{x.start + 1}\n"
                    )

        # now write the novel stuff, grouped by regions
        novel_region_index = 1
        for d1 in novel_ct.values():
            for ct in d1.values():
                gn = f"novel-{str(novel_region_index)}"
                for *_, _indices in ct.getregions():
                    v = [novel_list[ind] for ind in _indices]
                    v.sort(
                        key=lambda x: (x.start, x.end),
                        reverse=bool(v[0].strand == "-"),
                    )  # sort by start, then end
                    for i, x in enumerate(v):
                        f.write(
                            f"{x.record.record_line}\tSN:Z:{i + 1:06d}\tRG:Z:{'novel-unannotated'}\tgn:Z:{gn}\n"
                        )
                        if x.strand == "+":
                            f1.write(
                                f"{x.record.qID}\t{'novel-unannotated'}\t{gn}\t{i + 1:06d}\t{x.strand}\t{x.start + 1}\t{x.end}\n"
                            )
                        else:
                            f1.write(
                                f"{x.record.qID}\t{'novel-unannotated'}\t{gn}\t{i + 1:06d}\t{x.strand}\t{x.end}\t{x.start + 1}\n"
                            )
                    novel_region_index += 1

        logger.info(f"Output written to: {f.name}")
        logger.info(f"Output written to: {f1.name}")


@app.command(name="")
def main(
    gene_annotation_file: str = typer.Argument(..., help="Gene Annotation Text File"),
    input_fasta: str = typer.Argument(..., help="Input Fasta"),
    input_sam: str = typer.Argument(..., help="Input SAM"),
    output_prefix: str = typer.Argument(..., help="Output Prefix"),
    min_query_overlap: float = typer.Option(
        0.0,
        help="Minimum query overlap, in ratio",
    ),
    min_gene_overlap_bp: int = typer.Option(
        0,
        help="Minimum gene overlap, in bp",
    ),
    min_gene_overlap: float = typer.Option(
        0.5,
        help="Minimum gene overlap, in ratio",
    ),
):
    categorize_aln_by_annotation(
        gene_annotation_file,
        input_fasta,
        input_sam,
        output_prefix,
        min_gene_overlap_bp,
        min_query_overlap,
        min_gene_overlap,
    )


if __name__ == "__main__":
    typer.run(main)
