#!/usr/bin/env python
__author__ = "etseng@pacb.com"

"""
Identical to the collapse script provided in cDNA_primer (ToFU) GitHub.

Takes a SAM file and the input fasta/fastq used to produce the SAM file,
filter out alignments based on low coverage/identity and collapse/merge
any identical isoforms based on the aligned exonic structure.

Example:
collapse_isoforms_by_sam.py --input test.fq --fq -s test.fq.sorted.sam --dun-merge-5-shorter -o test

Suggested scripts to follow up with:
   get_abundance_post_collapse.py: create count (absolute and normalized) information post collapse
   filter_by_count.py: filter away based on FL count support
   filter_away_subset.py (if collapse is run with --dun-merge-5-shorter)
"""

from collections import defaultdict
from gzip import open as gzopen
from pathlib import Path
from typing import Optional, Union

import typer
from Bio import SeqIO
from bx.intervals import IntervalTree
from cupcake.logging import cupcake_logger as logger
from cupcake.sequence import GFF
from cupcake.tofu.branch import branch_simple2
from cupcake.tofu.compare_junctions import compare_junctions
from cupcake.tofu.utils import check_ids_unique

GFF_FIELDS = [
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute",
]


app = typer.Typer(name="collapse_isoforms_by_sam", add_completion=False)


def pick_rep(
    fa_fq_filename: str,
    gff_filename: Union[str, Path],
    group_filename: Union[str, Path],
    output_filename: Union[str, Path],
    is_fq: bool = False,
    pick_least_err_instead: bool = False,
    bad_gff_filename: Optional[str] = None,
) -> None:
    """
    For each group, select the representative record

    If is FASTA file (is_fa False) -- then always pick the longest one
    If is FASTQ file (is_fq True) -- then
          If pick_least_err_instead is True, pick the one w/ least number of expected base errors
          Else, pick the longest one
    """
    if Path(fa_fq_filename).suffix == ".gz":
        fd = SeqIO.to_dict(
            SeqIO.parse(gzopen(fa_fq_filename), "fastq" if is_fq else "fasta")
        )
    else:
        fd = SeqIO.to_dict(
            SeqIO.parse(open(fa_fq_filename), "fastq" if is_fq else "fasta")
        )

    coords = {}
    for line in open(gff_filename):
        #     0       1       2               3       4       5       6       7       8
        #     seqname source  feature         start   end     score   strand  frame   attribute
        # ex: chr1    PacBio  transcript      27567   29336   .       -       .       gene_id "PB.1"; transcript_id "PB.1.1";

        raw = {k: v for k, v in zip(GFF_FIELDS, line.strip().split("\t"))}
        if raw["feature"] == "transcript":
            tid = raw["attribute"].split("; ")[1].split()[1][1:-2]
            coords[
                tid
            ] = f'{raw["seqname"]}:{raw["start"]}-{raw["end"]}({raw["strand"]})'

    if bad_gff_filename is not None:
        for line in open(bad_gff_filename):
            raw = {k: v for k, v in zip(GFF_FIELDS, line.strip().split("\t"))}
            if raw["feature"] == "transcript":
                tid = raw["attribute"].split("; ")[1].split()[1][1:-2]
                coords[
                    tid
                ] = f'{raw["seqname"]}:{raw["start"]}-{raw["end"]}({raw["strand"]})'

    with open(output_filename, "w") as fout:
        for line in open(group_filename):
            pb_id, members = line.strip().split("\t")
            logger.info(f"Picking representative sequence for {pb_id}")
            best_rec = None
            # best_id = None
            # best_seq = None
            # best_qual = None
            best_err = 9999999
            err = 9999999
            max_len = 0
            for x in members.split(","):
                if is_fq and pick_least_err_instead:
                    err = sum(
                        i ** -(i / 10.0)
                        for i in fd[x].letter_annotations["phred_quality"]
                    )
                if (is_fq and pick_least_err_instead and err < best_err) or (
                    (not is_fq or not pick_least_err_instead)
                    and len(fd[x].seq) >= max_len
                ):
                    best_rec = fd[x]
                    # best_id = x
                    # best_seq = fd[x].seq
                    # if is_fq:
                    #    best_qual = fd[x].quality
                    #    best_err = err
                    max_len = len(fd[x].seq)

            _id_ = f"{pb_id}|{coords[pb_id]}|{best_rec.id}"
            best_rec.id = _id_
            SeqIO.write(best_rec, fout, "fastq" if is_fq else "fasta")


def collapse_fuzzy_junctions(
    gff_filename: Union[str, Path],
    group_filename: Union[str, Path],
    allow_extra_5exon: bool,
    internal_fuzzy_max_dist: int,
    max_5_diff: int,
    max_3_diff: int,
) -> defaultdict:
    def can_merge(m, r1, r2):
        if m == "exact":
            return True
        else:
            if not allow_extra_5exon:
                return False
        # below is continued only if (a) is 'subset' or 'super' AND (b) allow_extra_5exon is True
        if m == "subset":
            r1, r2 = r2, r1  # rotate so r1 is always the longer one
        if m == "super" or m == "subset":
            n2 = len(r2.ref_exons)
            # check that (a) r1 and r2 end on same 3' exon, that is the last acceptor site agrees
            # AND (b) the 5' start of r2 is sandwiched between the matching r1 exon coordinates
            if r1.strand == "+":
                return (
                    abs(r1.ref_exons[-1].start - r2.ref_exons[-1].start)
                    <= internal_fuzzy_max_dist
                    and r1.ref_exons[-n2].start
                    <= r2.ref_exons[0].start
                    < r1.ref_exons[-n2].end
                )
            else:
                return (
                    abs(r1.ref_exons[0].end - r2.ref_exons[0].end)
                    <= internal_fuzzy_max_dist
                    and r1.ref_exons[n2 - 1].start
                    <= r2.ref_exons[-1].end
                    < r1.ref_exons[n2].end
                )
        return False

    d = {}
    # chr --> strand --> tree
    recs = defaultdict(lambda: {"+": IntervalTree(), "-": IntervalTree()})
    fuzzy_match = defaultdict(lambda: [])
    for r in GFF.collapseGFFReader(gff_filename):
        d[r.seqid] = r
        has_match = False
        r.segments = r.ref_exons
        for r2 in recs[r.chr][r.strand].find(r.start, r.end):
            r2.segments = r2.ref_exons
            m = compare_junctions(
                r,
                r2,
                internal_fuzzy_max_dist=internal_fuzzy_max_dist,
                max_5_diff=max_5_diff,
                max_3_diff=max_3_diff,
            )
            if can_merge(m, r, r2):
                fuzzy_match[r2.seqid].append(r.seqid)
                has_match = True
                break
        if not has_match:
            recs[r.chr][r.strand].insert(r.start, r.end, r)
            fuzzy_match[r.seqid] = [r.seqid]

    group_info = {}
    with open(group_filename) as f:
        for line in f:
            pbid, members = line.strip().split("\t")
            group_info[pbid] = members.split(",")

    # pick for each fuzzy group the one that has the most exons
    keys = list(fuzzy_match.keys())
    keys.sort(key=lambda x: int(x.split(".")[1]))

    with open(f"{gff_filename}.fuzzy", "w") as f_gff, open(
        f"{group_filename}.fuzzy", "w"
    ) as f_group:
        for k in keys:
            all_members = []
            best_pbid, best_size, best_num_exons = (
                fuzzy_match[k][0],
                len(group_info[fuzzy_match[k][0]]),
                len(d[fuzzy_match[k][0]].ref_exons),
            )
            all_members += group_info[fuzzy_match[k][0]]
            for pbid in fuzzy_match[k][1:]:
                _num_exons = len(d[pbid].ref_exons)
                _size = len(group_info[pbid])
                all_members += group_info[pbid]
                if _num_exons > best_num_exons or (
                    _num_exons == best_num_exons and _size > best_size
                ):
                    best_pbid, best_size, best_num_exons = pbid, _size, _num_exons
            GFF.write_collapseGFF_format(f_gff, d[best_pbid])
            f_group.write(f'{best_pbid}\t{",".join(all_members)}\n')

    return fuzzy_match


@app.command()
def main(
    input_filename: str = typer.Option(..., "--input", help="Input FA/FQ filename"),
    sam: str = typer.Option(..., help="Sorted GMAP SAM filename"),
    fq: bool = typer.Option(False, "--fq", help="Input is a fastq file"),  # store_true
    prefix: str = typer.Option(..., "-p", "--prefix", help="Output filename prefix"),
    min_aln_coverage: float = typer.Option(
        0.99, "--min-coverage", "-c", help="Minimum alignment coverage"
    ),
    min_aln_identity: float = typer.Option(
        0.95, "--min-identity", "-i", help="Minimum alignment identity"
    ),
    max_fuzzy_junction: int = typer.Option(5, help="Max fuzzy junction dist"),
    max_5_diff: int = typer.Option(
        1000, help="Maximum allowed 5' difference if on same exon"
    ),
    max_3_diff: int = typer.Option(
        100, help="Maximum allowed 3' difference if on same exon"
    ),
    flnc_coverage: int = typer.Option(
        -1,
        help="Minimum # of FLNC reads, only use this for aligned FLNC reads, otherwise results undefined!",
    ),
    gen_mol_count: bool = typer.Option(
        False,
        help="Generate a .abundance.txt file based on the number of input sequences collapsed. Use only if input is FLNC or UMI-dedup output",
    ),  # store_true
    allow_extra_5exon: bool = typer.Option(
        True,
        "--dun-merge-5-shorter",
        help="Don't collapse shorter 5' transcripts (default: turned off)",
    ),  # store_false
) -> None:
    # sanity check that input file and input SAM exists
    if not Path(str(input_filename)).exists():
        raise FileNotFoundError(f"Input file {input_filename} does not exist. Abort.")

    if not Path(sam).exists():
        raise FileNotFoundError(f"SAM file {sam} does not exist. Abort.")

    # check for duplicate IDs
    check_ids_unique(input_filename, is_fq=fq)

    with open(f"{prefix}.ignored_ids.txt", "w") as ignored_fout:

        if flnc_coverage > 0:
            # keep these files closed *until* we need to write to them
            f_good = Path(f"{prefix}.collapsed.good.gff")
            f_bad = Path(f"{prefix}.collapsed.bad.gff")
            cov_threshold = flnc_coverage
        else:
            f_good = Path(f"{prefix}.collapsed.gff")
            f_bad = f_good
            cov_threshold = 1
        f_txt = Path(f"{prefix}.collapsed.group.txt")

        b = branch_simple2.BranchSimple(
            transfrag_filename=input_filename,
            cov_threshold=cov_threshold,
            min_aln_coverage=min_aln_coverage,
            min_aln_identity=min_aln_identity,
            is_fq=fq,
            max_5_diff=max_5_diff,
            max_3_diff=max_3_diff,
        )
        rec_iter = b.iter_gmap_sam(sam, ignored_fout)
        # recs is {'+': list of list of records, '-': list of list of records}
        for recs in rec_iter:
            for v in recs.values():
                for v2 in v:
                    if len(v2) > 0:
                        b.process_records(
                            records=v2,
                            allow_extra_5_exons=allow_extra_5exon,
                            skip_5_exon_alt=False,
                            f_good=f_good,
                            f_bad=f_bad,
                            f_group=f_txt,
                        )

    # need to further collapse those that have fuzzy junctions!
    if max_fuzzy_junction > 0:
        collapse_fuzzy_junctions(
            f_good,
            f_txt,
            allow_extra_5exon,
            internal_fuzzy_max_dist=max_fuzzy_junction,
            max_5_diff=max_5_diff,
            max_3_diff=max_3_diff,
        )
        Path(f_good.name).rename(f"{f_good.name}.unfuzzy")
        Path(f_txt.name).rename(f"{f_txt.name}.unfuzzy")
        Path(f"{f_good.name}.fuzzy").rename(f_good.name)
        Path(f"{f_txt.name}.fuzzy").rename(f_txt.name)

    if fq:
        outfile = f"{prefix}.collapsed.rep.fq"
    else:
        outfile = f"{prefix}.collapsed.rep.fa"
    if allow_extra_5exon:  # 5merge, pick longest
        pick_rep(
            fa_fq_filename=input_filename,
            gff_filename=f_good,
            group_filename=f_txt,
            output_filename=outfile,
            is_fq=fq,
            pick_least_err_instead=False,
            bad_gff_filename=f_bad.name,
        )
    else:
        pick_rep(
            fa_fq_filename=input_filename,
            gff_filename=f_good,
            group_filename=f_txt,
            output_filename=outfile,
            is_fq=fq,
            pick_least_err_instead=True,
            bad_gff_filename=f_bad.name,
        )

    if gen_mol_count:
        outfile = f"{prefix}.collapsed.abundance.txt"
        with open(outfile, "w") as f:
            f.write("pbid\tcount_fl\n")
            for line in open(f_txt.name):
                pbid, members = line.strip().split()
                f.write(f"{pbid}\t{members.count(',')+1}\n")

    logger.info(f"Ignored IDs written to: {ignored_fout.name}")
    logger.info(f"Output written to: {f_good.name}\n{f_txt.name}\n{outfile}\n")


if __name__ == "__main__":
    typer.run(main)
