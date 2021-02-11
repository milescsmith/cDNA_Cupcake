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

import sys
from gzip import open as gzopen
from collections import defaultdict
from pathlib import Path

from bx.intervals import IntervalTree
from Bio import SeqIO

from cupcake.tofu.utils import check_ids_unique
from cupcake.tofu.branch import branch_simple2
from cupcake.tofu.compare_junctions import compare_junctions
from cupcake.sequence import GFF

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


def pick_rep(
    fa_fq_filename,
    gff_filename,
    group_filename,
    output_filename,
    is_fq=False,
    pick_least_err_instead=False,
    bad_gff_filename=None,
):
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
    fout = open(output_filename, "w")

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

    for line in open(group_filename):
        pb_id, members = line.strip().split("\t")
        print(f"Picking representative sequence for {pb_id}", file=sys.stdout)
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
                    i ** -(i / 10.0) for i in fd[x].letter_annotations["phred_quality"]
                )
            if (is_fq and pick_least_err_instead and err < best_err) or (
                (not is_fq or not pick_least_err_instead) and len(fd[x].seq) >= max_len
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

    fout.close()


def collapse_fuzzy_junctions(
    gff_filename,
    group_filename,
    allow_extra_5exon,
    internal_fuzzy_max_dist,
    max_5_diff,
    max_3_diff,
):
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
    recs = defaultdict(
        lambda: {"+": IntervalTree(), "-": IntervalTree()}
    )  # chr --> strand --> tree
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
            group_info[pbid] = [x for x in members.split(",")]

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


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("--input", help="Input FA/FQ filename")
    parser.add_argument(
        "--fq",
        default=False,
        action="store_true",
        help="Input is a fastq file (default is fasta)",
    )
    parser.add_argument("-s", "--sam", required=True, help="Sorted GMAP SAM filename")
    parser.add_argument("-o", "--prefix", required=True, help="Output filename prefix")
    parser.add_argument(
        "-c",
        "--min-coverage",
        dest="min_aln_coverage",
        type=float,
        default=0.99,
        help="Minimum alignment coverage (default: 0.99)",
    )
    parser.add_argument(
        "-i",
        "--min-identity",
        dest="min_aln_identity",
        type=float,
        default=0.95,
        help="Minimum alignment identity (default: 0.95)",
    )
    parser.add_argument(
        "--max_fuzzy_junction",
        default=5,
        type=int,
        help="Max fuzzy junction dist (default: 5 bp)",
    )
    parser.add_argument(
        "--max_5_diff",
        default=1000,
        type=int,
        help="Maximum allowed 5' difference if on same exon (default: 1000 bp)",
    )
    parser.add_argument(
        "--max_3_diff",
        default=100,
        type=int,
        help="Maximum allowed 3' difference if on same exon (default: 100 bp)",
    )
    parser.add_argument(
        "--flnc_coverage",
        dest="flnc_coverage",
        type=int,
        default=-1,
        help="Minimum # of FLNC reads, only use this for aligned FLNC reads, otherwise results undefined!",
    )
    parser.add_argument(
        "--dun-merge-5-shorter",
        action="store_false",
        dest="allow_extra_5exon",
        default=True,
        help="Don't collapse shorter 5' transcripts (default: turned off)",
    )

    args = parser.parse_args()
    # sanity check that input file and input SAM exists
    if not Path(args.input).exists():
        print(f"Input file {args.input} does not exist. Abort.", file=sys.stderr)
        sys.exit(-1)

    if not Path(args.sam).exists:
        print(f"SAM file {args.sam} does not exist. Abort.", file=sys.stderr)
        sys.exit(-1)

    # check for duplicate IDs
    check_ids_unique(args.input, is_fq=args.fq)

    ignored_fout = open(f"{args.prefix}.ignored_ids.txt", "w")

    if args.flnc_coverage > 0:
        f_good = open(f"{args.prefix}.collapsed.good.gff", "w")
        f_bad = open(f"{args.prefix}.collapsed.bad.gff", "w")
        cov_threshold = args.flnc_coverage
    else:
        f_good = open(f"{args.prefix}.collapsed.gff", "w")
        f_bad = f_good
        cov_threshold = 1
    f_txt = open(f"{args.prefix}.collapsed.group.txt", "w")

    b = branch_simple2.BranchSimple(
        args.input,
        cov_threshold=cov_threshold,
        min_aln_coverage=args.min_aln_coverage,
        min_aln_identity=args.min_aln_identity,
        is_fq=args.fq,
        max_5_diff=args.max_5_diff,
        max_3_diff=args.max_3_diff,
    )
    iter = b.iter_gmap_sam(args.sam, ignored_fout)
    for (
        recs
    ) in iter:  # recs is {'+': list of list of records, '-': list of list of records}
        for v in recs.values():
            for v2 in v:
                if len(v2) > 0:
                    b.process_records(
                        v2, args.allow_extra_5exon, False, f_good, f_bad, f_txt
                    )

    ignored_fout.close()
    f_good.close()
    try:
        f_bad.close()
    except:
        pass
    f_txt.close()

    if (
        args.max_fuzzy_junction > 0
    ):  # need to further collapse those that have fuzzy junctions!
        collapse_fuzzy_junctions(
            f_good.name,
            f_txt.name,
            args.allow_extra_5exon,
            internal_fuzzy_max_dist=args.max_fuzzy_junction,
            max_5_diff=args.max_5_diff,
            max_3_diff=args.max_3_diff,
        )
        Path(f_good.name).rename(f"{f_good.name}.unfuzzy")
        Path(f_txt.name).rename(f"{f_txt.name}.unfuzzy")
        Path(f"{f_good.name}.fuzzy").rename(f_good.name)
        Path(f"{f_txt.name_}.fuzzy").rename(f_txt.name)

    if args.fq:
        outfile = f"{args.prefix}.collapsed.rep.fq"
    else:
        outfile = f"{args.prefix}.collapsed.rep.fa"
    if args.allow_extra_5exon:  # 5merge, pick longest
        pick_rep(
            args.input,
            f_good.name,
            f_txt.name,
            outfile,
            is_fq=args.fq,
            pick_least_err_instead=False,
            bad_gff_filename=f_bad.name,
        )
    else:
        pick_rep(
            args.input,
            f_good.name,
            f_txt.name,
            outfile,
            is_fq=args.fq,
            pick_least_err_instead=True,
            bad_gff_filename=f_bad.name,
        )

    print(f"Ignored IDs written to: {ignored_fout.name}", file=sys.stdout)
    newline = "\n"
    print(
        f"Output written to: {f_good.name}{newline}{f_txt.name}{newline}{outfile}{newline}{args}{newline}",
        file=sys.stdout,
    )


if __name__ == "__main__":
    main()