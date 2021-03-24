__author__ = "etseng@pacb.com"


import re
import sys
from collections import defaultdict, namedtuple
from csv import DictWriter
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from Bio import SeqIO
from bx.intervals import IntervalTree
from bx.intervals.cluster import ClusterTree
from cupcake.logging import cupcake_logger as logger
from cupcake.sequence import GFF
from cupcake.tofu import compare_junctions

seqid_rex = re.compile(r"(\\S+\\.\\d+)\\.(\\d+)")

MatchRecord = namedtuple(
    "MatchRecord", ["ref_id", "addon_id", "rec", "members", "seqrec"]
)


def find_representative_in_iso_list(records: List[GFF.gmapRecord]) -> GFF.gmapRecord:
    """
    :param records: list of GMAPRecord
    :return: representative record that is (a) the most number of exons or then (b) longest
    """
    rep = records[0]
    for r in records[1:]:
        if len(rep.ref_exons) < len(r.ref_exons) or (rep.end - rep.start) < (
            r.end - r.start
        ):
            rep = r
    return rep


def sanity_check_seqids(seqids: List[str]):
    for seqid in seqids:
        m = seqid_rex.match(seqid)
        if m is None:
            logger.error(
                f"Expected ID format (ex: PB.1.2) not followed by {seqid}! Abort!"
            )
            sys.exit(-1)


def get_fusion_id(seqid: str) -> str:
    m = seqid_rex.match(seqid)
    return m.group(1)


def write_reclist_to_gff_n_info(
    rec_list: Dict[str, Any],
    final_prefix: str,
    ref_name: str,
    addon_name: str,
    use_fq: bool = False,
) -> Dict[str, str]:
    # now go through the rec list and figure out in what order we are outputting the total records
    tree = defaultdict(lambda: {"+": ClusterTree(0, 0), "-": ClusterTree(0, 0)})
    tree_keys_numeric = set()
    tree_keys_alpha = set()
    for i, match_rec in enumerate(rec_list):
        tree[match_rec.rec.chr][match_rec.rec.strand].insert(
            match_rec.rec.start, match_rec.rec.end, i
        )

    for chrom in tree:
        try:
            k = int(chrom)
            tree_keys_numeric.add(k)
        except ValueError:
            tree_keys_alpha.add(chrom)
    tree_keys = sorted(tree_keys_numeric) + sorted(tree_keys_alpha)

    writer_info = DictWriter(
        Path(f"{final_prefix}.mega_info.txt".open("w")),
        fieldnames=["superPBID", ref_name, addon_name],
        delimiter="\t",
    )
    writer_info.writeheader()
    if use_fq:
        f_fq = Path(f"{final_prefix}.rep.fq")
    with open(f"{final_prefix}.gff", "w") as f_gff, open(
        f"{final_prefix}.group.txt", "w"
    ) as f_group:
        new_group_info = {}

        pb_i = 0
        for _chr in tree_keys:
            for _strand in ("+", "-"):
                for *_, _indices in tree[_chr][_strand].getregions():
                    # further sort these records by (start, end, num_exons)
                    _indices.sort(
                        key=lambda i: (
                            rec_list[i].rec.start,
                            rec_list[i].rec.end,
                            len(rec_list[i].rec.ref_exons),
                        )
                    )
                    pb_i += 1
                    for pb_j, recs_index in enumerate(_indices):
                        pbid = f"PB.{pb_i}.{pb_j + 1}"
                        match_rec = rec_list[recs_index]
                        new_group_info[pbid] = match_rec.members
                        match_rec.rec.seqid = pbid
                        GFF.write_collapseGFF_format(f_gff, match_rec.rec)
                        writer_info.writerow(
                            {
                                "superPBID": pbid,
                                ref_name: match_rec.ref_id,
                                addon_name: match_rec.addon_id,
                            }
                        )
                        f_group.write(f"{pbid}\t{','.join(match_rec.members)}\n")
                        if use_fq:
                            match_rec.seqrec.id = pbid
                            match_rec.seqrec.description = ""
                            SeqIO.write(match_rec.seqrec, f_fq, "fastq")

    return new_group_info


class MegaPBTree:
    """
    Structure for maintaining a non-redundant set of gene annotations
    Used to combine with different collapsed GFFs from different samples
    """

    def __init__(
        self,
        gff_filename: str,
        group_filename: str,
        internal_fuzzy_max_dist: int = 0,
        self_prefix: Optional[str] = None,
        allow_5merge: bool = False,
        fastq_filename: Optional[str] = None,
        max_3_diff: Optional[int] = None,
    ):
        self.gff_filename = gff_filename
        self.group_filename = group_filename
        self.self_prefix = self_prefix
        self.internal_fuzzy_max_dist = internal_fuzzy_max_dist
        self.max_3_diff = max_3_diff
        self.allow_5merge = allow_5merge
        self.record_d = {r.seqid: r for r in GFF.collapseGFFReader(gff_filename)}
        # sanity_check_seqids(self.record_d.keys()) # sanity check all IDs look like PB.1.2
        self.tree = defaultdict(
            lambda: {"+": IntervalTree(), "-": IntervalTree()}
        )  # chr --> strand --> tree
        self.fastq_dict = None
        if fastq_filename is not None:
            self.fastq_dict = MegaPBTree.read_fastq_to_dict(fastq_filename)

        # print >> sys.stderr, "self.internal_fuzzy_max_dist is", internal_fuzzy_max_dist
        # raw_input()
        self.read_gff_as_interval_tree()
        self.group_info = MegaPBTree.read_group(
            self.group_filename, self.self_prefix
        )  # ex: PB.1.1 --> [ RatHeart|i3_c123.... ]

    def read_gff_as_interval_tree(self):
        """
        Read a collapsed GFF file into an IntervalTree
        """
        for r in GFF.collapseGFFReader(self.gff_filename):
            self.tree[r.chr][r.strand].insert(r.start, r.end, r)

    @staticmethod
    def read_fastq_to_dict(
        fastq_filename: Union[str, Path]
    ) -> Dict[str, SeqIO.SeqRecord]:
        fastq_dict = {}
        for r in SeqIO.parse(open(fastq_filename), "fastq"):
            fastq_dict[r.id.split("|")[0]] = r
        return fastq_dict

    @staticmethod
    def read_group(
        group_filename: Union[str, Path], group_prefix: str
    ) -> Dict[str, List[str]]:
        group_info = {}
        with open(group_filename) as f:
            for line in f:
                pbid, members = line.strip().split("\t")
                if group_prefix is None:
                    group_info[pbid] = list(members.split(","))
                else:
                    group_info[pbid] = [
                        f"{group_prefix}|{x}" for x in members.split(",")
                    ]
        return group_info

    def match_record_to_tree(self, r: GFF.gmapRecord) -> GFF.gmapRecord:
        """
        r --- GMAPRecord
        tree --- dict of chromosome --> strand --> IntervalTree

        If exact match (every exon junction) or 5' truncated (allow_5merge is True), YIELD the matching GMAPRecord(s)
        *NOTE/UPDATE*: could have multiple matches! )
        """
        # if r.chr=='chr17' and r.start > 39604000:
        #    pdb.set_trace()
        matches = self.tree[r.chr][r.strand].find(r.start, r.end)
        for r2 in matches:
            r.segments = r.ref_exons
            r2.segments = r2.ref_exons
            n1 = len(r.segments)
            n2 = len(r2.segments)

            three_end_is_match = (
                self.max_3_diff is None
                or (r.strand == "+" and abs(r.end - r2.end) <= self.max_3_diff)
                or (r.strand == "-" and abs(r.start - r2.start) <= self.max_3_diff)
            )

            last_junction_match = False
            if n1 == 1:
                if n2 == 1:
                    last_junction_match = True
                else:
                    last_junction_match = False
            else:
                if n2 == 1:
                    last_junction_match = False
                else:
                    if r.strand == "+":
                        last_junction_match = (
                            abs(r.segments[-1].start - r2.segments[-1].start)
                            <= self.internal_fuzzy_max_dist
                        ) and (
                            abs(r.segments[0].end - r2.segments[0].end)
                            <= self.internal_fuzzy_max_dist
                        )
                    else:
                        last_junction_match = (
                            abs(r.segments[0].end - r2.segments[0].end)
                            <= self.internal_fuzzy_max_dist
                        ) and (
                            abs(r.segments[1].start - r2.segments[1].start)
                            <= self.internal_fuzzy_max_dist
                        )

            if (
                compare_junctions.compare_junctions(
                    r, r2, internal_fuzzy_max_dist=self.internal_fuzzy_max_dist
                )
                == "exact"
            ):  # is a match!
                if three_end_is_match:
                    yield r2
            elif (
                self.allow_5merge
            ):  # check if the shorter one is a subset of the longer one
                if len(r.segments) > len(r2.segments):
                    a, b = r, r2
                else:
                    a, b = r2, r
                # a is the longer one, b is the shorter one
                if (
                    compare_junctions.compare_junctions(
                        b, a, internal_fuzzy_max_dist=self.internal_fuzzy_max_dist
                    )
                    == "subset"
                ):
                    # we only know that a is a subset of b, verify that it is actually 5' truncated (strand-sensitive!)
                    # if + strand, last junction of (a,b) should match and 3' end not too diff
                    # if - strand, first exon of a should match first exon of b AND the next exon don't overlap
                    if three_end_is_match and last_junction_match:
                        yield r2

    def add_sample(
        self,
        gff_filename: Union[str, Path],
        group_filename: Union[str, Path],
        sample_prefix: str,
        output_prefix: str,
        fastq_filename: Union[str, Path] = None,
    ) -> None:
        combined = []  # list of (<matches to r2 or None>, r2)
        unmatched_recs = set(self.record_d.keys())

        for r in GFF.collapseGFFReader(gff_filename):
            match_rec_list = list(self.match_record_to_tree(r))
            if len(match_rec_list) > 0:  # found match(es)! put longer of r1/r2 in
                # if len(match_rec_list) > 1: pdb.set_trace()  #DEBUG
                combined.append((match_rec_list, r))
                for match_rec in match_rec_list:
                    try:
                        unmatched_recs.remove(match_rec.seqid)
                    except KeyError:
                        pass  # already deleted, OK, this can happen
            else:  # r is not present in current tree
                combined.append((None, r))
        # put whatever is left from the tree in
        for seqid in unmatched_recs:
            combined.append(([self.record_d[seqid]], None))

        # create a ClusterTree to re-calc the loci/transcripts
        final_tree = defaultdict(
            lambda: {"+": ClusterTree(0, 0), "-": ClusterTree(0, 0)}
        )
        for i, (r1s, r2) in enumerate(combined):
            if r1s is None:
                final_tree[r2.chr][r2.strand].insert(r2.start, r2.end, i)
            else:
                if r2 is not None:
                    rep = find_representative_in_iso_list(r1s + [r2])
                else:
                    rep = find_representative_in_iso_list(r1s)
                final_tree[rep.chr][rep.strand].insert(rep.start, rep.end, i)

        self.write_cluster_tree_as_gff(
            rec_list=combined,
            group_filename2=group_filename,
            sample_prefix2=sample_prefix,
            output_prefix=output_prefix,
            fastq_filename2=fastq_filename,
        )

    def write_cluster_tree_as_gff(
        self,
        rec_list: List[GFF.gmapRecord],
        group_filename2: Union[str, Path],
        sample_prefix2: str,
        output_prefix: str,
        fastq_filename2: Union[str, Path] = None,
    ) -> Dict[str, str]:
        """
        Write ClusterTree (chr --> dict --> (start, end, rec_list_index)) as collapsedGFF format
        Returns --- a new group_info!!!
        """
        use_fq = fastq_filename2 is not None and self.fastq_dict is not None
        if use_fq:
            fastq_dict2 = MegaPBTree.read_fastq_to_dict(fastq_filename2)
        group_info2 = MegaPBTree.read_group(group_filename2, sample_prefix2)

        # currently: rec_list is (r1s, r2) where r1s, r2 are records and could be None
        # make rec_list into list of MatchRec (ref_id, addon_id, representative rec, seqrec, group_info members)
        new_rec_list = []
        for r1s, r2 in rec_list:
            if r2 is None:
                for r1 in r1s:
                    new_rec_list.append(
                        MatchRecord(
                            ref_id=r1.seqid,
                            addon_id="NA",
                            rec=r1,
                            members=self.group_info[r1.seqid],
                            seqrec=self.fastq_dict[r1.seqid] if use_fq else None,
                        )
                    )
            elif r1s is None:
                new_rec_list.append(
                    MatchRecord(
                        ref_id="NA",
                        addon_id=r2.seqid,
                        rec=r2,
                        members=group_info2[r2.seqid],
                        seqrec=fastq_dict2[r2.seqid] if use_fq else None,
                    )
                )
            else:
                for r1 in r1s:
                    if len(r1s) > 1:
                        logger.info(f"matching {r1} to {r2}")
                    rep = find_representative_in_iso_list([r1, r2])
                    new_rec_list.append(
                        MatchRecord(
                            ref_id=r1.seqid,
                            addon_id=r2.seqid,
                            rec=rep,
                            members=self.group_info[r1.seqid] + group_info2[r2.seqid],
                            seqrec=self.fastq_dict[rep.seqid] if use_fq else None,
                        )
                    )
        new_group_info = write_reclist_to_gff_n_info(
            new_rec_list, output_prefix, self.self_prefix, sample_prefix2, use_fq
        )
        return new_group_info


class MegaPBTreeFusion(MegaPBTree):
    def __init__(
        self,
        gff_filename: Union[str, Path],
        group_filename: Union[str, Path],
        internal_fuzzy_max_dist: int = 0,
        self_prefix: str = None,
        fastq_filename: Union[str, Path] = None,
        fusion_max_dist: int = 10,
    ):
        """
        Differences with non-fusion MegaPBTree:

        1. allow_5merge is always FALSE. Not a parameter.
        2. fusion_max_dist --- maximum allowed distance on internal fusion sites to be called as equivalent fusions
        """
        super().__init__(
            gff_filename,
            group_filename,
            internal_fuzzy_max_dist,
            self_prefix,
            False,
            fastq_filename,
        )

        self.fusion_max_dist = fusion_max_dist

        # ex: PBfusion.1 -> [PBfusion.1.1, PBfusion.1.2]
        self.record_d_fusion = {
            fusion_id: records
            for fusion_id, records in GFF.collapseGFFFusionReader(gff_filename)
        }

    def junction_match_check_5(self, r1: GFF.gmapRecord, r2: GFF.gmapRecord) -> bool:
        if r1.strand == "+":
            return (
                abs(r1.ref_exons[0].start - r2.ref_exons[0].start)
                <= self.fusion_max_dist
            )
        else:
            return (
                abs(r1.ref_exons[-1].end - r2.ref_exons[-1].end) <= self.fusion_max_dist
            )

    def junction_match_check_3(self, r1: GFF.gmapRecord, r2: GFF.gmapRecord) -> bool:
        if r1.strand == "+":
            return (
                abs(r1.ref_exons[-1].end - r2.ref_exons[-1].end) <= self.fusion_max_dist
            )
        else:
            return (
                abs(r1.ref_exons[0].start - r2.ref_exons[0].start)
                <= self.fusion_max_dist
            )

    def match_record_to_tree(
        self, r: GFF.gmapRecord, check_5_dist: bool, check_3_dist: bool
    ) -> List[str]:
        """
        Matching a single record (locus).

        Major diff from non-fusion version:
        1. there could be multiple matches!
        2. no 5merge allowed
        3. additionally checks if the 5'/3' ends don't disagree too much (fusion_max_dist). this is used for fusion junctions.
        4. need to take care that fusions can be multi-chromosome! write output correctly!!!
        """
        matches = self.tree[r.chr][r.strand].find(r.start, r.end)
        result = []
        for r2 in matches:
            r.segments = r.ref_exons
            r2.segments = r2.ref_exons
            if (
                compare_junctions.compare_junctions(
                    r, r2, internal_fuzzy_max_dist=self.internal_fuzzy_max_dist
                )
                == "exact"
                and (not check_5_dist or self.junction_match_check_5(r, r2))
                and (not check_3_dist or self.junction_match_check_3(r, r2))
            ):  # is a match!
                result.append(r2.seqid)

        return result

    def check_records_match(
        self, records1: GFF.gmapRecord, records2: GFF.gmapRecord
    ) -> bool:
        """
        records1, records2 are two fusion records.
        They match iff:
        1. same number of records
        2. each record (a loci) matches
        """
        if len(records1) != len(records2):
            return False

        i = 0
        for r1, r2 in zip(records1, records2):
            # check: chr, strand, exons match
            if r1.chr != r2.chr or r1.strand != r2.strand:
                return False
            r1.segments = r1.ref_exons
            r2.segments = r2.ref_exons
            if (
                compare_junctions.compare_junctions(
                    r1, r2, internal_fuzzy_max_dist=self.internal_fuzzy_max_dist
                )
                != "exact"
            ):
                return False
            if i == 0:  # first record, only need 3' to agree
                if not self.junction_match_check_3(r1, r2):
                    return False
            elif i == len(records1) - 1:  # last record, only need 5' to agree
                if not self.junction_match_check_5(r1, r2):
                    return False
            else:
                if not self.junction_match_check_5(r1, r2):
                    return False
                if not self.junction_match_check_3(r1, r2):
                    return False
            i += 1

        return True

    def match_fusion_record(
        self, records: List[GFF.gmapRecord]
    ) -> Optional[GFF.gmapRecord]:
        """
        records --- in order, the records of a single fusion.
        """
        good = []
        # match the first record, requiring additionally that the precise 3' end matches
        cands = self.match_record_to_tree(
            records[0], check_5_dist=False, check_3_dist=True
        )
        # for each candidate (ex: PB.8.1, extract the full set of records and match them)
        for cand in cands:
            m = seqid_rex.match(cand)
            fusion_id = m.group(1)
            if self.check_records_match(records, self.record_d_fusion[fusion_id]):
                good.append(fusion_id)
        if len(good) == 0:
            return None
        elif len(good) == 1:
            return good[0]
        else:
            logger.error(
                "ERROR! more than one possible candidate in match_fusion_record! DEBUG."
            )
            logger.error(f"MATCHED: {good}")
            sys.exit(-1)

    def add_sample(
        self,
        gff_filename: Union[str, Path],
        group_filename: Union[str, Path],
        sample_prefix: str,
        output_prefix: str,
        fastq_filename: Optional[Union[str, Path]] = None,
    ) -> None:
        combined = (
            []
        )  # list of (r1 if r2 is None | r2 if r1 is None | longer of r1 or r2 if both not None)
        unmatched_recs = list(self.record_d_fusion.keys())

        for _, records in GFF.collapseGFFFusionReader(gff_filename):
            match_seqid = self.match_fusion_record(records)
            if match_seqid is not None:
                combined.append((self.record_d_fusion[match_seqid], records))
                try:
                    unmatched_recs.remove(match_seqid)
                except ValueError:
                    pass  # already deleted, OK, this happens for single-exon transcripts
            else:  # r is not present in current tree
                combined.append((None, records))
        # put whatever is left from the tree in
        for seqid in unmatched_recs:
            combined.append((self.record_d_fusion[seqid], None))

        # create a ClusterTree to re-calc the loci/transcripts
        final_tree = defaultdict(
            lambda: {"+": ClusterTree(0, 0), "-": ClusterTree(0, 0)}
        )
        for i, (r1s, r2s) in enumerate(combined):
            if r2s is None or (
                r1s is not None
                and r1s[0].end - r1s[0].start > r2s[0].end - r2s[0].start
            ):
                final_tree[r1s[0].chr][r1s[0].strand].insert(
                    r1s[0].start, r1s[0].end, i
                )
            else:
                final_tree[r2s[0].chr][r2s[0].strand].insert(
                    r2s[0].start, r2s[0].end, i
                )

        self.write_cluster_tree_as_gff(
            final_tree,
            combined,
            group_filename,
            sample_prefix,
            output_prefix,
            fastq_filename2=fastq_filename,
        )

    def write_cluster_tree_as_gff(
        self,
        cluster_tree: ClusterTree,
        rec_list: List[GFF.gmapRecord],
        group_filename2: Union[str, Path],
        sample_prefix2: str,
        output_prefix: str,
        fastq_filename2: Optional[Union[str, Path]] = None,
    ) -> Dict[str, str]:
        """
        Write ClusterTree (chr --> dict --> (start, end, rec_list_index)) as collapsedGFF format
        Returns --- a new group_info!!!
        """
        if fastq_filename2 is not None:
            fastq_dict2 = MegaPBTree.read_fastq_to_dict(fastq_filename2)
            f_fastq = Path(f"{output_prefix}.rep.fq")
        group_info2 = MegaPBTree.read_group(group_filename2, sample_prefix2)
        new_group_info = {}

        with open(f"{output_prefix}.mega_info.txt", "w") as f_mgroup:
            f_mgroup.write(f"pbid\t{self.self_prefix}\t{sample_prefix2}\n")
            fusion_index = 0
            chroms = list(cluster_tree.keys())
            chroms.sort()

            for (
                k
            ) in (
                chroms
            ):  # IMPORTANT: for fusion, this is *just* the chrom of the first record! Fusions can be multi-chrom
                for strand in ("+", "-"):
                    for *_, rec_indices in cluster_tree[k][strand].getregions():
                        for i in rec_indices:
                            fusion_index += 1
                            tID = f"PBfusion.{fusion_index}"
                            r1s, r2s = rec_list[i]
                            if r1s is None:  # r2s is not None
                                recs = r2s
                                r2_fusion_id = get_fusion_id(r2s[0].seqid)
                                new_group_info[tID] = group_info2[r2_fusion_id]
                                f_mgroup.write(f"{tID}\tNA\t{r2_fusion_id}\n")
                                if fastq_filename2 is not None:
                                    seqrec = fastq_dict2[r2_fusion_id]
                            elif r2s is None:  # r1 is not None
                                recs = r1s
                                r1_fusion_id = get_fusion_id(r1s[0].seqid)
                                new_group_info[tID] = self.group_info[r1_fusion_id]
                                f_mgroup.write(f"{tID}\t{r1_fusion_id}\tNA\n")
                                if fastq_filename2 is not None:
                                    seqrec = self.fastq_dict[r1_fusion_id]
                            else:  # both r1, r2 are not empty
                                r1_fusion_id = get_fusion_id(r1s[0].seqid)
                                r2_fusion_id = get_fusion_id(r2s[0].seqid)
                                r1_len = sum(x.end - x.start for x in r1s)
                                r2_len = sum(x.end - x.start for x in r2s)
                                if r1_len > r2_len:
                                    recs = r1s
                                    if fastq_filename2 is not None:
                                        seqrec = self.fastq_dict[r1_fusion_id]
                                else:
                                    recs = r2s
                                    if fastq_filename2 is not None:
                                        seqrec = fastq_dict2[r2_fusion_id]
                                new_group_info[tID] = (
                                    self.group_info[r1_fusion_id]
                                    + group_info2[r2_fusion_id]
                                )
                                f_mgroup.write(
                                    f"{tID}\t{r1_fusion_id}\t{r2_fusion_id}\n"
                                )

                            if fastq_filename2 is not None:
                                seqrec.id = tID
                                SeqIO.write(seqrec, f_fastq.open("w"), "fastq")

                            with open(f"{output_prefix}.group.txt", "w") as f_group:
                                f_group.write(
                                    f"{tID}\t{','.join(new_group_info[tID])}\n"
                                )

                            with open(f"{output_prefix}.gff", "w") as f_out:
                                # now write out the fusion transcript
                                for j, r in enumerate(recs):
                                    f_out.write(
                                        f'{r.chr}\tPacBio\ttranscript\t{r.start + 1}\t{r.end}\t.\t{strand}\t.\tgene_id "{tID}"; transcript_id "{tID}.{j + 1}";\n'
                                    )
                                    for exon in r.ref_exons:
                                        f_out.write(
                                            f'{r.chr}\tPacBio\texon\t{exon.start + 1}\t{exon.end}\t.\t{strand}\t.\tgene_id "{tID}"; transcript_id "{tID}.{j + 1}";\n'
                                        )
        return new_group_info
