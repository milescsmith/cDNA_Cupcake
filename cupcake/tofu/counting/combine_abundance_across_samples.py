__author__ = 'etseng@pacb.com'


import os, sys, re
from Bio import SeqIO
from collections import defaultdict
from cupcake.tofu import compare_junctions
from cupcake.io import GFF
from bx.intervals import IntervalTree
from bx.intervals.cluster import ClusterTree

seqid_rex = re.compile('(\S+\\.\d+)\\.(\d+)')

def sanity_check_seqids(seqids):
    for seqid in seqids:
        m = seqid_rex.match(seqid)
        if m is None:
            print >> sys.stderr, "Expected ID format (ex: PB.1.2) not followed by {0}! Abort!".format(seqid)
            sys.exit(-1)

def get_fusion_id(seqid):
    m = seqid_rex.match(seqid)
    return m.group(1)


class MegaPBTree(object):
    """
    Structure for maintaining a non-redundant set of gene annotations
    Used to combine with different collapsed GFFs from different samples
    """
    def __init__(self, gff_filename, group_filename, internal_fuzzy_max_dist=0, self_prefix=None, allow_5merge=False, fastq_filename=None):
        self.gff_filename = gff_filename
        self.group_filename = group_filename
        self.self_prefix = self_prefix
        self.internal_fuzzy_max_dist = internal_fuzzy_max_dist
        self.allow_5merge = allow_5merge
        self.record_d = dict((r.seqid, r) for r in GFF.collapseGFFReader(gff_filename))
        #sanity_check_seqids(self.record_d.keys()) # sanity check all IDs look like PB.1.2
        self.tree = defaultdict(lambda: {'+':IntervalTree(), '-':IntervalTree()}) # chr --> strand --> tree
        self.fastq_dict = None
        if fastq_filename is not None:
            self.fastq_dict = MegaPBTree.read_fastq_to_dict(fastq_filename)


        #print >> sys.stderr, "self.internal_fuzzy_max_dist is", internal_fuzzy_max_dist
        #raw_input()
        self.read_gff_as_interval_tree()
        self.group_info = MegaPBTree.read_group(self.group_filename, self.self_prefix) # ex: PB.1.1 --> [ RatHeart|i3_c123.... ]


    def read_gff_as_interval_tree(self):
        """
        Read a collapsed GFF file into an IntervalTree
        """
        for r in GFF.collapseGFFReader(self.gff_filename):
            self.tree[r.chr][r.strand].insert(r.start, r.end, r)

    @staticmethod
    def read_fastq_to_dict(fastq_filename):
        fastq_dict = {}
        for r in SeqIO.parse(open(fastq_filename), 'fastq'):
            fastq_dict[r.id.split('|')[0]] = r
        return fastq_dict

    @staticmethod
    def read_group(group_filename, group_prefix):
        group_info = {}
        with open(group_filename) as f:
            for line in f:
                pbid, members = line.strip().split('\t')
                if group_prefix is None:
                    group_info[pbid] = [x for x in members.split(',')]
                else:
                    group_info[pbid] = [group_prefix+'|'+x for x in members.split(',')]
        return group_info

    def match_record_to_tree(self, r):
        """
        r --- GMAPRecord
        tree --- dict of chromosome --> strand --> IntervalTree

        If exact match (every exon junction) or 5' truncated (allow_5merge is True), return the matching GMAPRecord
        Otherwise return None
        *NOTE*: the tree should be non-redundant so can return as soon as exact match is found!
        """
        matches = self.tree[r.chr][r.strand].find(r.start, r.end)
        for r2 in matches:
            r.segments = r.ref_exons
            r2.segments = r2.ref_exons
            if compare_junctions.compare_junctions(r, r2, internal_fuzzy_max_dist=self.internal_fuzzy_max_dist) == 'exact': # is a match!
                return r2
            elif self.allow_5merge: # check if the shorter one is a subset of the longer one
                if len(r.segments) > len(r2.segments):
                    a, b = r, r2
                else:
                    a, b = r2, r
                # a is the longer one, b is the shorter one
                if compare_junctions.compare_junctions(b, a, internal_fuzzy_max_dist=self.internal_fuzzy_max_dist) == 'subset':
                    # we only know that a is a subset of b, verify that it is actually 5' truncated (strand-sensitive!)
                    # if + strand, last exon of a should match last exon of b
                    # if - strand, first exon of a should match first exon of b
                    if (r.strand == '+' and compare_junctions.overlaps(a.segments[-1], b.segments[-1])) or \
                       (r.strand == '-' and compare_junctions.overlaps(a.segments[0], b.seq_exons[0])):
                        return r2

        return None

    def add_sample(self, gff_filename, group_filename, sample_prefix, output_prefix, fastq_filename=None):
        combined = [] # list of (r1 if r2 is None | r2 if r1 is None | longer of r1 or r2 if both not None)
        unmatched_recs = self.record_d.keys()

        for r in GFF.collapseGFFReader(gff_filename):
            match_rec = self.match_record_to_tree(r)
            if match_rec is not None:  # found a match! put longer of r1/r2 in
                combined.append((match_rec, r))
                try:
                    unmatched_recs.remove(match_rec.seqid)
                except ValueError:
                    pass # already deleted, OK, this happens for single-exon transcripts
            else:  # r is not present in current tree
                combined.append((None, r))
        # put whatever is left from the tree in
        for seqid in unmatched_recs:
            combined.append((self.record_d[seqid], None))

        # create a ClusterTree to re-calc the loci/transcripts
        final_tree = defaultdict(lambda: {'+': ClusterTree(0, 0), '-':ClusterTree(0, 0)})
        for i,(r1,r2) in enumerate(combined):
            if r2 is None or (r1 is not None and r1.end-r1.start > r2.end-r2.start):
                final_tree[r1.chr][r1.strand].insert(r1.start, r1.end, i)
            else:
                final_tree[r2.chr][r2.strand].insert(r2.start, r2.end, i)

        self.write_cluster_tree_as_gff(final_tree, combined, group_filename, sample_prefix, output_prefix, fastq_filename2=fastq_filename)


    def write_cluster_tree_as_gff(self, cluster_tree, rec_list, group_filename2, sample_prefix2, output_prefix, fastq_filename2=None):
        """
        Write ClusterTree (chr --> dict --> (start, end, rec_list_index)) as collapsedGFF format
        Returns --- a new group_info!!!
        """
        if fastq_filename2 is not None:
            fastq_dict2 = MegaPBTree.read_fastq_to_dict(fastq_filename2)
            f_fastq = open(output_prefix+'.rep.fq', 'w')
        group_info2 = MegaPBTree.read_group(group_filename2, sample_prefix2)
        new_group_info = {}
        f_out = open(output_prefix+'.gff', 'w')
        f_group = open(output_prefix+'.group.txt', 'w')
        f_mgroup = open(output_prefix + '.mega_info.txt', 'w')
        f_mgroup.write("pbid\t{0}\t{1}\n".format(self.self_prefix, sample_prefix2))
        loci_index = 0
        chroms = cluster_tree.keys()
        chroms.sort()
        for k in chroms:
            for strand in ('+', '-'):
                for _s, _e, rec_indices in cluster_tree[k][strand].getregions():
                    loci_index += 1
                    isoform_index = 0
                    for i in rec_indices:
                        isoform_index += 1
                        tID = "PB.{i}.{j}".format(i=loci_index,j=isoform_index)
                        r1, r2 = rec_list[i]
                        if r1 is None: # r2 is not None
                            r = r2
                            new_group_info[tID] = group_info2[r2.seqid]
                            f_mgroup.write("{tID}\tNA\t{group}\n".format(tID=tID, group=r2.seqid))
                            if fastq_filename2 is not None:
                                seqrec = fastq_dict2[r2.seqid]
                        elif r2 is None: # r1 is not None
                            r = r1
                            new_group_info[tID] = self.group_info[r1.seqid]
                            f_mgroup.write("{tID}\t{group}\tNA\n".format(tID=tID, group=r1.seqid))
                            if fastq_filename2 is not None:
                                seqrec = self.fastq_dict[r1.seqid]
                        else: # both r1, r2 are not empty
                            if (r1.end - r1.start > r2.end - r2.start):
                                r = r1
                                if fastq_filename2 is not None:
                                    seqrec = self.fastq_dict[r1.seqid]
                            else:
                                r = r2
                                if fastq_filename2 is not None:
                                    seqrec = fastq_dict2[r2.seqid]
                            new_group_info[tID] = self.group_info[r1.seqid] + group_info2[r2.seqid]
                            f_mgroup.write("{tID}\t{group1}\t{group2}\n".format(tID=tID, group1=r1.seqid, group2=r2.seqid))


                        if fastq_filename2 is not None:
                            seqrec.id = tID
                            SeqIO.write(seqrec, f_fastq, 'fastq')
                        f_group.write("{tID}\t{members}\n".format(tID=tID, members=",".join(new_group_info[tID])))
                        f_out.write("{chr}\tPacBio\ttranscript\t{s}\t{e}\t.\t{strand}\t.\tgene_id \"PB.{i}\"; transcript_id \"{tID}\";\n".format(\
                            chr=k, s=r.start+1, e=r.end, strand=strand, tID=tID, i=loci_index))
                        for exon in r.ref_exons:
                            f_out.write("{chr}\tPacBio\texon\t{s}\t{e}\t.\t{strand}\t.\tgene_id \"PB.{i}\"; transcript_id \"{tID}\";\n".format(\
                                chr=k, s=exon.start+1, e=exon.end, strand=strand, tID=tID, i=loci_index))
        f_out.close()
        f_group.close()
        f_mgroup.close()
        if fastq_filename2 is not None:
            f_fastq.close()
        return new_group_info


class MegaPBTreeFusion(MegaPBTree):

    def __init__(self, gff_filename, group_filename, internal_fuzzy_max_dist=0, self_prefix=None, fastq_filename=None, fusion_max_dist=10):
        """
        Differences with non-fusion MegaPBTree:

        1. allow_5merge is always FALSE. Not a parameter.
        2. fusion_max_dist --- maximum allowed distance on internal fusion sites to be called as equivalent fusions
        """
        super(MegaPBTreeFusion, self).__init__(gff_filename, group_filename, internal_fuzzy_max_dist, self_prefix, False, fastq_filename)

        self.fusion_max_dist = fusion_max_dist

        # ex: PBfusion.1 -> [PBfusion.1.1, PBfusion.1.2]
        self.record_d_fusion = dict((fusion_id, records) for fusion_id,records in GFF.collapseGFFFusionReader(gff_filename))

    def junction_match_check_5(self, r1, r2):
        if r1.strand == '+':
            return abs(r1.ref_exons[0].start-r2.ref_exons[0].start) <= self.fusion_max_dist
        else:
            return abs(r1.ref_exons[-1].end-r2.ref_exons[-1].end) <= self.fusion_max_dist

    def junction_match_check_3(self, r1, r2):
        if r1.strand == '+':
            return abs(r1.ref_exons[-1].end-r2.ref_exons[-1].end) <= self.fusion_max_dist
        else:
            return abs(r1.ref_exons[0].start-r2.ref_exons[0].start) <= self.fusion_max_dist

    def match_record_to_tree(self, r, check_5_dist, check_3_dist):
        """
        Matching a single record (locus).

        Major diff from non-fusion version:
        1. there could be multiple matches!
        2. no 5merge allowed
        3. additionally checks if the 5'/3' ends don't disagree too much (fusion_max_dist). this is used for fusion junctions.
        """
        matches = self.tree[r.chr][r.strand].find(r.start, r.end)
        result = []
        for r2 in matches:
            r.segments = r.ref_exons
            r2.segments = r2.ref_exons
            if compare_junctions.compare_junctions(r, r2, internal_fuzzy_max_dist=self.internal_fuzzy_max_dist) == 'exact' and \
                    (not check_5_dist or self.junction_match_check_5(r, r2)) and \
                    (not check_3_dist or self.junction_match_check_3(r, r2)): # is a match!
                result.append(r2.seqid)

        return result

    def check_records_match(self, records1, records2):
        """
        records1, records2 are two fusion records.
        They match iff:
        1. same number of records
        2. each record (a loci) matches
        """
        if len(records1)!=len(records2): return False

        i = 0
        for r1, r2 in zip(records1, records2):
            # check: chr, strand, exons match
            if r1.chr!=r2.chr or r1.strand!=r2.strand: return False
            r1.segments = r1.ref_exons
            r2.segments = r2.ref_exons
            if compare_junctions.compare_junctions(r1, r2, internal_fuzzy_max_dist=self.internal_fuzzy_max_dist)!='exact':
                return False
            if i == 0: # first record, only need 3' to agree
                if not self.junction_match_check_3(r1, r2): return False
            elif i == len(records1)-1: #last record, only need 5' to agree
                if not self.junction_match_check_5(r1, r2): return False
            else:
                if not self.junction_match_check_5(r1, r2): return False
                if not self.junction_match_check_3(r1, r2): return False
            i += 1

        return True

    def match_fusion_record(self, records):
        """
        records --- in order, the records of a single fusion.
        """
        good = []
        # match the first record, requiring additionally that the precise 3' end matches
        cands = self.match_record_to_tree(records[0], check_5_dist=False, check_3_dist=True)
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
            print >> sys.stderr, "ERROR! more than one possible candidate in match_fusion_record! DEBUG."
            print >> sys.stderr, "MATCHED:", good
            sys.exit(-1)


    def add_sample(self, gff_filename, group_filename, sample_prefix, output_prefix, fastq_filename=None):
        combined = [] # list of (r1 if r2 is None | r2 if r1 is None | longer of r1 or r2 if both not None)
        unmatched_recs = self.record_d_fusion.keys()

        for _id, records in GFF.collapseGFFFusionReader(gff_filename):
            match_seqid = self.match_fusion_record(records)
            if match_seqid is not None:
                combined.append((self.record_d_fusion[match_seqid], records))
                try:
                    unmatched_recs.remove(match_seqid)
                except ValueError:
                    pass # already deleted, OK, this happens for single-exon transcripts
            else:  # r is not present in current tree
                combined.append((None, records))
        # put whatever is left from the tree in
        for seqid in unmatched_recs:
            combined.append((self.record_d_fusion[seqid], None))

        #return combined

        # create a ClusterTree to re-calc the loci/transcripts
        final_tree = defaultdict(lambda: {'+': ClusterTree(0, 0), '-':ClusterTree(0, 0)})
        for i,(r1s,r2s) in enumerate(combined):
            if r2s is None or (r1s is not None and r1s[0].end-r1s[0].start > r2s[0].end-r2s[0].start):
                final_tree[r1s[0].chr][r1s[0].strand].insert(r1s[0].start, r1s[0].end, i)
            else:
                final_tree[r2s[0].chr][r2s[0].strand].insert(r2s[0].start, r2s[0].end, i)

        self.write_cluster_tree_as_gff(final_tree, combined, group_filename, sample_prefix, output_prefix, fastq_filename2=fastq_filename)


    def write_cluster_tree_as_gff(self, cluster_tree, rec_list, group_filename2, sample_prefix2, output_prefix, fastq_filename2=None):
        """
        Write ClusterTree (chr --> dict --> (start, end, rec_list_index)) as collapsedGFF format
        Returns --- a new group_info!!!
        """
        if fastq_filename2 is not None:
            fastq_dict2 = MegaPBTree.read_fastq_to_dict(fastq_filename2)
            f_fastq = open(output_prefix+'.rep.fq', 'w')
        group_info2 = MegaPBTree.read_group(group_filename2, sample_prefix2)
        new_group_info = {}
        f_out = open(output_prefix+'.gff', 'w')
        f_group = open(output_prefix+'.group.txt', 'w')
        f_mgroup = open(output_prefix + '.mega_info.txt', 'w')
        f_mgroup.write("pbid\t{0}\t{1}\n".format(self.self_prefix, sample_prefix2))
        fusion_index = 0
        chroms = cluster_tree.keys()
        chroms.sort()
        for k in chroms:
            for strand in ('+', '-'):
                for _s, _e, rec_indices in cluster_tree[k][strand].getregions():
                    for i in rec_indices:
                        fusion_index += 1
                        tID = "PBfusion.{i}".format(i=fusion_index)
                        r1s, r2s = rec_list[i]
                        if r1s is None: # r2s is not None
                            recs = r2s
                            r2_fusion_id = get_fusion_id(r2s[0].seqid)
                            new_group_info[tID] = group_info2[r2_fusion_id]
                            f_mgroup.write("{tID}\tNA\t{group}\n".format(tID=tID, group=r2_fusion_id))
                            if fastq_filename2 is not None:
                                seqrec = fastq_dict2[r2_fusion_id]
                        elif r2s is None: # r1 is not None
                            recs = r1s
                            r1_fusion_id = get_fusion_id(r1s[0].seqid)
                            new_group_info[tID] = self.group_info[r1_fusion_id]
                            f_mgroup.write("{tID}\t{group}\tNA\n".format(tID=tID, group=r1_fusion_id))
                            if fastq_filename2 is not None:
                                seqrec = self.fastq_dict[r1_fusion_id]
                        else: # both r1, r2 are not empty
                            r1_fusion_id = get_fusion_id(r1s[0].seqid)
                            r2_fusion_id = get_fusion_id(r2s[0].seqid)
                            r1_len = sum(x.end-x.start for x in r1s)
                            r2_len = sum(x.end-x.start for x in r2s)
                            if r1_len > r2_len:
                                recs = r1s
                                if fastq_filename2 is not None:
                                    seqrec = self.fastq_dict[r1_fusion_id]
                            else:
                                recs = r2s
                                if fastq_filename2 is not None:
                                    seqrec = fastq_dict2[r2_fusion_id]
                            new_group_info[tID] = self.group_info[r1_fusion_id] + group_info2[r2_fusion_id]
                            f_mgroup.write("{tID}\t{group1}\t{group2}\n".format(tID=tID, group1=r1_fusion_id, group2=r2_fusion_id))


                        if fastq_filename2 is not None:
                            seqrec.id = tID
                            SeqIO.write(seqrec, f_fastq, 'fastq')
                        f_group.write("{tID}\t{members}\n".format(tID=tID, members=",".join(new_group_info[tID])))

                        # now write out the fusion transcript
                        for j,r in enumerate(recs):
                            f_out.write("{chr}\tPacBio\ttranscript\t{s}\t{e}\t.\t{strand}\t.\tgene_id \"{gid}\"; transcript_id \"{gid}.{j}\";\n".format(\
                                chr=k, s=r.start+1, e=r.end, strand=strand, gid=tID, j=j+1))
                            for exon in r.ref_exons:
                                f_out.write("{chr}\tPacBio\texon\t{s}\t{e}\t.\t{strand}\t.\tgene_id \"{gid}\"; transcript_id \"{gid}.{j}\";\n".format(\
                                    chr=k, s=exon.start+1, e=exon.end, strand=strand, gid=tID, j=j+1))
        f_out.close()
        f_group.close()
        f_mgroup.close()
        if fastq_filename2 is not None:
            f_fastq.close()
        return new_group_info