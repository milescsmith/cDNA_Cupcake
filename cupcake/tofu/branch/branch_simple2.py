import os, sys
import numpy as np

import cupcake.io.BioReaders as BioReaders
import cupcake.tofu.branch.c_branch as c_branch
#from cupcake.tofu.branch.intersection_unique import IntervalTreeUnique, Interval, IntervalNodeUnique

from Bio import SeqIO

INTINF = 999999

class ContiVec:
    """
    Original struct: 'BC'

    """
    def __init__(self, size):
        self.baseC = np.zeros(size, dtype=np.int) # this was .B in the original code, base coverage
        self.altC_pos = np.zeros(size, dtype=np.int)  # this was .C in the original code, evidence for alternative junction
        self.altC_neg = np.zeros(size, dtype=np.int)
        self.transC = np.zeros(size, dtype=np.int)  # this is .isTrans in the original code

class BranchSimple:
    """
    BranchSimple is designed for just creating exons from PacBio's GMAP results
    Does not use Illumina
    """
    def __init__(self, transfrag_filename, cov_threshold=2, min_aln_coverage=.99, min_aln_identity=.85, is_fq=False):
        self.contiVec = None # current ContiVec object
        self.exons = None
        #self.MIN_EXON_SIZE = max_fuzzy_junction

        self.transfrag_filename = transfrag_filename
        self.transfrag_len_dict = dict((r.id.split()[0], len(r.seq)) for r in SeqIO.parse(open(transfrag_filename), 'fastq' if is_fq else 'fasta'))

        self.cov_threshold = cov_threshold # only output GTF records if >= this many GMAP records support it (this must be if I'm running non-clustered fasta on GMAP)

        self.min_aln_coverage = min_aln_coverage
        self.min_aln_identity = min_aln_identity

        self.cuff_index = 1


    def iter_gmap_sam(self, gmap_sam_filename, ignored_fout):
        """
        Iterate over a SORTED GMAP SAM file.
        Return a collection of records that overlap by at least 1 base.
        """
        def sep_by_strand(records):
            output = {'+':[], '-':[]}
            for r in records:
                output[r.flag.strand].append(r)
            return output

        records = None # holds the current set of records that overlap in coordinates
        iter = BioReaders.GMAPSAMReader(gmap_sam_filename, True, query_len_dict=self.transfrag_len_dict)
        for r in iter:
            if r.sID == '*':
                ignored_fout.write("{0}\tUnmapped.\n".format(r.qID))
            elif r.qCoverage < self.min_aln_coverage:
                ignored_fout.write("{0}\tCoverage {1:.3f} too low.\n".format(r.qID, r.qCoverage))
            elif r.identity < self.min_aln_identity:
                ignored_fout.write("{0}\tIdentity {1:.3f} too low.\n".format(r.qID, r.identity))
            else:
                break
        try:
            records = [r]
        except NameError:
            print >> sys.stderr, "No valid records from {0}!".format(gmap_sam_filename)
            return
        for r in iter:
            if r.sID == records[0].sID and r.sStart < records[-1].sStart:
                print >> sys.stderr, "SAM file is NOT sorted. ABORT!"
                sys.exit(-1)
            if r.qCoverage < self.min_aln_coverage:
                ignored_fout.write("{0}\tCoverage {1:.3f} too low.\n".format(r.qID, r.qCoverage))
            elif r.identity < self.min_aln_identity:
                ignored_fout.write("{0}\tIdentity {1:.3f} too low.\n".format(r.qID, r.identity))
            elif r.sID != records[0].sID or r.sStart > max(x.sEnd for x in records):
                yield sep_by_strand(records)
                records = [r]
            else:
                records.append(r)
        yield sep_by_strand(records)


    def parse_transfrag2contig(self, gmap_sam_records, skip_5_exon_alt=True):
        """
        GMAP SAM file MUST BE SORTED! (same criterion as cufflinks)

        Goes through a set of overlapping GMAP records (strand-aware)
        Fill in the baseC (base coverage)
        """
        records = gmap_sam_records
        # first figure out how long the "pseudo-chromosome" size is
        offset = records[0].sStart
        self.offset = offset
        self.chrom = records[0].sID
        self.strand = records[0].flag.strand
        # set of this fake "chromosome" to be as long as needed
        chrom_size = max(x.sEnd for x in records) - records[0].sStart
        self.contiVec = ContiVec(chrom_size)
        for r in records:
            for i,e in enumerate(r.segments):
                # fill base coverage
                self.contiVec.baseC[(e.start-offset):(e.end-offset)] += 1

                # in the original code, the mapped start altC was set to -MAX and end to MAX
                # add this alt. of beginning if 
                # (a) not first or last exon
                # (b) is first exon and strand - (so is 3')
                # (c) is first exon and strand + and not skip_5_exon_alt (so is 5')
                # (d) is last exon and strand - and not skip_5_exon_alt (so is 5')
                # (e) is last exon and strand + (so is 3')
                if (i != 0) or (i != len(r.segments)-1) or \
                    (i==0 and (self.strand=='-' or not skip_5_exon_alt)) or \
                    (i==len(r.segments)-1 and (self.strand=='+' or not skip_5_exon_alt)): 
                    self.contiVec.altC_neg[(e.start-offset)] -= INTINF
                # add this alt. of end if
                # (a) not first or last exon
                # (b) is last exon and strand + (so is 3')
                # (c) is last exon and strand - and not skip_5_exon_alt (so is 5')
                # (d) is first exon and strand + and not skip_5_exon_alt (so is 5')
                # (e) is first exon and strand - (so is 3')
                if (i != 0) or (i != len(r.segments)-1) or \
                    (i==len(r.segments)-1 and (self.strand=='-' or not skip_5_exon_alt)) or \
                    (i==0 and (self.strand=='+' or not skip_5_exon_alt)):
                    self.contiVec.altC_pos[(e.end-offset-1)] += INTINF  # adjust to 0-based
                # set .transC
                self.contiVec.transC[(e.start-offset):(e.end-offset)] += 1


    def exon_finding(self):
        """
        Go through contiVec to idfentify the exons using base coverage (.baseC) and alt junction evidence (.altC)
        """
        v = self.contiVec
        self.exons = c_branch.exon_finding(v.baseC, v.altC_neg, v.altC_pos, v.transC, \
                            len(v.baseC), 2, 0, offset=self.offset)


    def match_record(self, r, tolerate_middle=0, tolerate_end=10**4, ok_to_miss_matches=False, intervals_adjacent=True):
        """
        r --- a gmap record (GFF.gmapRecord) containing the full contig alignment of a FL long read
        """
        result = []
        num_exons = len(r.segments)
        for i,e in enumerate(r.segments):
            # allow the first and last exon to be longer or shorter in either the long read or the detected exons
            # however for middle exons the matching should be pretty precise
            # ToDo: parameterize this
            if i == 0 and num_exons > 1:
                tolerate_l = tolerate_end
                tolerate_r = tolerate_middle
            elif i == 0 and num_exons == 1:
                tolerate_l = tolerate_end
                tolerate_r = tolerate_end
            elif i == len(r.segments)-1 and num_exons > 1:
                tolerate_l = tolerate_middle
                tolerate_r = tolerate_end
            else:
                tolerate_l = tolerate_middle
                tolerate_r = tolerate_middle
            matches = exon_matching(self.exons, e, tolerate_l, tolerate_r, intervals_adjacent) # ToDo: CHANGE this back to c_branch later
            if matches is None:
                if not ok_to_miss_matches:
                    return None
            else:
                if (len(result) >= 1 and result[-1].value >= matches[0].value) and (not ok_to_miss_matches):
                    return None
                result += matches
        return result


##### THIS DOES NOT REALLY WORK.....use fuzzy collapse instead :((((
    # def correct_micro_exons(self):
    #     """
    #     ToDo: this should eventually be moved into c_branch.pyx, leaving here now so I don't have to compile Cython
    #           during heavy debugging phase.
    #
    #     The point here is after self.exon_finding(), correct those "micro exons" that are only 1 or 2 bps
    #     ex: exon 1 is 1 - 2
    #         exon 2 is 2-100
    #
    #     most likely this is a GMAP artifact. the answer is either 1-100 or 2-100
    #     there are two ways to correct it:
    #     <a> use the one that has higher baseC
    #     <b> if genome available, use the one that has canonical splice site (IMPLEMENT LATER)
    #
    #     Steps:
    #     --- identify all micro exons that are (i) X bp or shorter and (ii) adjacent to another exon (iii) directionality is the same
    #     --- use criterion above to decide what correct exon is
    #     """
    #     def correct_exon(_prev, _next):
    #         """
    #         note: _prev is always the microexon and _next is the downstream one
    #         ex: _prev is 1 - 2
    #             _next is 2 - 4
    #
    #         """
    #         print "baseC:", self.contiVec.baseC[_prev.start-self.offset-5:_prev.end-self.offset+5]
    #         print "altC pos:", self.contiVec.altC_pos[_prev.start-self.offset-5:_prev.end-self.offset+5]
    #         print "altC neg:", self.contiVec.altC_neg[_prev.start-self.offset-5:_prev.end-self.offset+5]
    #         if self.contiVec.baseC[_prev.start-self.offset] >= self.contiVec.baseC[_prev.end-self.offset]:
    #                 # prev has more evidence
    #             x = Interval(_prev.start, _next.end, _prev.interval.value)
    #             return IntervalNodeUnique(x.start, x.end, x)
    #         else: # next has more evidence
    #             x = Interval(_next.start, _next.end, _prev.interval.value)
    #             return IntervalNodeUnique(x.start, x.end, x)
    #
    #
    #     p = []
    #     self.exons.traverse(p.append)
    #
    #     i = 0
    #     while i < len(p) - 1:
    #         x = p[i]
    #         if x.end - x.start < self.MIN_EXON_SIZE:
    #             if x.end == p[i+1].start: # adjacent w/ downstream
    #                 new_e = correct_exon(x, p[i+1])
    #                 if new_e is not None:
    #                     print >> sys.stderr, "merging {0}, {1} --> {2}".format(x, p[i+1], new_e)
    #                     raw_input()
    #                     p.pop(i)
    #                     p[i] = new_e
    #                 else: i += 1 # nothing to do becuz directionality diff, advance
    #             else:
    #                 i += 1 # nothing to do, advance
    #         else:
    #             i += 1 # nothing to do, advance

    def process_records(self, records, allow_extra_5_exons, skip_5_exon_alt, f_good, f_bad, f_group, tolerate_end=100, starting_isoform_index=0, gene_prefix='PB'):
        """
        Given a set of records
        (1) process them by running through parse_transfrag2contig
        (2) call exons by exon_finding
        (3) go through each record, get the list of "nodes" they corresspond to
        (4) collapse identical records (53mergeing)

        Write out to GTF format
        """
        self.parse_transfrag2contig(records, skip_5_exon_alt)
        self.exon_finding()

        #self.correct_micro_exons() # DOES NOT WORK....

        p = []
        self.exons.traverse(p.append)
        node_d = dict((x.interval.value, x) for x in p)
        mat_size = max(x.interval.value for x in p) + 1
        result = []
        for r in records:
            stuff = self.match_record(r, tolerate_end=tolerate_end)#, tolerate_middle=self.MIN_EXON_SIZE)
            m = np.zeros((1, mat_size), dtype=np.int)
            for x in stuff: m[0, x.value]=1
            result.append((r.qID, r.flag.strand, m))

        result_merged = list(result)
        iterative_merge_transcripts(result_merged, node_d, allow_extra_5_exons)

        print >> sys.stderr, "merged {0} down to {1} transcripts".format(len(result), len(result_merged))

        self.isoform_index = starting_isoform_index
        # make the exon value --> interval dictionary
        a = []

        for ids, strand, m in result_merged:
            assert self.strand==strand
            if ids.count(',')+1 < self.cov_threshold:
                f_out = f_bad
            else:
                f_out = f_good
            self.isoform_index += 1
            segments = [node_d[x] for x in m.nonzero()[1]]
            f_group.write("{p}.{i}.{j}\t{ids}\n".format(ids=ids, i=self.cuff_index, j=self.isoform_index, p=gene_prefix))
            f_out.write("{chr}\tPacBio\ttranscript\t{s}\t{e}\t.\t{strand}\t.\tgene_id \"{p}.{i}\"; transcript_id \"{p}.{i}.{j}\";\n".format(\
                chr=self.chrom, s=segments[0].start+1, e=segments[-1].end, i=self.cuff_index, p=gene_prefix, j=self.isoform_index, strand=self.strand))

            i = 0
            j = 0
            for j in xrange(1, len(segments)):
                if segments[j].start != segments[j-1].end:
                    f_out.write("{chr}\tPacBio\texon\t{s}\t{e}\t.\t{strand}\t.\tgene_id \"{p}.{i}\"; transcript_id \"{p}.{i}.{j}\";\n".format(\
                    chr=self.chrom, s=segments[i].start+1, e=segments[j-1].end, p=gene_prefix,i=self.cuff_index, j=self.isoform_index, strand=self.strand))
                    i = j
            f_out.write("{chr}\tPacBio\texon\t{s}\t{e}\t.\t{strand}\t.\tgene_id \"{p}.{i}\"; transcript_id \"{p}.{i}.{j}\";\n".format(\
                    chr=self.chrom, s=segments[i].start+1, e=segments[j].end, p=gene_prefix, i=self.cuff_index, j=self.isoform_index, strand=self.strand))

        self.cuff_index += 1

        return result, result_merged

def iterative_merge_transcripts(result_list, node_d, merge5=True):
    """
    result_list --- list of (qID, strand, binary exon sparse matrix)
    """
    # sort by strand then starting position
    result_list.sort(key=lambda x: (x[1], x[2].nonzero()[1][0]))
    i = 0
    while i < len(result_list) - 1:
        j = i + 1
        while j < len(result_list):
            id1, strand1, m1 = result_list[i]
            id2, strand2, m2 = result_list[j]
            if (strand1 != strand2) or (m1.nonzero()[1][-1] < m2.nonzero()[1][0]):
                break
            else:
                flag, m3 = compare_exon_matrix(m1, m2, node_d, strand1, merge5)
                if flag:
                    result_list[i] = (id1+','+id2, strand1, m3)
                    result_list.pop(j)
                else:
                    j += 1
        i += 1

        
def compare_exon_matrix(m1, m2, node_d, strand, merge5=True):
    """
    m1, m2 are 1-d array where m1[0, i] is 1 if it uses the i-th exon, otherwise 0
    compare the two and merge them if they are compatible
    (i.e. only differ by first/last exon ends)

    merge5 -- if True, allow extra 5' exons as long as the rest is the same
              if False, then m1 and m2 must have the same first (5') exon and only allowed
                        if the difference is the very start

    return {True|False}, {merged array|None}
    """
    l1 = m1.nonzero()[1]
    l2 = m2.nonzero()[1]


    # let l1 be the one that has the earliest start
    if l2[0] < l1[0]: l1, l2 = l2, l1

    # does not intersect at all
    if l1[-1] < l2[0]: return False, None

    n1 = len(l1)
    n2 = len(l2)

    #  not ok if this is at the 3' end and does not share the last exon; if 5' end, ok to miss exons as long as rest agrees
    for i in xrange(n1):
        if l1[i] == l2[0]: break
        elif i > 0 and (strand=='-' and node_d[l1[i-1]].end!=node_d[l1[i]].start): return False, None # 3' end disagree
        elif i > 0 and (strand=='+' and not merge5 and node_d[l1[i-1]].end!=node_d[l1[i]].start): return False, None # 5' end disagree, in other words m1 has an extra 5' exon that m2 does not have and merge5 is no allowed
    # at this point: l1[i] == l2[0]

    for j in xrange(i, min(n1, n2+i)):
        # matching l1[j] with l2[j-i]
        if l1[j] != l2[j-i]: # they must not match
            return False, None

    # pre: l1 and l2 agree up to j, j-i
    if j == n1-1: # check that the remaining of l2 are adjacent
        if j-i == n2-1:
            return True, m1
        for k in xrange(j-i+1, n2):
            # case 1: this is the 3' end, check that there are no additional 3' exons
            if (strand=='+' and node_d[l2[k-1]].end!=node_d[l2[k]].start): return False, None
            # case 2: this is the 5' end, check that there are no additional 5' exons unless allowed
            if (strand=='-' and not merge5 and node_d[l2[k-1]].end!=node_d[l2[k]].start): return False, None
        m1[0, l2[j-i+1]:] = m1[0, l2[j-i+1]:] + m2[0, l2[j-i+1]:]
        return True, m1
    elif j-i == n2-1:
        for k in xrange(j+1, n1):
            # case 1, but for m1
            if (strand=='+' and node_d[l1[k-1]].end!=node_d[l1[k]].start): return False, None
            # case 2, but for m1
            if (strand=='-' and not merge5 and node_d[l1[k-1]].end!=node_d[l1[k]].start): return False, None
        return True, m1

    raise Exception, "Should not happen"


def trim_exon_left_to_right(m1, m2, node_d, max_distance):
    """

    """
    l1 = m1.nonzero()[1]
    l2 = m2.nonzero()[1]


    # let l1 be the one that has the earliest start
    if l2[0] < l1[0]: l1, l2 = l2, l1

    # does not intersect at all
    if l1[-1] < l2[0]: return False, None

    n1 = len(l1)
    n2 = len(l2)

    # find l1[i] == l2[0]
    # meanwhile, REJECT if either (a) a junction is encountered already; (b) the distance exceeds threshold
    for i in xrange(n1):
        if l1[i] == l2[0]: break
        elif i > 0 and node_d[l1[i-1]].end!=node_d[l1[i]].start: return False, None  # condition (a)
        elif node_d[l1[i]].end - node_d[l1[0]].start > max_distance: return False, None # condition (b)

    # walk down l1, l2 together until first junction is met
    for j in xrange(i, n1):
        if j-i > n2-1: return False, None # junction not yet encountered yet no more m2 (must be one-exon only), REJECT
        # matching l1[j] with l2[j-i]
        if l1[j] != l2[j-i]: # they must not match
            return False, None
        if node_d[l1[j]].start != node_d[l1[j-1]].end: # junction encountered
            break

    if node_d[l1[j]].start == node_d[l1[j-1]].end: return False, None # m1 is one-exon only

    # pre: l1 and l2 agree up to j, j-i
    #      l1[j-1], l1[j] is an exon junction
    # post: trim l1 down to begin at same place as l2
    m1[0, :l1[i]] = 0
    return True, m1


def exon_matching(exon_tree, ref_exon, match_extend_tolerate_left, match_extend_tolerate_right, intervals_adjacent=True):
    """
    exon_tree --- an IntervalTree made from .baseC/.altC using exon detection; probably only short read data
    ref_exon --- an Interval representing an exon; probably from PacBio
    match_extend_tolerate --- maximum difference between the matched start/end

    find a continuous exon path (consisting of 1 or more nodes for which the intervals must be adjacent)
    in exon_tree that matches to ref_exon
    """
    matches = exon_tree.find(ref_exon.start, ref_exon.end)
    if len(matches) == 0: # likely due to very low coverage on transcript
        return None
    # check that all matches are adjacent (no splicing! this just one integral exon)
    if (not intervals_adjacent) or c_branch.intervals_all_adjacent(matches):
        # check if the ends differ a little, if so, extend to min/max
        for i in xrange(len(matches)):
            d_start = abs(matches[i].start - ref_exon.start)
            #print "matching {0} to {1}".format(matches[i].start, ref_exon.start)
            #pdb.set_trace()
            if d_start <= match_extend_tolerate_left: # now find the furthest end that satisfies the results
                for j in xrange(len(matches)-1, i-1, -1):
                    if abs(matches[j].end - ref_exon.end) <= match_extend_tolerate_right:
                        return matches[i:(j+1)]
        return None
    else: # ack! could not find evidence for this :<
        return None



                
                
    
