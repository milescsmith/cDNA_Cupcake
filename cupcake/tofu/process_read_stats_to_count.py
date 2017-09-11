__author__ = 'etseng@pacb.com'


FIELD_EXPLAIN=\
"""#
# -----------------
# Field explanation
# -----------------
# count_fl: Number of associated FL reads
# count_nfl: Number of associated FL + unique nFL reads
# count_nfl_amb: Number of associated FL + unique nFL + weighted ambiguous nFL reads
# norm_fl: count_fl / total number of FL reads
# norm_nfl: count_nfl / total number of FL + unique nFL reads
# norm_nfl_amb: count_nfl_amb / total number of all reads
"""
# Total Number of FL reads: 28930
# Total Number of FL + unique nFL reads: 55937
# Total Number of all reads: 58666
#
#pbid    count_fl        count_nfl       count_nfl_amb   norm_fl norm_nfl        norm_nfl_amb

import os, sys
from collections import defaultdict
from csv import DictReader
from Bio import SeqIO


class ReadTally:
    def __init__(self):
        self.num_fl = 0  # FL is always unique
        self.num_nfl_unique = 0
        self.num_nfl_amb = 0
        self.tally = defaultdict(lambda: {'fl_only': 0, 'nfl_only': 0, 'nfl_amb_only': 0})
        self.nfl_amb_hits = defaultdict(lambda: [])

def gather_read_stat(read_stat_filename):
    o = ReadTally()
    for r in DictReader(open(read_stat_filename), delimiter='\t'):
        if r['id']=='id': continue # skip headers which probably is there before files were concated
        if r['is_fl']=='Y':
            o.num_fl += 1
            if r['stat']!='NA':
                o.tally[r['pbid']]['fl_only'] += 1
        else: # is nFL
            if r['stat'] == 'NA': # nFL unmapped
                o.num_nfl_amb += 1
            elif r['stat'] == 'unique':
                o.num_nfl_unique += 1
                o.tally[r['pbid']]['nfl_only'] += 1
            else: # not unique
                # for ambiguous, we will only know the weight after processing whole file
                assert r['stat'] == 'ambiguous'
                o.nfl_amb_hits[r['id']].append(r['pbid'])
    # now we can add the ambiguous stuff back
    for _id, _hits in o.nfl_amb_hits.iteritems():
        w = 1. / len(_hits)
        for pbid in _hits:
            o.tally[pbid]['nfl_amb_only'] += w
        o.num_nfl_amb += 1
    return o


def tally_count_file(fasta_filename, tally_obj, output_filename, is_pbid):
    # now properly form
    # nfl = fl_only + nfl_only
    # nfl_amb = fl_only + nfl_only + nfl_amb_only

    if is_pbid:
        allids = [r.id.split('|')[0] for r in SeqIO.parse(open(fasta_filename), 'fasta')]
    else:
        allids = [r.id for r in SeqIO.parse(open(fasta_filename), 'fasta')]

    o = tally_obj
    total_fl = o.num_fl
    total_nfl = total_fl + o.num_nfl_unique
    total_nfl_amb = total_nfl + o.num_nfl_amb
    with open(output_filename, 'w') as f:
        f.write(FIELD_EXPLAIN)
        f.write("# Total Number of FL reads: {0}\n".format(total_fl))
        f.write("# Total Number of FL + unique nFL reads: {0}\n".format(total_nfl))
        f.write("# Total Number of all reads: {0}\n".format(total_nfl_amb))
        f.write("#\n")
        f.write("pbid\tcount_fl\tcount_nfl\tcount_nfl_amb\tnorm_fl\tnorm_nfl\tnorm_nfl_amb\n")

        allids.sort(key=lambda x: map(int, x.split('.')[1:]))


        for pbid in allids:
            print "processing", pbid
            if pbid in o.tally:
                f.write("{0}\t0\t0\t0\t0\t0\t0\n".format(pbid))
            else:
                x = o.tally[pbid]
                count_fl = x['fl_only']
                count_nfl = x['fl_only'] + x['nfl_only']
                count_nfl_amb = x['fl_only'] + x['nfl_only'] + x['nfl_amb_only']
                f.write("{0}\t".format(pbid))
                f.write("{0}\t{1}\t{2}\t".format(count_fl, count_nfl, count_nfl_amb))
                f.write("{0}\t{1}\t{2}\n".format(count_fl*1./total_fl,\
                                             count_nfl*1./total_nfl,\
                                             count_nfl_amb*1./total_nfl_amb))

def main(fasta_filename, read_stat_filename, output_filename, is_pbid):
    print read_stat_filename, output_filename
    o = gather_read_stat(read_stat_filename)
    tally_count_file(fasta_filename, o, output_filename, is_pbid)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("fasta_filename", help="The reference fasta filename (ex: ref.fasta)")
    parser.add_argument("read_stat_filename", help="The processed read stat filename (ex: test.blasr.txt)")
    parser.add_argument("output_filename", help="Output count filename (ex: test.count.txt")
    parser.add_argument("--is_pbid", default=False, action="store_true", help="Ref fasta has PBIDs (ex: PB.1.2)")

    args = parser.parse_args()
    main(args.fasta_filename, args.read_stat_filename, args.output_filename, args.is_pbid)

