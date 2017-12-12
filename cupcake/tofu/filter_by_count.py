#!/usr/bin/env python
__author__ = 'etseng@pacb.com'

import os, sys
from cupcake.io import GFF
from csv import DictReader, DictWriter
from Bio import SeqIO

"""
Given the collapse script result, further filter by FL counts.

Input: <prefix> (must have .group.txt, .abundance.txt, .gff, .rep.fq)
Output: <output>.min_fl_<threshold>.{group|abundance|gff|rep fq}
"""

def filter_by_count(input_prefix, output_prefix, min_count, dun_use_group_count=False):

    group_filename = input_prefix + '.group.txt'
    count_filename = input_prefix + '.abundance.txt'
    gff_filename = input_prefix + '.gff'
    rep_filename = input_prefix + '.rep.fq'

    # read group
    group_max_count_fl = {}
    group_max_count_p = {}
    f = open(group_filename)
    for line in f:
        #ex: PB.1.1  i0HQ_54b0ca|c58773/f30p16/700
        pbid, members = line.strip().split('\t')
        group_max_count_fl[pbid] = 0
        group_max_count_p[pbid] = 0
        members = members.split(',')
        for m in members:
            i = m.find('|')
            if i > 0:
                tmp = m.split('|')[1].split('/')[1] #ex: tmp = f30p16
            else:
                tmp = m.split('/')[1]
            fl_count, p_count = tmp.split('p')
            fl_count = int(fl_count[1:])
            p_count = int(p_count)
            group_max_count_fl[pbid] = max(group_max_count_fl[pbid], fl_count)
            group_max_count_p[pbid] = max(group_max_count_p[pbid], p_count)
    f.close()

    # read abundance first
    f = open(count_filename)
    count_header = ''
    while True:
        cur_pos = f.tell()
        line = f.readline()
        if not line.startswith('#'):
            f.seek(cur_pos)
            break
        else:
            count_header += line
    d = dict((r['pbid'], r) for r in DictReader(f, delimiter='\t'))
    for k,v in d.iteritems():
        print k,v
    f.close()

    # group_max_count_p NOT used for now
    good = filter(lambda x: int(d[x]['count_fl']) >= min_count and (dun_use_group_count or group_max_count_fl[x] >= min_count) and group_max_count_p >= 0, d)

    # write output GFF
    f = open(output_prefix + '.gff', 'w')
    for r in GFF.collapseGFFReader(gff_filename):
        if r.seqid in good: GFF.write_collapseGFF_format(f, r)
    f.close()


    # write output rep.fq
    f = open(output_prefix + '.rep.fq', 'w')
    for r in SeqIO.parse(open(rep_filename), 'fastq'):
        if r.name.split('|')[0] in good:
           SeqIO.write(r, f, 'fastq')
    f.close()

    # write output to .abundance.txt
    f = open(output_prefix + '.abundance.txt', 'w')
    f.write(count_header)
    writer = DictWriter(f, fieldnames=['pbid','count_fl','count_nfl','count_nfl_amb','norm_fl','norm_nfl','norm_nfl_amb'], \
                        delimiter='\t', lineterminator='\n')
    writer.writeheader()
    for k in good:
        r = d[k]
        writer.writerow(r)
    f.close()

    print >> sys.stderr, "Output written to:", output_prefix + '.gff'
    print >> sys.stderr, "Output written to:", output_prefix + '.rep.fq'
    print >> sys.stderr, "Output written to:", output_prefix + '.abundance.txt'



if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("input_prefix")
    parser.add_argument("--min_count", type=int, default=2, help="Minimum FL count (default: 2)")
    parser.add_argument("--dun_use_group_count", action="store_true", default=False, help="Turn off more stringent min count (default: off)")

    args = parser.parse_args()
    output_prefix = "{i}.min_fl_{c}".format(i=args.input_prefix, c=args.min_count)
    filter_by_count(args.input_prefix, output_prefix, args.min_count, args.dun_use_group_count)