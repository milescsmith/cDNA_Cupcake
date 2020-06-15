import os, sys, pdb, subprocess, glob
from pickle import dump
from Bio import SeqIO

from cdna_cupcake.cupcake2.io.minimapIO import MiniReader
from cdna_cupcake.cupcake.io.SeqReaders import LazyFastaReader
from cdna_cupcake.cupcake2.tofu2.create_seed import create_seed_n_batch_files
import cupcake2.ice2.preClusterProcess as sp
import cupcake2.ice2.AlignerRunners as ar
import cupcake2.io.FileIO as FileIO
from collections import defaultdict


def detect_PCR_chimeras(orphans, fasta_d, window_size=20, min_T_count=16):
    """
    Remove from orphans, all sequences that have stretch of polyTs in the front
    Return: candidates that fit the PCR chimera criteria
    """
    result = []
    for _id in orphans:
        r = fasta_d[_id]
        if r.seq[:window_size].count("T") >= min_T_count:
            result.append(r.id)
    return result


def cleanup_precluster_intermediate_files(batch_index):
    """
    clean up seed<i>.S.fasta, seed<i>.orphans.fasta,
    batch<i>.fasta, batch<i>.remains.fasta, batch<i>.remains2.fasta
    """
    files = [
        "seed{}.S.fasta".format(batch_index),
        "seed{}.orphans.fasta".format(batch_index),
        "batch{}.fasta".format(batch_index),
        "batch{}.remains.fasta".format(batch_index),
        "batch{}.remains2.fasta".format(batch_index),
    ]

    files += glob.glob("batch{}*.minimap".format(batch_index))
    for file in files:
        try:
            os.remove(file)
        except:
            print("Failure to remove {}. Ignore.".format(file), file=sys.stderr)


def add_batch(batch_index, pCS, orphans, fasta_d, cpus, dun_use_partial):
    """
    1. align batch<i>.fasta against seed<i>.S.fasta, process -> write remains to batch<i>.remains.fasta
    2. align batch<i>.remains.fasta against seed<i>.orphans.fasta -> write remains to batch<i>.remains2.fasta
    3. self align batch<i>.remains2.fasta -> combine remains+orphans to new orphans
    4. write out seed<i+1>.S.fasta and seed<i+1>.orphans.fasta

    """
    cur_file = "batch{}.fasta".format(batch_index)
    seqids = {r.id for r in SeqIO.parse(open(cur_file), "fasta")}
    o = ar.run_minimap(cur_file, "seed{}.S.fasta".format(batch_index), cpus=cpus)
    print("processing", o, file=sys.stderr)
    pCS, remains = sp.process_align_to_pCS(
        o, seqids, pCS, MiniReader, dun_use_partial=dun_use_partial
    )
    print(
        "pCS: {}, tucked: {}, orphans: {}, remains: {}".format(
            len(pCS.S),
            sum(v == "T" for v in pCS.seq_stat.values()),
            len(orphans),
            len(remains),
        ),
        file=sys.stderr,
    )
    # write batch<i>.remains.fasta
    cur_file = "batch{}.remains.fasta".format(batch_index)
    FileIO.write_seqids_to_fasta(remains, cur_file, fasta_d)
    o = ar.run_minimap(cur_file, "seed{}.orphans.fasta".format(batch_index), cpus=cpus)
    print("processing", o, file=sys.stderr)
    pCS, orphans, remains = sp.process_align_to_orphan(
        o, remains, orphans, pCS, MiniReader, dun_use_partial=dun_use_partial
    )
    print(
        "pCS: {}, tucked: {}, orphans: {}, remains: {}".format(
            len(pCS.S),
            sum(v == "T" for v in pCS.seq_stat.values()),
            len(orphans),
            len(remains),
        ),
        file=sys.stderr,
    )
    # write batch<i>.remains2.fasta and self align
    cur_file = "batch{}.remains2.fasta".format(batch_index)
    FileIO.write_seqids_to_fasta(remains, cur_file, fasta_d)
    o = ar.run_minimap(cur_file, cur_file, cpus=cpus)
    print("processing", o, file=sys.stderr)
    pCS, remains = sp.process_self_align_into_seed(
        o, remains, MiniReader, pCS, dun_use_partial=dun_use_partial
    )
    print(
        "pCS: {}, tucked: {}, orphans: {}, remains: {}".format(
            len(pCS.S),
            sum(v == "T" for v in pCS.seq_stat.values()),
            len(orphans),
            len(remains),
        ),
        file=sys.stderr,
    )
    # combine remains+orphans to new orphans
    orphans = orphans.union(remains)
    FileIO.write_preClusterSet_to_fasta(
        pCS, "seed{}.S.fasta".format(batch_index + 1), fasta_d
    )
    FileIO.write_seqids_to_fasta(
        orphans, "seed{}.orphans.fasta".format(batch_index + 1), fasta_d
    )

    return pCS, orphans


def main(
    cpus,
    dun_make_bins=False,
    dun_use_partial=False,
    num_seqs_per_batch=100000,
    dun_cleanup_files=False,
):
    print("Indexing isoseq_flnc.fasta using LazyFastaReader...")
    d = LazyFastaReader("isoseq_flnc.fasta")

    print("Splitting input isoseq_flnc.fasta into seed/batches...")
    num_batchs = create_seed_n_batch_files(
        input="isoseq_flnc.fasta",
        fasta_d=d,
        seed_filename="seed0.fasta",
        batch_pre="batch",
        num_seqs_per_batch=num_seqs_per_batch,
    )

    # step1. run minimap of seed0 against itself and process
    o = ar.run_minimap("seed0.fasta", "seed0.fasta", cpus=cpus)
    seqids = {r.id for r in SeqIO.parse(open("seed0.fasta"), "fasta")}
    pCS, orphans = sp.process_self_align_into_seed(
        o, seqids, MiniReader, dun_use_partial=dun_use_partial
    )
    # keep stats
    size_S, size_tucked, size_orphans = (
        len(pCS.S),
        sum(v == "T" for v in pCS.seq_stat.values()),
        len(orphans),
    )
    print(
        "seed 0 initial: S {}, tucked {}, orphans {}".format(
            size_S, size_tucked, size_orphans
        )
    )

    # write out seed1.S.fasta and seed1.orphans.fasta
    FileIO.write_preClusterSet_to_fasta(pCS, "seed1.S.fasta", d)
    FileIO.write_seqids_to_fasta(orphans, "seed1.orphans.fasta", d)
    # step 2a. minimap batch1 against seed1.S and process

    for i in range(1, num_batchs):
        pCS, orphans = add_batch(
            i, pCS, orphans, d, cpus=cpus, dun_use_partial=dun_use_partial
        )
        cleanup_precluster_intermediate_files(i)

    # detect PCR chimeras from orphans
    chimeras = detect_PCR_chimeras(orphans, d)
    orphans = orphans.difference(chimeras)

    FileIO.write_seqids_to_fasta(orphans, "preCluster_out.orphans.fasta", d)
    FileIO.write_seqids_to_fasta(chimeras, "preCluster_out.chimeras.fasta", d)

    tucked_seqids = []
    # dump pCS, orphans, chimeras to a pickle
    # can't dump yet --- since pCS is an object
    # with open('preCluster.output.pickle', 'w') as f:
    #    dump({'pCS': pCS, 'chimeras': chimeras, 'orphans': orphans}, f)
    # write CSV file
    with open("preCluster.output.csv", "w") as f:
        f.write("seqid,stat\n")
        for x, stat in pCS.seq_stat.items():
            if stat == "T":
                f.write("{},tucked\n".format(x))
                tucked_seqids.append(x)
            elif stat == "M":
                f.write("{},{}\n".format(x, pCS.seq_map[x]))
        for x in orphans:
            f.write("{},orphan\n".format(x))
        for x in chimeras:
            f.write("{},chimera\n".format(x))

    # Liz: currently not using tucked...
    if len(tucked_seqids) > 0:
        FileIO.write_seqids_to_fasta(tucked_seqids, "preCluster_out.tucked.fasta", d)

    infof = open("preCluster.cluster_info.csv", "w")
    infof.write("cluster,size\n")
    # write out a directory per preCluster cid in preCluster_out/<cid>
    # Liz note: right now, write out even directories with just 1 sequence
    # (we know they have "tucked" support, so can run Partial/Arrow on it)
    # singlef = open("preCluster_out.singles.fasta", 'w')
    for cid in pCS.S:
        #    if pCS.S[cid].size == 1:
        #        r = d[pCS.S[cid].members[0]]
        #        singlef.write(">{0}\n{1}\n".format(r.id, r.seq))
        #    else:
        if True:
            if not dun_make_bins:
                dirname = os.path.join("preCluster_out", str(cid))
                os.makedirs(dirname)
                file = os.path.join(dirname, "isoseq_flnc.fasta")
                FileIO.write_seqids_to_fasta(pCS.S[cid].members, file, d)
            infof.write("{},{}\n".format(cid, len(pCS.S[cid].members)))
    # singlef.close()
    infof.close()

    if not dun_cleanup_files:  # clean up all seed* and batch* files
        for file in glob.glob("batch*fasta*"):
            os.remove(file)
        for file in glob.glob("seed*fasta*"):
            os.remove(file)


if __name__ == "__main__":
    import argparse
    from argparse import ArgumentParser

    parser = ArgumentParser("PreCluster processing of isoseq_flnc.fasta using minimap")
    parser.add_argument(
        "--cpus", default=12, type=int, help="Number of CPUs minimap uses (default: 12)"
    )
    parser.add_argument(
        "--dun_make_bins",
        default=False,
        action="store_true",
        help="Only write out CSV files, do not make the actual bins (default: OFF)",
    )
    parser.add_argument(
        "--dun_use_partial",
        default=False,
        action="store_true",
        help="Don't use partial hits (default: OFF)",
    )
    parser.add_argument(
        "--num_seqs_per_batch",
        default=100000,
        type=int,
        help="Number of seqs per batch (default: 100000, use less for targeted)",
    )
    parser.add_argument(
        "--dun_cleanup_files",
        action="store_true",
        default=False,
        help=argparse.SUPPRESS,
    )

    args = parser.parse_args()

    # basic sanity checking here
    if not os.path.exists("isoseq_flnc.fasta"):
        print(
            "Expects isoseq_flnc.fasta in local directory but failed! Abort.",
            file=sys.stderr,
        )
        sys.exit(-1)

    if os.path.exists("preCluster_out"):
        print("preCluster_out/ already exists! Abort.", file=sys.stderr)
        sys.exit(-1)

    main(
        args.cpus,
        args.dun_make_bins,
        args.dun_use_partial,
        args.num_seqs_per_batch,
        args.dun_cleanup_files,
    )
