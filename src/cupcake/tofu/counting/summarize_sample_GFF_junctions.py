#!/usr/bin/env python
__author__ = "etseng@pacb.com"

"""
Reporter for junction summary for one or more samples.
Recommended to run before scrubbing sample GFFs prior to chaining.

Suggested process is:
1. run collapse to get GFF for each sample
2. run this report script on all sample GFFs
3. run scrubber on all sample GFFs
4. run chaining using scrubbed (cleaned) sample GFFs
"""
import sys
from collections import defaultdict
from csv import DictWriter
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import typer
from Bio import SeqIO
from sklearn.cluster import Birch

import cupcake.sequence.GFF as GFF
from cupcake.logging import cupcake_logger as logger

app = typer.Typer(name="cupcake.tofu.counting.summarize_sample_GFF_junctions")


def sanity_check(
    sample_dirs: List[Dict[str, Union[str, Path]]],
    gff_filename: Union[str, Path],
    genome_filename: Optional[Union[str, Path]] = None,
    junction_filename: Optional[Union[str, Path]] = None,
) -> None:
    for d in sample_dirs.values():
        file = Path(d, gff_filename)
        if not file.exists():
            logger.error(f"Expected GFF file {file} does not exist. Abort!")
            sys.exit(-1)

    if genome_filename is not None and not Path(genome_filename).exists():
        logger.error(f"Genome file {genome_filename} given but does not exist. Abort!")
        sys.exit(-1)

    if junction_filename is not None and not Path(junction_filename).exists():
        logger.error(
            f"Junction file {junction_filename} given but does not exist. Abort!"
        )
        sys.exit(-1)


def read_config(filename: Path) -> Tuple[Dict[str, Path], List[str], Path, Path, Path]:
    """
    SAMPLE=<name>;<path>

    must also have
    GFF_FILENAME=

    optional:
    GENOME_FILENAME=
    JUNCTION_FILENAME=
    GROUP_FILENAME=

    Everything else will be ignored (so you can re-use sample.config for chain_samples.py)
    """
    sample_dirs: Dict[str, Path] = {}
    sample_names: List[str] = []
    gff_filename: Optional[Union[str, Path]] = None
    genome_filename: Optional[Union[str, Path]] = None
    junction_filename: Optional[Union[str, Path]] = None

    if not filename.exists():
        raise FileNotFoundError(f"The config file {filename} could not be found!")

    with open(filename) as f:
        for line in f:
            if line.startswith("tmpSAMPLE="):
                logger.error(
                    "Please only use SAMPLE=, not tmpSAMPLE= for junction reports!"
                )
                sys.exit(-1)
            elif line.startswith("SAMPLE="):
                name, path = line.strip()[len("SAMPLE=") :].split(";")
                if name.startswith("tmp_"):
                    logger.error(
                        f"Sample names are not allowed to start with tmp_! Please change {name} to something else."
                    )
                    sys.exit(-1)
                sample_dirs[name] = Path(path).resolve()
                sample_names.append(name)
            elif line.startswith("GFF_FILENAME="):
                gff_filename = Path(line.strip()[len("GFF_FILENAME=") :])
            elif line.startswith("GENOME_FILENAME="):
                genome_filename = Path(line.strip()[len("GENOME_FILENAME=") :])
            elif line.startswith("JUNCTION_FILENAME="):
                junction_filename = Path(line.strip()[len("JUNCTION_FILENAME=") :])

    if gff_filename is None:
        raise Exception(
            f"Expected GFF_FILENAME= but not in config file {filename}! Abort."
        )

    if len(sample_names) == 0:
        logger.error("No samples given. Exit.")
        sys.exit(-1)

    return sample_dirs, gff_filename, genome_filename, junction_filename


def read_annotation_junction_bed(junction_filename: Union[str, Path]) -> defaultdict:
    """
    junction.bed is in format:

    seqname, left (0-based), right (0-based), +/-

    following junction.bed format from TopHat
    http://ccb.jhu.edu/software/tophat/manual.shtml
    """
    junction = defaultdict(dict)  # (seqname, strand) --> (start, end)
    for line in open(junction_filename):
        chrom, left, right, strand = line.strip().split("\t")
        junction[chrom, strand][(int(left), int(right))] = 1
    return junction


def summarize_junctions(
    sample_dirs: Dict[str, Path],
    # sample_names: List[str],
    gff_filename: Union[str, Path],
    output_prefix: Union[str, Path],
    genome_d: Optional[Union[str, Path]] = None,
    junction_known: Optional[Union[str, Path]] = None,
) -> defaultdict:
    """
    1. for each sample, read all the GFF, store the junction information (both 0-based)

    """
    junc_by_chr_strand = defaultdict(
        lambda: defaultdict(list)
    )  # (seqname,strand) --> (donor,acceptor) --> samples it show up in (more than once possible)

    for sample_name, d in sample_dirs.items():
        for r in GFF.collapseGFFReader(Path(d, gff_filename)):
            n = len(r.ref_exons)
            if n == 1:
                continue  # ignore single exon transcripts
            for i in range(n - 1):
                donor = r.ref_exons[i].end - 1  # make it 0-based
                accep = r.ref_exons[i + 1].start  # start is already 0-based
                junc_by_chr_strand[r.seqname, r.strand][donor, accep].append(
                    sample_name
                )

    # write junction report
    with open(f"{output_prefix}.junction.bed", "w") as f1, open(
        f"{output_prefix}.junction_detail.txt", "w"
    ) as f:
        f1.write(f'track name=junctions description="{output_prefix}" useScore=1\n')

        JUNC_DETAIL_FIELDS = [
            "seqname",
            "left",
            "right",
            "strand",
            "num_transcript",
            "num_sample",
            "genome",
            "annotation",
            "label",
        ]

        writer = DictWriter(f, JUNC_DETAIL_FIELDS, delimiter="\t")
        writer.writeheader()
        keys = list(junc_by_chr_strand)
        keys.sort()
        for _seqname, _strand in keys:
            v = junc_by_chr_strand[_seqname, _strand]
            v_keys = list(v)
            v_keys.sort()
            labels = cluster_junctions(v_keys)
            for i, (_donor, _accep) in enumerate(v_keys):
                rec = {
                    "seqname": _seqname,
                    "left": _donor,
                    "right": _accep,
                    "strand": _strand,
                    "num_transcript": len(v[_donor, _accep]),
                    "num_sample": len(set(v[_donor, _accep])),
                }
                # f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t".format(_chr, _donor, _accep, _strand, len(v[_donor,_accep]), len(set(v[_donor,_accep]))))
                f1.write(
                    f"{_seqname}\t{_donor}\t{_accep + 1}\t{output_prefix}\t{len(v[_donor, _accep])}\t{_strand}\n"
                )
                # if genome is given, write acceptor-donor site
                if genome_d is None or _seqname not in genome_d:
                    rec["genome"] = "NA"
                    # f.write("NA\t")
                else:
                    up, down = (
                        genome_d[_seqname][(_donor + 1) : (_donor + 3)],
                        genome_d[_seqname][(_accep - 2) : _accep],
                    )
                    if _strand == "+":
                        rec["genome"] = f"{str(up.seq).upper()}-{str(down.seq).upper()}"
                        # f.write("{0}-{1}\t".format(str(up.seq).upper(), str(down.seq).upper()))
                    else:
                        rec[
                            "genome"
                        ] = f"{str(down.reverse_complement().seq).upper()}-{str(up.reverse_complement().seq).upper()}"
                        # f.write("{0}-{1}\t".format(str(down.reverse_complement().seq).upper(), str(up.reverse_complement().seq).upper()))
                # if annotation is given, check if matches with annotation
                if junction_known is None:
                    rec["annotation"] = "NA"
                    # f.write("NA\n")
                else:
                    if (_seqname, _strand) in junction_known and (
                        _donor,
                        _accep,
                    ) in junction_known[_seqname, _strand]:
                        rec["annotation"] = "Y"
                        # f.write("Y\t")
                    else:
                        rec["annotation"] = "N"
                        # f.write("N\t")
                rec["label"] = f"{_seqname}_{_strand}_{labels[i]}"
                writer.writerow(rec)
            # f.write("{c}_{s}_{lab}\n".format(c=_seqname, s=_strand, lab=labels[i]))

    return junc_by_chr_strand


def cluster_junctions(juncs: List[int]) -> np.ndarray:
    birch_model = Birch(threshold=3, n_clusters=None)
    X = np.array(juncs)
    birch_model.fit(X)

    return birch_model.labels_


@app.command(name="")
def main(
    config: Union[str, Path] = typer.Argument(..., help="Config filename"),
    output_prefix: str = typer.Argument(..., help="Output prefix"),
):

    try:
        (
            sample_dirs,
            gff_filename,
            genome_filename,
            junction_filename,
        ) = read_config(config)
    except FileNotFoundError as error:
        logger.error(error)

    sanity_check(sample_dirs, gff_filename, genome_filename, junction_filename)

    if genome_filename is not None:
        logger.info(f"Reading genome file {genome_filename}...")
        genome_d = SeqIO.to_dict(SeqIO.parse(open(genome_filename), "fasta"))
    else:
        logger.info("No genome file given. Ignore.")
        genome_d = None

    if junction_filename is not None:
        logger.info(f"Reading junction file {junction_filename}....")
        junction_bed = read_annotation_junction_bed(junction_filename)
    else:
        logger.info("No junction file given. Ignore.")
        junction_bed = None

    summarize_junctions(
        sample_dirs,
        gff_filename,
        output_prefix,
        genome_d,
        junction_bed,
    )


if __name__ == "__main__":
    typer.run(main)
