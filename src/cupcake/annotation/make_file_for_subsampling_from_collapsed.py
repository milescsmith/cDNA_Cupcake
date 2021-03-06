#!/usr/bin/env python
import sys
from csv import DictReader
from pathlib import Path
from typing import Dict, Optional, Tuple

import typer
from Bio import SeqIO

from cupcake import cupcake_logger as logger
from cupcake import set_verbosity, version_callback

app = typer.Typer(
    name="cupcake.annotation.make_file_for_subsampling_from_collapsed",
    help="Make subsample-ready file from Iso-Seq collapsed output",
)


def read_demux_fl_count_file(filename: str) -> Tuple[Dict[str, str], str]:
    d = {}
    reader = DictReader(open(filename), delimiter=",")
    assert "id" in reader.fieldnames
    for r in reader:
        d[r["id"]] = r
    samples = reader.fieldnames
    samples.remove("id")
    return d, samples


def make_file_for_subsample(
    input_prefix: str,
    output_prefix: str,
    demux_file=None,
    matchAnnot_parsed=None,
    sqanti_class=None,
    include_single_exons=False,
) -> None:
    """
    Two files must exist: .abundance.txt and .rep.fq so we can make the length
    """
    count_filename = f"{input_prefix}.abundance.txt"

    rep_filenames = [
        (f"{input_prefix}.rep.fq", "fastq"),
        (f"{input_prefix}.rep.fastq", "fastq"),
        (f"{input_prefix}.rep.fa", "fasta"),
        (f"{input_prefix}.rep.fasta", "fasta"),
    ]

    rep_filename = None
    rep_type = None
    for x, feature in rep_filenames:
        if Path(x).exists():
            rep_filename = x
            rep_type = feature

    if rep_filename is None:
        logger.error(
            "Expected to find input fasta or fastq files {input_prefix}.rep.fa or {input_prefix}.rep.fq. Not found. Abort!"
        )
        sys.exit(-1)

    if not include_single_exons:
        from cupcake.sequence.GFF import collapseGFFReader

        gff_filename = f"{input_prefix}.gff"
        logger.info(f"Reading {gff_filename} to exclude single exons...")
        # good_ids = []
        good_ids = [
            r.seqid for r in collapseGFFReader(gff_filename) if len(r.ref_exons) >= 2
        ]
        # for r in collapseGFFReader(gff_filename):
        #     if len(r.ref_exons) >= 2:
        #         good_ids.append(r.seqid)
    else:
        good_ids = []

    if demux_file is None and not Path(count_filename).exists():
        logger.error(f"Cannot find {count_filename}. Abort!")
        sys.exit(-1)

    if matchAnnot_parsed is not None and not Path(matchAnnot_parsed).exists():
        logger.error(f"Cannot find {matchAnnot_parsed}. Abort!")
        sys.exit(-1)

    if sqanti_class is not None and not Path(sqanti_class).exists():
        logger.error(f"Cannot find {sqanti_class}. Abort!")
        sys.exit(-1)

    if matchAnnot_parsed is not None:
        with open(matchAnnot_parsed) as ma:
            match_dict = {r["pbid"]: r for r in DictReader(ma, delimiter="\t")}
        for k in match_dict:
            match_dict[k]["category"] = match_dict[k]["score"]
    elif sqanti_class is not None:
        logger.info(f"Reading {sqanti_class} to get gene/isoform assignment...")
        match_dict = {}
        with open(sqanti_class) as sc:
            for r in DictReader(sc, delimiter="\t"):
                if r["associated_transcript"] == "novel":
                    refisoform = f"novel_{r['isoform']}"
                else:
                    refisoform = r["associated_transcript"]
                match_dict[r["isoform"]] = {
                    "refgene": r["associated_gene"],
                    "refisoform": refisoform,
                    "category": r["structural_category"],
                }
    else:
        match_dict = None
    with open(rep_filename) as rf:
        seqlen_dict = {
            r.id.split("|")[0]: len(r.seq) for r in SeqIO.parse(rf, rep_type)
        }

    to_write = {}
    if demux_file is None:
        to_write["all"] = {}
        with open(count_filename) as f:
            while True:
                cur = f.tell()
                if not f.readline().startswith("#"):
                    f.seek(cur)
                    break
            for r in DictReader(f, delimiter="\t"):
                if r["pbid"] in good_ids or include_single_exons:
                    to_write["all"][r["pbid"]] = r["count_fl"]
    else:
        d, samples = read_demux_fl_count_file(demux_file)
        for s in samples:
            to_write[s] = {}
        for pbid, d2 in d.items():
            for s in samples:
                if pbid in good_ids or include_single_exons:
                    to_write[s][pbid] = d2[s]

    for sample in to_write:
        with Path(f"{output_prefix}.{sample}.txt").open("a+") as h:
            if matchAnnot_parsed is None and sqanti_class is None:
                h.write("pbid\tpbgene\tlength\tfl_count\n")
            else:
                h.write(
                    "pbid\tpbgene\tlength\trefisoform\trefgene\tcategory\tfl_count\n"
                )
            for pbid in to_write[sample]:
                if matchAnnot_parsed is not None or sqanti_class is not None:
                    if pbid not in match_dict:
                        logger.warning(
                            f"Ignoring {pbid} because not in annotation (SQANTI/MatchAnnot) file."
                        )
                        continue
                    m = match_dict[pbid]
                    h.write(f"{pbid}\t{pbid.split('.')[1]}\t{seqlen_dict[pbid]}\t")
                    h.write(f'{m["refisoform"]}\t{m["refgene"]}\t{m["category"]}\t')
                else:
                    h.write(f'{pbid}\t{pbid.split(".")[1]}\t{seqlen_dict[pbid]}\t')
                h.write(f"{to_write[sample][pbid]}\n")
            logger.info(
                f"Output written to {Path(f'{output_prefix}.{sample}.txt').resolve()}."
            )


@set_verbosity
@app.command(name="")
def main(
    input_prefix: str = typer.Option(
        None,
        "--input",
        "-i",
        help="Collapsed prefix",
    ),
    output_prefix: str = typer.Option(
        "subsample", "--output", "-o", help="Output prefix"
    ),
    matchAnnot_parsed: Optional[str] = typer.Option(
        None, "--matchannot", "-m1", help="MatchAnnot parsed output"
    ),
    sqanti_class: Optional[str] = typer.Option(
        None, "--sqanti_class", "-m2", help="SQANTI classification file"
    ),
    demux: Optional[str] = typer.Option(
        None,
        help="Demuxed FL count file - if provided, will be used instead of the <input_prefix>.abundance.txt file",
    ),
    include_single_exons: bool = typer.Option(
        True, show_default=False, help="Include single exons [default: OFF]"
    ),
    verbose: int = typer.Option(3, "--verbose", "-v", count=True),
    version: bool = typer.Option(
        None,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Prints the version of the SQANTI3 package.",
    ),
) -> None:
    set_verbosity(verbose)
    make_file_for_subsample(
        input_prefix,
        output_prefix,
        demux,
        matchAnnot_parsed,
        sqanti_class,
        include_single_exons,
    )


if __name__ == "__main__":
    typer.run(main)
