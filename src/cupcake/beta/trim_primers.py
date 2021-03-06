"""
Experimemtal code for trimming primers & polyA tails from high error rate long reads
"""
from csv import DictWriter
from dataclasses import dataclass
from multiprocessing import Process
from pathlib import Path
from typing import List, Optional, Union

import parasail
import typer
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from cupcake import version_callback
from cupcake import cupcake_logger as logger


@dataclass
class ScoreTuple:
    score5: int
    end5: int
    score3: int
    end3: int
    endA: int


# for ONT using Clontech
SEQ_5P = "AAGCAGTGGTATCAACGCAGAGTACATGGGG"
SEQ_3P_REV = "GTATCAACGCAGAGTAC"

# SEQ_5P = 'GCAATGAAGTCGCAGGGTTGGGG'
# SEQ_5P = 'CAGGAAACAGCTATGACC'	#SEQ_5P = 'CAGGAAACAGCTATGACC'
# SEQ_5P = 'AAGCAGTGGTATCAACGCAGAGTACATGGGG'	#SEQ_3P_REV = 'AAGCAGTGGTATCAACGCAGAGTAC'
SEQ_3P_REV = "AAGCAGTGGTATCAACGCAGAGTAC"  # SEQ_3P_REV = 'ACTGGCCGTCGTTTTAC'

MINSCORE_5P = 20
MINSCORE_3P = 20
MIN_A_LEN = 20

SCOREMAT = parasail.matrix_create("ACGT", 2, -5)


app = typer.Typer(name="cupcake.beta.trim_primers")


def trim5p3p_helper(r: SeqRecord) -> ScoreTuple:
    """
    Search for 5' and 3' in the first and last 100 bp window
    """
    s1 = str(r.seq[:100])
    s2 = str(r.reverse_complement().seq[:100])

    o1 = parasail.sg_qx_trace(s1, SEQ_5P, 3, 1, SCOREMAT)
    o2 = parasail.sg_qe_db_trace(s2, SEQ_3P_REV, 3, 1, SCOREMAT)
    lenA = None
    if o2.score >= MINSCORE_3P:
        lenA = trimA(s2[o2.end_query + 1 :])

    if MIN_A_LEN == 0:
        end3 = len(r.seq) - o2.end_query - 1
        return ScoreTuple(
            score5=o1.score, end5=o1.end_query, score3=o2.score, end3=end3, endA=end3
        )
    elif lenA is not None:
        end3 = len(r.seq) - o2.end_query - 1
        endA = end3 - lenA + 1
        return ScoreTuple(
            score5=o1.score, end5=o1.end_query, score3=o2.score, end3=end3, endA=endA
        )
    else:
        end3 = len(r.seq) - o2.end_query - 1
        return ScoreTuple(
            score5=o1.score, end5=o1.end_query, score3=o2.score, end3=end3, endA=end3
        )


def trimA(rev_seq: str) -> Optional[int]:
    if len(rev_seq) == 0:
        return None
    n_rev_seq = len(rev_seq)
    mismatch = 0
    i = 0
    while mismatch < 2 and i < n_rev_seq:
        if rev_seq[i] != "T":
            mismatch += 1
        i += 1
    i -= 1
    if i >= MIN_A_LEN:
        return i
    else:
        return None


def trim5p3p_multithreaded(
    fastq_filename: Union[str, Path], output_prefix: str, chunks: int
) -> None:
    # first figure out how many records there are and record positions
    num_lines = 0
    for _ in open(fastq_filename, "r"):
        num_lines += 1
    num_records = num_lines // 4
    chunk_size = (num_records // chunks) + (num_records % chunks > 0)
    logger.info(f"{fastq_filename} has {num_records} records, {chunk_size} per chunk")

    pools = []
    records = []
    count = 0
    i = 1
    for r in SeqIO.parse(open(fastq_filename), "fastq"):
        count += 1
        records.append(r)
        if count >= chunk_size:
            p = Process(target=trim5p3p, args=(records, f"{output_prefix}.{str(i)}"))
            p.start()
            print(f"Starting worker {i}...")
            pools.append(p)
            records = []
            count = 0
            i += 1
    p = Process(target=trim5p3p, args=(records, f"{output_prefix}.{str(i)}"))
    p.start()
    logger.info(f"Starting worker {i}...")
    pools.append(p)

    for p in pools:
        p.join()


def trim5p3p(records: List[SeqRecord], output_prefix: str) -> None:
    with open(f"{output_prefix}.fl.fasta", "w") as f_FL, open(
        f"{output_prefix}.fl.clips", "w"
    ) as f_FL_clips, open(f"{output_prefix}.nfl.fasta", "w") as f_nFL, open(
        f"{output_prefix}.csv", "w"
    ) as f_csv:
        writer = DictWriter(f_csv, fieldnames=["id", "end5", "end3", "endA", "strand"])
        writer.writeheader()

        for r in records:
            r2 = r.reverse_complement()
            r2.id = r.id
            t1 = trim5p3p_helper(r)
            t2 = trim5p3p_helper(r2)

            is_fl_flag1 = (
                t1.score5 >= MINSCORE_5P
                and t1.score3 >= MINSCORE_3P
                and (MIN_A_LEN == 0 or t1.endA != t1.end3)
            )
            is_fl_flag2 = (
                t2.score5 >= MINSCORE_5P
                and t2.score3 >= MINSCORE_3P
                and (MIN_A_LEN == 0 or t2.endA != t2.end3)
            )

            if is_fl_flag1:
                if is_fl_flag2:
                    if t1.score5 + t1.score3 > t2.score5 + t2.score3:
                        strand = "+"
                    else:
                        strand = "-"
                else:  # pick t1
                    strand = "+"
            elif is_fl_flag2:
                strand = "-"
            else:
                strand = "NA"

            info = {
                "id": r.id,
                "end5": "NA",
                "end3": "NA",
                "endA": "NA",
                "strand": "NA",
            }

            if strand == "+":
                info["strand"] = "+"
                info["end5"] = t1.end5
                info["end3"] = t1.end3
                info["endA"] = t1.endA
                f_FL.write(f">{r.id}\n{r.seq[t1.end5:t1.endA]}\n")
                f_FL_clips.write(
                    f">{r.id}_5p strand:+ score:{t1.score5}\n{r.seq[:t1.end5]}\n"
                )
                f_FL_clips.write(
                    f">{r.id}_3p strand:+ score:{t1.score3}\n{r.seq[t1.endA:]}\n"
                )
            elif strand == "-":
                info["strand"] = "-"
                info["end5"] = t2.end5
                info["end3"] = t2.end3
                info["endA"] = t2.endA
                f_FL.write(f">{r2.id}\n{r2.seq[t2.end5:t2.endA]}\n")
                f_FL_clips.write(
                    f">{r.id}_5p strand:- score:{t2.score5}\n{r2.seq[:t2.end5]}\n"
                )
                f_FL_clips.write(
                    f">{r.id}_3p strand:- score:{t2.score3}\n{r2.seq[t2.endA:]}\n"
                )
            else:
                # non-fL, but we still wanna trim away the stuff
                if t1.score5 + t1.score3 > t2.score5 + t2.score3:
                    f_nFL.write(f">{r.id} strand:+?\n{r.seq[t1.end5:t1.endA]}\n")
                else:
                    f_nFL.write(f">{r.id} strand:-?\n{r2.seq[t2.end5:t2.endA]}\n")
            writer.writerow(info)


@app.command(name="")
def main(
    fastq_filename: str = typer.Argument(...),
    output_prefix: str = typer.Argument(...),
    chunks: int = typer.Option(
        10,
        "--chunks",
        "-n",
        help="Number of chunks (CPUs) to use, default 10",
    ),
    version: bool = typer.Option(
        None,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Prints the version of the SQANTI3 package.",
    ),
):
    trim5p3p_multithreaded(fastq_filename, output_prefix, chunks)


if __name__ == "__main__":
    typer.run(main)
