#!/usr/bin/env python
__author__ = "etseng@pacb.com"

import subprocess
import sys

from click.exceptions import FileError


def run_cmd(cmd):
    if subprocess.check_call(cmd, shell=True) != 0:
        raise Exception("Error cmd:").with_traceback(cmd)

input = sys.argv[1]  # ex: test.sam

if not input.endswith(".sam"):
    raise FileError("Only accepts files ending in .sam. Abort!")

prefix = input[:-4]

run_cmd("samtools view -bS {0}.sam > {0}.bam".format(prefix))
run_cmd("samtools sort {0}.bam > {0}.sorted.bam".format(prefix))
run_cmd(f"samtools index {prefix}.sorted.bam")
