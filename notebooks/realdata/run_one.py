#!/usr/bin/env python
"""Generate one fastQpick bootstrap replicate of the paired yeast library,
re-quantify it with kallisto, and write the per-transcript est_counts vector.

Usage: run_one.py <seed> <tmpdir> <out_counts.txt>
"""
import os
import sys
import shutil
import subprocess

import numpy as np
from fastQpick import fastQpick

HERE = os.path.dirname(os.path.abspath(__file__))
M1 = os.path.join(HERE, "data", "SRR453566_1.fastq.gz")
M2 = os.path.join(HERE, "data", "SRR453566_2.fastq.gz")
IDX = os.path.join(HERE, "yeast.idx")
KALLISTO = os.path.join(HERE, "kallisto", "kallisto")
N_READS = 5_725_730


def main():
    seed = int(sys.argv[1])
    tmpdir = sys.argv[2]
    out_counts = sys.argv[3]
    rep_dir = os.path.join(tmpdir, "rep")
    kq_dir = os.path.join(tmpdir, "kq")
    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)
    os.makedirs(tmpdir, exist_ok=True)

    # Full-size bootstrap replicate, with replacement, mates synchronized.
    fastQpick(
        input_files=[M1, M2],
        fraction=1.0,
        seeds=seed,
        output_dir=rep_dir,
        file_group_size=2,
        replacement=True,
        overwrite=True,
        verbose=False,
        fastq_to_length_dict={M1: N_READS, M2: N_READS},
    )
    r1 = os.path.join(rep_dir, "SRR453566_1.fastq")
    r2 = os.path.join(rep_dir, "SRR453566_2.fastq")

    subprocess.run(
        [KALLISTO, "quant", "-i", IDX, "-o", kq_dir, "-t", "4", r1, r2],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )

    # abundance.tsv columns: target_id length eff_length est_counts tpm
    est = np.loadtxt(os.path.join(kq_dir, "abundance.tsv"),
                     skiprows=1, usecols=3)
    np.savetxt(out_counts, est, fmt="%.6f")
    shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    main()
