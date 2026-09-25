#!/usr/bin/env python
"""Generate one fastQpick bootstrap replicate of the paired yeast library,
re-quantify it with kallisto, and write the per-transcript est_counts vector.

Usage: run_one.py <seed> <tmpdir> <out_counts.txt>

Written against the fastQpick 1.x API (seed=, without_replacement=, read_counts=,
threads=). Each replicate is a full-size sample with replacement (fraction=1.0) of
the two mate files, grouped so that mates stay synchronized. The replicate is written
uncompressed and deleted once kallisto has quantified it, to bound disk use.
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
N_READS = 5_725_730          # read pairs in SRR453566; skips fastQpick's counting pass
KALLISTO_THREADS = "4"


def main():
    seed = int(sys.argv[1])
    tmpdir = sys.argv[2]
    out_counts = sys.argv[3]
    rep_dir = os.path.join(tmpdir, "rep")
    kq_dir = os.path.join(tmpdir, "kq")
    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)
    os.makedirs(tmpdir, exist_ok=True)

    # Full-size bootstrap replicate, with replacement, mates synchronized. One thread per
    # replicate: drive.sh runs several replicates side by side instead.
    fastQpick(
        input_files=[M1, M2],
        fraction=1.0,
        seed=seed,
        output_dir=rep_dir,
        file_group_size=2,
        without_replacement=False,
        disable_gzip=True,
        read_counts=N_READS,
        threads=1,
        overwrite=True,
        verbose=False,
    )
    r1 = os.path.join(rep_dir, "SRR453566_1.fastq")
    r2 = os.path.join(rep_dir, "SRR453566_2.fastq")

    subprocess.run(
        [KALLISTO, "quant", "-i", IDX, "-o", kq_dir, "-t", KALLISTO_THREADS, r1, r2],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )

    # abundance.tsv columns: target_id length eff_length est_counts tpm
    est = np.loadtxt(os.path.join(kq_dir, "abundance.tsv"),
                     skiprows=1, usecols=3)
    np.savetxt(out_counts, est, fmt="%.6f")
    shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    main()
