# fastQpick

Fast and memory-efficient sampling of DNA-seq or RNA-seq FASTQ reads, with or without replacement. Sampling with replacement generates bootstrap replicates for uncertainty quantification in downstream analyses, and the same tool covers oversampling and subsampling (e.g., to equalize depth or to build smaller inputs for testing and benchmarking).

---

## Installation

### Install via PyPI
```bash
pip install fastQpick
```

### Install from Source Code

Using pip:
```bash
pip install git+https://github.com/pachterlab/fastQpick.git
```

---

## Usage

`fastQpick` runs from the command line or from Python, and both entry points accept the same options. Output files keep the input file names and are written to `--output-dir` (default `fastQpick_output/`), gzip-compressed by default. By default, sampling is with replacement at `--fraction 1` (a standard bootstrap replicate); `--without-replacement` gives subsampling. When sampling with replacement, repeated reads receive unique read names (`_1`, `_2`, ...) unless `--no-unique-headers` is set.

### Command-line examples

Generate one full-size bootstrap replicate (same number of reads as the input, sampled with replacement):
```bash
fastQpick -f 1 sample.fastq.gz
```

Generate 200 bootstrap replicates of a paired-end library. `-g 2` treats each consecutive pair of files as mates and keeps them synchronized. Replicates are written as `sample_R1.seed42.fastq.gz`, `sample_R1.seed43.fastq.gz`, ...:
```bash
fastQpick -f 1 -B 200 -g 2 -o bootstraps sample_R1.fastq.gz sample_R2.fastq.gz
```

Subsample 20% of the reads without replacement, with a fixed seed:
```bash
fastQpick -f 0.2 --without-replacement -s 7 -o subsampled sample.fastq.gz
```

Oversample to twice the original depth (any fraction of 1 or more samples with replacement):
```bash
fastQpick -f 2 sample.fastq.gz
```

Sample synchronized I1/R1/R2 triples (e.g., 10x Genomics), for every FASTQ file in a directory:
```bash
fastQpick -f 0.1 --without-replacement -g 3 -o subsampled fastq_dir/
```

Write plain (uncompressed) FASTQ:
```bash
fastQpick -f 1 --disable-gzip sample.fastq.gz
```

Also save the out-of-bag reads (the reads that were not drawn) of each replicate to `sample.oob.fastq.gz`, e.g., for cross-validation or the .632+ estimator. Without replacement, this yields complementary train/test splits:
```bash
fastQpick -f 1 --oob sample.fastq.gz
```

Write each sampled read once, with its multiplicity recorded in the header as `;size=<count>` (the USEARCH/VSEARCH abundance convention), instead of writing duplicate records. A full-size bootstrap replicate then contains about 63% as many records:
```bash
fastQpick -f 1 --collapse-duplicates sample.fastq.gz
```

The default mode first reads each file once to count its reads. If the counts are already known (e.g., from a QC report), pass them with `--read-counts`, comma-separated in input order, one per file or one per group with `-g`. fastQpick checks each count during the writing pass and stops with an error if it does not match the file:
```bash
fastQpick -f 1 -g 2 --read-counts 5725730 sample_R1.fastq.gz sample_R2.fastq.gz
```

Read counting, separate files, and separate replicates (`-n`, or a seed range) run in parallel. `-t/--threads` sets the total number of threads, shared between parallel jobs and gzip compression. It defaults to 4 (or to the number of available cores, if fewer). The output does not depend on the thread count:

```bash
fastQpick -f 1 -n 20 -t 8 sample.fastq.gz
```

### Choosing a sampling mode

| Mode | Flag | Output size | Peak memory (500M reads, `-f 1`, one replicate) | Reads a pipe | When to use |
|---|---|---|---|---|---|
| Default | (none) | exact | ~0.6 GB (~1.2 bytes/read) | no | An exact number of output reads is required, and the input is a file. |
| Single-pass | `-p` / `--single-pass` | exact in expectation | ~0.15 GB (constant) | yes | The input is a stream, or it is gzipped and you would rather not decompress it twice, or minimal memory matters more than an exact output size (relative standard deviation `1/sqrt(fraction * n)`). |

```bash
fastQpick -f 1 sample.fastq.gz                  # exact size, two passes
fastQpick -f 1 --single-pass sample.fastq.gz    # approximate size, constant memory, one pass
```

### Streaming from a pipe

The default mode reads the library twice, once to count the reads and once to
write the sample, so it needs a file they can re-open. The single-pass sampler never needs the
read count, so it can take the library on standard input as `-`:

```bash
zcat sample.fastq.gz | fastQpick -f 1 --single-pass -o out -    # from a pipe
fastq-dump --stdout SRR000001 | fastQpick -f 0.1 --single-pass -dr -o out -
```

The output is written to `out/stdin.fastq[.gz]`. A gzipped stream is detected and decompressed
automatically, so `cat sample.fastq.gz |` works as well as `zcat sample.fastq.gz |`.

Streaming is only available with `--single-pass`, only for a single input, and cannot be combined
with file grouping (`-g`), since each member of a group needs its own stream. fastQpick reports
an error rather than sampling incorrectly if any of these is violated.

On a gzipped library this is also the cheaper mode even when memory is plentiful: the two-pass
mode decompresses the whole file twice, which on a 500-million-read library costs about ten
minutes of CPU per extra pass.

### Python API

```python
from fastQpick import fastQpick

# 200 paired-end bootstrap replicates
fastQpick(
    input_files=["sample_R1.fastq.gz", "sample_R2.fastq.gz"],
    fraction=1.0,
    num_samples=200,
    file_group_size=2,
    output_dir="bootstraps",
)

# 20% subsample without replacement
fastQpick(
    input_files="sample.fastq.gz",
    fraction=0.2,
    without_replacement=True,
    seed=7,
    output_dir="subsampled",
)

# Skip the counting pass when the number of reads is already known
# (one count per file or per group, in input order, or a dict {path: count})
fastQpick(
    input_files="sample.fastq.gz",
    fraction=1.0,
    read_counts=5725730,
)
```

### Example: bootstrap standard errors for a quantification pipeline

```bash
fastQpick -f 1 -B 100 -g 2 -o boot sample_R1.fastq.gz sample_R2.fastq.gz
for seed in $(seq 42 141); do
    kallisto quant -i index.idx -o quant_${seed} boot/sample_R1.seed${seed}.fastq.gz boot/sample_R2.seed${seed}.fastq.gz
done
```
The spread of any downstream statistic across the `quant_*` runs estimates its sequencing sampling uncertainty. See the [tutorial notebooks](#tutorials) for a complete analysis.

---

## Documentation

- **Command-line Help**: Use the following command to see all available options:
  ```bash
  fastQpick --help
  ```

- **Python API Help**: Use the `help` function to explore the API:
  ```python
  help(fastQpick)
  ```


---

## Tutorials

Two Jupyter notebooks in [`notebooks/`](notebooks/) walk through `fastQpick` end-to-end:

- **[`intro.ipynb`](notebooks/intro.ipynb)** — Getting started on synthetic data. Simulates a small RNA-seq experiment with known transcript abundances, draws bootstrap replicates with replacement (`fraction=1.0`, `without_replacement=False`), and shows that the bootstrap standard errors recover the analytic multinomial sampling error.
- **[`yeast_example.ipynb`](notebooks/yeast_example.ipynb)** — Real-data application reproducing Figure 1 of the paper. Bootstraps a paired-end yeast RNA-seq dataset (SRA `SRR453566`), re-quantifies each replicate with `kallisto`, and characterizes the bootstrap distribution of the transcript abundance estimates.

---

## Features

- Time efficient - streams through the fastq and writes output in batches - generates a full-size (fraction=1, with replacement) bootstrap replicate of a 500M-read FASTQ in ~31 minutes in standard mode and ~22 minutes in single-pass mode (see [Benchmark](#benchmark) below).
- Memory efficient - the occurrence vector is sized to the largest per-read count actually drawn (one byte per read in the common case) and filled block by block, so neither the array of sampled indices nor a length-n counting temporary is ever materialized.
- Optional out-of-bag output (`--oob`) and multiplicity-tagged, duplicate-free output (`--collapse-duplicates`).
- Single-threaded by design where it matters: fastQpick never calls BLAS, so on import it pins the OpenBLAS pool bundled with numpy to one thread (`OPENBLAS_NUM_THREADS=1`, unless already set), which otherwise starts one busy-waiting thread per core in every process.
- Gzip-compressed output by default, using the ISA-L-accelerated [`isal`](https://github.com/pycompression/python-isal) library to keep compression from bottlenecking the write pass. Pass `--disable-gzip` (CLI) or `disable_gzip=True` (Python API) to write plain FASTQ instead.

---

## Benchmark

Table 1 of the manuscript is produced by `benchmarks/table1.py`. It times fastQpick (default two-pass mode and single-pass mode), seqtk, and seqkit on one synthetic 500-million-read library of 150 bp reads (158 GB uncompressed, 40 GB gzipped), generated by `benchmarks/make_inputs.sh`, in three conditions:

| Condition | Sample | Input | Output | Threads |
|---|---|---|---|---|
| bootstrap | full-size, with replacement (`-f 1`) | gzip | gzip | 4 |
| subsample20 | 20% without replacement (`-f 0.2 -dr`) | gzip | plain (seqtk cannot write gzip) | 4 |
| bootstrap10 | 10 full-size bootstrap replicates in one call (`-f 1 -B 10`) | gzip | gzip | 4 (fastQpick only) |

Every tool that exposes a thread count is given four (`--threads`, the fastQpick default; seqtk has none). With a single output file fastQpick uses the extra threads only for gzip compression; in the bootstrap10 condition it writes four replicates at a time. Every run starts from a cold page cache (the input is evicted with `posix_fadvise(DONTNEED)`), uses seed 42, and records wall time, the peak resident set of the largest process (`/usr/bin/time -v`), and the peak summed resident set of the whole process tree (sampled every 0.5 s, which is what matters for the multi-process bootstrap10 runs). `--replicates` (default 3) sets the number of timed runs per cell (the manuscript reports medians over at least three runs, except for bootstrap10, which is a single run); runs already recorded in the results TSV are skipped, so an interrupted benchmark can be resumed. The reported values are medians across replicates (`benchmarks/summarize_table1.py`). seqtk and seqkit sample each read independently with probability 0.2 (Bernoulli), so like fastQpick's single-pass mode their output size is exact only in expectation, and neither samples with replacement.

The hardware used in the manuscript: two Intel Xeon Gold 6152 CPUs (2.10 GHz, 44 cores), 754 GB of RAM, and a 12 TB 7200 rpm hard disk drive (ext4), running CentOS Stream 8 (kernel 4.18) and Python 3.10.

---

## License

fastQpick is licensed under the 2-clause BSD license. See the [LICENSE](LICENSE) file for details.

---

## Contributing

We welcome contributions! Please see the [CONTRIBUTING.md](CONTRIBUTING.md) file for guidelines on how to get involved.

---

## Manuscript

Read the manuscript describing fastQpick in the [bioRxiv preprint](https://www.biorxiv.org/content/10.64898/2026.06.23.734068v1) (DOI: 10.64898/2026.06.23.734068).

The figure and table of the manuscript are reproduced by two scripts, run with the current release (fastQpick v1.0.0):

- **Figure 1** (read-level bootstrap of the yeast library SRR453566): `notebooks/realdata/drive.sh [B] [P]` downloads kallisto, the reference, and the reads, quantifies the original library, generates and re-quantifies `B` (default 200) bootstrap replicates `P` (default 8) at a time, and renders `notebooks/figures/bootstrap_realdata.png` while printing the numbers quoted in the Application section. The same analysis is walked through in `notebooks/yeast_example.ipynb`.
- **Table 1** (runtime and memory on a 500-million-read library): `benchmarks/make_inputs.sh <dir>` generates the synthetic library and `benchmarks/table1.py --input <dir>/bench_500M.fastq.gz` times every cell of the table from a cold page cache (see [Benchmark](#benchmark)). `benchmarks/summarize_table1.py` reduces the runs to the medians reported in the table.

The `manuscript` git tag marks the code used for the original submission (v0.3.0); the revised manuscript was produced with v1.0.0.
