# fastQpick

Fast and memory-efficient sampling of DNA-seq or RNA-seq FASTQ data with replacement. Useful for generating bootstrap replicates to estimate technical variance in downstream analyses, and for subsampling large datasets for testing and benchmarking.

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

`fastQpick` runs from the command line or from Python, and both entry points accept the same options. Output files keep the input file names and are written to `--output-dir` (default `fastQpick_output/`), gzip-compressed by default. Sampling is with replacement unless `--without-replacement` is set.

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

### Choosing a sampling mode

| Mode | Flag | Output size | Peak memory (500M reads, `-f 1`) | When to use |
|---|---|---|---|---|
| Default | (none) | exact | ~9.4 GB (~20 bytes/read) | The machine has enough memory. Fastest exact mode. |
| Low-memory | `-l` / `--low-memory` | exact | ~1.4 GB (~3 bytes/read) | Memory is limiting and an exact number of output reads is required. |
| Single-pass | `-p` / `--one_pass` | exact in expectation | ~0.1 GB (constant) | Memory is limiting and speed or minimal memory matters more than an exact output size (relative standard deviation `1/sqrt(fraction * n)`). |

```bash
fastQpick -f 1 --low-memory sample.fastq.gz   # exact, low peak memory
fastQpick -f 1 --one_pass sample.fastq.gz     # approximate size, constant memory, one pass
```

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
fastQpick(
    input_files="sample.fastq.gz",
    fraction=1.0,
    fastq_to_length_dict={"sample.fastq.gz": 5725730},
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

- **[`intro.ipynb`](notebooks/intro.ipynb)** — Getting started on synthetic data. Simulates a small RNA-seq experiment with known transcript abundances, draws bootstrap replicates with replacement (`fraction=1.0`, `replacement=True`), and shows that the bootstrap standard errors recover the analytic multinomial sampling error.
- **[`yeast_example.ipynb`](notebooks/yeast_example.ipynb)** — Real-data application reproducing Figure 1 of the paper. Bootstraps a paired-end yeast RNA-seq dataset (SRA `SRR453566`), re-quantifies each replicate with `kallisto`, and characterizes the bootstrap distribution of the transcript abundance estimates.

To reproduce the figures exactly as they appear in the manuscript, check out the `manuscript` tag before running the notebooks:
```bash
git checkout manuscript
```

---

## Features

- Time efficient - streams through the fastq and writes output in batches - generates a full-size (fraction=1, with replacement) bootstrap replicate of a 500M-read FASTQ in ~26 minutes in standard mode, ~56 minutes in low-memory mode, and ~35 minutes in one-pass mode (see [Benchmark](#benchmark) below).
- Memory efficient - the occurrence vector is sized to the largest per-read count actually drawn (one byte per read in the common case), and low-memory mode further avoids materializing the array of sampled indices.
- Optional out-of-bag output (`--oob`) and multiplicity-tagged, duplicate-free output (`--collapse-duplicates`).
- Gzip-compressed output by default, using the ISA-L-accelerated [`isal`](https://github.com/pycompression/python-isal) library to keep compression from bottlenecking the write pass. Pass `--disable-gzip` (CLI) or `disable_gzip=True` (Python API) to write plain FASTQ instead.

---

## Benchmark

Benchmarks were run on a 500-million-read FASTQ file (143 GB uncompressed) with plain-text output, on a server with two Intel Xeon Gold 6152 CPUs (2.10 GHz, 44 cores / 88 threads in total), 754 GB of RAM, and a 12 TB 7200 rpm hard disk drive (HGST Ultrastar, ext4), running CentOS Stream 8 with Python 3.10. fastQpick uses a single process per output file. Each run started from a cold page cache.

| Mode | Runtime | Peak memory | Output size |
|---|---|---|---|
| Default | 26 min | 9.4 GB | exact |
| Low-memory (`-l`) | 56 min | 1.4 GB | exact |
| Single-pass (`-p`) | 35 min | 0.11 GB | exact in expectation |

Comparison with other subsampling tools on a 20% subsample without replacement (m = 100M reads) of the same file, each tool with its default thread count. "Cold" runs started with the input evicted from the page cache, "warm" runs with the input resident in the page cache. seqtk and seqkit sample each read with probability 0.2, so their output size (like that of fastQpick's single-pass mode) is exact in expectation. fastp has no random subsampling: `--reads_to_process` keeps the first m reads and stops reading, so it processes only 20% of the input. None of the three tools samples with replacement.

| Tool | Cold | Warm | Peak memory | Output reads |
|---|---|---|---|---|
| seqtk 1.5 (`sample -s 42 0.2`) | 16.8 min | 4.1 min | 3 MB | 99,979,934 |
| seqkit 2.13 (`sample -p 0.2`) | 17.9 min | 3.5 min | 35 MB | 100,009,778 |
| fastp 1.3 (`--reads_to_process`, first m reads) | 3.6 min | 2.0 min | 1.1 GB | 100,000,000 |
| fastQpick default | 23.1 min | 8.5 min | 5.2 GB | 100,000,000 |
| fastQpick low-memory (`-l`) | -- | 8.5 min | 0.56 GB | 100,000,000 |
| fastQpick single-pass (`-p`) | 17.6 min | 5.2 min | 67 MB | 99,995,639 |

---

## License

fastQpick is licensed under the 2-clause BSD license. See the [LICENSE](LICENSE) file for details.

---

## Contributing

We welcome contributions! Please see the [CONTRIBUTING.md](CONTRIBUTING.md) file for guidelines on how to get involved.

