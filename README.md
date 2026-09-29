# fastQpick

Fast and memory-efficient sampling of DNA-seq or RNA-seq FASTQ reads, with or without replacement. Enables bootstrap replicates for uncertainty quantification in downstream analyses, and oversampling or subsampling for depth equalization or testing and benchmarking.

## Installation

```bash
pip install fastQpick
```

## Usage

Generate one full-size bootstrap replicate (same number of reads as the input, sampled with replacement):
```bash
fastQpick -f 1 sample.fastq.gz
```

Generate 200 bootstrap replicates:
```bash
fastQpick -f 1 -B 200 sample.fastq.gz
```

Generate a 20% subsample without replacement:
```bash
fastQpick -f 0.2 --without-replacement sample.fastq.gz
```

Generate a 2x oversample:
```bash
fastQpick -f 2 sample.fastq.gz
```

Generate 200 bootstrap replicates of a paired-end library:
```bash
fastQpick -f 1 -B 200 -g 2 sample_R1.fastq.gz sample_R2.fastq.gz
```

Generate 200 bootstrap replicates of a single-cell RNA-seq library:
```bash
fastQpick -f 1 -B 200 -g 3 sample_I1.fastq.gz sample_R1.fastq.gz sample_R2.fastq.gz
```

Write plain (uncompressed) FASTQ:
```bash
fastQpick -f 1 --disable-gzip sample.fastq.gz
```

Stream from a pipe (single-pass mode only):
```bash
zcat sample.fastq.gz | fastQpick -f 1 --single-pass -o out -
```

Download a library from SRA and stream from a pipe (single-pass mode only):
```bash
fastq-dump --stdout SRR000001 | fastQpick -f 0.1 --single-pass -r -o out -
```

Also save the out-of-bag reads (the reads that were not drawn, about 36.8% of the library for a full-size bootstrap) of each replicate to `sample.oob.fastq.gz`, e.g., for cross-validation with the .632+ estimator. Without replacement, this yields complementary train/test splits:
```bash
fastQpick -f 1 --oob sample.fastq.gz
```

Write each sampled read once, with its multiplicity recorded in the header as `;size=<count>` (the USEARCH/VSEARCH abundance convention), instead of writing duplicate records. A full-size bootstrap replicate then contains about 63.2% as many records, and duplication-aware tools can process each distinct read once. By default, duplicates are written as ordinary records, so any downstream tool runs unchanged:
```bash
fastQpick -f 1 --collapse-duplicates sample.fastq.gz
```

## Example: bootstrap standard errors for a quantification pipeline

```bash
fastQpick -f 1 -B 100 -g 2 -o boot sample_R1.fastq.gz sample_R2.fastq.gz
for seed in $(seq 42 141); do
    kallisto quant -i index.idx -o quant_${seed} boot/sample_R1.seed${seed}.fastq.gz boot/sample_R2.seed${seed}.fastq.gz
done
```
The spread of any downstream statistic across the `quant_*` runs estimates its sequencing sampling uncertainty. See the [tutorial notebooks](#tutorials) for a complete analysis.

## Tutorials

Two Jupyter notebooks in [`notebooks/`](notebooks/) walk through `fastQpick` end-to-end:

- **[`intro.ipynb`](notebooks/intro.ipynb)** — Getting started on synthetic data. Simulates a small RNA-seq experiment with known transcript abundances, draws bootstrap replicates with replacement (`fraction=1.0`, `without_replacement=False`), and shows that the bootstrap standard errors recover the analytic multinomial sampling error.
- **[`yeast_example.ipynb`](notebooks/yeast_example.ipynb)** — Real-data application reproducing Figure 1 of the paper. Bootstraps a paired-end yeast RNA-seq dataset (SRA `SRR453566`), re-quantifies each replicate with `kallisto`, and characterizes the bootstrap distribution of the transcript abundance estimates.

## Choosing a sampling mode

| Mode | Flag | Output size | Reads a pipe | When to use |
|---|---|---|---|---|---|
| Default | (none) | exact (`floor(fraction * n)` reads) | no | Most cases. |
| Single-pass | `-p` / `--single-pass` | exact only in expectation | yes | At the cost of exactness of output size - (1) the input is a stream, (2) sampling a single sample quickly, or (3) write many samples with multithreading without linear memory growth.

## Python API

```python
from fastQpick import fastQpick
fastQpick(...)
```

The Python parameters are the same as the command-line flags, with underscores instead of hyphens.

## Documentation

```bash
fastQpick --help
```

```python
help(fastQpick)
```

## MCP server (LLM agents)

fastQpick ships a [Model Context Protocol](https://modelcontextprotocol.io) server so that LLM agents (Claude Code, Claude Desktop, Cursor, Codex, ...) can sample FASTQ files directly. It requires Python 3.10 or later:
```bash
pip install "fastQpick[mcp]"
```

The server exposes three tools: `sample_fastq` (the full set of options above), `count_fastq_reads` (read counts to pass back as `read_counts`), and `list_fastq_files` (the order in which a directory's files are grouped). Paths refer to the machine running the server. To register it with Claude Code:
```bash
claude mcp add fastqpick -- fastQpick-mcp
```
or, for any client that reads a JSON configuration:
```json
{"mcpServers": {"fastqpick": {"command": "fastQpick-mcp"}}}
```
Each sampling call runs the fastQpick command line in its own process and returns the output file paths and the end of the log. Streaming from standard input is not available through the server.

## License

fastQpick is licensed under the 2-clause BSD license. See the [LICENSE](LICENSE) file for details.

## Contributing

We welcome contributions! Please see the [CONTRIBUTING.md](CONTRIBUTING.md) file for guidelines on how to get involved.

## Manuscript

Read the manuscript describing fastQpick in the [bioRxiv preprint](https://www.biorxiv.org/content/10.64898/2026.06.23.734068v1) (DOI: 10.64898/2026.06.23.734068).

The `manuscript` git tag marks the code used for the original submission (v0.3.0); the revised manuscript was produced with v1.0.0.
