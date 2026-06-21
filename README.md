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

### Command-line Interface

Run `fastQpick` with a specified fraction and options:
```bash
fastQpick [OPTIONS] FASTQ_FILE1 FASTQ_FILE2 ...
```

### Python API

Use `fastQpick` in your Python code:
```python
from fastQpick import fastQpick

fastQpick(
    input_files=['FASTQ_FILE1', 'FASTQ_FILE2', ...],
    ...
)
```

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
- Gzip-compressed output by default, using the ISA-L-accelerated [`isal`](https://github.com/pycompression/python-isal) library to keep compression from bottlenecking the write pass. Pass `--disable-gzip` (CLI) or `disable_gzip=True` (Python API) to write plain FASTQ instead.

---

## License

fastQpick is licensed under the 2-clause BSD license. See the [LICENSE](LICENSE) file for details.

---

## Contributing

We welcome contributions! Please see the [CONTRIBUTING.md](CONTRIBUTING.md) file for guidelines on how to get involved.

