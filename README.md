# fastQpick

Fast and memory-efficient sampling of DNA-seq or RNA-seq FASTQ data with or without replacement.

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

Or clone the repository and build manually:
```bash
git clone https://github.com/pachterlab/fastQpick.git
cd fastQpick
python -m build
python -m pip install dist/fastQpick-x.x.x-py3-none-any.whl
```

---

## Usage

### Command-line Interface

Run `fastQpick` with a specified fraction and options:
```bash
fastQpick --fraction FRACTION [OPTIONS] FASTQ_FILE1 FASTQ_FILE2 ...
```

### Python API

Use `fastQpick` in your Python code:
```python
from fastQpick import fastQpick

fastQpick(
    input_file_list=['FASTQ_FILE1', 'FASTQ_FILE2', ...],
    fraction=FRACTION,
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

### Options
- input_files (str, list, or tuple)        List of input FASTQ files or directories containing FASTQ files. Required. Positional argument on command line.
-  fraction (int or float)                 The fraction of reads to sample, as a float greater than 0. Any value equal to or greater than 1 forces sampling with replacement.
-  seeds (int, str, or list)               Random seed(s). Provide a single int (e.g. 42), an inclusive dash-delimited range string (e.g. "1-300" for seeds 1 through 300), or a list mixing ints and range strings (e.g. [1, 2, "5-7"]). On the command line, pass space-separated values (e.g. -s 1 2 5-7). Default: 42
-  output_dir (str)                        Output directory. Default: ./fastQpick_output
-  gzip_output (bool)                      Whether or not to gzip the output. Default: False (uncompressed)
-  group_size (int)                        The size of grouped files. Provide each pair of files sequentially, separated by a space. E.g., I1, R1, R2 -  would have group_size=3. Default: 1 (unpaired)
-  disable_replacement (bool)              Sample without replacement. By default (flag omitted), sampling is done with replacement.
-  overwrite (bool)                        Overwrite existing output files. Default: False
-  low_memory (bool)                       Whether to use low memory mode (uses ~5x less memory than default, but adds marginal time to the data structure-generation preprocessing). Has no effect in one_pass mode. Default: False
-  one_pass (bool)                         Use the single-pass approximate sampler. Skips the read-counting pass and draws each read's output multiplicity independently (Poisson with replacement, Bernoulli without). Runs in O(1) memory with roughly half the read I/O; the output size equals fraction*n only in expectation. Default: False
-  verbose (bool)                          Whether to print progress information. Default: True

---

## Features

- Efficient sampling of large FASTQ files.
- Works with both single and paired-end sequencing data.
- Supports sampling with or without replacement.
- Command-line interface and Python API for seamless integration.
- Memory efficient - the occurrence vector is sized to the largest per-read count actually drawn (one byte per read in the common case), and low-memory mode further avoids materializing the array of sampled indices.
- Time efficient - streams through the fastq and writes output in batches - generates a full-size (fraction=1, with replacement) bootstrap replicate of a 500M-read FASTQ in ~26 minutes in standard mode, ~56 minutes in low-memory mode, and ~35 minutes in one-pass mode (see [Benchmark](#benchmark) below).

## Sampling modes: exact (two-pass) vs. one-pass

fastQpick offers two samplers that trade exactness against resource usage.

**Exact (default).** Two passes are made over each file: a counting pass to learn the number of reads `n`, then a writing pass. Given `n`, exactly `floor(fraction * n)` reads are sampled (with or without replacement). This is the right choice when the output read count must be exact.

**One-pass (`--one_pass` / `one_pass=True`).** A single pass is made. As each read streams by, its number of copies in the output is drawn independently:
- with replacement: `Poisson(fraction)`,
- without replacement: `Bernoulli(fraction)` (the read is kept with probability `fraction`).

Because the draw for a read depends on neither its position nor the total read count, no read is favored over another, so the sample is unbiased and `n` never needs to be known. This eliminates the counting pass (roughly halving read I/O) and runs in O(1) memory (no occurrence vector is built). The cost is that the output size is `fraction * n` only in expectation, with relative standard deviation `1 / sqrt(fraction * n)` — about 0.03% for a 100M-read library at `fraction=0.1`. Paired/grouped files stay synchronized because all members of a group draw from the same per-group sub-seed.

Choosing between the three:

| Property | Exact, default | Exact, low-memory (`-l`) | One-pass (`-p`) |
| --- | --- | --- | --- |
| Passes over each file | 2 | 2 | 1 |
| Peak memory | index array of `m` + length-`n` count temporary, then ~1 byte/read | ~1 byte/read occurrence vector only | O(1) |
| Output size | exactly `floor(fraction*n)` | exactly `floor(fraction*n)` | `fraction*n` in expectation |

Both exact modes write from the same occurrence vector (~1 byte per read once built), so they use similar memory during the long writing pass. They differ at the moment the vector is built: the default draws all `m` indices at once and counts them with `np.bincount` (fast, but the index array plus the length-`n` counting temporary raise the transient peak), whereas low-memory streams indices one at a time into the vector (no index array, lower peak, ~5× lower end-to-end on a 500M-read bootstrap, at the cost of some speed). `low_memory` has no effect in one-pass mode, which is already O(1) memory; setting both simply prints a note.

## Benchmark

Generating one full-size bootstrap replicate (`fraction=1`, with replacement) of a 500-million-read uncompressed FASTQ (143 GB, 150 bp reads), single-threaded on an Intel Xeon Gold 6152 (2.10 GHz). Each mode reads its own file with the page cache evicted immediately beforehand, so every run starts cold.

| Mode | Wall time | Peak memory |
| --- | --- | --- |
| Exact, default | ~26 min | ~9.4 GB |
| Exact, low-memory (`-l`) | ~56 min | ~1.4 GB |
| One-pass (`-p`) | ~35 min | ~0.11 GB |

Peak memory scales with read count, not read length or file size. The default mode is fastest in wall-clock time here because the 143 GB file fits in the machine's memory, so its second (writing) pass re-reads the file from the page cache rather than from disk; the one-pass mode instead reads and writes concurrently in a single pass, moving less total data but contending for the disk. The one-pass time advantage is expected when the library exceeds available memory or the sampled fraction is small, while its ~0.11 GB footprint is constant regardless of read count. The low-memory mode is the most CPU-intensive of the three because it draws each read's multiplicity through the standard-library random generators one read at a time.

## Low memory mode vs. standard
Low memory mode vs. standard, when fraction=1 (i.e., number of reads to sample is the same as the number of reads in the fastq):
- Adds an extra ~3.5 seconds per million reads per group_size (i.e., a 500M-read FASTQ took ~56 minutes in low-memory mode vs ~26 minutes in standard mode)
- Saves ~16MB RAM per million reads (i.e., a 500M-read FASTQ used ~1.4GB RAM in low-memory mode vs ~9.4GB RAM in standard mode)

---

## Examples

### 1. Sample 10% of reads with replacement from a FASTQ file:

**Command-line**
```bash
fastQpick --fraction 0.1 input.fastq
```

**Python**
```python
from fastQpick import fastQpick

fastQpick(
    input_files='input.fastq',
    fraction=0.1
)
```

Sampling is done with replacement by default. Pass `--disable_replacement` (CLI) or `replacement=False` (Python) to sample without replacement.

### 2. Sample 100% of reads with replacement from multiple paired FASTQ files (R1, R2) across three seeds (i.e., bootstrapping):

**Command-line**
```bash
fastQpick --fraction 1 -s 42 43 44 -g 2 input1_R1.fastq input1_R2.fastq
```

**Python**
```python
from fastQpick import fastQpick

fastQpick(
    input_files='input.fastq',
    fraction=1,
    seeds=[42, 43, 44],
    replacement=True,
    group_size=2,
)
```

Seeds can also be given as inclusive dash-delimited ranges, which is convenient for many bootstrap replicates. For example, `-s 1-300` (or `seeds="1-300"`) runs seeds 1 through 300, and values can be mixed: `-s 1 5 10-12` (or `seeds=[1, 5, "10-12"]`) runs seeds 1, 5, 10, 11, and 12.

### 3. Sample ~10% of reads in a single pass (approximate output size, O(1) memory):

**Command-line**
```bash
fastQpick --fraction 0.1 --one_pass input.fastq
```

**Python**
```python
from fastQpick import fastQpick

fastQpick(
    input_files='input.fastq',
    fraction=0.1,
    one_pass=True,
)
```

The one-pass sampler skips the counting pass and draws each read's multiplicity on the fly, so it never needs to know the read count. The output contains `fraction * n` reads in expectation. See [Sampling modes](#sampling-modes-exact-two-pass-vs-one-pass) for the exact/one-pass trade-off.

---

## License

fastQpick is licensed under the 2-clause BSD license. See the [LICENSE](LICENSE) file for details.

---

## Contributing

We welcome contributions! Please see the [CONTRIBUTING.md](CONTRIBUTING.md) file for guidelines on how to get involved.

