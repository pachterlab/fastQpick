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
-  low_memory (bool)                       Whether to use low memory mode (uses ~5.5x less memory than default, but adds marginal time to the data -  structure generation preprocessing). Default: False
-  verbose (bool)                          Whether to print progress information. Default: True

---

## Features

- Efficient sampling of large FASTQ files.
- Works with both single and paired-end sequencing data.
- Supports sampling with or without replacement.
- Command-line interface and Python API for seamless integration.
- Memory efficient - in low-memory mode, only uses as much memory as a list of (small) integers the length of the number of reads in the fastq file for each file.
- Time efficient - only passes through the fastq once and writes to output in batches - can process 600M reads in 10-15 minutes

## Low memory mode vs. standard
Low memory mode vs. standard, when fraction=1 (i.e., number of reads to sample is the same as the number of reads in the fastq):
- Adds an extra ~1-3 seconds per million reads per group_size (i.e., 500M reads would take 30 minutes instead of 20-25 minutes)
- Saves an extra ~40MB RAM per million reads (i.e., 500M reads would take 3.75GB RAM vs 20.6GB RAM)

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

---

## License

fastQpick is licensed under the 2-clause BSD license. See the [LICENSE](LICENSE) file for details.

---

## Contributing

We welcome contributions! Please see the [CONTRIBUTING.md](CONTRIBUTING.md) file for guidelines on how to get involved.

