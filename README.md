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
input_files (str, list, or tuple)       List of input FASTQ files or directories containing FASTQ files. Required. Positional argument on command line.
fraction (int or float)                 The fraction of reads to sample, as a float greater than 0. Any value equal to or greater than 1 will turn on the -r flag automatically.
seed (int or str)                       Random seed(s). Can provide multiple seeds separated by commas. Default: 42
output_dir (str)                        Output directory. Default: ./fastQpick_output
gzip_output (bool)                      Whether or not to gzip the output. Default: False (uncompressed)
group_size (int)                        The size of grouped files. Provide each pair of files sequentially, separated by a space. E.g., I1, R1, R2 would have group_size=3. Default: 1 (unpaired)
replacement (bool)                      Sample with replacement. Default: False (without replacement).
overwrite (bool)                        Overwrite existing output files. Default: False
low_memory (bool)                       Whether to use low memory mode (uses ~5.5x less memory than default, but adds marginal time to the data structure generation preprocessing). Default: False
verbose (bool)                          Whether to print progress information. Default: True

---

## Features

- Efficient sampling of large FASTQ files.
- Works with both single and paired-end sequencing data.
- Supports sampling with or without replacement.
- Command-line interface and Python API for seamless integration.
- Memory efficient - in low-memory mode, only uses as much memory as a list of integers the length of the number of reads in the fastq file for each file.
- Time efficient - only passes through the fastq once and writes to output in batches - can process 600M reads in 10-15 minutes

---

## Examples

1. Sample 10% of reads with replacement from a FASTQ file:

### Command-line Example
```bash
fastQpick --fraction 0.1 -r input.fastq
```

### Python Example
```python
from fastQpick import fastQpick

fastQpick(
    input_files='input.fastq',
    fraction=0.1,
    replacement=True
)
```

2. Sample 100% of reads with replacement from multiple paired FASTQ files (R1, R2) across three seeds (i.e., bootstrapping):

### Command-line Example
```bash
fastQpick --fraction 1 -s 42,43,44 -r -g 2 input1_R1.fastq input1_R2.fastq
```

### Python Example
```python
from fastQpick import fastQpick

fastQpick(
    input_files='input.fastq',
    fraction=1,
    seed="42,43,44",
    replacement=True,
    group_size=2,
)
```
---

## License

fastQpick is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Contributing

We welcome contributions! Please see the [CONTRIBUTING.md](CONTRIBUTING.md) file for guidelines on how to get involved.

