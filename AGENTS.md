# AGENTS.md

Guidance for AI coding agents working in this repository, and for agents that use fastQpick as a tool.

## Project

fastQpick is a Python CLI and library for fast, memory-efficient sampling of FASTQ reads, with replacement (bootstrap replicates, oversampling) or without (subsampling). It is designed for libraries of hundreds of millions of reads: it makes at most two streaming passes over each file and picks the smallest data structure that fits the sampling task.

## Setup and commands

```bash
pip install -e .              # editable install
pip install -e ".[mcp]"       # plus the MCP server (Python >= 3.10)
pytest tests/                 # full test suite (fast; fixtures are small)
pytest tests/test_fastQpick.py::test_single_file   # one test
python -m build               # wheel + sdist
fastQpick -f 0.1 -dr input.fastq                   # CLI
fastQpick-mcp                 # MCP server on stdio
```

There is no linter or formatter configured. Match the style of the surrounding code, including its comment density: comments explain why, not what.

## Layout

| Path | Contents |
|---|---|
| `fastQpick/main.py` | The whole pipeline: `fastQpick()` (single entry point for the CLI and the Python API), both samplers, the writers, and the argparse `main()`. |
| `fastQpick/utils.py` | Stateless helpers (read counting, seed parsing, file grouping, config snapshot). |
| `fastQpick/mcp_server.py` | MCP server (`fastQpick-mcp`), a thin wrapper that runs the CLI in a subprocess. |
| `fastQpick/__init__.py` | Package logger; pins `OPENBLAS_NUM_THREADS=1` before numpy is imported. |
| `tests/` | End-to-end tests on temporary FASTQ fixtures. |
| `benchmarks/` | Table 1 of the manuscript (`table1.py`, `summarize_table1.py`, `make_inputs.sh`). |
| `notebooks/` | Tutorials and Figure 1 (`notebooks/realdata/drive.sh`). |

## Architecture

Two samplers, selected by `single_pass`:

- Default (two-pass, exact). A counting pass learns the number of reads `n` (skippable with `read_counts`), then a writing pass emits exactly `floor(fraction * n)` reads. `make_occurrence_list` builds, for every read index, the number of times it is written; `write_fastq` streams the file with `pyfastx` and emits each read that many times.
- Single-pass (approximate). Each read's multiplicity is drawn independently as it streams, `Poisson(fraction)` with replacement or `Bernoulli(fraction)` without. Constant memory, one pass, can read standard input (`-`); the output size is exact only in expectation.

Parallelism: occurrence lists are built serially in the parent process, and only the writes are dispatched to a `ProcessPoolExecutor`. `threads` is a total budget split between worker processes and gzip deflate threads by `split_thread_budget`. Output never depends on `threads`.

## Invariants to preserve

Tests guard most of these. Read the relevant code before changing sampling logic.

- Grouped files (`file_group_size > 1`, e.g. R1/R2 or I1/R1/R2) share one occurrence list (two-pass) or one spawned `SeedSequence` sub-seed (single-pass), so mates stay synchronized. Grouped files must have equal read counts.
- `fraction >= 1` forces sampling with replacement.
- Reproducibility: output is byte-identical for a fixed seed. Two-pass mode uses one `np.random.default_rng(seed)` per seed, consumed in a fixed serial order. Do not use the global `np.random` or the stdlib `random` module.
- `make_occurrence_list` memory behavior:
  - The dense occurrence vector is filled block by block (`occurrence_block_size`), with block totals from `rng.multinomial` (with replacement) or `rng.hypergeometric` (without). Do not materialize the length-`m` index array or a length-`n` bincount temporary, and do not use `np.add.at` or `random.sample`.
  - A `Counter` is used instead of a dense array only when the sample is very sparse (`m < n / counter_sparsity_threshold`).
  - The dtype starts at `uint8` and widens only if a realized count exceeds it. Do not size it by `m`.
- User-supplied read counts are verified during the writing pass (`check_record_count`); a mismatch raises rather than silently biasing the sample.
- `save_params_to_config_file` inspects the caller's stack frame (`levels_up`). Adding or removing a wrapper frame around it silently snapshots the wrong arguments.
- `fastq_to_length_dict` is a module global that persists across in-process calls. This is why the MCP server runs each request in a subprocess.
- `fastQpick()` refuses a non-empty output directory unless `overwrite=True`.

## Tests

`tests/test_fastQpick.py` runs `fastQpick()` end to end and checks format, exact counts, uniqueness, mate synchronization, determinism, and occurrence-list dtype and uniformity. Single-pass tests check that the output size lies within 6 standard deviations of its mean rather than asserting an exact count. `tests/test_mcp_server.py` is skipped if `mcp` is not installed. Add a test for any behavior change, and run the full suite before finishing.

## Using fastQpick from an agent

With the MCP server registered (`claude mcp add fastqpick -- fastQpick-mcp`, or `{"mcpServers": {"fastqpick": {"command": "fastQpick-mcp"}}}`), three tools are available:

- `sample_fastq`: sample with or without replacement. Returns the output directory, the absolute paths of the written FASTQ files, the config file, and the end of the log.
- `count_fastq_reads`: read counts per file. Pass them to `sample_fastq` as `read_counts` when the same library is sampled more than once.
- `list_fastq_files`: the files fastQpick reads from a directory, in the order used for grouping.

Common requests and the corresponding arguments:

| Request | Arguments |
|---|---|
| One bootstrap replicate | `fraction=1` |
| 100 paired-end bootstrap replicates | `fraction=1, num_samples=100, file_group_size=2`, files ordered R1, R2 |
| 10% subsample | `fraction=0.1, without_replacement=True` |
| Train/test split | `fraction=0.8, without_replacement=True, oob=True` |
| Lowest memory, or a very large gzipped single sample | `single_pass=True` (output size approximate) |

Paths are resolved on the machine running the server, relative to its working directory. A full-size replicate of a 500-million-read gzipped library takes about 30 minutes, so warn the user before launching large jobs, and do not set `overwrite=True` unless the user asked to replace existing output.
