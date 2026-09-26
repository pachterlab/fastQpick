"""Model Context Protocol (MCP) server exposing fastQpick to LLM agents.

Run with ``fastQpick-mcp`` (stdio transport) after ``pip install "fastQpick[mcp]"``.

Each sampling request runs the fastQpick CLI in a fresh subprocess rather than calling
``fastQpick()`` in-process, for three reasons: the stdio transport owns this process's stdout,
which must carry only protocol messages; ``fastQpick()`` keeps read counts in the module-global
``fastq_to_length_dict``, which would otherwise persist between requests; and a long run does not
block the server's event loop.
"""
import asyncio
import multiprocessing
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from typing import List, Optional, Union

from fastQpick._version import __version__

try:  # mcp >= 2
    from mcp.server.mcpserver import MCPServer
    server_version_kwargs = {"version": __version__}
except ImportError:  # mcp 1.x, whose FastMCP does not take a version
    from mcp.server.fastmcp import FastMCP as MCPServer
    server_version_kwargs = {}

from fastQpick.main import valid_fastq_extensions, STDIN_SENTINEL
from fastQpick.utils import count_reads as _count_reads_in_file

log_tail_lines = 40  # lines of fastQpick's stderr returned with each result

server = MCPServer(
    name="fastQpick",
    instructions=(
        "fastQpick samples reads from FASTQ files, with replacement (bootstrap replicates, "
        "oversampling) or without replacement (subsampling). Paths are resolved on the machine "
        "running this server. For paired-end or multi-read libraries (e.g. R1/R2 or I1/R1/R2), "
        "list the mates consecutively and set file_group_size to the number of files per library "
        "so that mates stay synchronized; grouped files must have equal read counts. Use "
        "count_fastq_reads first when the same library will be sampled repeatedly, and pass the "
        "counts to sample_fastq as read_counts to skip its counting pass. Large libraries can take "
        "tens of minutes per replicate."
    ),
    **server_version_kwargs,
)


def _check_input_paths(input_files):
    for path in input_files:
        if path == STDIN_SENTINEL:
            raise ValueError("Standard input ('-') is not available through the MCP server; pass a file path.")
        if not os.path.exists(path):
            raise FileNotFoundError(f"File or directory '{path}' not found.")


def _list_output_files(output_dir):
    output_files = []
    for root, _, files in os.walk(output_dir):
        for name in sorted(files):
            if name.endswith(valid_fastq_extensions):
                output_files.append(os.path.abspath(os.path.join(root, name)))
    return sorted(output_files)


def build_cli_args(input_files, fraction, output_dir, num_samples, seed, without_replacement,
                   file_group_size, single_pass, disable_gzip, collapse_duplicates, oob,
                   no_unique_headers, read_counts, threads, overwrite):
    args = ["-f", str(fraction), "-o", output_dir, "-n", str(num_samples), "-g", str(file_group_size)]
    if seed is not None:
        args += ["-s", str(seed)]
    if without_replacement:
        args.append("-dr")
    if single_pass:
        args.append("-p")
    if disable_gzip:
        args.append("-z")
    if collapse_duplicates:
        args.append("-c")
    if oob:
        args.append("--oob")
    if no_unique_headers:
        args.append("--no-unique-headers")
    if read_counts:
        args += ["--read-counts", ",".join(str(int(count)) for count in read_counts)]
    if threads is not None:
        args += ["-t", str(threads)]
    if overwrite:
        args.append("-w")
    # "--" keeps input paths that start with "-" from being parsed as options
    return args + ["--"] + list(input_files)


@server.tool()
async def sample_fastq(
    input_files: List[str],
    fraction: float = 1.0,
    output_dir: str = "fastQpick_output",
    num_samples: int = 1,
    seed: Union[int, str] = 42,
    without_replacement: bool = False,
    file_group_size: int = 1,
    single_pass: bool = False,
    disable_gzip: bool = False,
    collapse_duplicates: bool = False,
    oob: bool = False,
    no_unique_headers: bool = False,
    read_counts: Optional[List[int]] = None,
    threads: Optional[int] = None,
    overwrite: bool = False,
) -> dict:
    """Sample reads from FASTQ files with or without replacement and write the samples to output_dir.

    Args:
        input_files: FASTQ files (.fastq, .fq, optionally .gz) or directories containing them.
        fraction: Fraction of reads to sample (> 0). Values >= 1 always sample with replacement;
            fraction=1 with replacement is a standard bootstrap replicate.
        output_dir: Output directory. Must be empty or absent unless overwrite is true.
        num_samples: Number of independent replicates. Replicate i uses seed + i and is written as
            <name>.seed<seed>.fastq[.gz] when more than one replicate is requested.
        seed: Random seed, or a comma-separated list / dash range string (e.g. "1-10"), which
            then sets the number of replicates. Output is reproducible for a fixed seed.
        without_replacement: Subsample without replacement (ignored when fraction >= 1).
        file_group_size: Number of consecutive input files that form one library (2 for R1/R2,
            3 for I1/R1/R2). Mates are sampled at the same read indices.
        single_pass: Read each input once with constant memory; the output size is then exact
            only in expectation.
        disable_gzip: Write plain FASTQ instead of gzip.
        collapse_duplicates: Write each sampled read once with ";size=<count>" in its header.
        oob: Also write the reads not drawn into <name>.oob.fastq[.gz].
        no_unique_headers: Keep original headers for repeated reads instead of adding _1, _2, ...
        read_counts: Known read count per file (or per group), in input order, to skip the counting
            pass. A wrong count raises an error rather than biasing the sample.
        threads: Total thread budget (default: 4, or fewer if fewer cores are available).
        overwrite: Allow writing into a non-empty output_dir.

    Returns:
        The output directory, the sampled FASTQ files written there, and the tail of the fastQpick log.
    """
    _check_input_paths(input_files)
    if fraction <= 0:
        raise ValueError(f"fraction must be greater than 0, got {fraction}.")
    cli_args = build_cli_args(input_files, fraction, output_dir, num_samples, seed, without_replacement,
                              file_group_size, single_pass, disable_gzip, collapse_duplicates, oob,
                              no_unique_headers, read_counts, threads, overwrite)
    process = await asyncio.create_subprocess_exec(
        sys.executable, "-c", "from fastQpick.main import main; main()", *cli_args,
        stdin=asyncio.subprocess.DEVNULL,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
    )
    try:
        stdout, stderr = await process.communicate()
    except asyncio.CancelledError:  # the client cancelled the request: do not leave the job running
        process.kill()
        await process.wait()
        raise
    log = (stdout + stderr).decode(errors="replace").replace("\r", "\n").splitlines()
    log_tail = "\n".join(line for line in log if line.strip())
    log_tail = "\n".join(log_tail.splitlines()[-log_tail_lines:])
    if process.returncode != 0:
        raise RuntimeError(f"fastQpick exited with status {process.returncode}:\n{log_tail}")
    return {
        "output_dir": os.path.abspath(output_dir),
        "output_files": _list_output_files(output_dir),
        "config_file": os.path.join(os.path.abspath(output_dir), "fastQpick_config.json"),
        "log_tail": log_tail,
    }


@server.tool()
async def count_fastq_reads(input_files: List[str], threads: int = 4) -> dict:
    """Count the reads in each FASTQ file, in parallel.

    The counts can be passed to sample_fastq as read_counts (in the same order) to skip its
    counting pass, which saves one full read of each file per call.

    Args:
        input_files: FASTQ files (.fastq, .fq, optionally .gz).
        threads: Number of files counted at once.

    Returns:
        A mapping from each file path to its number of reads.
    """
    _check_input_paths(input_files)
    for path in input_files:
        if not os.path.isfile(path) or not path.endswith(valid_fastq_extensions):
            raise ValueError(f"'{path}' is not a FASTQ file (expected one of {', '.join(valid_fastq_extensions)}).")
    loop = asyncio.get_running_loop()
    # "spawn" rather than fork: forking a process that runs an event loop can copy held locks
    with ProcessPoolExecutor(max_workers=max(1, min(threads, len(input_files))),
                             mp_context=multiprocessing.get_context("spawn")) as executor:
        counts = await asyncio.gather(*(loop.run_in_executor(executor, _count_reads_in_file, path) for path in input_files))
    return dict(zip(input_files, counts))


@server.tool()
def list_fastq_files(directory: str) -> List[str]:
    """List the FASTQ files that fastQpick would read from a directory, in the order it reads them.

    This is the order used for file grouping, so check it before setting file_group_size on a
    directory input: mates must be adjacent (e.g. sample_R1, sample_R2, next_R1, next_R2).

    Args:
        directory: Directory to scan (not recursive).
    """
    if not os.path.isdir(directory):
        raise NotADirectoryError(f"'{directory}' is not a directory.")
    return [os.path.join(directory, name) for name in sorted(os.listdir(directory))
            if os.path.isfile(os.path.join(directory, name)) and name.endswith(valid_fastq_extensions)]


def main():
    server.run()


if __name__ == "__main__":
    main()
