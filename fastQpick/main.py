import argparse
import gzip
import os
import random
import numpy as np
from tqdm import tqdm
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, FIRST_COMPLETED, wait
from pydantic import ConfigDict, Field, validate_call
from typing import Union
import pyfastx  # to loop through fastq (faster than custom python code)

from fastQpick import logger
from fastQpick._version import __version__
from fastQpick.utils import save_params_to_config_file, is_directory_effectively_empty, group_items, count_reads, parse_seed

# Global variables
valid_fastq_extensions = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
batch_size = 200000  # for buffer
fastq_to_length_dict = {}  # set to empty, and the user can provide otherwise it will be calculated

def write_fastq(input_fastq, output_path, occurrence_list, total_reads, gzip_output, seed = None, unique_headers = False, verbose = True):
    if gzip_output:
        open_func = gzip.open
        write_mode = "wt"
    else:
        open_func = open
        write_mode = "w"
    
    buffer = []  # Temporary storage for the batch

    input_fastq_read_only = pyfastx.Fastx(input_fastq)

    # use tqdm if verbose else silently loop
    iterator = (
        tqdm(input_fastq_read_only, desc=f"Iterating through seed {seed}, file {input_fastq}", unit="read", total=total_reads)
        if verbose else input_fastq_read_only
    )
    
    with open_func(output_path, write_mode) as f:
        if not unique_headers:  # original (non-unique) headers
            for i, (name, seq, qual) in enumerate(iterator):
                # Add the FASTQ entry to the buffer
                buffer.extend([f"@{name}\n{seq}\n+\n{qual}\n"] * occurrence_list[i])
                
                # If the buffer reaches the batch size, write all at once and clear the buffer
                if (i + 1) % batch_size == 0:
                    f.writelines(buffer)
                    buffer.clear()  # Clear the buffer after writing
        else:  # unique headers
            for i, (name, seq, qual) in enumerate(iterator):
                if occurrence_list[i] > 0:  # not strictly necessary for coding logic, but saves time (if 0 > 0 is faster than saying for j in range(0))
                    buffer.extend([f"@{name}_{j}\n{seq}\n+\n{qual}\n" for j in range(1, occurrence_list[i]+1)])
                    
                # If the buffer reaches the batch size, write all at once and clear the buffer
                if (i + 1) % batch_size == 0:
                    f.writelines(buffer)
                    buffer.clear()
            
        # Write any remaining entries in the buffer
        if buffer:
            f.writelines(buffer)
            buffer.clear()

def occurrence_chunk_stream(rng, fraction, replacement, chunk_size):
    # One-pass sampling: lazily yield the number of times each successive read should appear in
    # the output, drawing variates in vectorized chunks of chunk_size for speed. The stream is
    # infinite and is consumed only as far as there are reads, so the file length never needs to
    # be known in advance. With replacement each read's multiplicity is Poisson(fraction); without
    # replacement each read is kept independently with probability fraction (a Bernoulli draw),
    # which can only appear 0 or 1 time. Because the draws depend on neither the read position nor
    # the total read count, no read is favored over another and the sample is unbiased.
    while True:
        if replacement:
            chunk = rng.poisson(fraction, size=chunk_size)
        else:
            chunk = (rng.random(chunk_size) < fraction).astype(np.uint8)
        for count in chunk:
            yield count

def write_fastq_one_pass(input_fastq, output_path, fraction, replacement, child_seed, gzip_output, seed=None, unique_headers=False, verbose=True):
    # Single-pass writer. A read's output multiplicity is drawn as the file streams by, so the read
    # count is never needed and peak memory is constant (only the flush buffer). All members of a
    # group are passed the same child_seed, so re-seeding here reproduces the identical multiplicity
    # sequence for every member and keeps mate pairs synchronized without a shared occurrence vector.
    if gzip_output:
        open_func = gzip.open
        write_mode = "wt"
    else:
        open_func = open
        write_mode = "w"

    rng = np.random.default_rng(child_seed)
    occurrence_stream = occurrence_chunk_stream(rng, fraction, replacement, batch_size)

    buffer = []  # Temporary storage for the batch

    input_fastq_read_only = pyfastx.Fastx(input_fastq)

    # total is unknown without a counting pass, so tqdm shows throughput rather than a percentage
    iterator = (
        tqdm(input_fastq_read_only, desc=f"Iterating through seed {seed}, file {input_fastq}", unit="read")
        if verbose else input_fastq_read_only
    )

    with open_func(output_path, write_mode) as f:
        if not unique_headers:  # original (non-unique) headers
            for i, ((name, seq, qual), count) in enumerate(zip(iterator, occurrence_stream)):
                if count:
                    buffer.append(f"@{name}\n{seq}\n+\n{qual}\n" * count)

                if (i + 1) % batch_size == 0:
                    f.writelines(buffer)
                    buffer.clear()
        else:  # unique headers
            for i, ((name, seq, qual), count) in enumerate(zip(iterator, occurrence_stream)):
                if count > 0:
                    buffer.extend([f"@{name}_{j}\n{seq}\n+\n{qual}\n" for j in range(1, count+1)])

                if (i + 1) % batch_size == 0:
                    f.writelines(buffer)
                    buffer.clear()

        # Write any remaining entries in the buffer
        if buffer:
            f.writelines(buffer)
            buffer.clear()

def smallest_uint_dtype(max_value):
    # Smallest unsigned numpy integer type that can hold max_value. Occurrence counts are tiny
    # (~Poisson(fraction)), so this is almost always uint8, which is what makes the dense occurrence
    # vector cost ~1 byte per read rather than 4-8.
    for dtype in (np.uint8, np.uint16, np.uint32):
        if max_value <= np.iinfo(dtype).max:
            return dtype
    return np.uint64

# A Counter entry costs on the order of 100 bytes (dict slot plus two Python int objects), whereas a
# dense occurrence vector costs total_reads * dtype_bytes and the per-read count almost always fits
# in one byte. The sparse representation therefore only saves memory when the number of sampled reads
# is a small fraction of the file; below this ratio the Counter (which also avoids allocating the
# length-n bincount temporary) wins, above it the dense array wins.
counter_sparsity_threshold = 100

def make_occurrence_list(file, seed, total_reads, number_of_reads_to_sample, replacement, low_memory, rng=None, verbose=True):
    if verbose:
        logger.info(f"Calculating total reads and determining random indices for seed {seed}, file {file}")

    use_counter = number_of_reads_to_sample < (total_reads / counter_sparsity_threshold)

    if low_memory:
        if replacement:
            random_indices = (random.choice(range(total_reads)) for _ in range(number_of_reads_to_sample))
        else:
            random_indices = (index for index in random.sample(range(total_reads), k=number_of_reads_to_sample))

        # Count occurrences. random_indices is a lazy generator here, so it cannot be
        # handed to np.bincount (which needs a materialized array); accumulate instead.
        if use_counter:
            occurrence_list = Counter(random_indices)
        else:
            # Without replacement each index is distinct, so the count is exactly one and uint8 is
            # provably sufficient. With replacement the count is ~Poisson(fraction); uint16 leaves an
            # enormous safety margin (a single read would need to be drawn 65,535 times to overflow).
            dtype_occurences_list = np.uint16 if replacement else np.uint8
            occurrence_list = np.zeros(total_reads, dtype=dtype_occurences_list)
            for index in random_indices:
                occurrence_list[index] += 1
            # Downcast to the realized maximum (almost always uint8) so the vector that persists
            # through the write pass is as small as in the default path.
            if replacement:
                realized_dtype = smallest_uint_dtype(int(occurrence_list.max()) if occurrence_list.size else 0)
                if realized_dtype != occurrence_list.dtype:
                    occurrence_list = occurrence_list.astype(realized_dtype)
    else:
        # rng is a seeded numpy Generator threaded down from sample_multiple_files so that the default
        # path is reproducible (the legacy global np.random it previously used was never seeded).
        if replacement:
            # Drawing integers directly avoids materializing a length-n arange just to index into it.
            dtype_random_indices = np.uint32 if total_reads <= np.iinfo(np.uint32).max else np.uint64
            random_indices = rng.integers(0, total_reads, size=number_of_reads_to_sample, dtype=dtype_random_indices)
        else:
            # Sampling without replacement requires a permutation, which numpy materializes internally.
            random_indices = rng.choice(total_reads, size=number_of_reads_to_sample, replace=False)

        # Count occurrences
        if use_counter:
            occurrence_list = Counter(random_indices.tolist())
        else:
            counts = np.bincount(random_indices, minlength=total_reads)
            # Size the occurrence vector to the realized maximum count rather than to the sample size,
            # which over-provisions by 4-8x; the maximum is almost always small enough for uint8.
            dtype_occurences_list = smallest_uint_dtype(int(counts.max()) if counts.size else 0)
            occurrence_list = counts.astype(dtype_occurences_list)
            del counts

    del random_indices

    return occurrence_list

def insert_seed_suffix(filename, seed):
    # Insert a per-seed marker before the FASTQ extension, e.g. "R1.fastq.gz" -> "R1.seed3.fastq.gz",
    # so that distinct seeds write to distinct output files instead of overwriting one another.
    for ext in sorted(valid_fastq_extensions, key=len, reverse=True):
        if filename.endswith(ext):
            return f"{filename[:-len(ext)]}.seed{seed}{ext}"
    return f"{filename}.seed{seed}"

def resolve_output_path(file, output_directory, gzip_output, multiple_seeds, seed):
    # Build the per-file output path, disambiguating by seed when needed and forcing the
    # extension to match the requested gzip setting.
    output_basename = os.path.basename(file)
    if multiple_seeds:  # disambiguate output files when more than one seed is sampled
        output_basename = insert_seed_suffix(output_basename, seed)
    output_path = os.path.join(output_directory, output_basename)
    if output_directory:
        os.makedirs(output_directory, exist_ok=True)

    if gzip_output and not output_path.endswith(".gz"):
        output_path += ".gz"
    elif not gzip_output and output_path.endswith(".gz"):
        output_path = output_path[:-3]
    return output_path

def two_pass_write_jobs(file_list, fraction, seed, output, gzip_output, replacement, low_memory, unique_headers, multiple_seeds, rng, verbose):
    # Lazily yield (callable, kwargs) write jobs for one seed. Each group's occurrence list is built
    # HERE, in the parent process, so the shared `rng` is consumed in exactly the same order as the
    # single-process path; only the I/O-bound write_fastq call is handed off as a job. This keeps the
    # default-path output byte-identical to the serial version while letting the writes run
    # concurrently across files. Members of a group share one occurrence list, computed once.
    for file in file_list:
        files_total = (file, ) if isinstance(file, str) else file

        total_reads = fastq_to_length_dict[files_total[0]]
        number_of_reads_to_sample = int(fraction * total_reads)

        occurrence_list = make_occurrence_list(file=files_total[0], seed=seed, total_reads=total_reads, number_of_reads_to_sample=number_of_reads_to_sample, replacement=replacement, low_memory=low_memory, rng=rng, verbose=verbose)

        for member in files_total:
            output_path = resolve_output_path(member, output, gzip_output, multiple_seeds, seed)
            yield write_fastq, dict(input_fastq=member, output_path=output_path, occurrence_list=occurrence_list, total_reads=total_reads, gzip_output=gzip_output, seed=seed, unique_headers=unique_headers, verbose=verbose)

def bootstrap_single_file_one_pass(files_total = None, child_seed = None, gzip_output = None, output_directory = None, seed = None, fraction = None, replacement = None, unique_headers = False, multiple_seeds = False, verbose=True):
    # One-pass counterpart of bootstrap_single_file. No occurrence vector is materialized and the
    # file length is never counted; every member of the group is written from the same child_seed
    # so their sampled multiplicities match read-for-read, keeping mate pairs synchronized.
    if isinstance(files_total, str):
        files_total = (files_total, )

    for file in files_total:
        output_path = resolve_output_path(file, output_directory, gzip_output, multiple_seeds, seed)

        write_fastq_one_pass(input_fastq = file, output_path = output_path, fraction = fraction, replacement = replacement, child_seed = child_seed, gzip_output = gzip_output, seed = seed, unique_headers = unique_headers, verbose = verbose)

def run_write_jobs(jobs, max_workers):
    # Execute (callable, kwargs) write jobs. The write stage is dominated by I/O (read, optional
    # compression, write), so with more than one file the jobs are dispatched to a process pool to
    # run concurrently; with a single file they run inline to avoid pool overhead. A sliding window
    # caps the number of in-flight jobs at max_workers, so the parent never materializes more than
    # max_workers occurrence lists at once and peak memory stays bounded.
    if max_workers <= 1:
        for func, kwargs in jobs:
            func(**kwargs)
        return

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        in_flight = set()
        for func, kwargs in jobs:
            while len(in_flight) >= max_workers:
                done, in_flight = wait(in_flight, return_when=FIRST_COMPLETED)
                for future in done:
                    future.result()  # surface any worker exception
            in_flight.add(executor.submit(func, **kwargs))
        for future in in_flight:
            future.result()

def sample_multiple_files(file_list, fraction, seed_list, output, gzip_output, replacement, low_memory, unique_headers, one_pass, verbose):
    multiple_seeds = len(seed_list) > 1
    cpu_count = os.cpu_count() or 1
    # Number of individual output files (groups expanded); this is the most write jobs that can run
    # at once and so caps the worker count.
    num_individual_files = sum(len(file) if isinstance(file, tuple) else 1 for file in file_list)
    for seed in seed_list:
        if one_pass:
            # Derive one independent sub-seed per group from the master seed. Members of a group
            # share their sub-seed (handled inside bootstrap_single_file_one_pass) so mates stay
            # synchronized, while different groups draw independently, mirroring the two-pass path
            # where each file consumes fresh random state. Each group is fully independent, so a whole
            # group is one parallel write job.
            child_seeds = np.random.SeedSequence(seed).spawn(len(file_list))
            jobs = (
                (bootstrap_single_file_one_pass,
                 dict(files_total=file, child_seed=child_seed, gzip_output=gzip_output, output_directory=output, seed=seed, fraction=fraction, replacement=replacement, unique_headers=unique_headers, multiple_seeds=multiple_seeds, verbose=verbose))
                for file, child_seed in zip(file_list, child_seeds)
            )
            max_workers = min(len(file_list), cpu_count)
        else:
            # Seed both RNGs once per seed: stdlib random drives the low-memory path, and a numpy
            # Generator drives the default path. A single Generator is shared across the seed's files
            # so that successive files draw independent samples while remaining reproducible. The
            # occurrence lists are built serially in this process (inside two_pass_write_jobs), so the
            # shared RNG is consumed in order and the output is independent of how the writes are
            # scheduled across processes.
            random.seed(seed)
            rng = np.random.default_rng(seed)
            jobs = two_pass_write_jobs(file_list, fraction, seed, output, gzip_output, replacement, low_memory, unique_headers, multiple_seeds, rng, verbose)
            max_workers = min(num_individual_files, cpu_count)

        run_write_jobs(jobs, max_workers=max_workers)
    
def make_fastq_to_length_dict(file_list, verbose=True):
    global fastq_to_length_dict
    for file in file_list:
        if isinstance(file, tuple):
            if all(specific_file in fastq_to_length_dict for specific_file in file):
                continue
            if verbose:
                logger.info(f"Counting {file[0]}")
            count = count_reads(file[0])
            for i in range(len(file)):
                fastq_to_length_dict[file[i]] = count
        elif isinstance(file, str):
            if file in fastq_to_length_dict:
                continue
            if verbose:
                logger.info(f"Counting {file}")
            count = count_reads(file)
            fastq_to_length_dict[file] = count
    if verbose:
        logger.info(f"fastq_to_length_dict: {fastq_to_length_dict}")

@validate_call
def fastQpick(
    input_files: str | list | tuple,
    fraction: float = 1.0,
    seeds: int | str | list = 42,
    output_dir: str = "fastQpick_output",
    gzip_output: bool = False,
    file_group_size: int = 1,
    replacement: bool = True,
    overwrite: bool = False,
    low_memory: bool = False,
    unique_headers: Union[bool, None] = None,
    one_pass: bool = False,
    verbose: bool = True,
    **kwargs
):
    """
    Fast and memory-efficient sampling of DNA-Seq or RNA-seq fastq data with or without replacement.

    Parameters
    ----------
    input_files (str, list, or tuple)       Input FASTQ files or directories containing FASTQ files.
    fraction (int or float)                 The fraction of reads to sample, as a float greater than 0. Any value equal to or greater than 1 turns on sampling with replacement automatically.
    seeds (int, str, or list)               Random seed(s). Provide a single int (e.g. 42), an inclusive dash-delimited range string (e.g. "1-300" for seeds 1 through 300), or a list mixing ints and range strings (e.g. [1, 2, "5-7"]).
    output_dir (str)                        Output directory.
    gzip_output (bool)                      Gzip the output files.
    file_group_size (int)                        The size of grouped files. Provide each pair of files sequentially, separated by a space. E.g., I1, R1, R2 would have file_group_size=3.
    replacement (bool)                      Sample with replacement. Automatically enabled if fraction >= 1.
    overwrite (bool)                        Overwrite existing output files.
    low_memory (bool)                       Activate low memory mode. Disabled if one_pass is True.
    unique_headers (bool)                   Add a unique identifier to the header names of the output files. Default True if replacement is True, False if replacement is False.
    one_pass (bool)                         Use the single-pass approximate sampler. Uses low memory and runs faster in some cases (when the fraction is small or the input is large compared to available memory).
    verbose (bool)                          Whether to print progress information.

    kwargs
    ------
    fastq_to_length_dict (dict)             Dictionary of FASTQ file paths to number of reads in each file. If not provided, will be calculated.
    """
    # check if fastq_to_length_dict is in kwargs
    if "fastq_to_length_dict" in kwargs and isinstance(kwargs["fastq_to_length_dict"], dict):
        global fastq_to_length_dict
        fastq_to_length_dict = kwargs["fastq_to_length_dict"]

    # Check overwrite
    if not overwrite:
        if os.path.exists(output_dir) and not is_directory_effectively_empty(output_dir):  # check if dir exists and is not empty
            raise FileExistsError(f"Output directory '{output_dir}' already exists. Please specify a different output directory or set the overwrite flag to True.")

    seeds = parse_seed(seeds)  # normalize int / range string / iterable into a flat list of integer seeds (also keeps the config snapshot JSON-serializable)

    # Save arguments to a config file
    os.makedirs(output_dir, exist_ok=True)
    config_file = os.path.join(output_dir, "fastQpick_config.json")
    save_params_to_config_file(config_file)

    # type checking
    # if fraction >= 1, set replacement to True
    if float(fraction) >= 1.0:
        replacement = True

    # low_memory only governs the two-pass occurrence-vector construction; the one-pass sampler
    # never builds that vector, so the flag has nothing to act on there.
    if one_pass and low_memory and verbose:
        logger.warning("low_memory has no effect in one_pass mode, which already runs in O(1) memory.")

    # go through files, and only keep those that are valid fastq files or that are a folder containing valid fastq files in the direct subdirectory
    if isinstance(input_files, str):
        input_files = [input_files]
    elif not isinstance(input_files, (tuple, list)):
        raise ValueError("Input file list must be a string, tuple of strings, or list of strings.")

    input_files_parsed = []
    for path in input_files:
        if not isinstance(path, str):
            raise ValueError("Input file list must be a string, tuple of strings, or list of strings.")
        if not os.path.exists(path):
            raise FileNotFoundError(f"File or directory '{path}' not found.")
        elif os.path.isdir(path):
            input_files_before_path = len(input_files_parsed)
            # sorted() makes directory expansion deterministic, which matters for grouped/paired files
            for subpath in sorted(os.listdir(path)):
                full_subpath = os.path.join(path, subpath)
                if os.path.isfile(full_subpath) and subpath.endswith(tuple(valid_fastq_extensions)):
                    input_files_parsed.append(full_subpath)
            if len(input_files_parsed) == input_files_before_path:
                raise ValueError(f"No valid FASTQ files found in directory '{path}'.")
        elif os.path.isfile(path) and not path.endswith(tuple(valid_fastq_extensions)):
            raise ValueError(f"File '{path}' is not a valid FASTQ file.")
        elif os.path.isfile(path) and path.endswith(tuple(valid_fastq_extensions)):
            input_files_parsed.append(path)

    file_group_size = int(file_group_size)  # make sure file_group_size is an int (not a string)
    fraction = float(fraction)  # make sure fraction is a float (not a string)

    if file_group_size > 1:
        input_files_parsed = group_items(input_files_parsed, group_size=file_group_size)
    
    if unique_headers is None:
        unique_headers = replacement  # default to True if replacement is True, False if replacement is False
    if replacement and not unique_headers:
        logger.warning(f"unique_headers is {unique_headers} but replacement is {replacement}. This may lead to duplicate header names in the output files, which can cause issues for downstream tools. Consider setting unique_headers to {replacement} to match the replacement setting.")
    if not replacement and unique_headers:
        logger.warning(f"unique_headers is {unique_headers} but replacement is {replacement}. This may lead to unnecessarily verbose header names in the output files. Consider setting unique_headers to {replacement} to match the replacement setting.")
    
    # Count reads in each file and store in a dictionary. The one-pass sampler does not need the
    # counts, so this pass is skipped entirely in that mode.
    if not one_pass:
        make_fastq_to_length_dict(input_files_parsed, verbose=verbose)

    # Do the sampling
    sample_multiple_files(file_list=input_files_parsed, fraction=fraction, seed_list=seeds, output=output_dir, gzip_output=gzip_output, replacement=replacement, low_memory=low_memory, unique_headers=unique_headers, one_pass=one_pass, verbose=verbose)

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description="Fast and memory-efficient sampling of DNA-Seq or RNA-seq fastq data with or without replacement.")
    parser.add_argument("-f", "--fraction", required=True, default=False, help="The fraction of reads to sample, as a float greater than 0. Any value equal to or greater than 1 turns on sampling with replacement automatically.")
    parser.add_argument("-s", "--seeds", required=False, default=42, nargs="+", help='Random seed(s). Provide a single int (e.g. 42), an inclusive dash-delimited range string (e.g. "1-300" for seeds 1 through 300), or a list mixing ints and range strings (e.g. [1, 2, "5-7"]).')
    parser.add_argument("-o", "--output-dir", required=False, type=str, default="fastQpick_output", help="Output directory.")
    parser.add_argument("-z", "--gzip-output", required=False, default=False, help="Gzip the output files.")
    parser.add_argument("-g", "--file-group-size", required=False, default=1, help="The size of grouped files. Provide each pair of files sequentially, separated by a space. E.g., I1, R1, R2 would have file_group_size=3.")
    parser.add_argument("-dr", "--disable-replacement", action="store_true", help="Sample with replacement. Automatically enabled if fraction >= 1.")
    parser.add_argument("-w", "--overwrite", action="store_true", help="Overwrite existing output files.")
    parser.add_argument("-l", "--low-memory", action="store_true", help="Activate low memory mode. Disabled if one_pass is True.")
    parser.add_argument("-u", "--unique-headers", action="store_true", default=None, help="Add a unique identifier to the header names of the output files. Defaults to True when sampling with replacement, False otherwise.")
    parser.add_argument("-p", "--one_pass", action="store_true", help="Use the single-pass approximate sampler. Uses low memory and runs faster in some cases (when the fraction is small or the input is large compared to available memory).")
    parser.add_argument("-q", "--quiet", action="store_false", help="Whether to print progress information.")
    parser.add_argument("-v", "--version", action="version", version=f"fastQpick {__version__}", help="Show program's version number and exit")

    # Positional argument for input files (indefinite number)
    parser.add_argument("input_files", nargs="+", help="Input FASTQ files or directories containing FASTQ files.")

    # Parse arguments
    args = parser.parse_args()
            
    fastQpick(input_files=args.input_files,
              fraction=args.fraction,
              seeds=args.seeds,
              output_dir=args.output_dir,
              gzip_output=args.gzip_output,
              file_group_size=args.file_group_size,
              replacement=not args.disable_replacement,
              overwrite=args.overwrite,
              low_memory=args.low_memory,
              unique_headers=args.unique_headers,
              one_pass=args.one_pass,
              verbose=args.quiet)
