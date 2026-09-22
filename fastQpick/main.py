import argparse
import contextlib
import io
import itertools
import os
import random
import sys
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

try:
    from isal import igzip as gzip
    gzip_compresslevel = 1
except ImportError:
    import gzip
    gzip_compresslevel = 6

try:
    # Multithreaded deflate. Writing a full-size replicate pushes the whole library back
    # through the compressor, and a single ISA-L thread tops out near 400 MB/s, which makes
    # compression the bottleneck rather than the sampling. Four threads reach ~1.6 GB/s,
    # comfortably above the rate at which the reader can supply records.
    from isal import igzip_threaded
except ImportError:
    igzip_threaded = None

# Global variables
valid_fastq_extensions = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
batch_size = 200000  # for buffer
fastq_to_length_dict = {}  # set to empty, and the user can provide otherwise it will be calculated
gzip_output_threads = 4  # deflate threads when igzip_threaded is available
STDIN_SENTINEL = "-"  # read the library from standard input instead of from a file

def open_output(path, gzip_output):
    # Text-mode output handle using the fastest gzip writer available.
    if not gzip_output:
        return open(path, "w")
    if igzip_threaded is not None:
        return igzip_threaded.open(path, "wt", compresslevel=gzip_compresslevel,
                                   threads=gzip_output_threads)
    return gzip.open(path, "wt", compresslevel=gzip_compresslevel)

def stream_fastq_records(fh, chunk_size=1 << 22):
    # Yield (name, seq, qual) from a non-seekable stream, matching what pyfastx.Fastx yields
    # for a file. pyfastx is not used here because it silently drops the final record when
    # its input is a pipe rather than a file. A gzip stream is decompressed on the fly, so
    # both "zcat f.gz |" and "cat f.gz |" work.
    head = fh.read(2)
    if head[:2] == b"\x1f\x8b":
        fh = gzip.open(_Prepended(head, fh), "rb")
        head = b""
    pending = []  # complete lines not yet consumed as a record
    tail = head   # partial trailing line carried between chunks
    while True:
        chunk = fh.read(chunk_size)
        if not chunk:
            break
        lines = (tail + chunk).split(b"\n")
        tail = lines.pop()
        pending.extend(lines)
        complete = len(pending) - len(pending) % 4
        for i in range(0, complete, 4):
            yield (pending[i][1:].decode(), pending[i + 1].decode(), pending[i + 3].decode())
        del pending[:complete]
    if tail:
        pending.append(tail)
    if len(pending) % 4:
        raise ValueError(f"Truncated FASTQ record on the input stream: {len(pending)} "
                         "trailing lines do not form a complete four-line record.")
    for i in range(0, len(pending), 4):
        yield (pending[i][1:].decode(), pending[i + 1].decode(), pending[i + 3].decode())

class _Prepended(io.RawIOBase):
    # Pushes already-consumed magic bytes back in front of a non-seekable stream, so that the
    # gzip reader sees the header it needs without the stream having to be rewound.
    def __init__(self, head, stream):
        self._head, self._stream = head, stream
    def readable(self):
        return True
    def readinto(self, buf):
        if self._head:
            n = min(len(buf), len(self._head))
            buf[:n], self._head = self._head[:n], self._head[n:]
            return n
        data = self._stream.read(len(buf))
        buf[:len(data)] = data
        return len(data)

def open_input_stream(input_fastq):
    # Records from a file (indexed reader) or from standard input (streaming reader).
    if input_fastq == STDIN_SENTINEL:
        return stream_fastq_records(sys.stdin.buffer)
    return pyfastx.Fastx(input_fastq)

def write_reads_tagged(read_count_pairs, f, f_oob, unique_headers, collapse_duplicates):
    # General writer used when collapse_duplicates and/or oob is requested (the plain writers below
    # keep their specialized loops, since they are the hot path). read_count_pairs yields
    # ((name, seq, qual), count). With collapse_duplicates each sampled read is written once and its
    # multiplicity is recorded in the header as ";size=<count>" (the USEARCH/VSEARCH abundance
    # convention), so the multiset is preserved without physically duplicating records. With f_oob,
    # the out-of-bag reads (count == 0) are written, unmodified, to that second handle.
    buffer = []
    oob_buffer = []
    i = -1  # stays -1 for an empty input
    for i, ((name, seq, qual), count) in enumerate(read_count_pairs):
        if count:
            if collapse_duplicates:
                buffer.append(f"@{name};size={count}\n{seq}\n+\n{qual}\n")
            elif unique_headers:
                buffer.extend([f"@{name}_{j}\n{seq}\n+\n{qual}\n" for j in range(1, count+1)])
            else:
                buffer.append(f"@{name}\n{seq}\n+\n{qual}\n" * int(count))
        elif f_oob is not None:
            oob_buffer.append(f"@{name}\n{seq}\n+\n{qual}\n")

        if (i + 1) % batch_size == 0:
            f.writelines(buffer)
            buffer.clear()
            if oob_buffer:
                f_oob.writelines(oob_buffer)
                oob_buffer.clear()

    f.writelines(buffer)
    if oob_buffer:
        f_oob.writelines(oob_buffer)
    return i + 1

def write_fastq(input_fastq, output_path, occurrence_list, total_reads, gzip_output, seed = None, unique_headers = False, collapse_duplicates = False, oob_path = None, verbose = True):
    open_func = lambda path, mode=None: open_output(path, gzip_output)
    write_mode = None
    
    buffer = []  # Temporary storage for the batch

    input_fastq_read_only = pyfastx.Fastx(input_fastq)

    # use tqdm if verbose else silently loop
    iterator = (
        tqdm(input_fastq_read_only, desc=f"Iterating through seed {seed}, file {input_fastq}", unit="read", total=total_reads)
        if verbose else input_fastq_read_only
    )
    
    if collapse_duplicates or oob_path:
        # occurrence_list is indexed by read position (dense array or Counter), so pair each read with its count lazily
        counts = map(occurrence_list.__getitem__, itertools.count())
        with contextlib.ExitStack() as stack:
            f = stack.enter_context(open_func(output_path, write_mode))
            f_oob = stack.enter_context(open_func(oob_path, write_mode)) if oob_path else None
            try:
                num_records = write_reads_tagged(zip(iterator, counts), f, f_oob, unique_headers, collapse_duplicates)
            except IndexError:
                num_records = None  # the file holds more records than total_reads
        check_record_count(input_fastq, num_records, total_reads)
        return

    i = -1  # stays -1 for an empty file
    with open_func(output_path, write_mode) as f:
        try:
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
        except IndexError:
            i = None  # the file holds more records than total_reads (dense occurrence array overran)

        # Write any remaining entries in the buffer
        if buffer:
            f.writelines(buffer)
            buffer.clear()

    check_record_count(input_fastq, None if i is None else i + 1, total_reads)

def check_record_count(input_fastq, num_records, total_reads):
    # The occurrence list is built for total_reads records, which comes from the counting pass or
    # from a user-supplied read count. A wrong supplied count would silently bias the sample (reads
    # past the end are never drawn, or drawn indices point past the last read), so the writing pass
    # confirms the length it actually saw. num_records is None when the file overran total_reads.
    if num_records != total_reads:
        seen = f"more than {total_reads}" if num_records is None or num_records > total_reads else str(num_records)
        raise ValueError(f"'{input_fastq}' contains {seen} reads, but the read count used for sampling was "
                         f"{total_reads}. Check the value passed to read_counts (--read-counts), or omit it "
                         "so that the reads are counted.")

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

def write_fastq_one_pass(input_fastq, output_path, fraction, replacement, child_seed, gzip_output, seed=None, unique_headers=False, collapse_duplicates=False, oob_path=None, verbose=True):
    # Single-pass writer. A read's output multiplicity is drawn as the file streams by, so the read
    # count is never needed and peak memory is constant (only the flush buffer). All members of a
    # group are passed the same child_seed, so re-seeding here reproduces the identical multiplicity
    # sequence for every member and keeps mate pairs synchronized without a shared occurrence vector.
    open_func = lambda path, mode=None: open_output(path, gzip_output)
    write_mode = None

    rng = np.random.default_rng(child_seed)
    occurrence_stream = occurrence_chunk_stream(rng, fraction, replacement, batch_size)

    buffer = []  # Temporary storage for the batch

    input_fastq_read_only = open_input_stream(input_fastq)

    # total is unknown without a counting pass, so tqdm shows throughput rather than a percentage
    iterator = (
        tqdm(input_fastq_read_only, desc=f"Iterating through seed {seed}, file {input_fastq}", unit="read")
        if verbose else input_fastq_read_only
    )

    if collapse_duplicates or oob_path:
        with contextlib.ExitStack() as stack:
            f = stack.enter_context(open_func(output_path, write_mode))
            f_oob = stack.enter_context(open_func(oob_path, write_mode)) if oob_path else None
            write_reads_tagged(zip(iterator, occurrence_stream), f, f_oob, unique_headers, collapse_duplicates)
        return

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

# Chunk length for the dense low-memory samplers (see make_occurrence_list).
low_memory_chunk_size = 1_000_000

def make_occurrence_list(file, seed, total_reads, number_of_reads_to_sample, replacement, low_memory, rng=None, verbose=True):
    if verbose:
        logger.info(f"Calculating total reads and determining random indices for seed {seed}, file {file}")

    use_counter = number_of_reads_to_sample < (total_reads / counter_sparsity_threshold)

    if low_memory and not use_counter:
        # Dense low-memory path. Draws are processed in fixed-size chunks and accumulated straight
        # into the occurrence array, so neither the full array of m sampled indices nor the length-n
        # bincount temporary of the default path is ever materialized; memory is the occurrence array
        # plus one chunk. With replacement, each chunk of indices is drawn uniformly and added in place.
        # Without replacement, random.sample would hold a set of m indices as Python ints, which costs
        # far more than the occurrence array itself, so instead the reads are visited in chunks: the
        # number drawn from each chunk is Hypergeometric in the reads and draws still remaining, and
        # those reads are chosen uniformly within the chunk. Both are exactly uniform sampling.
        if replacement:
            # Counts are ~Poisson(fraction); uint16 leaves an enormous safety margin (a single read would
            # need to be drawn 65,535 times to overflow) and is downcast to the realized maximum below.
            occurrence_list = np.zeros(total_reads, dtype=np.uint16)
            dtype_random_indices = np.uint32 if total_reads <= np.iinfo(np.uint32).max else np.uint64
            remaining_draws = number_of_reads_to_sample
            while remaining_draws > 0:
                k = min(low_memory_chunk_size, remaining_draws)
                np.add.at(occurrence_list, rng.integers(0, total_reads, size=k, dtype=dtype_random_indices), 1)
                remaining_draws -= k
            realized_dtype = smallest_uint_dtype(int(occurrence_list.max()) if occurrence_list.size else 0)
            if realized_dtype != occurrence_list.dtype:
                occurrence_list = occurrence_list.astype(realized_dtype)
        else:
            # Without replacement each index is distinct, so the count is exactly one and uint8 is provably sufficient.
            occurrence_list = np.zeros(total_reads, dtype=np.uint8)
            remaining_draws = number_of_reads_to_sample
            for start in range(0, total_reads, low_memory_chunk_size):
                chunk_len = min(low_memory_chunk_size, total_reads - start)
                remaining_reads = total_reads - start
                if remaining_draws <= 0:
                    break
                if remaining_draws >= remaining_reads:
                    k = chunk_len
                else:
                    k = int(rng.hypergeometric(remaining_draws, remaining_reads - remaining_draws, chunk_len))
                if k:
                    occurrence_list[start + rng.choice(chunk_len, size=k, replace=False)] = 1
                    remaining_draws -= k
        return occurrence_list

    if low_memory:
        # Sparse case: the sample is small, so a lazy stdlib generator feeding a Counter is cheapest.
        if replacement:
            random_indices = (random.choice(range(total_reads)) for _ in range(number_of_reads_to_sample))
        else:
            random_indices = (index for index in random.sample(range(total_reads), k=number_of_reads_to_sample))
        occurrence_list = Counter(random_indices)
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

def insert_suffix(filename, suffix):
    # Insert a marker before the FASTQ extension, e.g. "R1.fastq.gz" -> "R1.seed3.fastq.gz".
    for ext in sorted(valid_fastq_extensions, key=len, reverse=True):
        if filename.endswith(ext):
            return f"{filename[:-len(ext)]}.{suffix}{ext}"
    return f"{filename}.{suffix}"

def insert_seed_suffix(filename, seed):
    # Per-seed marker, so that distinct seeds write to distinct output files instead of overwriting one another.
    return insert_suffix(filename, f"seed{seed}")

def resolve_output_path(file, output_directory, gzip_output, multiple_seeds, seed):
    # Build the per-file output path, disambiguating by seed when needed and forcing the
    # extension to match the requested gzip setting.
    output_basename = "stdin.fastq" if file == STDIN_SENTINEL else os.path.basename(file)
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

def two_pass_write_jobs(file_list, fraction, seed, output, gzip_output, replacement, low_memory, unique_headers, collapse_duplicates, oob, multiple_seeds, rng, verbose):
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
            oob_path = insert_suffix(output_path, "oob") if oob else None
            yield write_fastq, dict(input_fastq=member, output_path=output_path, occurrence_list=occurrence_list, total_reads=total_reads, gzip_output=gzip_output, seed=seed, unique_headers=unique_headers, collapse_duplicates=collapse_duplicates, oob_path=oob_path, verbose=verbose)

def bootstrap_single_file_one_pass(files_total = None, child_seed = None, gzip_output = None, output_directory = None, seed = None, fraction = None, replacement = None, unique_headers = False, collapse_duplicates = False, oob = False, multiple_seeds = False, verbose=True):
    # One-pass counterpart of bootstrap_single_file. No occurrence vector is materialized and the
    # file length is never counted; every member of the group is written from the same child_seed
    # so their sampled multiplicities match read-for-read, keeping mate pairs synchronized.
    if isinstance(files_total, str):
        files_total = (files_total, )

    for file in files_total:
        output_path = resolve_output_path(file, output_directory, gzip_output, multiple_seeds, seed)
        oob_path = insert_suffix(output_path, "oob") if oob else None

        write_fastq_one_pass(input_fastq = file, output_path = output_path, fraction = fraction, replacement = replacement, child_seed = child_seed, gzip_output = gzip_output, seed = seed, unique_headers = unique_headers, collapse_duplicates = collapse_duplicates, oob_path = oob_path, verbose = verbose)

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

def sample_multiple_files(file_list, fraction, seed_list, output, gzip_output, replacement, low_memory, unique_headers, one_pass, verbose, collapse_duplicates=False, oob=False):
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
                 dict(files_total=file, child_seed=child_seed, gzip_output=gzip_output, output_directory=output, seed=seed, fraction=fraction, replacement=replacement, unique_headers=unique_headers, collapse_duplicates=collapse_duplicates, oob=oob, multiple_seeds=multiple_seeds, verbose=verbose))
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
            jobs = two_pass_write_jobs(file_list, fraction, seed, output, gzip_output, replacement, low_memory, unique_headers, collapse_duplicates, oob, multiple_seeds, rng, verbose)
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

def apply_read_counts(read_counts, file_list):
    # Record user-supplied read counts in fastq_to_length_dict so that make_fastq_to_length_dict
    # skips those files. read_counts is a dict {path: count}, or an int / list of ints given in
    # input order, either one per group (every member of a group gets its group's count) or one per
    # individual file (members of a group must then agree, since they share one occurrence list).
    if isinstance(read_counts, dict):
        pairs = list(read_counts.items())
    else:
        counts = [read_counts] if isinstance(read_counts, int) else list(read_counts)
        groups = [file if isinstance(file, tuple) else (file,) for file in file_list]
        members = [member for group in groups for member in group]
        if len(counts) == len(groups):
            pairs = [(member, count) for group, count in zip(groups, counts) for member in group]
        elif len(counts) == len(members):
            pairs = list(zip(members, counts))
            count_of = dict(pairs)
            for group in groups:
                if len({count_of[member] for member in group}) > 1:
                    raise ValueError(f"The read counts given for the grouped files {group} differ "
                                     f"({[count_of[member] for member in group]}). Files in a group must "
                                     "have the same number of reads.")
        else:
            raise ValueError(f"read_counts has {len(counts)} value(s), but there are {len(members)} input "
                             f"file(s) in {len(groups)} group(s). Give one count per file or one per group.")

    for path, count in pairs:
        if isinstance(count, bool) or not isinstance(count, int) or count < 0:
            raise ValueError(f"Read count for '{path}' must be a non-negative integer, got {count!r}.")
        fastq_to_length_dict[path] = count

@validate_call
def fastQpick(
    input_files: str | list | tuple,
    fraction: float = 1.0,
    seed: int | str | list = 42,
    num_samples: int = 1,
    output_dir: str = "fastQpick_output",
    disable_gzip: bool = False,
    file_group_size: int = 1,
    without_replacement: bool = False,
    overwrite: bool = False,
    low_memory: bool = False,
    unique_headers: Union[bool, None] = None,
    one_pass: bool = False,
    collapse_duplicates: bool = False,
    oob: bool = False,
    read_counts: Union[int, list, dict, None] = None,
    verbose: bool = True,
    **kwargs
):
    """
    Fast and memory-efficient sampling of DNA-Seq or RNA-seq fastq data with or without replacement.

    Parameters
    ----------
    input_files (str, list, or tuple)       Input FASTQ files or directories containing FASTQ files.
    fraction (int or float)                 The fraction of reads to sample, as a float greater than 0. Any value equal to or greater than 1 turns on sampling with replacement automatically.
    seed (int)                              Random seed.
    num_samples (int)                       Number of independent samples (replicates) to generate.
    output_dir (str)                        Output directory.
    disable_gzip (bool)                     Write plain (uncompressed) FASTQ. Output is gzip-compressed by default.
    file_group_size (int)                   The size of grouped files. Provide each pair of files sequentially, separated by a space. E.g., I1, R1, R2 would have file_group_size=3.
    without_replacement (bool)              Sample without replacement. Automatically disabled if fraction >= 1.
    overwrite (bool)                        Overwrite existing output files.
    low_memory (bool)                       Activate low memory mode. Disabled if one_pass is True.
    unique_headers (bool)                   Add a unique identifier to the header names of the output files. Default False if without_replacement is True, True if without_replacement is False.
    one_pass (bool)                         Use the single-pass approximate sampler. Uses low memory and runs faster in some cases (when the fraction is small or the input is large compared to available memory).
    collapse_duplicates (bool)              Write each sampled read once and record its multiplicity in the header as ";size=<count>" instead of writing duplicate records. Reduces output size when sampling with replacement. Overrides unique_headers.
    oob (bool)                              Also write the out-of-bag reads (reads not selected in a sample) to a separate "<name>.oob.fastq[.gz]" file for each output file.
    read_counts (int, list, or dict)        Number of reads in each input file, to skip the counting pass of the two-pass modes. Either a dict mapping each file path to its read count, or an int / list of ints in input order (after directory expansion), with one count per file or one per group. Ignored when one_pass is True. The writing pass checks each count and raises an error if it does not match the file.
    verbose (bool)                          Whether to print progress information.

    kwargs
    ------
    fastq_to_length_dict (dict)             Dictionary of FASTQ file paths to number of reads in each file. If not provided, will be calculated.
    """
    replacement = not without_replacement
    gzip_output = not disable_gzip  # output is gzip-compressed by default

    # check if fastq_to_length_dict is in kwargs
    if "fastq_to_length_dict" in kwargs and isinstance(kwargs["fastq_to_length_dict"], dict):
        global fastq_to_length_dict
        fastq_to_length_dict = kwargs["fastq_to_length_dict"]

    # Check overwrite
    if not overwrite:
        if os.path.exists(output_dir) and not is_directory_effectively_empty(output_dir):  # check if dir exists and is not empty
            raise FileExistsError(f"Output directory '{output_dir}' already exists. Please specify a different output directory or set the overwrite flag to True.")

    # Normalize the seed specification into a flat list of integer seeds. The seed argument doubles as a
    # (hidden) way to request multiple samples for backwards compatibility: if it expands to more than one
    # seed (a list or a dash-delimited range string), each seed produces one sample and num_samples is
    # overridden to match. For a single seed, num_samples consecutive seeds are derived from it so that
    # num_samples controls the number of independent replicates. With the defaults (seed=42, num_samples=1)
    # this yields the single seed 42. The user need not be aware of this seed/num_samples interplay.
    seeds = parse_seed(seed)
    if len(seeds) > 1:
        num_samples = len(seeds)
    else:
        seeds = list(range(seeds[0], seeds[0] + num_samples))

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
        if path == STDIN_SENTINEL:
            input_files_parsed.append(path)
            continue
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

    # Standard input is not seekable and its length is not known in advance, so it is only
    # compatible with the one-pass sampler reading a single, ungrouped stream.
    if STDIN_SENTINEL in input_files_parsed:
        if len(input_files_parsed) > 1:
            raise ValueError("Reading from standard input ('-') cannot be combined with other input files.")
        if not one_pass:
            raise ValueError("Reading from standard input ('-') requires one_pass=True (--one-pass on the "
                             "command line): the default and low-memory modes must read the library twice, "
                             "which a stream does not allow.")
        if file_group_size > 1:
            raise ValueError("Reading from standard input ('-') cannot be combined with file grouping "
                             "(file_group_size > 1), which needs one stream per group member.")

    if file_group_size > 1:
        input_files_parsed = group_items(input_files_parsed, group_size=file_group_size)
    
    if collapse_duplicates:
        unique_headers = False  # each read is written once, so its header is already unique
    if unique_headers is None:
        unique_headers = replacement  # default to True if replacement is True, False if replacement is False
    if replacement and not unique_headers and not collapse_duplicates:
        logger.warning(f"unique_headers is {unique_headers} but replacement is {replacement}. This may lead to duplicate header names in the output files, which can cause issues for downstream tools. Consider setting unique_headers to {replacement} to match the replacement setting.")
    if not replacement and unique_headers:
        logger.warning(f"unique_headers is {unique_headers} but replacement is {replacement}. This may lead to unnecessarily verbose header names in the output files. Consider setting unique_headers to {replacement} to match the replacement setting.")
    
    # Count reads in each file and store in a dictionary. The one-pass sampler does not need the
    # counts, so this pass is skipped entirely in that mode. Files with a user-supplied count are
    # not counted.
    if read_counts is not None and one_pass:
        if verbose:
            logger.warning("read_counts has no effect in one_pass mode, which never counts the reads.")
    elif read_counts is not None:
        apply_read_counts(read_counts, input_files_parsed)
    if not one_pass:
        make_fastq_to_length_dict(input_files_parsed, verbose=verbose)

    # Do the sampling
    sample_multiple_files(file_list=input_files_parsed, fraction=fraction, seed_list=seeds, output=output_dir, gzip_output=gzip_output, replacement=replacement, low_memory=low_memory, unique_headers=unique_headers, one_pass=one_pass, verbose=verbose, collapse_duplicates=collapse_duplicates, oob=oob)

def parse_read_counts(value):
    # argparse type for --read-counts: a comma-separated list of integers. A single token is used
    # rather than nargs="+" so that the flag cannot swallow the positional input files after it.
    try:
        return [int(token) for token in value.split(",") if token.strip()]
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid read counts '{value}'. Expected comma-separated integers, e.g. 1000,1000.")

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description="Fast and memory-efficient sampling of DNA-Seq or RNA-seq fastq data with or without replacement.")
    parser.add_argument("-f", "--fraction", required=True, default=False, help="The fraction of reads to sample, as a float greater than 0. Any value equal to or greater than 1 turns on sampling with replacement automatically.")
    parser.add_argument("-s", "--seed", required=False, default=42, nargs="+", help='Random seed.')
    parser.add_argument("-B", "-n", "--num-samples", required=False, type=int, default=1, help="Number of independent samples (replicates) to generate")
    parser.add_argument("-o", "--output-dir", required=False, type=str, default="fastQpick_output", help="Output directory.")
    parser.add_argument("-z", "--disable-gzip", action="store_true", help="Write plain (uncompressed) FASTQ. Output is gzip-compressed by default.")
    parser.add_argument("-g", "--file-group-size", required=False, default=1, help="The size of grouped files. Provide each pair of files sequentially, separated by a space. E.g., I1, R1, R2 would have file_group_size=3.")
    parser.add_argument("-dr", "--without-replacement", action="store_true", help="Sample without replacement. Automatically disabled if fraction >= 1.")
    parser.add_argument("-w", "--overwrite", action="store_true", help="Overwrite existing output files.")
    parser.add_argument("-l", "--low-memory", action="store_true", help="Activate low memory mode. Disabled if one_pass is True.")
    parser.add_argument("-u", "--unique-headers", action="store_true", default=None, help="Add a unique identifier to the header names of the output files. Defaults to True when sampling with replacement, False otherwise.")
    parser.add_argument("-p", "--one_pass", action="store_true", help="Use the single-pass approximate sampler. Uses low memory and runs faster in some cases (when the fraction is small or the input is large compared to available memory).")
    parser.add_argument("-c", "--collapse-duplicates", action="store_true", help='Write each sampled read once and record its multiplicity in the header as ";size=<count>" instead of writing duplicate records. Overrides --unique-headers.')
    parser.add_argument("--oob", action="store_true", help='Also write the out-of-bag reads (reads not selected in a sample) to a separate "<name>.oob.fastq[.gz]" file.')
    parser.add_argument("--read-counts", required=False, type=parse_read_counts, default=None, help="Comma-separated number of reads in each input file, in input order (one per file, or one per group with -g), e.g. --read-counts 5725730,5725730. Skips the counting pass. Ignored with --one_pass. An error is raised if a count does not match its file.")
    parser.add_argument("-q", "--quiet", action="store_false", help="Whether to print progress information.")
    parser.add_argument("-v", "--version", action="version", version=f"fastQpick {__version__}", help="Show program's version number and exit")

    # Positional argument for input files (indefinite number)
    parser.add_argument("input_files", nargs="+", help="Input FASTQ files or directories containing FASTQ files.")

    # Parse arguments
    args = parser.parse_args()
            
    fastQpick(input_files=args.input_files,
              fraction=args.fraction,
              seed=args.seed,
              num_samples=args.num_samples,
              output_dir=args.output_dir,
              disable_gzip=args.disable_gzip,
              file_group_size=args.file_group_size,
              without_replacement=args.without_replacement,
              overwrite=args.overwrite,
              low_memory=args.low_memory,
              unique_headers=args.unique_headers,
              one_pass=args.one_pass,
              collapse_duplicates=args.collapse_duplicates,
              oob=args.oob,
              read_counts=args.read_counts,
              verbose=args.quiet)
