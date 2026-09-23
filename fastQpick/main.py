import argparse
import contextlib
import io
import itertools
import os
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
from fastQpick.utils import save_params_to_config_file, is_directory_effectively_empty, group_items, count_reads, parse_seed, available_cpus

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
gzip_output_threads = 4  # most deflate threads per output file when igzip_threaded is available
default_threads = 4  # total thread budget when the caller does not set one
STDIN_SENTINEL = "-"  # read the library from standard input instead of from a file

def open_output(path, gzip_output, gzip_threads=gzip_output_threads):
    # Text-mode output handle using the fastest gzip writer available. gzip_threads is the number of
    # background deflate threads; 0 compresses in the calling thread.
    if not gzip_output:
        return open(path, "w")
    if igzip_threaded is not None:
        return igzip_threaded.open(path, "wt", compresslevel=gzip_compresslevel,
                                   threads=gzip_threads)
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

def write_fastq(input_fastq, output_path, occurrence_list, total_reads, gzip_output, seed = None, unique_headers = False, collapse_duplicates = False, oob_path = None, verbose = True, gzip_threads = gzip_output_threads):
    open_func = lambda path, mode=None: open_output(path, gzip_output, gzip_threads)
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
    # Single-pass sampling: lazily yield the number of times each successive read should appear in
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

def write_fastq_single_pass(input_fastq, output_path, fraction, replacement, child_seed, gzip_output, seed=None, unique_headers=False, collapse_duplicates=False, oob_path=None, verbose=True, gzip_threads=gzip_output_threads):
    # Single-pass writer. A read's output multiplicity is drawn as the file streams by, so the read
    # count is never needed and peak memory is constant (only the flush buffer). All members of a
    # group are passed the same child_seed, so re-seeding here reproduces the identical multiplicity
    # sequence for every member and keeps mate pairs synchronized without a shared occurrence vector.
    open_func = lambda path, mode=None: open_output(path, gzip_output, gzip_threads)
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
# is a small fraction of the file; below this ratio the Counter wins, above it the dense array wins.
counter_sparsity_threshold = 100

# Number of reads per block when the dense occurrence vector is filled block by block (see
# make_occurrence_list). Working memory beyond the occurrence vector itself is O(block size).
occurrence_block_size = 1_000_000

def make_occurrence_list(file, seed, total_reads, number_of_reads_to_sample, replacement, rng, verbose=True):
    # Return how many times each read (by position) is drawn in an exact uniform sample of
    # number_of_reads_to_sample reads. rng is a seeded numpy Generator threaded down from
    # sample_multiple_files so that the output is reproducible.
    if verbose:
        logger.info(f"Calculating total reads and determining random indices for seed {seed}, file {file}")

    n, m = total_reads, number_of_reads_to_sample

    if m < n / counter_sparsity_threshold:
        # Sparse case: the m sampled indices are few, so they are drawn at once and counted in a Counter.
        if replacement:
            dtype_random_indices = np.uint32 if n <= np.iinfo(np.uint32).max else np.uint64
            random_indices = rng.integers(0, n, size=m, dtype=dtype_random_indices)
        else:
            random_indices = rng.choice(n, size=m, replace=False)
        return Counter(random_indices.tolist())

    # Dense case: the reads are split into blocks and the occurrence vector is filled one block at a
    # time, so neither the m sampled indices nor a length-n counting temporary is ever materialized.
    # The number of draws landing in each block is drawn first (Multinomial with replacement,
    # Hypergeometric without), and the draws are then placed uniformly within the block, which is
    # exactly uniform sampling of m reads from the whole file.
    occurrence_list = np.zeros(n, dtype=np.uint8)
    remaining_draws = m
    if replacement:
        block_starts = np.arange(0, n, occurrence_block_size)
        block_lengths = np.minimum(occurrence_block_size, n - block_starts)
        draws_per_block = rng.multinomial(m, block_lengths / n)
        for start, block_length, k in zip(block_starts, block_lengths, draws_per_block):
            if not k:
                continue
            counts = np.bincount(rng.integers(0, block_length, size=k, dtype=np.uint32), minlength=block_length)
            # Counts are ~Poisson(fraction), so uint8 almost always suffices; widen only if a count
            # actually exceeds it, keeping the dtype at the smallest one that holds the realized maximum.
            max_count = int(counts.max())
            if max_count > np.iinfo(occurrence_list.dtype).max:
                occurrence_list = occurrence_list.astype(smallest_uint_dtype(max_count))
            occurrence_list[start:start + block_length] = counts
    else:
        # Without replacement each read is drawn at most once, so uint8 is provably sufficient.
        for start in range(0, n, occurrence_block_size):
            if remaining_draws <= 0:
                break
            block_length = min(occurrence_block_size, n - start)
            remaining_reads = n - start
            if remaining_draws >= remaining_reads:
                k = block_length
            else:
                k = int(rng.hypergeometric(remaining_draws, remaining_reads - remaining_draws, block_length))
            if k:
                occurrence_list[start + rng.choice(block_length, size=k, replace=False)] = 1
                remaining_draws -= k
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

def two_pass_write_jobs(file_list, fraction, seed, output, gzip_output, replacement, unique_headers, collapse_duplicates, oob, multiple_seeds, rng, verbose, gzip_threads=gzip_output_threads):
    # Lazily yield (callable, kwargs) write jobs for one seed. Each group's occurrence list is built
    # HERE, in the parent process, so the shared `rng` is consumed in exactly the same order as the
    # single-process path; only the I/O-bound write_fastq call is handed off as a job. This keeps the
    # default-path output byte-identical to the serial version while letting the writes run
    # concurrently across files. Members of a group share one occurrence list, computed once.
    for file in file_list:
        files_total = (file, ) if isinstance(file, str) else file

        total_reads = fastq_to_length_dict[files_total[0]]
        number_of_reads_to_sample = int(fraction * total_reads)

        occurrence_list = make_occurrence_list(file=files_total[0], seed=seed, total_reads=total_reads, number_of_reads_to_sample=number_of_reads_to_sample, replacement=replacement, rng=rng, verbose=verbose)

        for member in files_total:
            output_path = resolve_output_path(member, output, gzip_output, multiple_seeds, seed)
            oob_path = insert_suffix(output_path, "oob") if oob else None
            yield write_fastq, dict(input_fastq=member, output_path=output_path, occurrence_list=occurrence_list, total_reads=total_reads, gzip_output=gzip_output, seed=seed, unique_headers=unique_headers, collapse_duplicates=collapse_duplicates, oob_path=oob_path, verbose=verbose, gzip_threads=gzip_threads)

def bootstrap_single_file_single_pass(files_total = None, child_seed = None, gzip_output = None, output_directory = None, seed = None, fraction = None, replacement = None, unique_headers = False, collapse_duplicates = False, oob = False, multiple_seeds = False, verbose=True, gzip_threads=gzip_output_threads):
    # Single-pass counterpart of bootstrap_single_file. No occurrence vector is materialized and the
    # file length is never counted; every member of the group is written from the same child_seed
    # so their sampled multiplicities match read-for-read, keeping mate pairs synchronized.
    if isinstance(files_total, str):
        files_total = (files_total, )

    for file in files_total:
        output_path = resolve_output_path(file, output_directory, gzip_output, multiple_seeds, seed)
        oob_path = insert_suffix(output_path, "oob") if oob else None

        write_fastq_single_pass(input_fastq = file, output_path = output_path, fraction = fraction, replacement = replacement, child_seed = child_seed, gzip_output = gzip_output, seed = seed, unique_headers = unique_headers, collapse_duplicates = collapse_duplicates, oob_path = oob_path, verbose = verbose, gzip_threads = gzip_threads)

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

def resolve_threads(threads):
    # An explicit thread count is honored as given. The default is a small fixed budget, as in most
    # bioinformatics tools, so that an unconfigured run neither monopolizes a shared node nor varies
    # with the machine; it is lowered on machines with fewer cores than that.
    if threads is None:
        return min(default_threads, available_cpus())
    return threads

def split_thread_budget(threads, num_jobs):
    # Divide a total thread budget between parallel write jobs and the deflate threads inside each
    # job. Each worker spends one thread on its record loop (which holds the GIL), and whatever is
    # left of its share goes to background compression, up to gzip_output_threads. With a single job
    # on a machine with 5+ cores this reproduces the previous behavior (one process, four deflate
    # threads); with many jobs every thread becomes a worker and compression runs inline.
    workers = max(1, min(num_jobs, threads))
    gzip_threads = min(gzip_output_threads, threads // workers - 1)
    return workers, gzip_threads

def single_pass_jobs_all_seeds(file_list, seed_list, output, gzip_output, fraction, replacement, unique_headers, collapse_duplicates, oob, multiple_seeds, verbose, gzip_threads):
    for seed in seed_list:
        # Derive one independent sub-seed per group from the master seed. Members of a group share
        # their sub-seed (handled inside bootstrap_single_file_single_pass) so mates stay synchronized,
        # while different groups draw independently, mirroring the two-pass path where each file
        # consumes fresh random state. Each (seed, group) pair is fully independent, so it is one
        # parallel write job.
        child_seeds = np.random.SeedSequence(seed).spawn(len(file_list))
        for file, child_seed in zip(file_list, child_seeds):
            yield bootstrap_single_file_single_pass, dict(files_total=file, child_seed=child_seed, gzip_output=gzip_output, output_directory=output, seed=seed, fraction=fraction, replacement=replacement, unique_headers=unique_headers, collapse_duplicates=collapse_duplicates, oob=oob, multiple_seeds=multiple_seeds, verbose=verbose, gzip_threads=gzip_threads)

def two_pass_jobs_all_seeds(file_list, seed_list, fraction, output, gzip_output, replacement, unique_headers, collapse_duplicates, oob, multiple_seeds, verbose, gzip_threads):
    for seed in seed_list:
        # Seed a numpy Generator once per seed. A single Generator is shared across the seed's files so
        # that successive files draw independent samples while remaining reproducible. This generator
        # is consumed lazily in the parent, so a seed's RNG is (re)seeded only after every
        # occurrence list of the previous seed has been built; the RNG state each list sees is
        # therefore the same as in a serial run, however the writes are scheduled across processes.
        rng = np.random.default_rng(seed)
        yield from two_pass_write_jobs(file_list, fraction, seed, output, gzip_output, replacement, unique_headers, collapse_duplicates, oob, multiple_seeds, rng, verbose, gzip_threads=gzip_threads)

def sample_multiple_files(file_list, fraction, seed_list, output, gzip_output, replacement, unique_headers, single_pass, verbose, collapse_duplicates=False, oob=False, threads=None):
    multiple_seeds = len(seed_list) > 1
    threads = resolve_threads(threads)
    # Jobs from every seed share one pool, so replicates of a single file run in parallel as well as
    # distinct files. A job is a whole group in single-pass mode (its members re-derive the same
    # multiplicities from a shared sub-seed) and an individual output file in two-pass mode (group
    # members share an occurrence list built in the parent).
    if single_pass:
        num_jobs = len(file_list) * len(seed_list)
    else:
        num_individual_files = sum(len(file) if isinstance(file, tuple) else 1 for file in file_list)
        num_jobs = num_individual_files * len(seed_list)
    if STDIN_SENTINEL in file_list:
        num_jobs = 1  # a single stream cannot be read by several workers at once
    max_workers, gzip_threads = split_thread_budget(threads, num_jobs)

    common = dict(file_list=file_list, seed_list=seed_list, output=output, gzip_output=gzip_output, fraction=fraction, replacement=replacement, unique_headers=unique_headers, collapse_duplicates=collapse_duplicates, oob=oob, multiple_seeds=multiple_seeds, verbose=verbose, gzip_threads=gzip_threads)
    if single_pass:
        jobs = single_pass_jobs_all_seeds(**common)
    else:
        jobs = two_pass_jobs_all_seeds(**common)
    run_write_jobs(jobs, max_workers=max_workers)

def make_fastq_to_length_dict(file_list, verbose=True, threads=None):
    # Count the reads of every file (or of the first member of every group; members share its count)
    # that does not already have a count. Counting is I/O- and decompression-bound and each file is
    # independent, so with more than one file to count the files are counted in parallel.
    global fastq_to_length_dict
    to_count = []  # (file counted, keys that receive its count)
    for file in file_list:
        members = file if isinstance(file, tuple) else (file,)
        if all(member in fastq_to_length_dict for member in members):
            continue
        to_count.append((members[0], members))

    if to_count:
        if verbose:
            for counted, _ in to_count:
                logger.info(f"Counting {counted}")
        paths = [counted for counted, _ in to_count]
        max_workers = min(len(paths), resolve_threads(threads))
        if max_workers <= 1:
            counts = [count_reads(path) for path in paths]
        else:
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                counts = list(executor.map(count_reads, paths))
        for (_, members), count in zip(to_count, counts):
            for member in members:
                fastq_to_length_dict[member] = count
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
    unique_headers: Union[bool, None] = None,
    single_pass: bool = False,
    collapse_duplicates: bool = False,
    oob: bool = False,
    read_counts: Union[int, list, dict, None] = None,
    threads: Union[int, None] = None,
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
    unique_headers (bool)                   Add a unique identifier to the header names of the output files. Default False if without_replacement is True, True if without_replacement is False.
    single_pass (bool)                      Read the input once instead of twice, using constant memory. The output size is exact only in expectation (relative standard deviation 1/sqrt(fraction * n)). Required to read from standard input.
    collapse_duplicates (bool)              Write each sampled read once and record its multiplicity in the header as ";size=<count>" instead of writing duplicate records. Reduces output size when sampling with replacement. Overrides unique_headers.
    oob (bool)                              Also write the out-of-bag reads (reads not selected in a sample) to a separate "<name>.oob.fastq[.gz]" file for each output file.
    read_counts (int, list, or dict)        Number of reads in each input file, to skip the counting pass of the two-pass modes. Either a dict mapping each file path to its read count, or an int / list of ints in input order (after directory expansion), with one count per file or one per group. Ignored when single_pass is True. The writing pass checks each count and raises an error if it does not match the file.
    threads (int)                           Total number of threads (CPU cores) to use, shared between parallel file/replicate jobs and gzip compression. Defaults to 4, or to the number of available cores if fewer. Output does not depend on this value.
    verbose (bool)                          Whether to print progress information.

    kwargs
    ------
    fastq_to_length_dict (dict)             Dictionary of FASTQ file paths to number of reads in each file. If not provided, will be calculated.
    """
    if "one_pass" in kwargs:  # deprecated name of single_pass
        logger.warning("one_pass is deprecated; use single_pass instead.")
        single_pass = single_pass or bool(kwargs.pop("one_pass"))
    replacement = not without_replacement
    if threads is not None and threads < 1:
        raise ValueError(f"threads must be a positive integer, got {threads}.")
    threads = resolve_threads(threads)
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
    # compatible with the single-pass sampler reading a single, ungrouped stream.
    if STDIN_SENTINEL in input_files_parsed:
        if len(input_files_parsed) > 1:
            raise ValueError("Reading from standard input ('-') cannot be combined with other input files.")
        if not single_pass:
            raise ValueError("Reading from standard input ('-') requires single_pass=True (--single-pass on the "
                             "command line): the default (exact) mode must read the library twice, "
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
    
    # Count reads in each file and store in a dictionary. The single-pass sampler does not need the
    # counts, so this pass is skipped entirely in that mode. Files with a user-supplied count are
    # not counted.
    if read_counts is not None and single_pass:
        if verbose:
            logger.warning("read_counts has no effect in single_pass mode, which never counts the reads.")
    elif read_counts is not None:
        apply_read_counts(read_counts, input_files_parsed)
    if not single_pass:
        make_fastq_to_length_dict(input_files_parsed, verbose=verbose, threads=threads)

    # Do the sampling
    sample_multiple_files(file_list=input_files_parsed, fraction=fraction, seed_list=seeds, output=output_dir, gzip_output=gzip_output, replacement=replacement, unique_headers=unique_headers, single_pass=single_pass, verbose=verbose, collapse_duplicates=collapse_duplicates, oob=oob, threads=threads)

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
    parser.add_argument("-u", "--unique-headers", action="store_true", default=None, help="Add a unique identifier to the header names of the output files. Defaults to True when sampling with replacement, False otherwise.")
    parser.add_argument("-p", "--single-pass", action="store_true", help="Read the input once instead of twice, using constant memory. The output size is exact only in expectation (relative standard deviation 1/sqrt(fraction * n)). Required to read from standard input.")
    parser.add_argument("--one-pass", "--one_pass", dest="single_pass", action="store_true", help=argparse.SUPPRESS)  # deprecated alias
    parser.add_argument("-c", "--collapse-duplicates", action="store_true", help='Write each sampled read once and record its multiplicity in the header as ";size=<count>" instead of writing duplicate records. Overrides --unique-headers.')
    parser.add_argument("--oob", action="store_true", help='Also write the out-of-bag reads (reads not selected in a sample) to a separate "<name>.oob.fastq[.gz]" file.')
    parser.add_argument("--read-counts", required=False, type=parse_read_counts, default=None, help="Comma-separated number of reads in each input file, in input order (one per file, or one per group with -g), e.g. --read-counts 5725730,5725730. Skips the counting pass. Ignored with --single-pass. An error is raised if a count does not match its file.")
    parser.add_argument("-t", "--threads", required=False, type=int, default=None, help="Total number of threads to use, shared between parallel file/replicate jobs and gzip compression. Default: 4, or the number of available cores if fewer. Output does not depend on this value.")
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
              unique_headers=args.unique_headers,
              single_pass=args.single_pass,
              collapse_duplicates=args.collapse_duplicates,
              oob=args.oob,
              read_counts=args.read_counts,
              threads=args.threads,
              verbose=args.quiet)
