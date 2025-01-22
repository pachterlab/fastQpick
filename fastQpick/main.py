import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import os
import random
from tqdm import tqdm
import pyfastx  # to loop through fastq (faster than custom python code)
from pysam import BGZFile  # to write to BGZF gzip files

from fastQpick.utils import save_params_to_config_file, is_directory_effectively_empty, pair_items, count_reads

# Global variables
valid_fastq_extensions = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
use_buffer = True
batch_size = 100000
use_bgzf = True  # False to use regular gzip
fastq_to_length_dict = {}  # set to empty, and the user can provide otherwise it will be calculated

def write_bgzf_without_buffer(f, name, seq, qual, occurrences):
    """Writes to a BGZF file."""
    # Join and encode the multiplied string
    f.write((f"@{name}\n{seq}\n+\n{qual}\n" * occurrences).encode("utf-8"))

def write_plain_without_buffer(f, name, seq, qual, occurrences):
    """Writes to a plain text file."""
    # Write directly using writelines
    f.writelines(f"@{name}\n{seq}\n+\n{qual}\n" * occurrences)

def write_bgzf(f, buffer):
    """Writes the buffer to a BGZF file."""
    f.write("".join(buffer).encode("utf-8"))

def write_plain(f, buffer):
    """Writes the buffer to a plain text file."""
    f.writelines(buffer)

def write_fastq(input_fastq, output_path, occurrence_list, total_reads, gzip_output, seed = None, verbose = True):
    if gzip_output:
        global use_bgzf
        open_func = BGZFile if use_bgzf else gzip.open
        write_mode = "wb" if use_bgzf else "wt"
        if use_buffer:
            write_func = write_bgzf if use_bgzf else write_plain
        else:
            write_func = write_bgzf_without_buffer if use_bgzf else write_plain_without_buffer
    else:
        open_func = open
        write_mode = "w"
        write_func = write_plain if use_buffer else write_plain_without_buffer
    
    buffer = []  # Temporary storage for the batch

    input_fastq_read_only = pyfastx.Fastx(input_fastq)

    # use tqdm if verbose else silently loop
    iterator = (
        tqdm(input_fastq_read_only, desc=f"Iterating through seed {seed}, file {input_fastq}", unit="read", total=total_reads)
        if verbose else input_fastq_read_only
    )

    if use_buffer:
        with open_func(output_path, write_mode) as f:
            for i, (name, seq, qual) in enumerate(iterator):
                # Add the FASTQ entry to the buffer
                buffer.extend([f"@{name}\n{seq}\n+\n{qual}\n"] * occurrence_list[i])
                
                # If the buffer reaches the batch size, write all at once and clear the buffer
                if (i + 1) % batch_size == 0:
                    write_func(f, buffer)
                    buffer.clear()  # Clear the buffer after writing
            
            # Write any remaining entries in the buffer
            if buffer:
                f.writelines(buffer)
                buffer.clear()
    else:
        with open_func(output_path, write_mode) as f:
            for i, (name, seq, qual) in enumerate(iterator):
                write_func(f, name, seq, qual, occurrence_list[i])

def make_occurrence_list(file, seed, total_reads, number_of_reads_to_sample, replacement, verbose):
    if replacement:
        random_indices = random.choices(range(total_reads), k=number_of_reads_to_sample)  # with replacement
    else:
        random_indices = random.sample(range(total_reads), k=number_of_reads_to_sample)  # without replacement

    # Initialize a list with zeros
    occurrence_list = [0] * total_reads

    # use tqdm if verbose, else just silently loop through
    iterator = (
        tqdm(random_indices, desc=f"Counting occurrences for seed {seed}, file {file}", unit="read", total=number_of_reads_to_sample)
        if verbose else random_indices
    )

    # Count occurrences (I don't use a counter in order to save memory, as a counter is essentially a dictionary)
    for index in iterator:
        occurrence_list[index] += 1

    del random_indices

    return occurrence_list

def bootstrap_single_file(file = None, file1 = None, file2 = None, gzip_output = None, output_path = None, output_path1 = None, output_path2 = None, seed = None, fraction = None, replacement = None, verbose=True):
    # Create output directory if it doesn't exist
    output_path_args = [output_path, output_path1, output_path2]
    for output_path_arg in output_path_args:
        if output_path_arg and os.path.dirname(output_path_arg):
            os.makedirs(os.path.dirname(output_path_arg), exist_ok=True)

    if file:
        if not output_path:
            output_path = file.replace(".fastq", f"_bootstrapped_seed{seed}.fastq").replace(".fq", f"_bootstrapped_seed{seed}.fq")
        if gzip_output and not output_path.endswith(".gz"):
            output_path += ".gz"

        if verbose:
            print(f"Calculating total reads and determining random indices for seed {seed}, file {file}")
        total_reads = fastq_to_length_dict[file]
        number_of_reads_to_sample = int(fraction * total_reads)

        occurrence_list = make_occurrence_list(file=file, seed=seed, total_reads=total_reads, number_of_reads_to_sample=number_of_reads_to_sample, replacement=replacement, verbose=verbose)

        # write fastq
        write_fastq(input_fastq = file, output_path = output_path, occurrence_list = occurrence_list, total_reads = total_reads, gzip_output = gzip_output, seed = seed, verbose = verbose)

    elif file1 and file2:
        if not output_path1:
            output_path1 = file1.replace(".fastq", f"_bootstrapped_seed{seed}.fastq").replace(".fq", f"_bootstrapped_seed{seed}.fq")
        if gzip_output and not output_path1.endswith(".gz"):
            output_path1 += ".gz"
        if not output_path2:
            output_path2 = file.replace(".fastq", f"_bootstrapped_seed{seed}.fastq").replace(".fq", f"_bootstrapped_seed{seed}.fq")
        if gzip_output and not output_path2.endswith(".gz"):
            output_path2 += ".gz"
                
        # Paired-end FASTQ files
        if verbose:
            print(f"Calculating total reads and determining random indices for seed {seed}, file {file1}")
        total_reads = fastq_to_length_dict[file1]
        number_of_reads_to_sample = int(fraction * total_reads)

        occurrence_list = make_occurrence_list(file=file1, seed=seed, total_reads=total_reads, number_of_reads_to_sample=number_of_reads_to_sample, replacement=replacement, verbose=verbose)

        # TODO: multithread this
        # write fastqs
        write_fastq(input_fastq = file1, output_path = output_path1, occurrence_list = occurrence_list, total_reads = total_reads, gzip_output = gzip_output, seed = seed, verbose = verbose)
        write_fastq(input_fastq = file2, output_path = output_path2, occurrence_list = occurrence_list, total_reads = total_reads, gzip_output = gzip_output, seed = seed, verbose = verbose)

    else:
        raise ValueError("You must provide either a single FASTQ file or paired-end FASTQ files.")


def process_seed_and_file(seed, file, fraction, gzip_output, replacement, output_directory, verbose):
    random.seed(seed)
    # print(f"seed {seed}, file {file}")
    if isinstance(file, tuple):
        assert len(file) == 2, "Paired-end FASTQ files must be a 2-tuple."
        output_path1 = os.path.join(output_directory, os.path.basename(file[0]))
        output_path2 = os.path.join(output_directory, os.path.basename(file[1]))
        bootstrap_single_file(file1 = file[0], file2 = file[1], gzip_output = gzip_output, output_path = None, output_path1 = output_path1, output_path2 = output_path2, seed = seed, fraction = fraction, replacement = replacement, verbose = verbose)
    elif isinstance(file, str):
        output_path = os.path.join(output_directory, os.path.basename(file))
        bootstrap_single_file(file = file, gzip_output = gzip_output, output_path = output_path, seed = seed, fraction = fraction, replacement = replacement, verbose = verbose)

def sample_multiple_files(file_list, fraction, seed_list, output, threads, gzip_output, replacement, verbose):
    # Use ThreadPoolExecutor to process seeds in parallel
    with ThreadPoolExecutor(max_workers=threads) as executor:
        # Submit tasks for all combinations of seeds and files
        futures = [
            executor.submit(
                process_seed_and_file,
                seed,
                file,
                fraction,
                gzip_output,
                replacement,
                output,
                verbose
            )
            for seed in seed_list
            for file in file_list
        ]

# TODO: multithread across file list    
def make_fastq_to_length_dict(file_list, verbose=True):
    global fastq_to_length_dict
    for file in file_list:
        if isinstance(file, tuple):
            assert len(file) == 2, "Paired-end FASTQ files must be a 2-tuple."
            if file[0] in fastq_to_length_dict and file[1] in fastq_to_length_dict:
                continue
            if verbose:
                print(f"Counting {file[0]}")
            count = count_reads(file[0])
            fastq_to_length_dict[file[0]] = count
            fastq_to_length_dict[file[1]] = count
        elif isinstance(file, str):
            if file in fastq_to_length_dict:
                continue
            if verbose:
                print(f"Counting {file}")
            count = count_reads(file)
            fastq_to_length_dict[file] = count
    if verbose:
        print("fastq_to_length_dict:", fastq_to_length_dict)

def fastQpick(input_file_list, fraction, seed=42, output_dir="fastQpick_output", threads=1, gzip_output=False, paired=False, replacement=False, overwrite=False, verbose=True, **kwargs):
    """
    Fast and memory-efficient sampling of DNA-Seq or RNA-seq fastq data with or without replacement.

    Parameters
    ----------
    input_file_list (list)      List of input FASTQ files or directories containing FASTQ files.
    fraction (int or float)     The fraction of reads to sample, as a float greater than 0. Any value equal to or greater than 1 will turn on the -r flag automatically.
    seed (int or str)           Random seed(s). Can provide multiple seeds separated by commas. Default: 42
    output_dir (str)            Output directory. Default: ./fastQpick_output
    threads (int)               Number of threads. Default: 2
    gzip_output (bool)          Whether or not to gzip the output. Default: False (uncompressed)
    paired (bool)               Whether or not the fastq files are paired-end. If paired-end, provide each pair of files sequentially, separated by a space. Default: False
    replacement (bool)          Sample with replacement. Default: False (without replacement).
    overwrite (bool)            Overwrite existing output files. Default: False
    verbose (bool)              Whether to print progress information. Default: True

    kwargs
    ------
    fastq_to_length_dict (dict) Dictionary of FASTQ file paths to number of reads in each file. If not provided, will be calculated.
    use_bgzf (bool)             Use BGZF compression for gzip output. Default: True
    """
    # check if fastq_to_length_dict is in kwargs
    if "fastq_to_length_dict" in kwargs and isinstance(kwargs["fastq_to_length_dict"], dict):
        global fastq_to_length_dict
        fastq_to_length_dict = kwargs["fastq_to_length_dict"]

    # Check overwrite
    if not overwrite:
        if os.path.exists(output_dir) and not is_directory_effectively_empty(output_dir):  # check if dir exists and is not empty
            raise FileExistsError(f"Output directory '{output_dir}' already exists. Please specify a different output directory or set the overwrite flag to True.")

    # Save arguments to a config file
    os.makedirs(output_dir, exist_ok=True)
    config_file = os.path.join(output_dir, "fastQpick_config.json")
    save_params_to_config_file(config_file)

    # type checking
    # if fraction >= 1, set replacement to True
    if float(fraction) >= 1.0:
        replacement = True

    # go through files, and only keep those that are valid fastq files or that are a folder containing valid fastq files in the direct subdirectory
    input_file_list_parsed = []
    if isinstance(input_file_list, str):
        input_file_list_parsed = [input_file_list]
    elif isinstance(input_file_list, tuple) or isinstance(input_file_list, list):
        for path in input_file_list:
            if not isinstance(path, str):
                raise ValueError("Input file list must be a string, tuple of strings, or list of strings.")
            if not os.path.exists(path):
                raise FileNotFoundError(f"File or directory '{path}' not found.")
            elif os.path.isfile(path) and not path.endswith(tuple(valid_fastq_extensions)):
                raise ValueError(f"File '{path}' is not a valid FASTQ file.")
            elif os.path.isdir(path):
                input_files_before_path = input_file_list_parsed.copy()
                for subpath in os.listdir(path):
                    if os.path.isfile(subpath) and subpath.endswith(tuple(valid_fastq_extensions)):
                        input_file_list_parsed.append(subpath)
                if input_files_before_path == input_file_list_parsed:
                    raise ValueError(f"No valid FASTQ files found in directory '{path}'.")
            elif os.path.isfile(path) and path.endswith(tuple(valid_fastq_extensions)):
                input_file_list_parsed.append(path)
    else:
        raise ValueError("Input file list must be a string, tuple of strings, or list of strings.")

    if isinstance(seed, int):  # if a single int is passed as a seed
        seed = [seed]
    elif isinstance(seed, str):  # if a string of comma-separated ints is passed as a seed (like on the command line)
        seed = [int(specific_seed) for specific_seed in seed.split(",")]

    if paired:
        input_file_list_parsed = pair_items(input_file_list_parsed)
    
    # Count reads in each file and store in a dictionary
    make_fastq_to_length_dict(input_file_list_parsed, verbose=verbose)

    # Do the sampling
    sample_multiple_files(file_list=input_file_list_parsed, fraction=fraction, seed_list=seed, output=output_dir, threads=threads, gzip_output=gzip_output, replacement=replacement, verbose=verbose)

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description="Fast and memory-efficient sampling of DNA-Seq or RNA-seq fastq data with or without replacement.")
    parser.add_argument("-f", "--fraction", required=True, default=False, help="The fraction of reads to sample, as a float greater than 0. Any value equal to or greater than 1 will turn on the -r flag automatically.")
    parser.add_argument("-s", "--seed", required=False, default=42, help="Random seed(s). Can provide multiple seeds separated by commas. Default: 42")
    parser.add_argument("-o", "--output_dir", required=False, type=str, default="fastQpick_output", help="Output file path. Default: ./fastQpick_output")
    parser.add_argument("-t", "--threads", required=False, type=int, default=1, help="Number of threads. Default: 2")
    parser.add_argument("-z", "--gzip_output", required=False, default=False, help="Whether or not to gzip the output. Default: False (uncompressed)")
    parser.add_argument("-p", "--paired", required=False, default=False, help="Whether or not the fastq files are paired-end. If paired-end, provide each pair of files sequentially, separated by a space. Default: False")
    parser.add_argument("-r", "--replacement", required=False, default=False, help="Sample with replacement. Default: False (without replacement).")
    parser.add_argument("-w", "--overwrite", required=False, default=False, help="Overwrite existing output files. Default: False")
    parser.add_argument("-q", "--quiet", required=False, default=False, help="Turn off verbose output. Default: False")

    # Positional argument for input files (indefinite number)
    parser.add_argument("input_file_list", nargs="+", help="Input FASTQ file(s) (one after the other, space-separated) or FASTQ folder(s)")

    # Parse arguments
    args = parser.parse_args()
    verbose = not args.quiet
            
    fastQpick(input_file_list=args,
              fraction=args.fraction,
              seed=args.seed,
              output=args.output_dir,
              threads=args.threads,
              gzip_output=args.gzip_output,
              paired=args.paired,
              replacement=args.replacement,
              verbose=verbose)
