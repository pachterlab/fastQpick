import argparse
import gzip
import os
import random
from tqdm import tqdm
import pyfastx  # to loop through fastq (faster than custom python code)

from fastQpick._version import __version__
from fastQpick.utils import save_params_to_config_file, is_directory_effectively_empty, group_items, count_reads

# Global variables
valid_fastq_extensions = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
batch_size = 200000  # for buffer
fastq_to_length_dict = {}  # set to empty, and the user can provide otherwise it will be calculated

def write_fastq(input_fastq, output_path, occurrence_list, total_reads, gzip_output, seed = None, verbose = True):
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
        for i, (name, seq, qual) in enumerate(iterator):
            # Add the FASTQ entry to the buffer
            buffer.extend([f"@{name}\n{seq}\n+\n{qual}\n"] * occurrence_list[i])
            
            # If the buffer reaches the batch size, write all at once and clear the buffer
            if (i + 1) % batch_size == 0:
                f.writelines(buffer)
                buffer.clear()  # Clear the buffer after writing
        
        # Write any remaining entries in the buffer
        if buffer:
            f.writelines(buffer)
            buffer.clear()

def make_occurrence_list(file, seed, total_reads, number_of_reads_to_sample, replacement, low_memory, verbose):
    if verbose:
        print(f"Calculating total reads and determining random indices for seed {seed}, file {file}")
    if replacement:
        if low_memory:
            random_indices = (random.choice(range(total_reads)) for _ in range(number_of_reads_to_sample))
        else:
            random_indices = tuple(random.choices(range(total_reads), k=number_of_reads_to_sample))  # with replacement
    else:
        if low_memory:
            random_indices = (index for index in random.sample(range(total_reads), k=number_of_reads_to_sample))
        else:
            random_indices = tuple(random.sample(range(total_reads), k=number_of_reads_to_sample))  # without replacement

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

def bootstrap_single_file(files_total = None, gzip_output = None, output_directory = None, seed = None, fraction = None, replacement = None, low_memory = False, verbose=True):
    if isinstance(files_total, str):
        files_total = (files_total, )

    total_reads = fastq_to_length_dict[files_total[0]]
    number_of_reads_to_sample = int(fraction * total_reads)

    occurrence_list = make_occurrence_list(file=files_total[0], seed=seed, total_reads=total_reads, number_of_reads_to_sample=number_of_reads_to_sample, replacement=replacement, low_memory=low_memory, verbose=verbose)
    
    for file in files_total:
        # Create output directory if it doesn't exist
        output_path = os.path.join(output_directory, os.path.basename(file))
        if output_directory:
            os.makedirs(output_directory, exist_ok=True)

        if gzip_output and not output_path.endswith(".gz"):
            output_path += ".gz"
        elif not gzip_output and output_path.endswith(".gz"):
            output_path = output_path[:-3]

        # write fastq
        write_fastq(input_fastq = file, output_path = output_path, occurrence_list = occurrence_list, total_reads = total_reads, gzip_output = gzip_output, seed = seed, verbose = verbose)

def sample_multiple_files(file_list, fraction, seed_list, output, gzip_output, replacement, low_memory, verbose):
    for seed in seed_list:
        random.seed(seed)
        for file in file_list:
            bootstrap_single_file(files_total = file, gzip_output = gzip_output, output_directory = output, seed = seed, fraction = fraction, replacement = replacement, low_memory = low_memory, verbose = verbose)
    
def make_fastq_to_length_dict(file_list, verbose=True):
    global fastq_to_length_dict
    for file in file_list:
        if isinstance(file, tuple):
            if all(specific_file in fastq_to_length_dict for specific_file in file):
                continue
            if verbose:
                print(f"Counting {file[0]}")
            count = count_reads(file[0])
            for i in range(len(file)):
                fastq_to_length_dict[file[i]] = count
        elif isinstance(file, str):
            if file in fastq_to_length_dict:
                continue
            if verbose:
                print(f"Counting {file}")
            count = count_reads(file)
            fastq_to_length_dict[file] = count
    if verbose:
        print("fastq_to_length_dict:", fastq_to_length_dict)

def fastQpick(input_files, fraction, seed=42, output_dir="fastQpick_output", gzip_output=False, group_size=1, replacement=False, overwrite=False, low_memory=False, verbose=True, **kwargs):
    """
    Fast and memory-efficient sampling of DNA-Seq or RNA-seq fastq data with or without replacement.

    Parameters
    ----------
    input_files (str, list, or tuple)       List of input FASTQ files or directories containing FASTQ files.
    fraction (int or float)                 The fraction of reads to sample, as a float greater than 0. Any value equal to or greater than 1 will turn on the -r flag automatically.
    seed (int or str)                       Random seed(s). Can provide multiple seeds separated by commas. Default: 42
    output_dir (str)                        Output directory. Default: ./fastQpick_output
    gzip_output (bool)                      Whether or not to gzip the output. Default: False (uncompressed)
    group_size (int)                        The size of grouped files. Provide each pair of files sequentially, separated by a space. E.g., I1, R1, R2 would have group_size=3. Default: 1 (unpaired)
    replacement (bool)                      Sample with replacement. Default: False (without replacement).
    overwrite (bool)                        Overwrite existing output files. Default: False
    low_memory (bool)                       Whether to use low memory mode (uses ~5.5x less memory than default, but adds marginal time to the data structure generation preprocessing). Default: False
    verbose (bool)                          Whether to print progress information. Default: True

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

    # Save arguments to a config file
    os.makedirs(output_dir, exist_ok=True)
    config_file = os.path.join(output_dir, "fastQpick_config.json")
    save_params_to_config_file(config_file)

    # type checking
    # if fraction >= 1, set replacement to True
    if float(fraction) >= 1.0:
        replacement = True

    # go through files, and only keep those that are valid fastq files or that are a folder containing valid fastq files in the direct subdirectory
    input_files_parsed = []
    if isinstance(input_files, str):
        input_files_parsed = [input_files]
    elif isinstance(input_files, tuple) or isinstance(input_files, list):
        for path in input_files:
            if not isinstance(path, str):
                raise ValueError("Input file list must be a string, tuple of strings, or list of strings.")
            if not os.path.exists(path):
                raise FileNotFoundError(f"File or directory '{path}' not found.")
            elif os.path.isfile(path) and not path.endswith(tuple(valid_fastq_extensions)):
                raise ValueError(f"File '{path}' is not a valid FASTQ file.")
            elif os.path.isdir(path):
                input_files_before_path = input_files_parsed.copy()
                for subpath in os.listdir(path):
                    if os.path.isfile(subpath) and subpath.endswith(tuple(valid_fastq_extensions)):
                        input_files_parsed.append(subpath)
                if input_files_before_path == input_files_parsed:
                    raise ValueError(f"No valid FASTQ files found in directory '{path}'.")
            elif os.path.isfile(path) and path.endswith(tuple(valid_fastq_extensions)):
                input_files_parsed.append(path)
    else:
        raise ValueError("Input file list must be a string, tuple of strings, or list of strings.")

    if isinstance(seed, int):  # if a single int is passed as a seed
        seed = [seed]
    elif isinstance(seed, str):  # if a string of comma-separated ints is passed as a seed (like on the command line)
        seed = [int(specific_seed) for specific_seed in seed.split(",")]

    group_size = int(group_size)  # make sure group_size is an int (not a string)
    fraction = float(fraction)  # make sure fraction is a float (not a string)

    if group_size > 1:
        input_files_parsed = group_items(input_files_parsed, group_size=group_size)
    
    # Count reads in each file and store in a dictionary
    make_fastq_to_length_dict(input_files_parsed, verbose=verbose)

    # Do the sampling
    sample_multiple_files(file_list=input_files_parsed, fraction=fraction, seed_list=seed, output=output_dir, gzip_output=gzip_output, replacement=replacement, low_memory=low_memory, verbose=verbose)

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description="Fast and memory-efficient sampling of DNA-Seq or RNA-seq fastq data with or without replacement.")
    parser.add_argument("-f", "--fraction", required=True, default=False, help="The fraction of reads to sample, as a float greater than 0. Any value equal to or greater than 1 will turn on the -r flag automatically.")
    parser.add_argument("-s", "--seed", required=False, default=42, help="Random seed(s). Can provide multiple seeds separated by commas. Default: 42")
    parser.add_argument("-o", "--output_dir", required=False, type=str, default="fastQpick_output", help="Output file path. Default: ./fastQpick_output")
    parser.add_argument("-z", "--gzip_output", required=False, default=False, help="Whether or not to gzip the output. Default: False (uncompressed)")
    parser.add_argument("-g", "--group_size", required=False, default=1, help="The size of grouped files. Provide each pair of files sequentially, separated by a space. E.g., I1, R1, R2 would have group_size=3. Default: 1 (unpaired)")
    parser.add_argument("-r", "--replacement", action="store_true", help="Sample with replacement. Default: False (without replacement).")
    parser.add_argument("-w", "--overwrite", action="store_true", help="Overwrite existing output files. Default: False")
    parser.add_argument("-l", "--low_memory", action="store_true", help="Whether to use low memory mode (uses ~5.5x less memory than default, but adds marginal time to the data structure generation preprocessing). Default: False")
    parser.add_argument("-q", "--quiet", action="store_false", help="Turn off verbose output. Default: False")
    parser.add_argument("-v", "--version", action="version", version=f"fastQpick {__version__}", help="Show program's version number and exit")

    # Positional argument for input files (indefinite number)
    parser.add_argument("input_files", nargs="+", help="Input FASTQ file(s) (one after the other, space-separated) or FASTQ folder(s)")

    # Parse arguments
    args = parser.parse_args()
            
    fastQpick(input_files=args.input_files,
              fraction=args.fraction,
              seed=args.seed,
              output_dir=args.output_dir,
              gzip_output=args.gzip_output,
              group_size=args.group_size,
              replacement=args.replacement,
              overwrite=args.overwrite,
              low_memory=args.low_memory,
              verbose=args.quiet)
