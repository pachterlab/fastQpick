import gzip
import inspect
import json
import os
from collections import OrderedDict
import pyfastx

def count_reads(filepath):
    fastq_file = pyfastx.Fastx(filepath)
    num_reads = sum(1 for _ in fastq_file)
    return num_reads

def read_fastq(fastq_file, include_plus_line=False):
    is_gzipped = fastq_file.endswith(".gz")
    open_func = gzip.open if is_gzipped else open
    open_mode = "rt" if is_gzipped else "r"

    try:
        if include_plus_line:
            with open_func(fastq_file, open_mode) as file:
                while True:
                    header = file.readline().strip()
                    sequence = file.readline().strip()
                    plus_line = file.readline().strip()
                    quality = file.readline().strip()

                    if not header:
                        break

                    yield header, sequence, plus_line, quality
        else:  # copy-paste the code so that it doesn't have to check the conditional every iteration
            with open_func(fastq_file, open_mode) as file:
                while True:
                    header = file.readline().strip()
                    sequence = file.readline().strip()
                    plus_line = file.readline().strip()
                    quality = file.readline().strip()

                    if not header:
                        break

                    yield header, sequence, quality
    except Exception as e:
        raise RuntimeError(f"Error reading FASTQ file '{fastq_file}': {e}")

def make_function_parameter_to_value_dict(levels_up = 1):
    # Collect parameters in a dictionary
    params = OrderedDict()

    # Get the caller's frame (one level up in the stack)
    frame = inspect.currentframe()

    for _ in range(levels_up):
        if frame is None:
            break
        frame = frame.f_back

    function_args, varargs, varkw, values = inspect.getargvalues(frame)

    # handle explicit function arguments
    for arg in function_args:
        params[arg] = values[arg]
    
    # handle *args
    if varargs:
        params["*args"] = values[varargs]
    
    # handle **kwargs
    if varkw:
        for key, value in values[varkw].items():
            params[key] = value
    
    return params


def save_params_to_config_file(out_file="run_config.json"):
    out_file_directory = os.path.dirname(out_file)
    if not out_file_directory:
        out_file_directory = "."
    else:
        os.makedirs(out_file_directory, exist_ok=True)

    # Collect parameters in a dictionary
    params = make_function_parameter_to_value_dict(levels_up = 2)

    # Write to JSON
    with open(out_file, "w") as file:
        json.dump(params, file, indent=4)


def is_directory_effectively_empty(directory_path):
    # Get all non-hidden entries, excluding system files like `.DS_Store`
    entries = [
        entry for entry in os.listdir(directory_path)
        if entry not in {".DS_Store"} and not entry.startswith(".")
    ]
    return len(entries) == 0

def group_items(file_list, group_size=2):
    if len(file_list) % group_size != 0:
        raise ValueError(f"The list length must be divisible by {group_size} to form groups.")
    return [tuple(file_list[i:i + group_size]) for i in range(0, len(file_list), group_size)]