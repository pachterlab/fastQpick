import os
import tempfile
import pytest
from fastQpick import fastQpick
from fastQpick.utils import read_fastq, count_reads
from pdb import set_trace as st

@pytest.fixture
def temp_fastq_file():
    content = """@Header1
AAAAAAAAAAAAAAAAAAAAA
+
IIIIIIIIIIIIIIIIIIIII
@Header2
CCCCCCCCCCCCCCCCCCCCC
+
IIIIIIIIIIIIIIIIIIIII
@Header3
GGGGGGGGGGGGGGGGGGGGG
+
IIIIIIIIIIIIIIIIIIIII
@Header4
TTTTTTTTTTTTTTTTTTTTT
+
IIIIIIIIIIIIIIIIIIIII
@Header5
AAAAAAAAAAAAACCCCCCCC
+
IIIIIIIIIIIIIIIIIIIII
"""
    # Create a temporary file
    with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".fastq") as temp_file:
        temp_file.write(content)
        temp_file.seek(0)  # Move to the start of the file
        yield temp_file.name  # Provide the file path to the test

    # Cleanup after the test
    os.remove(temp_file.name)

# Fixture to create two temporary FASTQ files
@pytest.fixture
def temp_paired_fastq_files():
    content_1 = """@Header1_1
AAAAAAAAAAAAAAAAAAAAA
+
IIIIIIIIIIIIIIIIIIIII
@Header2_1
CCCCCCCCCCCCCCCCCCCCC
+
IIIIIIIIIIIIIIIIIIIII
@Header3_1
GGGGGGGGGGGGGGGGGGGGG
+
IIIIIIIIIIIIIIIIIIIII
@Header4_1
TTTTTTTTTTTTTTTTTTTTT
+
IIIIIIIIIIIIIIIIIIIII
"""
    
    content_2 = """@Header1_2
AAAAAAAACCCCCCCCCCCCC
+
IIIIIIIIIIIIIIIIIIIII
@Header2_2
AAAAAAAGGGGGGGGGGGGGG
+
IIIIIIIIIIIIIIIIIIIII
@Header3_2
AAAAAAATTTTTTTTTTTTTT
+
IIIIIIIIIIIIIIIIIIIII
@Header4_2
CCCCCCCCAAAAAAAAAAAAA
+
IIIIIIIIIIIIIIIIIIIII
"""

    # Create two temporary files
    with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".fastq") as temp_file1, \
         tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".fastq") as temp_file2:
        temp_file1.write(content_1)
        temp_file2.write(content_2)
        temp_file1.seek(0)
        temp_file2.seek(0)
        yield [temp_file1.name, temp_file2.name]  # Yield the paths of both files

    # Cleanup after the test
    os.remove(temp_file1.name)
    os.remove(temp_file2.name)




def is_gzipped(file_path):
    with open(file_path, "rb") as f:
        magic_number = f.read(2)
        return magic_number == b"\x1f\x8b"

def validate_fastq_format(file_path, ground_truth=None):
    for header, seq, plus_line, qual in read_fastq(file_path, include_plus_line=True):
        assert header.startswith("@"), f"Header does not start with '@': {header}"
        assert len(seq) == len(qual), f"Sequence and quality lengths do not match: {seq} {qual}"
        assert plus_line.startswith("+"), f"Plus line does not start with '+': {plus_line}"

        if ground_truth:
            assert header in ground_truth, f"Header not found in ground truth: {header}"
            assert seq == ground_truth[header]["sequence"], f"Sequence mismatch - expected: {seq}; got: {ground_truth[header]['sequence']}"
            assert plus_line == ground_truth[header]["plus_line"], f"Plus line mismatch - expected: {plus_line}; got: {ground_truth[header]['plus_line']}"
            assert qual == ground_truth[header]["quality"], f"Quality mismatch - expected: {qual}; got: {ground_truth[header]['quality']}"

def count_number_of_unique_headers(file_path):
    headers = set()
    for header, _, _, _ in read_fastq(file_path, include_plus_line=True):
        headers.add(header)
    return len(headers)

def make_fastq_dict(file_path):
    fastq_dict = {}
    for header, seq, plus_line, qual in read_fastq(file_path, include_plus_line=True):
        fastq_dict[header] = {}
        fastq_dict[header]["sequence"] = seq
        fastq_dict[header]["plus_line"] = plus_line
        fastq_dict[header]["quality"] = qual
    return fastq_dict

        
def check_pairwise_agreement(temp_paired_fastq_files, temp_output_dir, gzip_output):
    file1_base_name = os.path.basename(temp_paired_fastq_files[0])
    file2_base_name = os.path.basename(temp_paired_fastq_files[1])
    
    output_fastq_file1 = os.path.join(temp_output_dir, file1_base_name)
    output_fastq_file2 = os.path.join(temp_output_dir, file2_base_name)

    if gzip_output:
        output_fastq_file1 += ".gz"
        output_fastq_file2 += ".gz"

    for (header1, seq1, plus_line1, qual1), (header2, seq2, plus_line2, qual2) in zip(
        read_fastq(output_fastq_file1, include_plus_line=True), 
        read_fastq(output_fastq_file2, include_plus_line=True)
    ):
        # Split headers up to the last underscore
        split_header1 = header1.rsplit('_', 1)[0]
        split_header2 = header2.rsplit('_', 1)[0]

        # Assert that the two headers are equal
        assert split_header1 == split_header2, f"Headers do not match: {split_header1} != {split_header2}"

def run_all_single_file_tests(temp_output_dir, temp_fastq_file, gzip_output, fraction, replacement):
    # Assert that the output directory exists
        assert os.path.exists(temp_output_dir), "Output directory does not exist!"

        # Optionally, verify the output files
        output_files = os.listdir(temp_output_dir)
        assert len(output_files) > 0, "No output files were created!"

        file_base_name = os.path.basename(temp_fastq_file)
        output_fastq_file = os.path.join(temp_output_dir, file_base_name)

        if gzip_output:
            output_fastq_file += ".gz"

        input_fastq_dict = make_fastq_dict(temp_fastq_file)
        validate_fastq_format(output_fastq_file, ground_truth=input_fastq_dict)

        output_is_gzipped = is_gzipped(output_fastq_file)
        assert output_is_gzipped == gzip_output, f"Gzipped output - expected: {gzip_output}; got: {output_is_gzipped}"

        num_reads_truth = count_reads(temp_fastq_file)
        num_reads_output = count_reads(output_fastq_file)

        assert num_reads_output == num_reads_truth * fraction, f"Number of reads mismatch - expected: {num_reads_truth * fraction}; got: {num_reads_output}"

        num_unique_reads = count_number_of_unique_headers(output_fastq_file)

        if not replacement:
            assert num_unique_reads == num_reads_output, f"Number of unique reads mismatch - expected: {num_reads_output}; got: {num_unique_reads}"

        if replacement and fraction > 1:
            assert num_unique_reads < num_reads_output, f"Number of unique reads mismatch - expected: less than {num_reads_output}; got: {num_unique_reads}"

def test_single_file(temp_fastq_file):
    fraction = 0.6
    seed = 42
    gzip_output = False
    paired = False
    replacement = False
    
    with tempfile.TemporaryDirectory() as temp_output_dir:
        fastQpick(input_file_list=temp_fastq_file,
                fraction=fraction,
                seed=seed,
                output_dir=temp_output_dir,
                threads=1,
                gzip_output=gzip_output,
                paired=paired,
                replacement=replacement,
                overwrite=True
                )
        
        run_all_single_file_tests(temp_output_dir=temp_output_dir, temp_fastq_file=temp_fastq_file, gzip_output=gzip_output, fraction=fraction, replacement=replacement)

        st()

def test_single_file_bootstrapped(temp_fastq_file):
    fraction = 1
    seed = 42
    gzip_output = False
    paired = False
    replacement = True
    
    with tempfile.TemporaryDirectory() as temp_output_dir:
        fastQpick(input_file_list=temp_fastq_file,
                fraction=fraction,
                seed=seed,
                output_dir=temp_output_dir,
                threads=1,
                gzip_output=gzip_output,
                paired=paired,
                replacement=replacement,
                overwrite=True
                )
        
        run_all_single_file_tests(temp_output_dir=temp_output_dir, temp_fastq_file=temp_fastq_file, gzip_output=gzip_output, fraction=fraction, replacement=replacement)

        st()

def test_single_file_oversampled(temp_fastq_file):
    fraction = 3
    seed = 42
    gzip_output = False
    paired = False
    replacement = True
    
    with tempfile.TemporaryDirectory() as temp_output_dir:
        fastQpick(input_file_list=temp_fastq_file,
                fraction=fraction,
                seed=seed,
                output_dir=temp_output_dir,
                threads=1,
                gzip_output=gzip_output,
                paired=paired,
                replacement=replacement,
                overwrite=True
                )
        
        run_all_single_file_tests(temp_output_dir=temp_output_dir, temp_fastq_file=temp_fastq_file, gzip_output=gzip_output, fraction=fraction, replacement=replacement)

        st()
        
def test_single_gzipped(temp_fastq_file):
    fraction = 0.6
    seed = 42
    gzip_output = True
    paired = False
    replacement = False
    
    with tempfile.TemporaryDirectory() as temp_output_dir:
        fastQpick(input_file_list=temp_fastq_file,
                fraction=fraction,
                seed=seed,
                output_dir=temp_output_dir,
                threads=1,
                gzip_output=gzip_output,
                paired=paired,
                replacement=replacement,
                overwrite=True
                )
        
        run_all_single_file_tests(temp_output_dir=temp_output_dir, temp_fastq_file=temp_fastq_file, gzip_output=gzip_output, fraction=fraction, replacement=replacement)

        st()


def test_paired_files(temp_paired_fastq_files):
    fraction = 0.75
    seed = 42
    gzip_output = False
    paired = False
    replacement = False
    
    with tempfile.TemporaryDirectory() as temp_output_dir:
        fastQpick(input_file_list=temp_paired_fastq_files,
                fraction=fraction,
                seed=seed,
                output_dir=temp_output_dir,
                threads=1,
                gzip_output=gzip_output,
                paired=paired,
                replacement=replacement,
                overwrite=True
                )
        
        for fastq_file in temp_paired_fastq_files:
            run_all_single_file_tests(temp_output_dir=temp_output_dir, temp_fastq_file=fastq_file, gzip_output=gzip_output, fraction=fraction, replacement=replacement)

        check_pairwise_agreement(temp_paired_fastq_files=temp_paired_fastq_files, temp_output_dir=temp_output_dir, gzip_output=gzip_output)

        st()

def test_paired_files_bootstrapped(temp_paired_fastq_files):
    fraction = 1
    seed = 42
    gzip_output = False
    paired = False
    replacement = True
    
    with tempfile.TemporaryDirectory() as temp_output_dir:
        fastQpick(input_file_list=temp_paired_fastq_files,
                fraction=fraction,
                seed=seed,
                output_dir=temp_output_dir,
                threads=1,
                gzip_output=gzip_output,
                paired=paired,
                replacement=replacement,
                overwrite=True
                )
        
        for fastq_file in temp_paired_fastq_files:
            run_all_single_file_tests(temp_output_dir=temp_output_dir, temp_fastq_file=fastq_file, gzip_output=gzip_output, fraction=fraction, replacement=replacement)

        check_pairwise_agreement(temp_paired_fastq_files=temp_paired_fastq_files, temp_output_dir=temp_output_dir, gzip_output=gzip_output)

        st()