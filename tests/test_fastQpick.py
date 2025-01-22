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
@Header5_2
CCCCCCGGGGGGGGGGGGGGG
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
        yield temp_file1.name, temp_file2.name  # Yield the paths of both files

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

def check_for_read_uniqueness(file_path):
    headers = set()
    for header, _, _, _ in read_fastq(file_path, include_plus_line=True):
        assert header not in headers, f"Non-unique header found: {header}"
        headers.add(header)

def make_fastq_dict(file_path):
    fastq_dict = {}
    for header, seq, plus_line, qual in read_fastq(file_path, include_plus_line=True):
        fastq_dict[header] = {}
        fastq_dict[header]["sequence"] = seq
        fastq_dict[header]["plus_line"] = plus_line
        fastq_dict[header]["quality"] = qual
    return fastq_dict

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
        
        # Assert that the output directory exists
        assert os.path.exists(temp_output_dir), "Output directory does not exist!"

        # Optionally, verify the output files
        output_files = os.listdir(temp_output_dir)
        assert len(output_files) > 0, "No output files were created!"

        file_base_name = os.path.basename(temp_fastq_file)
        output_fastq_file = os.path.join(temp_output_dir, file_base_name)

        input_fastq_dict = make_fastq_dict(temp_fastq_file)
        validate_fastq_format(output_fastq_file, ground_truth=input_fastq_dict)

        output_is_gzipped = is_gzipped(output_fastq_file)
        assert output_is_gzipped == gzip_output, f"Gzipped output - expected: {gzip_output}; got: {output_is_gzipped}"

        num_reads_truth = count_reads(temp_fastq_file)
        num_reads_output = count_reads(output_fastq_file)

        assert num_reads_output == num_reads_truth * fraction, f"Number of reads mismatch - expected: {num_reads_truth * fraction}; got: {num_reads_output}"

        if not replacement:
            check_for_read_uniqueness(output_fastq_file)

        st()
        

# Example test using the paired FASTQ files
def test_paired_files(temp_paired_fastq_files):
    pass


