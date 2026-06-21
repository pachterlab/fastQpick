import os
import tempfile
import pytest
from fastQpick import fastQpick
from fastQpick.utils import read_fastq, count_reads, parse_seed
from fastQpick.main import insert_seed_suffix
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




# A larger file is used for the one-pass tests so that the (random) output size concentrates
# tightly enough around its expectation for a loose statistical check to be stable.
@pytest.fixture
def temp_large_fastq_file():
    n = 20000
    lines = []
    for i in range(n):
        lines.append(f"@Header{i}\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n")
    with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".fastq") as temp_file:
        temp_file.write("".join(lines))
        temp_file.seek(0)
        yield temp_file.name
    os.remove(temp_file.name)


@pytest.fixture
def temp_large_paired_fastq_files():
    n = 20000
    lines_1 = "".join(f"@Header{i}_1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n" for i in range(n))
    lines_2 = "".join(f"@Header{i}_2\nTTTTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n" for i in range(n))
    with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".fastq") as temp_file1, \
         tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".fastq") as temp_file2:
        temp_file1.write(lines_1)
        temp_file2.write(lines_2)
        temp_file1.seek(0)
        temp_file2.seek(0)
        yield [temp_file1.name, temp_file2.name]
    os.remove(temp_file1.name)
    os.remove(temp_file2.name)


def test_one_pass_bernoulli_without_replacement(temp_large_fastq_file):
    # Without replacement the one-pass sampler keeps each read with probability `fraction`
    # (a Bernoulli draw), so every output read is unique and the count is Binomial(n, fraction).
    fraction = 0.5
    n = count_reads(temp_large_fastq_file)
    with tempfile.TemporaryDirectory() as temp_output_dir:
        fastQpick(input_files=temp_large_fastq_file, fraction=fraction, seed=42,
                  output_dir=temp_output_dir, without_replacement=True, one_pass=True, overwrite=True, verbose=False,
                  disable_gzip=True)
        output_fastq_file = os.path.join(temp_output_dir, os.path.basename(temp_large_fastq_file))

        input_fastq_dict = make_fastq_dict(temp_large_fastq_file)
        validate_fastq_format(output_fastq_file, ground_truth=input_fastq_dict)

        num_out = count_reads(output_fastq_file)
        num_unique = count_number_of_unique_headers(output_fastq_file)
        assert num_unique == num_out, "without replacement all output reads must be unique"

        # Loose check: |count - mean| within 6 standard deviations of Binomial(n, f).
        mean = fraction * n
        std = (n * fraction * (1 - fraction)) ** 0.5
        assert abs(num_out - mean) < 6 * std, f"output size {num_out} far from expected {mean:.0f}"


def test_one_pass_poisson_with_replacement(temp_large_fastq_file):
    # With replacement each read's multiplicity is Poisson(fraction); oversampling (fraction > 1)
    # must therefore produce duplicate reads and an output larger than the input.
    fraction = 2.0
    n = count_reads(temp_large_fastq_file)
    with tempfile.TemporaryDirectory() as temp_output_dir:
        fastQpick(input_files=temp_large_fastq_file, fraction=fraction, seed=42,
                  output_dir=temp_output_dir, without_replacement=False, unique_headers=False,
                  one_pass=True, overwrite=True, verbose=False, disable_gzip=True)
        output_fastq_file = os.path.join(temp_output_dir, os.path.basename(temp_large_fastq_file))

        num_out = count_reads(output_fastq_file)
        num_unique = count_number_of_unique_headers(output_fastq_file)
        assert num_unique < num_out, "oversampling with replacement must yield duplicate reads"

        mean = fraction * n
        std = (fraction * n) ** 0.5  # Poisson variance == mean
        assert abs(num_out - mean) < 6 * std, f"output size {num_out} far from expected {mean:.0f}"


def test_one_pass_pairwise_agreement(temp_large_paired_fastq_files):
    # Grouped files must remain synchronized in one-pass mode: both mates draw identical
    # multiplicities from a shared per-group sub-seed.
    fraction = 0.6
    with tempfile.TemporaryDirectory() as temp_output_dir:
        fastQpick(input_files=temp_large_paired_fastq_files, fraction=fraction, seed=42,
                  output_dir=temp_output_dir, file_group_size=2, without_replacement=True, one_pass=True,
                  overwrite=True, verbose=False, disable_gzip=True)

        out1 = count_reads(os.path.join(temp_output_dir, os.path.basename(temp_large_paired_fastq_files[0])))
        out2 = count_reads(os.path.join(temp_output_dir, os.path.basename(temp_large_paired_fastq_files[1])))
        assert out1 == out2, "grouped files must have equal output sizes when synchronized"

        check_pairwise_agreement(temp_paired_fastq_files=temp_large_paired_fastq_files,
                                 temp_output_dir=temp_output_dir, gzip_output=False)


def test_one_pass_is_deterministic(temp_large_fastq_file):
    # The same seed must reproduce byte-identical output across runs.
    fraction = 0.4
    with tempfile.TemporaryDirectory() as dir1, tempfile.TemporaryDirectory() as dir2:
        for out_dir in (dir1, dir2):
            fastQpick(input_files=temp_large_fastq_file, fraction=fraction, seed=123,
                      output_dir=out_dir, without_replacement=False, one_pass=True, overwrite=True, verbose=False,
                      disable_gzip=True)
        base = os.path.basename(temp_large_fastq_file)
        with open(os.path.join(dir1, base)) as f1, open(os.path.join(dir2, base)) as f2:
            assert f1.read() == f2.read(), "one-pass output must be deterministic for a fixed seed"


def test_occurrence_dtype_matches_data():
    # The occurrence vector must be sized to the realized maximum count, not the sample size.
    # For these moderate fractions the maximum count is small, so the dtype should be uint8.
    import numpy as np
    from fastQpick.main import make_occurrence_list, smallest_uint_dtype

    assert smallest_uint_dtype(0) == np.uint8
    assert smallest_uint_dtype(255) == np.uint8
    assert smallest_uint_dtype(256) == np.uint16
    assert smallest_uint_dtype(70000) == np.uint32

    n, m = 200000, 200000  # fraction = 1, dense regime
    for low_memory in (False, True):
        occ = make_occurrence_list("f", 0, n, m, replacement=True, low_memory=low_memory,
                                   rng=np.random.default_rng(0), verbose=False)
        assert occ.dtype == np.uint8, f"expected uint8, got {occ.dtype} (low_memory={low_memory})"
        # dtype must actually hold the data, and the total must equal the number sampled
        assert occ.max() <= np.iinfo(occ.dtype).max
        assert int(occ.sum()) == m, f"occurrence total {int(occ.sum())} != sampled {m}"

    # Without replacement every count is 0 or 1.
    occ = make_occurrence_list("f", 0, n, m // 2, replacement=False, low_memory=False,
                               rng=np.random.default_rng(0), verbose=False)
    assert occ.max() == 1 and int(occ.sum()) == m // 2


def test_default_mode_is_reproducible(temp_large_fastq_file):
    # The default (two-pass, numpy) path must be reproducible for a fixed seed.
    fraction = 0.5
    with tempfile.TemporaryDirectory() as dir1, tempfile.TemporaryDirectory() as dir2:
        for out_dir in (dir1, dir2):
            fastQpick(input_files=temp_large_fastq_file, fraction=fraction, seed=7,
                      output_dir=out_dir, without_replacement=False, overwrite=True, verbose=False,
                      disable_gzip=True)
        base = os.path.basename(temp_large_fastq_file)
        with open(os.path.join(dir1, base)) as f1, open(os.path.join(dir2, base)) as f2:
            assert f1.read() == f2.read(), "default-mode output must be deterministic for a fixed seed"


def test_parse_seed():
    # single int / single token string
    assert parse_seed(42) == [42]
    assert parse_seed("42") == [42]
    assert parse_seed("1-5") == [1, 2, 3, 4, 5]
    assert parse_seed("7-7") == [7]
    assert parse_seed(" 3 - 4 ") == [3, 4]

    # iterables of mixed ints and range strings
    assert parse_seed([42, 43, 44]) == [42, 43, 44]
    assert parse_seed(["42", "43"]) == [42, 43]
    assert parse_seed([1, 2, "5-7"]) == [1, 2, 5, 6, 7]
    assert parse_seed((1, "3-4")) == [1, 3, 4]
    assert parse_seed(range(1, 4)) == [1, 2, 3]

    with pytest.raises(ValueError):
        parse_seed("5-1")  # end less than start
    with pytest.raises(ValueError):
        parse_seed("a-b")  # non-integer range
    with pytest.raises(ValueError):
        parse_seed("foo")  # non-integer seed
    with pytest.raises(ValueError):
        parse_seed("42,43,44")  # comma syntax no longer supported

def test_seed_range_produces_multiple_outputs(temp_fastq_file):
    fraction = 0.6
    seed = "1-3"
    gzip_output = False

    with tempfile.TemporaryDirectory() as temp_output_dir:
        fastQpick(input_files=temp_fastq_file,
                fraction=fraction,
                seed=seed,
                output_dir=temp_output_dir,
                disable_gzip=not gzip_output,
                file_group_size=1,
                without_replacement=True,
                overwrite=True
                )

        # One distinct output file per seed should be present, suffixed with the seed
        output_files = sorted(f for f in os.listdir(temp_output_dir) if f.endswith(".fastq"))
        expected = sorted(insert_seed_suffix(os.path.basename(temp_fastq_file), s) for s in (1, 2, 3))
        assert output_files == expected, f"Expected {expected}, got {output_files}"


def test_num_samples_produces_multiple_outputs(temp_fastq_file):
    # num_samples derives consecutive seeds from a single base seed, producing one output per replicate.
    fraction = 0.6
    seed = 5
    num_samples = 3

    with tempfile.TemporaryDirectory() as temp_output_dir:
        fastQpick(input_files=temp_fastq_file,
                fraction=fraction,
                seed=seed,
                num_samples=num_samples,
                output_dir=temp_output_dir,
                disable_gzip=True,
                file_group_size=1,
                without_replacement=True,
                overwrite=True
                )

        # num_samples consecutive seeds (5, 6, 7) each yield one suffixed output file
        output_files = sorted(f for f in os.listdir(temp_output_dir) if f.endswith(".fastq"))
        expected = sorted(insert_seed_suffix(os.path.basename(temp_fastq_file), s) for s in (5, 6, 7))
        assert output_files == expected, f"Expected {expected}, got {output_files}"

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
    group_size = 1
    replacement = False
    
    with tempfile.TemporaryDirectory() as temp_output_dir:
        fastQpick(input_files=temp_fastq_file,
                fraction=fraction,
                seed=seed,
                output_dir=temp_output_dir,
                disable_gzip=not gzip_output,
                file_group_size=group_size,
                without_replacement=not replacement,
                unique_headers=False,
                overwrite=True
                )
        
        run_all_single_file_tests(temp_output_dir=temp_output_dir, temp_fastq_file=temp_fastq_file, gzip_output=gzip_output, fraction=fraction, replacement=replacement)

def test_single_file_bootstrapped(temp_fastq_file):
    fraction = 1
    seed = 42
    gzip_output = False
    group_size = 1
    replacement = True
    
    with tempfile.TemporaryDirectory() as temp_output_dir:
        fastQpick(input_files=temp_fastq_file,
                fraction=fraction,
                seed=seed,
                output_dir=temp_output_dir,
                disable_gzip=not gzip_output,
                file_group_size=group_size,
                without_replacement=not replacement,
                unique_headers=False,
                overwrite=True
                )
        
        run_all_single_file_tests(temp_output_dir=temp_output_dir, temp_fastq_file=temp_fastq_file, gzip_output=gzip_output, fraction=fraction, replacement=replacement)

        # st()

def test_single_file_oversampled(temp_fastq_file):
    fraction = 3
    seed = 42
    gzip_output = False
    group_size = 1
    replacement = True
    
    with tempfile.TemporaryDirectory() as temp_output_dir:
        fastQpick(input_files=temp_fastq_file,
                fraction=fraction,
                seed=seed,
                output_dir=temp_output_dir,
                disable_gzip=not gzip_output,
                file_group_size=group_size,
                without_replacement=not replacement,
                unique_headers=False,
                overwrite=True
                )
        
        run_all_single_file_tests(temp_output_dir=temp_output_dir, temp_fastq_file=temp_fastq_file, gzip_output=gzip_output, fraction=fraction, replacement=replacement)

        # st()
        
def test_single_gzipped(temp_fastq_file):
    fraction = 0.6
    seed = 42
    gzip_output = True
    group_size = 1
    replacement = False
    
    with tempfile.TemporaryDirectory() as temp_output_dir:
        fastQpick(input_files=temp_fastq_file,
                fraction=fraction,
                seed=seed,
                output_dir=temp_output_dir,
                disable_gzip=not gzip_output,
                file_group_size=group_size,
                without_replacement=not replacement,
                unique_headers=False,
                overwrite=True
                )
        
        run_all_single_file_tests(temp_output_dir=temp_output_dir, temp_fastq_file=temp_fastq_file, gzip_output=gzip_output, fraction=fraction, replacement=replacement)

        # st()


def test_paired_files(temp_paired_fastq_files):
    fraction = 0.75
    seed = 42
    gzip_output = False
    group_size = 2
    replacement = False
    
    with tempfile.TemporaryDirectory() as temp_output_dir:
        fastQpick(input_files=temp_paired_fastq_files,
                fraction=fraction,
                seed=seed,
                output_dir=temp_output_dir,
                disable_gzip=not gzip_output,
                file_group_size=group_size,
                without_replacement=not replacement,
                unique_headers=False,
                overwrite=True
                )
        
        for fastq_file in temp_paired_fastq_files:
            run_all_single_file_tests(temp_output_dir=temp_output_dir, temp_fastq_file=fastq_file, gzip_output=gzip_output, fraction=fraction, replacement=replacement)

        check_pairwise_agreement(temp_paired_fastq_files=temp_paired_fastq_files, temp_output_dir=temp_output_dir, gzip_output=gzip_output)

        # st()

def test_paired_files_bootstrapped(temp_paired_fastq_files):
    fraction = 1
    seed = 42
    gzip_output = False
    group_size = 2
    replacement = True
    
    with tempfile.TemporaryDirectory() as temp_output_dir:
        fastQpick(input_files=temp_paired_fastq_files,
                fraction=fraction,
                seed=seed,
                output_dir=temp_output_dir,
                disable_gzip=not gzip_output,
                file_group_size=group_size,
                without_replacement=not replacement,
                unique_headers=False,
                overwrite=True
                )
        
        for fastq_file in temp_paired_fastq_files:
            run_all_single_file_tests(temp_output_dir=temp_output_dir, temp_fastq_file=fastq_file, gzip_output=gzip_output, fraction=fraction, replacement=replacement)

        check_pairwise_agreement(temp_paired_fastq_files=temp_paired_fastq_files, temp_output_dir=temp_output_dir, gzip_output=gzip_output)

        # st()