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

def read_headers(file_path):
    return [header for header, _, _, _ in read_fastq(file_path, include_plus_line=True)]


@pytest.mark.parametrize("mode", [dict(), dict(low_memory=True), dict(one_pass=True)])
def test_collapse_duplicates_matches_expanded_output(temp_large_fastq_file, mode):
    # The collapsed output must encode exactly the same multiset of reads as the ordinary output
    # for the same seed: each sampled read once, with its multiplicity in a ";size=<count>" tag.
    base = os.path.basename(temp_large_fastq_file)
    with tempfile.TemporaryDirectory() as dir_expanded, tempfile.TemporaryDirectory() as dir_collapsed:
        common = dict(input_files=temp_large_fastq_file, fraction=1.0, seed=7, overwrite=True, verbose=False, disable_gzip=True, **mode)
        fastQpick(output_dir=dir_expanded, unique_headers=False, **common)
        fastQpick(output_dir=dir_collapsed, collapse_duplicates=True, **common)

        expanded_counts = {}
        for header in read_headers(os.path.join(dir_expanded, base)):
            expanded_counts[header] = expanded_counts.get(header, 0) + 1

        collapsed_counts = {}
        for header in read_headers(os.path.join(dir_collapsed, base)):
            name, size = header.rsplit(";size=", 1)
            assert name not in collapsed_counts, "a collapsed read must be written only once"
            collapsed_counts[name] = int(size)

        assert collapsed_counts == expanded_counts
        assert max(collapsed_counts.values()) > 1, "a full-size bootstrap should contain duplicated reads"


@pytest.mark.parametrize("mode", [dict(), dict(low_memory=True), dict(one_pass=True)])
@pytest.mark.parametrize("without_replacement,fraction", [(False, 1.0), (True, 0.3)])
def test_oob_is_complement_of_sample(temp_large_fastq_file, mode, without_replacement, fraction):
    # The out-of-bag file must contain exactly the input reads that are absent from the sample.
    base = os.path.basename(temp_large_fastq_file)
    with tempfile.TemporaryDirectory() as temp_output_dir:
        fastQpick(input_files=temp_large_fastq_file, fraction=fraction, seed=11, output_dir=temp_output_dir,
                  without_replacement=without_replacement, unique_headers=False, oob=True, overwrite=True,
                  verbose=False, disable_gzip=True, **mode)
        oob_file = os.path.join(temp_output_dir, base.replace(".fastq", ".oob.fastq"))
        validate_fastq_format(oob_file, ground_truth=make_fastq_dict(temp_large_fastq_file))

        sampled = set(read_headers(os.path.join(temp_output_dir, base)))
        oob_headers = read_headers(oob_file)
        all_headers = set(read_headers(temp_large_fastq_file))

        assert len(oob_headers) == len(set(oob_headers)), "out-of-bag reads must be written once"
        assert sampled.isdisjoint(oob_headers)
        assert sampled | set(oob_headers) == all_headers
        if not without_replacement:
            # a full-size bootstrap leaves ~1/e of the reads out of bag
            assert abs(len(oob_headers) / len(all_headers) - 0.3679) < 0.02


def test_oob_pairwise_agreement(temp_large_paired_fastq_files):
    # Out-of-bag files of grouped inputs must stay synchronized, like the samples themselves.
    with tempfile.TemporaryDirectory() as temp_output_dir:
        fastQpick(input_files=temp_large_paired_fastq_files, fraction=1.0, seed=3, output_dir=temp_output_dir,
                  file_group_size=2, oob=True, collapse_duplicates=True, overwrite=True, verbose=False, disable_gzip=True)
        outputs = [os.path.join(temp_output_dir, os.path.basename(f)) for f in temp_large_paired_fastq_files]
        for suffix in (".fastq", ".oob.fastq"):
            headers_1 = read_headers(outputs[0].replace(".fastq", suffix))
            headers_2 = read_headers(outputs[1].replace(".fastq", suffix))
            assert [h.replace("_1", "", 1) for h in headers_1] == [h.replace("_2", "", 1) for h in headers_2]


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

# --- streaming input (reading the library from a pipe) ------------------------------

def _records_bytes(n, rl=4):
    return b"".join(f"@r{i} desc\n{'ACGT' * rl}\n+\n{'I' * (4 * rl)}\n".encode() for i in range(n))


@pytest.mark.parametrize("chunk_size", [1 << 22, 64, 7])
def test_stream_fastq_records_keeps_final_record(chunk_size):
    # Regression test for the reason this reader exists: pyfastx.Fastx silently drops the
    # last record when its input is a pipe. Small chunk sizes exercise the carry logic that
    # splices records straddling a chunk boundary.
    import io
    from fastQpick.main import stream_fastq_records

    raw = _records_bytes(500)
    records = list(stream_fastq_records(io.BytesIO(raw), chunk_size=chunk_size))
    assert len(records) == 500
    assert records[0] == ("r0 desc", "ACGTACGTACGTACGT", "I" * 16)
    assert records[-1] == ("r499 desc", "ACGTACGTACGTACGT", "I" * 16)


def test_stream_fastq_records_without_trailing_newline():
    import io
    from fastQpick.main import stream_fastq_records

    raw = _records_bytes(10).rstrip(b"\n")
    assert len(list(stream_fastq_records(io.BytesIO(raw)))) == 10


def test_stream_fastq_records_rejects_truncated_record():
    import io
    from fastQpick.main import stream_fastq_records

    raw = _records_bytes(10) + b"@r10 desc\nACGT\n"
    with pytest.raises(ValueError, match="Truncated FASTQ record"):
        list(stream_fastq_records(io.BytesIO(raw)))


def test_stream_fastq_records_reads_gzip_stream():
    import gzip as _gzip
    import io
    from fastQpick.main import stream_fastq_records

    raw = _gzip.compress(_records_bytes(300))
    assert len(list(stream_fastq_records(io.BytesIO(raw)))) == 300


@pytest.mark.parametrize("compress_stdin", [False, True])
def test_stdin_matches_file_input(tmp_path, temp_large_fastq_file, compress_stdin):
    # A library piped in must produce exactly the output the same library on disk produces.
    import gzip as _gzip
    import shutil
    import subprocess

    exe = shutil.which("fastQpick")
    if exe is None:
        pytest.skip("fastQpick console script is not on PATH")

    from_file = tmp_path / "from_file"
    from_pipe = tmp_path / "from_pipe"
    common = ["-f", "1", "-s", "42", "-p", "-z", "-w", "-q", "-o"]

    subprocess.run([exe, *common, str(from_file), temp_large_fastq_file], check=True)

    payload = open(temp_large_fastq_file, "rb").read()
    if compress_stdin:
        payload = _gzip.compress(payload)
    subprocess.run([exe, *common, str(from_pipe), "-"], input=payload, check=True)

    expected = (from_file / os.path.basename(temp_large_fastq_file)).read_bytes()
    assert (from_pipe / "stdin.fastq").read_bytes() == expected


def test_stdin_requires_one_pass(tmp_path, temp_large_fastq_file):
    # The two-pass modes need to read the library twice, which a stream cannot do.
    with pytest.raises(ValueError, match="requires one_pass"):
        fastQpick(input_files="-", fraction=1.0, seed=42, output_dir=str(tmp_path / "out"),
                  one_pass=False, disable_gzip=True, overwrite=True, verbose=False)


def test_stdin_rejects_additional_inputs(tmp_path, temp_large_fastq_file):
    with pytest.raises(ValueError, match="cannot be combined with other input files"):
        fastQpick(input_files=["-", temp_large_fastq_file], fraction=1.0, seed=42,
                  output_dir=str(tmp_path / "out"), one_pass=True, disable_gzip=True,
                  overwrite=True, verbose=False)


def test_gzip_output_round_trips(tmp_path, temp_large_fastq_file):
    # The output writer uses multithreaded ISA-L deflate when available; the bytes it emits
    # must still be ordinary gzip that any reader can open.
    import gzip as _gzip

    gz_dir = tmp_path / "gz"
    plain_dir = tmp_path / "plain"
    for out, disable in ((gz_dir, False), (plain_dir, True)):
        fastQpick(input_files=temp_large_fastq_file, fraction=1.0, seed=42, output_dir=str(out),
                  one_pass=True, disable_gzip=disable, overwrite=True, verbose=False)

    base = os.path.basename(temp_large_fastq_file)
    assert _gzip.open(gz_dir / f"{base}.gz", "rt").read() == (plain_dir / base).read_text()


# --- user-supplied read counts (skipping the counting pass) --------------------------

@pytest.fixture
def fresh_length_dict(monkeypatch):
    # fastq_to_length_dict is a module global that caches counts across calls; isolate each test.
    import fastQpick.main as main_module
    monkeypatch.setattr(main_module, "fastq_to_length_dict", {})
    return main_module


@pytest.mark.parametrize("read_counts", [20000, [20000], "dict"])
def test_read_counts_skips_counting_and_matches(tmp_path, temp_large_fastq_file, fresh_length_dict, monkeypatch, read_counts):
    common = dict(input_files=temp_large_fastq_file, fraction=1.0, seed=42, disable_gzip=True, overwrite=True, verbose=False)
    fastQpick(output_dir=str(tmp_path / "counted"), **common)

    fresh_length_dict.fastq_to_length_dict.clear()
    def fail(*args, **kwargs):
        raise AssertionError("count_reads must not run when read_counts is given")
    monkeypatch.setattr(fresh_length_dict, "count_reads", fail)
    if read_counts == "dict":
        read_counts = {temp_large_fastq_file: 20000}
    fastQpick(output_dir=str(tmp_path / "given"), read_counts=read_counts, **common)

    base = os.path.basename(temp_large_fastq_file)
    assert (tmp_path / "given" / base).read_bytes() == (tmp_path / "counted" / base).read_bytes()


@pytest.mark.parametrize("wrong_count", [19999, 20001])
@pytest.mark.parametrize("fraction, low_memory, oob", [(1.0, False, False), (1.0, True, False), (1.0, False, True), (0.001, False, False)])
def test_wrong_read_count_raises(tmp_path, temp_large_fastq_file, fresh_length_dict, wrong_count, fraction, low_memory, oob):
    # Covers the dense-array and Counter (fraction=0.001) occurrence lists and the tagged writer (oob).
    with pytest.raises(ValueError, match="read count used for sampling"):
        fastQpick(input_files=temp_large_fastq_file, fraction=fraction, seed=42, output_dir=str(tmp_path / "out"),
                  low_memory=low_memory, oob=oob, read_counts=wrong_count, disable_gzip=True, overwrite=True, verbose=False)


@pytest.mark.parametrize("read_counts", [[20000], [20000, 20000]])
def test_read_counts_per_group_or_per_file(tmp_path, temp_large_paired_fastq_files, fresh_length_dict, read_counts):
    fastQpick(input_files=temp_large_paired_fastq_files, fraction=1.0, seed=42, output_dir=str(tmp_path / "out"),
              file_group_size=2, read_counts=read_counts, unique_headers=False, disable_gzip=True, overwrite=True, verbose=False)
    check_pairwise_agreement(temp_paired_fastq_files=temp_large_paired_fastq_files, temp_output_dir=str(tmp_path / "out"), gzip_output=False)


@pytest.mark.parametrize("read_counts, message", [
    ([20000, 19999], "differ"),
    ([20000, 20000, 20000], "one count per file or one per group"),
    ([-1], "non-negative integer"),
])
def test_invalid_read_counts_rejected(tmp_path, temp_large_paired_fastq_files, fresh_length_dict, read_counts, message):
    with pytest.raises(ValueError, match=message):
        fastQpick(input_files=temp_large_paired_fastq_files, fraction=1.0, seed=42, output_dir=str(tmp_path / "out"),
                  file_group_size=2, read_counts=read_counts, disable_gzip=True, overwrite=True, verbose=False)


def test_read_counts_cli(tmp_path, temp_large_paired_fastq_files):
    import shutil
    import subprocess

    exe = shutil.which("fastQpick")
    if exe is None:
        pytest.skip("fastQpick console script is not on PATH")

    common = ["-f", "1", "-s", "42", "-g", "2", "-z", "-w", "-q"]
    subprocess.run([exe, *common, "-o", str(tmp_path / "counted"), *temp_large_paired_fastq_files], check=True)
    subprocess.run([exe, *common, "--read-counts", "20000,20000", "-o", str(tmp_path / "given"), *temp_large_paired_fastq_files], check=True)
    for path in temp_large_paired_fastq_files:
        base = os.path.basename(path)
        assert (tmp_path / "given" / base).read_bytes() == (tmp_path / "counted" / base).read_bytes()

    wrong = subprocess.run([exe, *common, "--read-counts", "123", "-o", str(tmp_path / "wrong"), *temp_large_paired_fastq_files], capture_output=True, text=True)
    assert wrong.returncode != 0 and "read count used for sampling" in wrong.stderr


# --- threads (parallel counting, files, and replicates) --------------------------------

def test_split_thread_budget():
    from fastQpick.main import split_thread_budget
    assert split_thread_budget(88, 1) == (1, 4)   # one job: all four deflate threads
    assert split_thread_budget(8, 2) == (2, 3)    # 4 threads per worker: 1 record loop + 3 deflate
    assert split_thread_budget(8, 8) == (8, 0)    # one thread per worker: compress inline
    assert split_thread_budget(4, 100) == (4, 0)  # never more workers than threads
    assert split_thread_budget(1, 1) == (1, 0)


@pytest.mark.parametrize("mode", [dict(), dict(low_memory=True), dict(one_pass=True)])
@pytest.mark.parametrize("disable_gzip", [True, False])
def test_output_independent_of_threads(tmp_path, temp_large_paired_fastq_files, fresh_length_dict, mode, disable_gzip):
    # Parallel counting, and writing several files and seeds at once, must not change the output.
    import gzip as _gzip
    outputs = {}
    for threads in (1, 3, 8):
        fresh_length_dict.fastq_to_length_dict.clear()
        out_dir = tmp_path / f"t{threads}"
        fastQpick(input_files=temp_large_paired_fastq_files, fraction=1.0, seed="1-3", file_group_size=2,
                  output_dir=str(out_dir), threads=threads, disable_gzip=disable_gzip, overwrite=True, verbose=False, **mode)
        opener = open if disable_gzip else _gzip.open
        outputs[threads] = {name: opener(out_dir / name, "rb").read() for name in sorted(os.listdir(out_dir)) if ".fastq" in name}
    assert len(outputs[1]) == 6  # 2 mates x 3 seeds
    assert outputs[1] == outputs[3] == outputs[8]


def test_resolve_threads(monkeypatch):
    import fastQpick.main as main_module
    monkeypatch.setattr(main_module, "available_cpus", lambda: 88)
    assert main_module.resolve_threads(None) == 4
    assert main_module.resolve_threads(16) == 16  # explicit values are not capped
    monkeypatch.setattr(main_module, "available_cpus", lambda: 2)
    assert main_module.resolve_threads(None) == 2


def test_invalid_threads_rejected(tmp_path, temp_large_fastq_file):
    with pytest.raises(ValueError, match="threads"):
        fastQpick(input_files=temp_large_fastq_file, fraction=0.5, output_dir=str(tmp_path / "out"), threads=0, overwrite=True, verbose=False)
