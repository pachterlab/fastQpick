import asyncio
import gzip
import os
import pytest

pytest.importorskip("mcp")

from fastQpick.mcp_server import count_fastq_reads, list_fastq_files, sample_fastq
from fastQpick.utils import count_reads


def write_fastq(path, num_reads, prefix):
    with open(path, "w") as f:
        for i in range(num_reads):
            f.write(f"@{prefix}{i}\nACGTACGTAC\n+\nIIIIIIIIII\n")


@pytest.fixture
def paired_fastq_dir(tmp_path):
    write_fastq(tmp_path / "sample_R1.fastq", 200, "read")
    write_fastq(tmp_path / "sample_R2.fastq", 200, "read")
    return tmp_path


def test_list_and_count(paired_fastq_dir):
    files = list_fastq_files(str(paired_fastq_dir))
    assert [os.path.basename(f) for f in files] == ["sample_R1.fastq", "sample_R2.fastq"]
    counts = asyncio.run(count_fastq_reads(files))
    assert counts == {files[0]: 200, files[1]: 200}


def test_sample_fastq_paired_subsample(paired_fastq_dir, tmp_path):
    files = list_fastq_files(str(paired_fastq_dir))
    out = tmp_path / "out"
    result = asyncio.run(sample_fastq(files, fraction=0.25, output_dir=str(out), without_replacement=True,
                                      file_group_size=2, read_counts=[200], num_samples=2))
    assert len(result["output_files"]) == 4
    assert os.path.isfile(result["config_file"])
    for path in result["output_files"]:
        assert count_reads(path) == 50
    r1, r2 = (sorted(p for p in result["output_files"] if mate in p) for mate in ("_R1", "_R2"))
    for a, b in zip(r1, r2):  # mates stay synchronized
        with gzip.open(a, "rt") as fa, gzip.open(b, "rt") as fb:
            assert fa.readlines()[::4] == fb.readlines()[::4]


def test_sample_fastq_reports_errors(paired_fastq_dir, tmp_path):
    files = list_fastq_files(str(paired_fastq_dir))
    with pytest.raises(RuntimeError, match="fastQpick exited"):  # wrong read count
        asyncio.run(sample_fastq(files[:1], output_dir=str(tmp_path / "bad"), read_counts=[199]))
    with pytest.raises(ValueError):
        asyncio.run(sample_fastq(["-"], output_dir=str(tmp_path / "stdin")))
