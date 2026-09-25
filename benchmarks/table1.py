#!/usr/bin/env python3
"""Benchmark driver for Table 1 of the fastQpick manuscript.

Runs every (condition, tool) cell of the table on one gzipped FASTQ library and appends
one row per run to a TSV. Every run starts from a cold page cache: the input is evicted
with posix_fadvise(DONTNEED) immediately before the timer starts.

Conditions (all with seed 42):
  bootstrap    full-size sample with replacement (f = 1), gzipped output. fastQpick only:
               seqtk and seqkit do not sample with replacement.
  subsample20  20% sample without replacement (m = 0.2 n), gzipped input and plain output
               for every tool, since seqtk can only write uncompressed. seqtk and seqkit
               keep each read with probability 0.2 (Bernoulli), like fastQpick single-pass.
  bootstrap10  ten full-size bootstrap replicates in one invocation (fastQpick -B 10) with
               a budget of --multi-threads threads, gzipped output. This is where the
               two-pass mode amortizes its single counting pass over the replicates.

Tools: fastQpick default (two-pass, exact), fastQpick single-pass, seqtk, seqkit. Every
tool that exposes a thread count is given --threads (default 4, the fastQpick default;
seqtk has no thread option). With a single output file fastQpick uses the extra threads
only for gzip compression; in the bootstrap10 condition it writes four replicates at a
time. The ordering is replicate-major, so slow drift shows up as
between-replicate variance rather than as a tool effect.

Peak memory is recorded two ways: the largest resident set of any single process, as
reported by /usr/bin/time -v, and the peak of the summed resident set of the whole process
tree, sampled every 0.5 s (needed for the bootstrap10 runs, which use several worker
processes; for a single-process run the two agree).

Restartable: a (rep, condition, tool) row already present in the TSV is skipped.

Example (one replicate, as used to check the values in the paper):
    python benchmarks/table1.py --input trash/bench_500M_v2.fastq.gz --replicates 1
Summarize with benchmarks/summarize_table1.py.
"""
import argparse
import datetime
import os
import platform
import shutil
import subprocess
import sys
import threading
import time

import psutil

HERE = os.path.dirname(os.path.abspath(__file__))
SEED = 42
SUBSAMPLE_FRACTION = 0.2
MULTI_B = 10

TIME_FIELDS = {
    "wall_s": "Elapsed (wall clock) time",
    "maxrss_kb": "Maximum resident set size",
    "cpu_pct": "Percent of CPU this job got",
}
# fastQpick >= 1.0 pins numpy's OpenBLAS pool itself (it never calls BLAS, and the pool
# otherwise spins one thread per core); the driver pins it too so that older versions are
# benchmarked under the same conditions.
RUN_ENV = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
COLUMNS = ["rep", "condition", "tool", "wall_s", "maxrss_kb", "tree_rss_kb", "cpu_pct",
           "out_bytes", "out_reads", "rc", "started"]


def log(logfile, msg):
    line = f"[{datetime.datetime.now():%F %T}] {msg}"
    print(line, flush=True)
    with open(logfile, "a") as fh:
        fh.write(line + "\n")


def evict(path):
    """Drop a file from the page cache so the next read is a cold read from disk."""
    subprocess.run(["sync"], check=True)
    fd = os.open(path, os.O_RDONLY)
    try:
        os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
    finally:
        os.close(fd)


def parse_time_v(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            if ": " not in line:
                continue
            key, value = line.strip().rsplit(": ", 1)
            for name, prefix in TIME_FIELDS.items():
                if key.startswith(prefix):
                    out[name] = value
    if "wall_s" in out:  # h:mm:ss or m:ss -> seconds
        secs = 0.0
        for part in out["wall_s"].split(":"):
            secs = secs * 60 + float(part)
        out["wall_s"] = f"{secs:.1f}"
    out["cpu_pct"] = out.get("cpu_pct", "").rstrip("%")
    return out


def watch_tree(root_pid, stop, peak):
    """Record the peak summed RSS of a process and all of its descendants."""
    try:
        root = psutil.Process(root_pid)
    except psutil.NoSuchProcess:
        return
    while not stop.is_set():
        total = 0
        try:
            procs = [root] + root.children(recursive=True)
        except psutil.NoSuchProcess:
            break
        for proc in procs:
            try:
                total += proc.memory_info().rss
            except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
                pass
        peak[0] = max(peak[0], total)
        stop.wait(0.5)


def count_reads(outdir, gz):
    """Number of reads in the first FASTQ written to outdir (all files are counted when
    the output is a single file; for a multi-replicate run this is one replicate)."""
    files = sorted(f for f in os.listdir(outdir) if f.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")))
    if not files:
        return ""
    path = os.path.join(outdir, files[0])
    if gz:
        cmd = f"pigz -dc -p 8 {path} | wc -l"
    else:
        cmd = f"wc -l < {path}"
    lines = int(subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True).stdout)
    return str(lines // 4)


def dir_bytes(path):
    return sum(os.path.getsize(os.path.join(root, f)) for root, _, files in os.walk(path) for f in files)


def commands(args, outdir):
    """(condition, tool) -> (shell command, gzipped output?)."""
    fq = args.fastqpick
    inp = args.input
    t1, tm = args.threads, args.multi_threads
    return {
        ("bootstrap", "fastQpick_default"):
            (f"{fq} -f 1 -s {SEED} -t {t1} -o {outdir} -w -q {inp}", True),
        ("bootstrap", "fastQpick_single_pass"):
            (f"{fq} -f 1 -p -s {SEED} -t {t1} -o {outdir} -w -q {inp}", True),
        ("subsample20", "fastQpick_default"):
            (f"{fq} -f {SUBSAMPLE_FRACTION} -dr -z -s {SEED} -t {t1} -o {outdir} -w -q {inp}", False),
        ("subsample20", "fastQpick_single_pass"):
            (f"{fq} -f {SUBSAMPLE_FRACTION} -dr -p -z -s {SEED} -t {t1} -o {outdir} -w -q {inp}", False),
        ("subsample20", "seqtk"):
            (f"{args.seqtk} sample -s {SEED} {inp} {SUBSAMPLE_FRACTION} > {outdir}/out.fastq", False),
        ("subsample20", "seqkit"):
            (f"{args.seqkit} sample -p {SUBSAMPLE_FRACTION} -s {SEED} -j {t1} -o {outdir}/out.fastq {inp}", False),
        ("bootstrap10", "fastQpick_default"):
            (f"{fq} -f 1 -B {MULTI_B} -s {SEED} -t {tm} -o {outdir} -w -q {inp}", True),
        ("bootstrap10", "fastQpick_single_pass"):
            (f"{fq} -f 1 -B {MULTI_B} -p -s {SEED} -t {tm} -o {outdir} -w -q {inp}", True),
    }


def existing_rows(tsv):
    done = set()
    if os.path.exists(tsv):
        with open(tsv) as fh:
            next(fh, None)
            for line in fh:
                rep, condition, tool = line.split("\t")[:3]
                done.add((int(rep), condition, tool))
    return done


def version_of(cmd):
    try:
        out = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=30)
        return (out.stdout + out.stderr).strip().splitlines()[0]
    except Exception as exc:  # noqa: BLE001
        return f"unavailable ({exc})"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="gzipped FASTQ library (see make_inputs.sh)")
    ap.add_argument("--replicates", type=int, default=3, help="timed runs per cell (default 3)")
    ap.add_argument("--threads", type=int, default=4, help="threads per tool for the single-sample conditions (default 4, the fastQpick default)")
    ap.add_argument("--multi-threads", type=int, default=4, help="fastQpick thread budget for the bootstrap10 condition (default 4)")
    ap.add_argument("--results", default=os.path.join(HERE, "results", "table1.tsv"), help="TSV of runs (appended; restartable)")
    ap.add_argument("--workdir", default=os.path.join(HERE, "work"), help="scratch directory for tool output (deleted after each run)")
    ap.add_argument("--conditions", nargs="+", default=["bootstrap", "subsample20", "bootstrap10"])
    ap.add_argument("--tools", nargs="+", default=["fastQpick_default", "fastQpick_single_pass", "seqtk", "seqkit"])
    ap.add_argument("--fastqpick", default="fastQpick")
    ap.add_argument("--seqtk", default=os.environ.get("SEQTK", "seqtk"))
    ap.add_argument("--seqkit", default=os.environ.get("SEQKIT", "seqkit"))
    ap.add_argument("--no-count", action="store_true", help="skip counting the reads in the output")
    args = ap.parse_args()

    args.input = os.path.abspath(args.input)
    results_dir = os.path.dirname(os.path.abspath(args.results))
    times_dir = os.path.join(results_dir, "times")
    os.makedirs(times_dir, exist_ok=True)
    os.makedirs(args.workdir, exist_ok=True)
    logfile = os.path.splitext(os.path.abspath(args.results))[0] + ".log"
    outdir = os.path.join(os.path.abspath(args.workdir), "out")

    if not os.path.exists(args.results):
        with open(args.results, "w") as fh:
            fh.write("\t".join(COLUMNS) + "\n")
    done = existing_rows(args.results)

    repo = os.path.dirname(HERE)
    git = subprocess.run(["git", "-C", repo, "describe", "--always", "--dirty"], capture_output=True, text=True).stdout.strip()
    mem_gb = psutil.virtual_memory().total / 1e9
    log(logfile, f"=== Table 1 benchmark: input={args.input} ({os.path.getsize(args.input):,} bytes) "
                 f"replicates={args.replicates} threads={args.threads} multi_threads={args.multi_threads}")
    log(logfile, f"host={platform.node()} cpus={os.cpu_count()} ram={mem_gb:.0f}GB python={platform.python_version()} "
                 f"kernel={platform.release()} git={git}")
    log(logfile, "environment: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 (BLAS is never used; fastQpick >= 1.0 pins this itself)")
    log(logfile, f"versions: {version_of(args.fastqpick + ' --version')} | seqtk {version_of(args.seqtk + ' 2>&1 | grep -i version')} | "
                 f"{version_of(args.seqkit + ' version')}")

    plan = commands(args, outdir)
    for rep in range(1, args.replicates + 1):
        log(logfile, f"=========== REP {rep} ===========")
        for condition in args.conditions:
            for tool in args.tools:
                if (condition, tool) not in plan:
                    continue  # e.g. seqtk cannot bootstrap
                if (rep, condition, tool) in done:
                    log(logfile, f"skip {condition}/{tool} rep {rep} (already recorded)")
                    continue
                cmd, gz = plan[(condition, tool)]
                tag = f"r{rep}_{condition}_{tool}"
                timefile = os.path.join(times_dir, tag + ".time")
                shutil.rmtree(outdir, ignore_errors=True)
                os.makedirs(outdir)
                evict(args.input)
                started = f"{datetime.datetime.now():%F %T}"
                log(logfile, f">>> {tag}: {cmd}")
                proc = subprocess.Popen(["/usr/bin/time", "-v", "-o", timefile, "bash", "-c", cmd], env=RUN_ENV)
                stop, peak = threading.Event(), [0]
                watcher = threading.Thread(target=watch_tree, args=(proc.pid, stop, peak), daemon=True)
                watcher.start()
                rc = proc.wait()
                stop.set()
                watcher.join()
                stats = parse_time_v(timefile)
                out_bytes = dir_bytes(outdir) if os.path.isdir(outdir) else 0
                out_reads = "" if (args.no_count or rc != 0) else count_reads(outdir, gz)
                row = dict(rep=rep, condition=condition, tool=tool, wall_s=stats.get("wall_s", ""),
                           maxrss_kb=stats.get("maxrss_kb", ""), tree_rss_kb=peak[0] // 1024,
                           cpu_pct=stats.get("cpu_pct", ""), out_bytes=out_bytes, out_reads=out_reads,
                           rc=rc, started=started)
                with open(args.results, "a") as fh:
                    fh.write("\t".join(str(row[c]) for c in COLUMNS) + "\n")
                done.add((rep, condition, tool))
                log(logfile, f"<<< {tag} rc={rc} wall={row['wall_s']}s maxRSS={row['maxrss_kb']}KB "
                             f"treeRSS={row['tree_rss_kb']}KB cpu={row['cpu_pct']}% out={out_bytes:,}B reads={out_reads or 'n/a'}")
                shutil.rmtree(outdir, ignore_errors=True)
    log(logfile, "=== ALL DONE ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
