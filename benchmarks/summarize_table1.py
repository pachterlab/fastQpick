#!/usr/bin/env python3
"""Summarize the runs recorded by table1.py into the cells of Table 1.

For every (condition, tool) the median wall time (min) and median peak memory (GB) over
the recorded replicates are printed, together with the individual values. Peak memory is
the largest single-process resident set (/usr/bin/time) for the single-process conditions
and the peak summed resident set of the process tree for bootstrap10, which runs several
worker processes.

Usage: python benchmarks/summarize_table1.py [results/table1.tsv]
"""
import csv
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = ["fastQpick_default", "fastQpick_single_pass", "seqtk", "seqkit"]
CONDITIONS = ["bootstrap", "subsample20", "bootstrap10"]
TREE_RSS_CONDITIONS = {"bootstrap10"}


def main():
    tsv = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "results", "table1.tsv")
    with open(tsv) as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t") if r["rc"] == "0"]

    cells = {}
    for r in rows:
        key = (r["condition"], r["tool"])
        mem_kb = r["tree_rss_kb"] if r["condition"] in TREE_RSS_CONDITIONS else r["maxrss_kb"]
        cells.setdefault(key, []).append((float(r["wall_s"]) / 60, float(mem_kb) / 1e6, r["out_reads"]))

    print(f"{'condition':<12} {'tool':<24} {'n':>2} {'runtime (min)':>14} {'memory (GB)':>12}   runs (min) / memory (GB) / output reads")
    for condition in CONDITIONS:
        for tool in TOOLS:
            runs = cells.get((condition, tool))
            if not runs:
                continue
            med_t = statistics.median(t for t, _, _ in runs)
            med_m = statistics.median(m for _, m, _ in runs)
            detail = "; ".join(f"{t:.1f} / {m:.3g} / {n or 'n/a'}" for t, m, n in runs)
            print(f"{condition:<12} {tool:<24} {len(runs):>2} {med_t:>14.1f} {med_m:>12.3g}   {detail}")

    # Table 1 layout of the manuscript: rows = tools, column pairs = conditions.
    print("\nTable 1 layout (median runtime min | median memory GB):")
    header = f"{'Tool':<24}" + "".join(f"{c:>26}" for c in CONDITIONS)
    print(header)
    for tool in TOOLS:
        line = f"{tool:<24}"
        for condition in CONDITIONS:
            runs = cells.get((condition, tool))
            if runs:
                line += f"{statistics.median(t for t, _, _ in runs):>13.1f}{statistics.median(m for _, m, _ in runs):>13.3g}"
            else:
                line += f"{'--':>13}{'--':>13}"
        print(line)


if __name__ == "__main__":
    main()
