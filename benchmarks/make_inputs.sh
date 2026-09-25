#!/bin/bash
# Build the synthetic 500-million-read library used in Table 1 (about 158 GB uncompressed
# and 40 GB gzipped; allow a few hours). Records are fixed-width (150 bp, 316 bytes), with
# uniform random bases and a position-dependent binned quality profile so that the file
# compresses about 4x, as real data does. gen_fastq.py is deterministic for a given seed.
# Usage: make_inputs.sh <output directory>
set -euo pipefail
DIR=${1:?output directory}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$DIR"
PLAIN=$DIR/bench_500M.fastq
GZ=$PLAIN.gz
[ -s "$PLAIN" ] || python "$HERE/gen_fastq.py" --reads 500000000 --seed 1 --out "$PLAIN"
[ -s "$GZ" ] || pigz -6 -p 32 -c "$PLAIN" > "$GZ"
ls -l "$PLAIN" "$GZ"
