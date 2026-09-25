#!/bin/bash
# Reproduce Figure 1 (read-level bootstrap of the yeast RNA-seq library SRR453566).
#
#   1. Downloads kallisto, the Ensembl R64-1-1 yeast cDNA, and the paired reads if absent,
#      builds the index, and quantifies the original library (counts/baseline.txt).
#   2. Generates B full-size fastQpick bootstrap replicates (seeds 1..B), re-quantifies
#      each with kallisto (run_one.py), P replicates at a time.
#   3. Renders notebooks/figures/bootstrap_realdata.png and prints the numbers quoted in
#      the Application section (analyze.py).
#
# Usage: drive.sh [B] [P]     (default B=200 replicates, P=8 in parallel)
# Restartable: replicates whose counts file already exists are skipped.
set -euo pipefail
cd "$(dirname "$0")"
B=${1:-200}
P=${2:-8}
ACCESSION=SRR453566

log() { echo "[$(date '+%F %T')] $*" | tee -a progress.log; }

# ---- inputs -------------------------------------------------------------------
if [ ! -x kallisto/kallisto ]; then
  log "downloading kallisto 0.51.1"
  wget -q https://github.com/pachterlab/kallisto/releases/download/v0.51.1/kallisto_linux-v0.51.1.tar.gz -O kallisto.tar.gz
  tar xzf kallisto.tar.gz
fi
if [ ! -s cdna.fa.gz ]; then
  log "downloading yeast cDNA (Ensembl release 110, R64-1-1)"
  wget -q http://ftp.ensembl.org/pub/release-110/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz -O cdna.fa.gz
fi
[ -s yeast.idx ] || { log "building kallisto index"; kallisto/kallisto index -i yeast.idx cdna.fa.gz; }
mkdir -p data
for r in 1 2; do
  f=data/${ACCESSION}_${r}.fastq.gz
  [ -s "$f" ] || { log "downloading $f"; wget -q "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR453/${ACCESSION}/${ACCESSION}_${r}.fastq.gz" -O "$f"; }
done

# ---- baseline quantification ------------------------------------------------------
mkdir -p counts
if [ ! -s counts/baseline.txt ]; then
  log "quantifying the original library"
  kallisto/kallisto quant -i yeast.idx -o kout_baseline -t 4 data/${ACCESSION}_1.fastq.gz data/${ACCESSION}_2.fastq.gz
  awk 'NR>1{print $1}' kout_baseline/abundance.tsv > targets.txt
  awk 'NR>1{printf "%.6f\n", $4}' kout_baseline/abundance.tsv > counts/baseline.txt
fi

# ---- bootstrap replicates ----------------------------------------------------------
log "fastQpick $(fastQpick --version | awk '{print $2}'), git $(git -C ../.. rev-parse --short HEAD)$(git -C ../.. diff --quiet || echo '-dirty')"
log "generating $B bootstrap replicates, $P at a time"
seq 1 "$B" | xargs -P "$P" -I{} bash -c '
  s=$(printf "%04d" {})
  out=counts/rep_${s}.txt
  [ -s "$out" ] && exit 0
  python run_one.py {} tmp/w{} "$out" 2>/dev/null && echo "[$(date "+%F %T")] done {}" >> progress.log
'
n=$(ls counts/rep_*.txt | wc -l)
[ "$n" -eq "$B" ] || { log "only $n of $B replicates finished"; exit 1; }
rm -rf tmp
log "ALL_DONE ($n replicates)"

# ---- figure ------------------------------------------------------------------------
python analyze.py | tee -a progress.log
