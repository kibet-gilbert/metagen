#!/usr/bin/env bash
# =============================================================================
#  run_bracken.sh
#  SLURM wrapper to launch the BRACKEN_REESTIMATE standalone subworkflow.
#
#  Submit with:   sbatch run_bracken.sh
#  Or dry-run:    bash run_bracken.sh
# =============================================================================
#SBATCH -p batch
#SBATCH -J bracken_reestimate
#SBATCH --qos vip
#SBATCH -n 12
#SBATCH --mem 240G
#SBATCH --time 48:00:00
#SBATCH -o logs/bracken_%j.out
#SBATCH -e logs/bracken_%j.err

set -uEo pipefail
trap 'echo "Error $? on line $LINENO" && exit 1' ERR

mkdir -p logs

# ---------------------------------------------------------------------------
# Load Nextflow
# ---------------------------------------------------------------------------
module load nextflow/25.04

export NXF_OPTS='-Xms1g -Xmx4g'
export NXF_ANSI_LOG=false

# ---------------------------------------------------------------------------
# Paths — edit these before running
# ---------------------------------------------------------------------------
# INDIR="$PWD"                                # Analysis root (where data lives)
# INDIR="/home/gkibet/bioinformatics/github/metagenomics/data/visualization_WWP2/"
INDIR="/home/gkibet/bioinformatics/github/metagenomics/data/visualization/"
PIPELINE_DIR="${INDIR}../../scripts/metagen/" # project root (where main.nf is)
KREPORTS="${INDIR}/kreports/wastewater/kraken2/*.kraken2_report.txt"
DB="/export/data/ilri/sarscov2/databases/metagenomicsDBs/20231129_k2_nt_20231129"  # <-- set your DB path

OUTDIR="${INDIR}/kreports/wastewater"
WORKDIR="/var/scratch/${USER}/bracken_work"

mkdir -p "$WORKDIR" "$OUTDIR"

# ---------------------------------------------------------------------------
# Bracken settings — adjust as needed
# ---------------------------------------------------------------------------
READLEN=150         # must match the --read_len used during kraken2-build
LEVEL="S"           # S=species  G=genus  F=family  O=order  C=class  P=phylum
THRESHOLD=10        # minimum read count to re-estimate a taxon
BUILD_DB=false      # set to true if your DB lacks a .kmer_distrib file

# ---------------------------------------------------------------------------
# Resume token (comment out on first run; set after a failed run to resume)
# ---------------------------------------------------------------------------
# RESUME="-resume <session_id>"
RESUME="-resume"

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
echo "=========================================================="
echo "  BRACKEN REESTIMATE  — $(date)"
echo "  kreports : $KREPORTS"
echo "  database : $DB"
echo "  outdir   : $OUTDIR"
echo "  workdir  : $WORKDIR"
echo "  node     : ${SLURMD_NODENAME:-local}"
echo "=========================================================="

nextflow run \
    "${PIPELINE_DIR}/subworkflows/local/bracken_reestimate/main.nf" \
    -entry BRACKEN_REESTIMATE_STANDALONE \
    -c "${PIPELINE_DIR}/conf/nextflow.bracken.config" \
    -w "${WORKDIR}" \
    --kreports  "${KREPORTS}" \
    --db        "${DB}" \
    --outdir    "${OUTDIR}" \
    --bracken_readlen   ${READLEN} \
    --bracken_level     "${LEVEL}" \
    --bracken_threshold ${THRESHOLD} \
    --bracken_build     ${BUILD_DB} \
    -profile singularity \
    ${RESUME} \
    |& tee "${INDIR}/logs/bracken_run.log"

COMBINED_REPORT="${INDIR}/kreports/bracken/combined/all_samples_${LEVEL}.tsv"
mv "${INDIR}/kreports/bracken/combined/all_samples.txt" ${COMBINED_REPORT}

echo "Done — $(date)"
