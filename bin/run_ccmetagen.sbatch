#!/usr/bin/env bash
# =============================================================================
#  run_ccmetagen.sh
#  SLURM wrapper to launch the CCMETAGEN_TAXONOMY standalone subworkflow.
#
#  Submit with:   sbatch run_ccmetagen.sh
#  Dry-run:       bash run_ccmetagen.sh
# =============================================================================
#SBATCH -p batch
#SBATCH -J ccmetagen
#SBATCH --qos vip
#SBATCH -n 16
#SBATCH --mem 240G
#SBATCH --time 72:00:00
#SBATCH -o logs/ccmetagen_%j.out
#SBATCH -e logs/ccmetagen_%j.err

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
# Project root and paths
# ---------------------------------------------------------------------------
INDIR="/home/gkibet/bioinformatics/github/metagenomics/data/20251201_sludge24_NextSeq2k/"
PIPELINE_DIR="${INDIR}../../scripts/metagen/"
OUTDIR="${INDIR}/results_arg/ccmetagen"
WORKDIR="/var/scratch/${USER}/ccmetagen_work"
mkdir -p "${WORKDIR}" "${OUTDIR}" logs

# ---------------------------------------------------------------------------
# DB MODE — uncomment the block that matches your situation
# ---------------------------------------------------------------------------

# ── MODE 1: Use a pre-built KMA DB directory (FASTEST — recommended) ────────
#DB_MODE="prebuilt"
#DB_PATH="/export/data/ilri/sarscov2/databases/ccmetagen/RefSeq_kma"
# Set DB_PATH to your KMA index directory containing *.name *.seq *.len files.

# ── MODE 2: Export from a local BLAST DB then clean + KMA-index ─────────────
# DB_MODE="blast_export"
# BLAST_PREFIX="/export/data/ncbi/blast/db/v5/core_nt"   # path WITHOUT extension
# BLAST_NAME="core_nt"

# ── MODE 3: Build KMA DB from FASTA (NCBI_nt download) ──────────────────────
DB_MODE="build"
DB_KEYWORD="core_nt"     # core_nt | nt | refseq
# TAXON="bacteria"         # only used for targetted taxon (refseq)

# ── MODE 4: Download a pre-built archive (Unimelb mediaflux) ─────────────────
# DB_MODE="download"
# DB_KEYWORD="RefSeq"      # RefSeq | NCBI_nt | NCBI_nt_no_env

# ---------------------------------------------------------------------------
# Reads — edit the glob to match your FASTQ location and naming convention
# Paired-end example : "fastq/*_{R1,R2}.fastq.gz"
# Single-end example : "fastq/*.fastq.gz"  (also set SINGLE_END=true)
# ---------------------------------------------------------------------------
READS="fastq/*_{R1,R2}.fastq.gz"
SINGLE_END=false

# ---------------------------------------------------------------------------
# Resume (uncomment after a failed run; Nextflow picks up the last session)
# ---------------------------------------------------------------------------
# RESUME="-resume"
RESUME="-resume"

# ---------------------------------------------------------------------------
# Build nextflow arguments based on DB_MODE
# ---------------------------------------------------------------------------
DB_ARGS=""

case "$DB_MODE" in
    prebuilt)
        DB_ARGS="--ccmetagen_db_mode prebuilt --db '${DB_PATH}'"
        ;;
    blast_export)
        DB_ARGS="--ccmetagen_db_mode blast_export \
                 --blast_db_prefix '${BLAST_PREFIX:-}' \
                 --blast_db_name '${BLAST_NAME:-core_nt}'"
        ;;
    build)
        DB_ARGS="--ccmetagen_db_mode build \
                 --ccmetagen_keyword '${DB_KEYWORD:-core_nt}' \
                 --ccmetagen_taxon '${TAXON:-bacteria}'"
        ;;
    download)
        DB_ARGS="--ccmetagen_db_mode download \
                 --ccmetagen_keyword '${DB_KEYWORD:-NCBI_nt}'"
        ;;
    *)
        echo "ERROR: unknown DB_MODE '${DB_MODE}'"
        exit 1
        ;;
esac

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
echo "============================================================"
echo "  CCMETAGEN TAXONOMY — $(date)"
echo "  reads    : $READS"
echo "  db mode  : $DB_MODE"
echo "  outdir   : $OUTDIR"
echo "  workdir  : $WORKDIR"
echo "  node     : ${SLURMD_NODENAME:-local}"
echo "============================================================"

nextflow run \
    "${PIPELINE_DIR}/subworkflows/local/ccmetagen_taxonomy/main.nf" \
    -entry  CCMETAGEN_TAXONOMY_STANDALONE \
    -c      "${PIPELINE_DIR}/conf/nextflow.ccmetagen.config" \
    -w      "${WORKDIR}" \
    --reads "${READS}" \
    --outdir "${OUTDIR}" \
    --single_end ${SINGLE_END} \
    ${DB_ARGS} \
    -profile singularity \
    ${RESUME} \
    |& tee "${INDIR}/logs/ccmetagen_run.log"

echo "Done — $(date)"
