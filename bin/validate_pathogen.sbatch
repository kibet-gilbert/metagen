#!/usr/bin/env bash
# =============================================================================
# validate_pathogen.sbatch
#
# SLURM submission script for the standalone pathogen validation pipeline.
# Mirrors the structure of ccmetagen.batch.sbatch in the metagen repo.
#
# Submits the Nextflow orchestrator as a single long-running job;
# Nextflow then submits individual SLURM jobs for each process.
#
# Usage (from the metagenomics project directory):
#   sbatch scripts/validate_pathogen.sbatch \
#     --samplesheet data/samplesheet_validation.csv \
#     --outdir      /var/scratch/gkibet/results/validation/run01 \
#     --taxids      "666,562,573,470,287,485" \
#     --blast_db    /export/data/bio/ncbi/blast/db/v5 \
#     --db_name     nt
#
# Or set all options inside the script and run:
#   sbatch scripts/validate_pathogen.sbatch
# =============================================================================

#SBATCH --job-name=validaPathogen_nf
#SBATCH --output=logs/validate_pathogen_%j.log
#SBATCH --error=logs/validate_pathogen_%j.err
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=800G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=gkibet@ilri.org
# #SBATCH --time=72:00:00

set -Eeuo pipefail

# =============================================================================
# CONFIGURATION — edit these or override via command-line arguments
# =============================================================================

# ── Input / output ────────────────────────────────────────────────────────────
#Setting work Directory in HPC
INDIR=$PWD
FASTQDIR=${INDIR}/fastq
src_dir=${INDIR##*/}
DATADIR=/export/data/ilri/sarscov2/databases/
SCRIPTDIR=~/bioinformatics/github/metagenomics/scripts/
WORKDIR=/var/scratch/$USER/${src_dir}
src_dir1=${WORKDIR##*/}

# Nextflow
WKDIR=${WORKDIR}/work
NXF_SINGULARITY_CACHEDIR=${WKDIR}/singularity
export NXF_SINGULARITY_CACHEDIR=${WKDIR}/singularity
NXF_HOME=${WKDIR}
mkdir -p $WORKDIR
# cd $WORKDIR

# ── Project and pipeline paths ────────────────────────────────────────────────
PIPELINE_DIR="${SCRIPTDIR}/metagen"
NXF_SCRIPT="${PIPELINE_DIR}/subworkflows/local/validate_pathogen_detection/main.nf"
NXF_CONFIG="${PIPELINE_DIR}/conf/nextflow.validate_pathogen.config"
NXF_PROFILE="ilri_hpc"
SAMPLESHEET="${INDIR}/samplesheet_validation.csv"
OUTDIR="${WORKDIR}/results/pathogen_validation"

# ── Target pathogen taxids (NCBI) ─────────────────────────────────────────────
# Vibrio cholerae  : 666     Salmonella Typhi       : 90370
# Salmonella Typhim: 90371   Escherichia coli       : 562
# Klebsiella pneu  : 573     Acinetobacter baum     : 470
# Pseudomonas aerug: 287     Neisseria gonorrhoeae  : 485
# Haemophilus inf (S): 727     Staphylococcus aureus  (S) : 1280
# Mononegavirales (F): 11157   Bunyavirales           () : 1980410
# Flaviviridae    (F): 11050
#TARGET_TAXIDS="11157,11050,1980410,666,90370,90371,562,573,470,287,485,727,1280"
TARGET_TAXIDS="11157,11050,1980410,573,485"

# ── BLAST database ────────────────────────────────────────────────────────────
BLAST_DB_DIR="/export/data/bio/ncbi/blast/db/v5"
BLAST_DB_NAME="nt"
BLAST_DB_MODE="local"   # local | custom | remote

# ── Extraction options ────────────────────────────────────────────────────────
INCLUDE_CHILDREN=true
MAX_READS_PER_TAXID=500

# ── BLAST options ─────────────────────────────────────────────────────────────
BLAST_TASK="megablast"
BLAST_EVALUE="1e-20"
BLAST_PERC_ID=95
BLAST_QCOV=80
BLAST_MAX_SEQS=5

# ── Validation thresholds ─────────────────────────────────────────────────────
MIN_BLAST_HITS=3
MIN_PIDENT=95.0
MIN_QCOVS=80.0
MIN_BITSCORE=100.0

# ── Skip flags ────────────────────────────────────────────────────────────────
SKIP_BLAST=false          # true = extract reads only, skip BLAST

# =============================================================================
# OVERRIDE WITH COMMAND-LINE ARGUMENTS
# Supports:  --samplesheet  --outdir  --taxids  --blast_db  --db_name
# =============================================================================
while [[ $# -gt 0 ]]; do
  case "$1" in
    --samplesheet)  SAMPLESHEET="$2";  shift 2 ;;
    --outdir)       OUTDIR="$2";       shift 2 ;;
    --taxids)       TARGET_TAXIDS="$2"; shift 2 ;;
    --blast_db)     BLAST_DB_DIR="$2"; shift 2 ;;
    --db_name)      BLAST_DB_NAME="$2"; shift 2 ;;
    --db_mode)      BLAST_DB_MODE="$2"; shift 2 ;;
    --skip_blast)   SKIP_BLAST=true;   shift   ;;
    --resume)       RESUME="-resume";  shift   ;;
    *) echo "[WARN] Unknown argument: $1"; shift ;;
  esac
done
RESUME="${RESUME:-}"

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================
echo "======================================================================"
echo "Pathogen Detection Validation Pipeline"
echo "ILRI Genomics Platform"
echo "Started     : $(date)"
echo "Job ID      : ${SLURM_JOB_ID}"
echo "Node        : $(hostname)"
echo "Samplesheet : ${SAMPLESHEET}"
echo "PWD         : ${PWD}"
echo "Outdir      : ${OUTDIR}"
echo "Target taxa : ${TARGET_TAXIDS}"
echo "BLAST DB    : ${BLAST_DB_DIR}/${BLAST_DB_NAME} (${BLAST_DB_MODE})"
echo "Skip BLAST  : ${SKIP_BLAST}"
echo "Resume      : ${RESUME}"
echo "======================================================================"

# Load required modules
module purge
module load nextflow/26.04  2>/dev/null || \
  module load nextflow         2>/dev/null || \
  { echo "[ERROR] Nextflow not available as a module — checking PATH"; }

nextflow -version || { echo "[ERROR] Nextflow not found"; exit 1; }
export NXF_OPTS='-Xms1g -Xmx4g'
export NXF_ANSI_LOG=false

# ── Working directory setup ───────────────────────────────────────────────────
WORK_DIR="${WORKDIR}/validate_pathogen"
mkdir -p "${OUTDIR}" "${WORK_DIR}" logs
echo "[INFO] Nextflow work dir: ${WORK_DIR}"

# ── Apptainer cache ───────────────────────────────────────────────────────────
export APPTAINER_CACHEDIR="${WORK_DIR}/.apptainer_cache"
export SINGULARITY_CACHEDIR="${APPTAINER_CACHEDIR}"
mkdir -p "${APPTAINER_CACHEDIR}"

# =============================================================================
# RUN NEXTFLOW
# =============================================================================

nextflow run "${NXF_SCRIPT}" \
  -c              "${NXF_CONFIG}" \
  -profile        "${NXF_PROFILE}" \
  -work-dir       "${WORK_DIR}" \
  -ansi-log       false \
  ${RESUME} \
  \
  --samplesheet       "${SAMPLESHEET}" \
  --outdir            "${OUTDIR}" \
  --target_taxids     "${TARGET_TAXIDS}" \
  \
  --blast_db_dir      "${BLAST_DB_DIR}" \
  --blast_db_name     "${BLAST_DB_NAME}" \
  --blast_db_mode     "${BLAST_DB_MODE}" \
  \
  --include_children  "${INCLUDE_CHILDREN}" \
  --max_reads_per_taxid "${MAX_READS_PER_TAXID}" \
  \
  --blast_task        "${BLAST_TASK}" \
  --blast_evalue      "${BLAST_EVALUE}" \
  --blast_perc_id     "${BLAST_PERC_ID}" \
  --blast_qcov        "${BLAST_QCOV}" \
  --blast_max_seqs    "${BLAST_MAX_SEQS}" \
  \
  --min_blast_hits    "${MIN_BLAST_HITS}" \
  --min_pident        "${MIN_PIDENT}" \
  --min_qcovs         "${MIN_QCOVS}" \
  --min_bitscore      "${MIN_BITSCORE}" \
  \
  --skip_blast        "${SKIP_BLAST}" \
  --publish_dir_mode  "copy" \
  \
  2>&1 | tee "${INDIR}/nextflow_run_${SLURM_JOB_ID}.log"

EXIT_CODE=${PIPESTATUS[0]}

# =============================================================================
# POST-RUN SUMMARY
# =============================================================================
echo ""
echo "======================================================================"
echo "Pipeline exit code: ${EXIT_CODE}"
echo "Finished    : $(date)"
echo ""

if [[ ${EXIT_CODE} -eq 0 ]]; then
  echo "[SUCCESS] Validation complete."
  echo ""
  echo "Key output files:"
  echo "  ${OUTDIR}/validation_summary/all_samples_validation.tsv"
  echo "  ${OUTDIR}/validation/          (per-sample TSV + JSON)"
  echo "  ${OUTDIR}/blast/               (raw BLAST hits)"
  echo "  ${OUTDIR}/extracted_reads/     (KrakenTools FASTA outputs)"
  echo "  ${OUTDIR}/pipeline_info/       (timeline, report, trace)"
  echo ""

  # Print a quick summary of validation outcomes
  SUMMARY="${OUTDIR}/validation_summary/all_samples_validation.tsv"
  if [[ -f "${SUMMARY}" ]]; then
    echo "Validation outcome summary:"
    awk -F'\t' 'NR>1{counts[$NF]++} END{for(s in counts) printf "  %-20s %d\n", s, counts[s]}' \
      "${SUMMARY}"
  fi

  # Optional: clean Nextflow work directory after successful run
  # Uncomment to free scratch space (disable if you need -resume capability)
  # echo ""
  # echo "[INFO] Cleaning work directory: ${WORK_DIR}"
  # rm -rf "${WORK_DIR}"

else
  echo "[FAILED] Pipeline exited with code ${EXIT_CODE}."
  echo ""
  echo "Troubleshooting:"
  echo "  1. Check the log: ${OUTDIR}/pipeline_info/nextflow_run_${SLURM_JOB_ID}.log"
  echo "  2. Check per-process errors in: ${WORK_DIR}/"
  echo "  3. Re-run with -resume after fixing the issue:"
  echo ""
  echo "     sbatch scripts/validate_pathogen.sbatch --resume \\"
  echo "       --samplesheet ${SAMPLESHEET} \\"
  echo "       --outdir ${OUTDIR}"
  echo ""
fi

echo "======================================================================"
exit ${EXIT_CODE}
