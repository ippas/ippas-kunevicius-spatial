#!/usr/bin/env bash
# ============================================================
# Space Ranger count analysis for autismFMT Visium data
#
# Uses:
# - Space Ranger 4.0.1
# - refdata-gex-GRCm39-2024-A
# - metadata_autismFMT.tsv
# - automatic alignment without Loupe JSON files
#
# Output directory is provided at runtime with --output-dir.
# Missing raw sample folders are skipped and reported.
# ============================================================

set -euo pipefail

PROJECT_DIR="/home/mateusz/projects/ippas-kunevicius-spatial"

SPACERANGER="${PROJECT_DIR}/tools/spaceranger-4.0.1/spaceranger"
METADATA_FILE="${PROJECT_DIR}/data/metadata_autismFMT.tsv"

DATA_FOLDER=""
OUTPUT_DIR=""
SAMPLES=""
LOCAL_CORES=""
LOCAL_MEM=""
WITHOUT_JSON="false"

usage() {
    cat <<EOF

Usage:

bash preprocessing/spaceranger_analysis/spaceranger-count-analysis.sh \\
  --data-folder raw \\
  --output-dir data/spaceranger_count_raw_without_json \\
  --samples all \\
  --without-json \\
  --localcores 10 \\
  --localmem 50

EOF
}

error_exit() {
    echo "ERROR: $1" >&2
    exit 1
}

message() {
    echo "[$(date '+%F %T')] $1"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --data-folder)
            DATA_FOLDER="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --samples)
            SAMPLES="$2"
            shift 2
            ;;
        --without-json)
            WITHOUT_JSON="true"
            shift
            ;;
        --localcores)
            LOCAL_CORES="$2"
            shift 2
            ;;
        --localmem)
            LOCAL_MEM="$2"
            shift 2
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            error_exit "Unknown argument: $1"
            ;;
    esac
done

[[ -n "$DATA_FOLDER" ]] || error_exit "Missing --data-folder"
[[ -n "$OUTPUT_DIR" ]] || error_exit "Missing --output-dir"
[[ -n "$SAMPLES" ]] || error_exit "Missing --samples"
[[ -n "$LOCAL_CORES" ]] || error_exit "Missing --localcores"
[[ -n "$LOCAL_MEM" ]] || error_exit "Missing --localmem"
[[ "$WITHOUT_JSON" == "true" ]] || error_exit "Add --without-json"

[[ -x "$SPACERANGER" ]] || error_exit "Space Ranger not found: $SPACERANGER"
[[ -f "$METADATA_FILE" ]] || error_exit "Metadata file not found: $METADATA_FILE"

if [[ "$DATA_FOLDER" = /* ]]; then
    DATA_DIR="$DATA_FOLDER"
else
    DATA_DIR="${PROJECT_DIR}/${DATA_FOLDER}"
fi

if [[ "$OUTPUT_DIR" = /* ]]; then
    OUTPUT_DIR="$OUTPUT_DIR"
else
    OUTPUT_DIR="${PROJECT_DIR}/${OUTPUT_DIR}"
fi

TRANSCRIPTOME="${DATA_DIR}/reference/refdata-gex-GRCm39-2024-A"

[[ -d "$DATA_DIR" ]] || error_exit "Data directory not found: $DATA_DIR"
[[ -d "$TRANSCRIPTOME" ]] || error_exit "Reference not found: $TRANSCRIPTOME"

mkdir -p "$OUTPUT_DIR"
LOG_DIR="${OUTPUT_DIR}/logs"
mkdir -p "$LOG_DIR"

message "Space Ranger: $SPACERANGER"
message "Data directory: $DATA_DIR"
message "Reference: $TRANSCRIPTOME"
message "Output directory: $OUTPUT_DIR"
message "Samples: $SAMPLES"
message "Loupe JSON alignment: disabled"
message "Resources per sample: ${LOCAL_CORES} cores, ${LOCAL_MEM} GB RAM"
echo

if [[ "$SAMPLES" == "all" ]]; then
    mapfile -t SAMPLE_ARRAY < <(
        tail -n +2 "$METADATA_FILE" | cut -f1 | tr -d '\r'
    )
else
    IFS=',' read -r -a SAMPLE_ARRAY <<< "$SAMPLES"
fi

MISSING_SAMPLES=()
COMPLETED_SAMPLES=()
SKIPPED_SAMPLES=()

for SAMPLE_ID in "${SAMPLE_ARRAY[@]}"; do

    SAMPLE_ID="$(echo "$SAMPLE_ID" | xargs)"
    [[ -n "$SAMPLE_ID" ]] || continue

    METADATA_ROW="$(
        awk -F '\t' -v sample="$SAMPLE_ID" \
        'NR > 1 && $1 == sample {print; exit}' \
        "$METADATA_FILE" | tr -d '\r'
    )"

    [[ -n "$METADATA_ROW" ]] || error_exit "Sample not found in metadata: $SAMPLE_ID"

    IFS=$'\t' read -r \
        METADATA_SAMPLE_ID \
        SLIDE_ID \
        SLIDE_AREA \
        _ <<< "$METADATA_ROW"

    SAMPLE_DIR="${DATA_DIR}/${SAMPLE_ID}"

    if [[ ! -d "$SAMPLE_DIR" ]]; then
        message "SKIP: raw data directory missing for ${SAMPLE_ID}"
        MISSING_SAMPLES+=("$SAMPLE_ID")
        continue
    fi

    mapfile -t TIFF_FILES < <(
        find "$SAMPLE_DIR" \
            -maxdepth 1 \
            -type f \
            \( -iname "*.tif" -o -iname "*.tiff" \) \
            -print | sort
    )

    [[ ${#TIFF_FILES[@]} -eq 1 ]] || \
        error_exit "Expected one TIFF image for ${SAMPLE_ID}; found ${#TIFF_FILES[@]}"

    IMAGE_FILE="${TIFF_FILES[0]}"

    R1_COUNT="$(find "$SAMPLE_DIR" -maxdepth 1 -type f -name "*_R1_*.fastq.gz" | wc -l)"
    R2_COUNT="$(find "$SAMPLE_DIR" -maxdepth 1 -type f -name "*_R2_*.fastq.gz" | wc -l)"

    [[ "$R1_COUNT" -gt 0 ]] || error_exit "No R1 FASTQ found for ${SAMPLE_ID}"
    [[ "$R2_COUNT" -gt 0 ]] || error_exit "No R2 FASTQ found for ${SAMPLE_ID}"

    RUN_ID="$SAMPLE_ID"
    RUN_OUTPUT="${OUTPUT_DIR}/${RUN_ID}"

    if [[ -d "${RUN_OUTPUT}/outs" ]]; then
        message "SKIP: completed output already exists for ${SAMPLE_ID}"
        SKIPPED_SAMPLES+=("$SAMPLE_ID")
        continue
    fi

    if [[ -e "$RUN_OUTPUT" ]]; then
        error_exit "Incomplete output exists: ${RUN_OUTPUT}. Remove it before rerunning."
    fi

    COMMAND=(
        "$SPACERANGER" count
        "--id=${RUN_ID}"
        "--description=autismFMT_${SAMPLE_ID}_without_json"
        "--transcriptome=${TRANSCRIPTOME}"
        "--fastqs=${SAMPLE_DIR}"
        "--sample=${SAMPLE_ID}"
        "--image=${IMAGE_FILE}"
        "--slide=${SLIDE_ID}"
        "--area=${SLIDE_AREA}"
        "--create-bam=false"
        "--localcores=${LOCAL_CORES}"
        "--localmem=${LOCAL_MEM}"
    )

    echo
    message "START: ${SAMPLE_ID}"
    message "Slide: ${SLIDE_ID}; area: ${SLIDE_AREA}"
    message "Image: ${IMAGE_FILE}"
    message "Output: ${RUN_OUTPUT}"
    echo

    LOG_FILE="${LOG_DIR}/${SAMPLE_ID}_$(date '+%Y%m%d_%H%M%S').log"

    set +e
    (
        cd "$OUTPUT_DIR"
        "${COMMAND[@]}"
    ) 2>&1 | tee "$LOG_FILE"

    SPACERANGER_STATUS="${PIPESTATUS[0]}"
    set -e

    [[ "$SPACERANGER_STATUS" -eq 0 ]] || \
        error_exit "Space Ranger failed for ${SAMPLE_ID}. Log: ${LOG_FILE}"

    [[ -d "${RUN_OUTPUT}/outs" ]] || \
        error_exit "Run completed but outs/ folder is missing for ${SAMPLE_ID}"

    COMPLETED_SAMPLES+=("$SAMPLE_ID")

    message "DONE: ${SAMPLE_ID}"

done

echo
message "Finished."
message "Completed: ${#COMPLETED_SAMPLES[@]}"
message "Already completed and skipped: ${#SKIPPED_SAMPLES[@]}"
message "Missing raw directories: ${#MISSING_SAMPLES[@]}"

if [[ ${#MISSING_SAMPLES[@]} -gt 0 ]]; then
    echo "Missing samples: ${MISSING_SAMPLES[*]}"
fi
