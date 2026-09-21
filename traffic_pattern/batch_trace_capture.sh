#!/usr/bin/env bash
# Batch DUMPI Trace Capture for LAMMPS Communication Analysis
# Runs all combinations of cases and rank counts
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
INSTALL_DIR="${PROJECT_DIR}/install"
LMP_BIN="${INSTALL_DIR}/bin/lmp"
LIBDUMPI="${INSTALL_DIR}/lib/libdumpi.so"
DUMPI2CCDG="${PROJECT_DIR}/dumpi2ccdg/dumpi2ccdg"

# Define experiment matrix
declare -a CASES=(
    "lj_0.1k"
    "lj_1k"
    "lj_10k"
    "lj_mix_1k"
)
declare -a RANKS=(4 8 16 32)

# Summary file (local to traffic_pattern)
SUMMARY="${SCRIPT_DIR}/EXPERIMENT_SUMMARY.md"
echo "# Experiment Summary" > "${SUMMARY}"
echo "" >> "${SUMMARY}"
echo "| # | Case | Atoms | Ranks | Run Dir | Status | Nodes | Comm Nodes | Cross Edges | Total Time (s) |" >> "${SUMMARY}"
echo "|---|------|-------|-------|---------|--------|-------|------------|-------------|----------------|" >> "${SUMMARY}"

EXP_NUM=0
TOTAL_EXP=$((${#CASES[@]} * ${#RANKS[@]}))

for CASE in "${CASES[@]}"; do
    CASE_DIR="${SCRIPT_DIR}/cases/${CASE}"
    if [ ! -f "${CASE_DIR}/in.lammps" ]; then
        echo "WARNING: Case dir ${CASE_DIR} has no in.lammps, skipping"
        continue
    fi

    for N_RANKS in "${RANKS[@]}"; do
        EXP_NUM=$((EXP_NUM + 1))
        TIMESTAMP=$(date +%Y%m%d_%H%M%S)
        RUN_DIR="${PROJECT_DIR}/runs/${CASE}_${N_RANKS}r_${TIMESTAMP}"
        CCDG_FILE="${RUN_DIR}/${CASE}_global.ccdg"

        echo ""
        echo "=========================================="
        echo "[${EXP_NUM}/${TOTAL_EXP}] ${CASE} x ${N_RANKS} ranks"
        echo "  Output: ${RUN_DIR}"
        echo "=========================================="

        # --- Step 1: Trace Capture ---
        rm -rf "${RUN_DIR}"
        mkdir -p "${RUN_DIR}"
        cp -r "${CASE_DIR}"/* "${RUN_DIR}/"
        cd "${RUN_DIR}"

        export LD_PRELOAD="${LIBDUMPI}"
        export LD_LIBRARY_PATH="${INSTALL_DIR}/lib:${LD_LIBRARY_PATH:-}"
        export DUMPI_OUTDIR="${RUN_DIR}"

        echo "  [Step 1/3] Running LAMMPS (${N_RANKS} ranks)..."
        if ! mpirun -np "${N_RANKS}" --allow-run-as-root \
            "${LMP_BIN}" -in in.lammps > lammps.log 2>&1; then
            echo "  FAILED: LAMMPS run failed"
            echo "| ${EXP_NUM} | ${CASE} | - | ${N_RANKS} | ${RUN_DIR} | LAMMPS_FAIL | - | - | - | - |" >> "${SUMMARY}"
            continue
        fi

        # Check trace files
        META_COUNT=$(ls dumpi-*.meta 2>/dev/null | wc -l)
        BIN_COUNT=$(ls dumpi-*.bin 2>/dev/null | wc -l)
        echo "  Trace files: ${META_COUNT} meta, ${BIN_COUNT} bin"

        if [ "${META_COUNT}" -eq 0 ]; then
            echo "  FAILED: No trace files generated"
            echo "| ${EXP_NUM} | ${CASE} | - | ${N_RANKS} | ${RUN_DIR} | NO_TRACE | - | - | - | - |" >> "${SUMMARY}"
            continue
        fi

        cd "${SCRIPT_DIR}"

        # --- Step 2: CCDG Conversion ---
        echo "  [Step 2/3] Converting to CCDG..."
        if [ -x "${DUMPI2CCDG}" ]; then
            cd "${PROJECT_DIR}/dumpi2ccdg"
            if ! ./dumpi2ccdg "${RUN_DIR}" > "${CCDG_FILE}" 2>/dev/null; then
                echo "  WARNING: CCDG conversion may have failed"
            fi
            cd "${SCRIPT_DIR}"
        else
            echo "  SKIP: dumpi2ccdg not found at ${DUMPI2CCDG}"
        fi

        # --- Step 3: Extract stats from CCDG ---
        echo "  [Step 3/3] Extracting statistics..."
        STATUS="OK"
        NODES="?"
        COMM_NODES="?"
        CROSS_EDGES="?"
        TOTAL_TIME="?"

        if [ -f "${CCDG_FILE}" ]; then
            # Extract stats from CCDG header (lines before JSON)
            NODES=$(grep -oP 'Total nodes:\s+\K\d+' "${CCDG_FILE}" | head -1)
            COMM_NODES=$(grep -oP 'Communication nodes:\s+\K\d+' "${CCDG_FILE}" | head -1)
            TOTAL_TIME=$(grep -oP 'Total compute time:\s+\K[\d.]+' "${CCDG_FILE}" | head -1)
            CROSS_EDGES=$(grep -oP 'Cross-rank edges:\s+\K\d+' "${CCDG_FILE}" | head -1)
            STATUS="OK"
        fi

        echo "  Result: ${STATUS} | Nodes=${NODES} | Comm=${COMM_NODES} | CrossEdges=${CROSS_EDGES} | Time=${TOTAL_TIME}s"
        echo "| ${EXP_NUM} | ${CASE} | - | ${N_RANKS} | ${RUN_DIR##*/} | ${STATUS} | ${NODES} | ${COMM_NODES} | ${CROSS_EDGES} | ${TOTAL_TIME} |" >> "${SUMMARY}"
    done
done

echo ""
echo "=========================================="
echo "Batch capture complete! ${EXP_NUM} experiments."
echo "Summary: ${SUMMARY}"
echo "=========================================="
cat "${SUMMARY}"
