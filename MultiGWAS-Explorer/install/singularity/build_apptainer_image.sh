#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd -P)"
DEF_FILE="${SCRIPT_DIR}/MultiGWAS-Explorer_pipeline.def"
OUTPUT_IMAGE="${1:-${REPO_ROOT}/MultiGWAS-Explorer_pipeline.sif}"
APPTAINER_BIN="${APPTAINER_BIN:-apptainer}"

cd "${REPO_ROOT}"
"${APPTAINER_BIN}" build "${OUTPUT_IMAGE}" "${DEF_FILE}"
echo "Built Apptainer image: ${OUTPUT_IMAGE}"
