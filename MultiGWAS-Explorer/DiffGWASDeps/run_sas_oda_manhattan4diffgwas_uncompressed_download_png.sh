#!/usr/bin/env bash
set -euo pipefail

DEPS_DIR="${DEPS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)}"
WORKDIR="${WORKDIR:-$(cd "${DEPS_DIR}/.." && pwd -P)}"
ODA_HELPER_SCRIPT="${ODA_HELPER_SCRIPT:-${DEPS_DIR}/run_sas_codes_or_script_in_ODA.pl}"
if [[ -f "${ODA_HELPER_SCRIPT}" ]]; then
  ODA_PERL=(perl "${ODA_HELPER_SCRIPT}")
else
  ODA_PERL=(perl -S run_sas_codes_or_script_in_ODA.pl)
fi
RUNNER_CONFIG_JSON="${RUNNER_CONFIG_JSON:-}"
if [[ -n "${RUNNER_CONFIG_JSON}" ]]; then
  eval "$("perl" "${DEPS_DIR}/emit_diff_gwas_runner_env.pl" --config "${RUNNER_CONFIG_JSON}")"
fi
PROJECT_TAG="${PROJECT_TAG:-PGC_SCZ}"
DEFAULT_DATA_TSV="/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects.stdized.wide_beta_se_p_p_lt_0p05.final.tsv"
DATA_TSV="${DATA_TSV:-${DEFAULT_DATA_TSV}}"
DEFAULT_REMOTE_DATA_BASENAME="PGC_SCZ_female_vs_male_diff_effects.stdized.wide_beta_se_p_p_lt_0p05.final.tsv"
# Override DATA_TSV to point at a different local uncompressed subset. By default
# the uploaded remote filename matches the local basename.
REMOTE_DATA_BASENAME="${REMOTE_DATA_BASENAME:-$(basename "${DATA_TSV}")}"
MACRO_SAS="${DEPS_DIR}/Manhattan4DiffGWASs_png.sas"
RUN_SAS_TEMPLATE="${DEPS_DIR}/run_sas_oda_manhattan4diffgwas_uncompressed.sas"
SCHEMA_CONFIG_JSON="${SCHEMA_CONFIG_JSON:-${EXTRACTOR_CONFIG_JSON:-${WORKDIR}/configs/preset_pgc_scz_sex_diff.json}}"
SCHEMA_INCLUDE_HELPER="${SCHEMA_INCLUDE_HELPER:-${DEPS_DIR}/generate_sas_wide_import_include.pl}"
RENDER_SAS_HELPER="${RENDER_SAS_HELPER:-${DEPS_DIR}/render_sas_template.pl}"
SESSION_ID="${SESSION_ID:-mysession}"
CLEAN_ODA_INPUT="${CLEAN_ODA_INPUT:-1}"
MANHATTAN_GWAS_MODE="${MANHATTAN_GWAS_MODE:-multi}"
DEFAULT_MULTI_P_VAR="ALL_STD_P"
DEFAULT_MULTI_OTHER_P_VARS="ASN_STD_P EUR_STD_P ALL_DIFF_P ASN_DIFF_P EUR_DIFF_P"
DEFAULT_MULTI_GWAS_LABEL_NAMES="All standardized P|Asian standardized P|European standardized P|All female-vs-male diff P|Asian female-vs-male diff P|European female-vs-male diff P"
DEFAULT_SINGLE_P_VAR="ALL_STD_P"
DEFAULT_SINGLE_OTHER_P_VARS=""
DEFAULT_SINGLE_GWAS_LABEL_NAMES="All standardized P"
if [[ "${MANHATTAN_GWAS_MODE}" == "single" ]]; then
  DEFAULT_MANHATTAN_P_VAR="${DEFAULT_SINGLE_P_VAR}"
  DEFAULT_MANHATTAN_OTHER_P_VARS="${DEFAULT_SINGLE_OTHER_P_VARS}"
  DEFAULT_MANHATTAN_GWAS_LABEL_NAMES="${DEFAULT_SINGLE_GWAS_LABEL_NAMES}"
  DEFAULT_MANHATTAN_FIG_HEIGHT="420"
else
  DEFAULT_MANHATTAN_P_VAR="${DEFAULT_MULTI_P_VAR}"
  DEFAULT_MANHATTAN_OTHER_P_VARS="${DEFAULT_MULTI_OTHER_P_VARS}"
  DEFAULT_MANHATTAN_GWAS_LABEL_NAMES="${DEFAULT_MULTI_GWAS_LABEL_NAMES}"
  DEFAULT_MANHATTAN_FIG_HEIGHT="820"
fi
MANHATTAN_P_VAR="${MANHATTAN_P_VAR:-${DEFAULT_MANHATTAN_P_VAR}}"
MANHATTAN_OTHER_P_VARS="${MANHATTAN_OTHER_P_VARS:-${DEFAULT_MANHATTAN_OTHER_P_VARS}}"
MANHATTAN_GWAS_LABEL_NAMES="${MANHATTAN_GWAS_LABEL_NAMES:-${DEFAULT_MANHATTAN_GWAS_LABEL_NAMES}}"
MANHATTAN_FIG_HEIGHT="${MANHATTAN_FIG_HEIGHT:-${DEFAULT_MANHATTAN_FIG_HEIGHT}}"
OUTPUT_PREFIX="${OUTPUT_PREFIX:-${PROJECT_TAG}_SAS_manhattan}"
HTML_TITLE="${HTML_TITLE:-${PROJECT_TAG} SAS Manhattan Plot}"
OPEN_RESULT="${OPEN_RESULT:-1}"

open_html_result() {
  local target="${1:-}"
  [[ -n "${target}" && -e "${target}" ]] || return 1
  if command -v cygstart >/dev/null 2>&1; then
    echo "Opening HTML result with cygstart..."
    cygstart "${target}" >/dev/null 2>&1 &
    return 0
  fi
  if [[ "$(uname -s 2>/dev/null || printf '')" == "Darwin" ]] && command -v open >/dev/null 2>&1; then
    echo "Opening HTML result with open..."
    open "${target}" >/dev/null 2>&1 &
    return 0
  fi
  if command -v xdg-open >/dev/null 2>&1; then
    echo "Opening HTML result with xdg-open..."
    xdg-open "${target}" >/dev/null 2>&1 &
    return 0
  fi
  if command -v python3 >/dev/null 2>&1; then
    echo "Opening HTML result with python3 webbrowser..."
    python3 -c 'import pathlib, sys, webbrowser; webbrowser.open(pathlib.Path(sys.argv[1]).resolve().as_uri())' "${target}" >/dev/null 2>&1 &
    return 0
  fi
  if command -v python >/dev/null 2>&1; then
    echo "Opening HTML result with python webbrowser..."
    python -c 'import pathlib, sys, webbrowser; webbrowser.open(pathlib.Path(sys.argv[1]).resolve().as_uri())' "${target}" >/dev/null 2>&1 &
    return 0
  fi
  echo "WARN: Could not auto-open ${target}. Install xdg-open, use macOS open, or set OPEN_RESULT=0." >&2
  return 1
}

cd "${WORKDIR}"

stamp="$(date +%Y%m%d_%H%M%S)"
PNG_OUT="${WORKDIR}/${OUTPUT_PREFIX}.png"
HTML_OUT="${WORKDIR}/${OUTPUT_PREFIX}_png.html"
RUN_SAS_RENDERED="${WORKDIR}/run_sas_oda_manhattan4diffgwas_uncompressed.${stamp}.sas"
IMPORT_BLOCK_RENDERED="${WORKDIR}/auto_wide_import_manhattan_uncompressed.${stamp}.sas"

trap 'rm -f "${RUN_SAS_RENDERED}" "${IMPORT_BLOCK_RENDERED}"' EXIT

perl "${SCHEMA_INCLUDE_HELPER}" \
  --config "${SCHEMA_CONFIG_JSON}" \
  --dataset scz_mh \
  --source-type plain \
  --remote-basename "${REMOTE_DATA_BASENAME}" > "${IMPORT_BLOCK_RENDERED}"

perl "${RENDER_SAS_HELPER}" \
  --template "${RUN_SAS_TEMPLATE}" \
  --output "${RUN_SAS_RENDERED}" \
  --replace "MANHATTAN_P_VAR=${MANHATTAN_P_VAR}" \
  --replace "MANHATTAN_OTHER_P_VARS=${MANHATTAN_OTHER_P_VARS}" \
  --replace "MANHATTAN_GWAS_LABEL_NAMES=${MANHATTAN_GWAS_LABEL_NAMES}" \
  --replace "MANHATTAN_FIG_HEIGHT=${MANHATTAN_FIG_HEIGHT}" \
  --replace "OUTPUT_PREFIX=${OUTPUT_PREFIX}" \
  --replace "HTML_TITLE=${HTML_TITLE}" \
  --replace-file "WIDE_IMPORT_BLOCK=${IMPORT_BLOCK_RENDERED}"

rm -f "${PNG_OUT}" "${HTML_OUT}"

echo "[1/5] Uploading PNG Manhattan macro to SAS ODA..."
"${ODA_PERL[@]}" \
  --upload-file "${MACRO_SAS}" \
  --persistent \
  --session-id "${SESSION_ID}" \
  --output-prefix "upload_manhattan_png_macro_${stamp}"

echo "[2/5] Uploading uncompressed Manhattan subset to SAS ODA..."
"${ODA_PERL[@]}" \
  --upload-file "${DATA_TSV}" \
  --persistent \
  --session-id "${SESSION_ID}" \
  --output-prefix "upload_manhattan_subset_tsv_${stamp}"

echo "[3/5] Running SAS Manhattan plot..."
"${ODA_PERL[@]}" \
  --file "${RUN_SAS_RENDERED}" \
  --persistent \
  --session-id "${SESSION_ID}" \
  --output-prefix "run_manhattan_png_tsv_${stamp}"

echo "[4/5] Downloading PNG and small HTML wrapper..."
"${ODA_PERL[@]}" \
  --download-file "~/${OUTPUT_PREFIX}.png" \
  --download-local-path "${PNG_OUT}" \
  --persistent \
  --session-id "${SESSION_ID}" \
  --output-prefix "download_manhattan_png_tsv_${stamp}"

"${ODA_PERL[@]}" \
  --download-file "~/${OUTPUT_PREFIX}_png.html" \
  --download-local-path "${HTML_OUT}" \
  --persistent \
  --session-id "${SESSION_ID}" \
  --output-prefix "download_manhattan_html_tsv_${stamp}"

if [[ ! -s "${PNG_OUT}" ]]; then
  echo "ERROR: Expected downloaded PNG was not created or is empty: ${PNG_OUT}" >&2
  exit 1
fi

if [[ ! -s "${HTML_OUT}" ]]; then
  echo "ERROR: Expected downloaded HTML was not created or is empty: ${HTML_OUT}" >&2
  exit 1
fi

echo "Verified PNG:  ${PNG_OUT} ($(wc -c < "${PNG_OUT}") bytes)"
echo "Verified HTML: ${HTML_OUT} ($(wc -c < "${HTML_OUT}") bytes)"

if [[ "${CLEAN_ODA_INPUT}" == "1" ]]; then
  echo "[5/5] Removing uploaded TSV input from SAS ODA to save space..."
  "${ODA_PERL[@]}" \
    --delete-file "${REMOTE_DATA_BASENAME}" \
    --persistent \
    --session-id "${SESSION_ID}" \
    --output-prefix "cleanup_manhattan_input_tsv_${stamp}"
else
  echo "[5/5] Keeping uploaded TSV input in SAS ODA because CLEAN_ODA_INPUT=${CLEAN_ODA_INPUT}."
fi

echo "Done."
echo "Downloaded PNG:  ${PNG_OUT}"
echo "Downloaded HTML: ${HTML_OUT}"

if [[ "${OPEN_RESULT}" == "1" ]]; then
  open_html_result "${HTML_OUT}" || true
else
  echo "Not opening result because OPEN_RESULT=${OPEN_RESULT}."
fi
