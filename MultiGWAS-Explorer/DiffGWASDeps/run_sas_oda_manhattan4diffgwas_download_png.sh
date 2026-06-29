#!/usr/bin/env bash
set -euo pipefail

DEPS_DIR="${DEPS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)}"
WORKDIR="${WORKDIR:-$(cd "${DEPS_DIR}/.." && pwd -P)}"
RUNNER_CONFIG_JSON="${RUNNER_CONFIG_JSON:-}"
if [[ -n "${RUNNER_CONFIG_JSON}" ]]; then
  eval "$("perl" "${DEPS_DIR}/emit_diff_gwas_runner_env.pl" --config "${RUNNER_CONFIG_JSON}")"
fi
PROJECT_TAG="${PROJECT_TAG:-PGC_SCZ}"
DEFAULT_DATA_GZ="/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects.stdized.wide_beta_se_p_p_lt_0p05.final.tsv.gz"
DATA_GZ="${DATA_GZ:-${DEFAULT_DATA_GZ}}"
DEFAULT_REMOTE_DATA_BASENAME="PGC_SCZ_female_vs_male_diff_effects.stdized.wide_beta_se_p_p_lt_0p05.final.tsv.gz"
# Override DATA_GZ to point at a different local gz subset. By default the
# uploaded remote filename matches the local basename.
REMOTE_DATA_BASENAME="${REMOTE_DATA_BASENAME:-$(basename "${DATA_GZ}")}"
MACRO_SAS="${DEPS_DIR}/Manhattan4DiffGWASs_png.sas"
RUN_SAS_TEMPLATE="${DEPS_DIR}/run_sas_oda_manhattan4diffgwas.sas"
SCHEMA_CONFIG_JSON="${SCHEMA_CONFIG_JSON:-${EXTRACTOR_CONFIG_JSON:-${WORKDIR}/configs/preset_pgc_scz_sex_diff.json}}"
SCHEMA_INCLUDE_HELPER="${SCHEMA_INCLUDE_HELPER:-${DEPS_DIR}/generate_sas_wide_import_include.pl}"
RENDER_SAS_HELPER="${RENDER_SAS_HELPER:-${DEPS_DIR}/render_sas_template.pl}"
LOCAL_SAS_DEBUG_EMITTER="${LOCAL_SAS_DEBUG_EMITTER:-${DEPS_DIR}/emit_local_sas_debug_script.pl}"
SESSION_ID="${SESSION_ID:-mysession}"
USE_PERSISTENT_SESSION="${USE_PERSISTENT_SESSION:-0}"
CLEAN_ODA_INPUT="${CLEAN_ODA_INPUT:-1}"
SKIP_DATA_UPLOAD="${SKIP_DATA_UPLOAD:-0}"
KEEP_REMOTE_PLOT_DATA="${KEEP_REMOTE_PLOT_DATA:-0}"
EMIT_LOCAL_SAS_DEBUG="${EMIT_LOCAL_SAS_DEBUG:-0}"
LOCAL_SAS_DEBUG_ONLY="${LOCAL_SAS_DEBUG_ONLY:-0}"
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
DEFAULT_MANHATTAN_FIG_WIDTH="${DEFAULT_MANHATTAN_FIG_WIDTH:-1800}"
MANHATTAN_FIG_WIDTH="${MANHATTAN_FIG_WIDTH:-${DEFAULT_MANHATTAN_FIG_WIDTH}}"
MANHATTAN_FIG_HEIGHT="${MANHATTAN_FIG_HEIGHT:-${DEFAULT_MANHATTAN_FIG_HEIGHT}}"
MANHATTAN_FONTSIZE="${MANHATTAN_FONTSIZE:-2.4}"
MANHATTAN_Y_AXIS_LABEL_SIZE="${MANHATTAN_Y_AXIS_LABEL_SIZE:-2.4}"
MANHATTAN_Y_AXIS_VALUE_SIZE="${MANHATTAN_Y_AXIS_VALUE_SIZE:-2.2}"
MANHATTAN_GWAS_LABEL_SIZE="${MANHATTAN_GWAS_LABEL_SIZE:-2.4}"
MANHATTAN_GWAS_LABEL_HALO_SIZE="${MANHATTAN_GWAS_LABEL_HALO_SIZE:-2.4}"
OUTPUT_PREFIX="${OUTPUT_PREFIX:-${PROJECT_TAG}_SAS_manhattan}"
HTML_TITLE="${HTML_TITLE:-${PROJECT_TAG} SAS Manhattan Plot}"
# Open the downloaded HTML by default. Set OPEN_RESULT=0 for non-interactive
# validation runs when you only want to download and verify the files.
OPEN_RESULT="${OPEN_RESULT:-1}"

cd "${WORKDIR}"

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

ODA_HELPER_SCRIPT="${ODA_HELPER_SCRIPT:-${DEPS_DIR}/run_sas_codes_or_script_in_ODA.pl}"
if [[ ! -f "${ODA_HELPER_SCRIPT}" ]]; then
  ODA_HELPER_SCRIPT="${WORKDIR}/run_sas_codes_or_script_in_ODA.pl"
fi
if [[ -f "${ODA_HELPER_SCRIPT}" ]]; then
  ODA_PERL_BASE=(perl "${ODA_HELPER_SCRIPT}")
else
  ODA_PERL_BASE=(perl -S run_sas_codes_or_script_in_ODA.pl)
fi

if [[ "${USE_PERSISTENT_SESSION}" == "1" ]]; then
  ODA_PERL=("${ODA_PERL_BASE[@]}" --persistent --session-id "${SESSION_ID}")
else
  ODA_PERL=("${ODA_PERL_BASE[@]}")
fi

resolve_python_bin() {
  local cand
  local record_file="${WORKDIR}/.venv-pipeline/.python-bin"
  if [[ -f "${record_file}" ]]; then
    cand="$(head -n 1 "${record_file}" 2>/dev/null || true)"
    if [[ -n "${cand}" ]]; then
      if command -v "${cand}" >/dev/null 2>&1 || [[ -x "${cand}" ]]; then
        printf '%s\n' "${cand}"
        return 0
      fi
    fi
  fi
  for cand in \
    "${PIPELINE_PYTHON_BIN:-}" \
    "${WORKDIR}/.venv-pipeline/bin/python" \
    "${WORKDIR}/.venv-pipeline/bin/python3" \
    "${WORKDIR}/.venv-pipeline/Scripts/python.exe" \
    "${WORKDIR}/.venv-pipeline/Scripts/python" \
    python3 \
    python; do
    [[ -n "${cand}" ]] || continue
    if command -v "${cand}" >/dev/null 2>&1; then
      printf '%s\n' "${cand}"
      return 0
    fi
    if [[ -x "${cand}" ]]; then
      printf '%s\n' "${cand}"
      return 0
    fi
  done
  return 1
}

PYTHON_BIN="$(resolve_python_bin || true)"
if [[ -z "${PYTHON_BIN}" ]]; then
  echo "ERROR: Could not resolve a Python interpreter for PNG post-processing." >&2
  exit 1
fi

stamp="$(date +%Y%m%d_%H%M%S)"
PNG_OUT="${WORKDIR}/${OUTPUT_PREFIX}.png"
HTML_OUT="${WORKDIR}/${OUTPUT_PREFIX}_png.html"
RUN_SAS_RENDERED="${WORKDIR}/run_sas_oda_manhattan4diffgwas.${stamp}.sas"
LOCAL_DEBUG_SAS_RENDERED="${WORKDIR}/run_sas_local_debug_manhattan4diffgwas.${stamp}.sas"
IMPORT_BLOCK_RENDERED="${WORKDIR}/auto_wide_import_manhattan.${stamp}.sas"

trap 'rm -f "${RUN_SAS_RENDERED}" "${IMPORT_BLOCK_RENDERED}"' EXIT

oda_download_many() {
  local output_prefix="$1"
  shift
  "${ODA_PERL[@]}" "$@" --output-prefix "${output_prefix}"
}

oda_delete_many() {
  local output_prefix="$1"
  shift
  "${ODA_PERL[@]}" "$@" --output-prefix "${output_prefix}"
}

recenter_png_horizontally() {
  local png_path="$1"
  "${PYTHON_BIN}" - "$png_path" <<'PY'
from PIL import Image
import sys

png_path = sys.argv[1]
img = Image.open(png_path).convert("RGBA")
w, h = img.size
rgb = img.convert("RGB")

left = None
right = None
for x in range(w):
    column_has_content = False
    for y in range(h):
        r, g, b = rgb.getpixel((x, y))
        if min(r, g, b) < 250:
            column_has_content = True
            break
    if column_has_content:
        if left is None:
            left = x
        right = x

if left is None or right is None:
    sys.exit(0)

right_margin = (w - 1) - right
if abs(left - right_margin) <= 4:
    sys.exit(0)

content = img.crop((left, 0, right + 1, h))
canvas = Image.new("RGBA", (w, h), (255, 255, 255, 255))
new_left = max(0, (w - content.width) // 2)
canvas.paste(content, (new_left, 0))
canvas.save(png_path)
PY
}

remote_data_exists() {
  local check_output
  check_output="$(
    "${ODA_PERL[@]}" \
      --dir4listing '~' \
      --output-prefix "check_remote_manhattan_data_${stamp}" 2>&1 || true
  )"
  printf '%s\n' "${check_output}" | grep -Fxq "${REMOTE_DATA_BASENAME}"
}

remote_data_size_bytes() {
  local info_output
  info_output="$(
    "${ODA_PERL[@]}" \
      --file-info "~/${REMOTE_DATA_BASENAME}" \
      --output-prefix "check_remote_manhattan_size_${stamp}" 2>&1 || true
  )"
  printf '%s\n' "${info_output}" | awk -F '\t' '$1=="SIZE"{print $2}' | tail -n 1
}

remote_data_matches_local_size() {
  [[ -s "${DATA_GZ}" ]] || return 1
  local local_size remote_size
  local_size="$(wc -c < "${DATA_GZ}" | tr -d '[:space:]')"
  remote_size="$(remote_data_size_bytes | tr -d '[:space:]')"
  [[ -n "${local_size}" && -n "${remote_size}" && "${local_size}" == "${remote_size}" ]]
}

delete_partial_remote_data() {
  "${ODA_PERL[@]}" \
    --delete-file "${REMOTE_DATA_BASENAME}" \
    --output-prefix "delete_partial_manhattan_subset_${stamp}" >/dev/null 2>&1 || true
}

upload_data_with_integrity_check() {
  local upload_cmd local_size remote_size
  local_size="$(wc -c < "${DATA_GZ}" | tr -d '[:space:]')"
  upload_cmd=("${ODA_PERL[@]}" --upload-file "${DATA_GZ}" --output-prefix "upload_manhattan_subset_${stamp}")
  if command -v timeout >/dev/null 2>&1; then
    timeout "${ODA_UPLOAD_TIMEOUT_SECONDS:-7200}" "${upload_cmd[@]}"
  else
    "${upload_cmd[@]}"
  fi
  remote_size="$(remote_data_size_bytes | tr -d '[:space:]')"
  if [[ -z "${remote_size}" || "${remote_size}" != "${local_size}" ]]; then
    echo "ERROR: Remote upload size mismatch for ${REMOTE_DATA_BASENAME}. local=${local_size} remote=${remote_size:-missing}" >&2
    delete_partial_remote_data
    return 1
  fi
}

perl "${SCHEMA_INCLUDE_HELPER}" \
  --config "${SCHEMA_CONFIG_JSON}" \
  --dataset scz_mh \
  --source-type gzip \
  --remote-basename "${REMOTE_DATA_BASENAME}" > "${IMPORT_BLOCK_RENDERED}"

perl "${RENDER_SAS_HELPER}" \
  --template "${RUN_SAS_TEMPLATE}" \
  --output "${RUN_SAS_RENDERED}" \
  --replace "MANHATTAN_P_VAR=${MANHATTAN_P_VAR}" \
  --replace "MANHATTAN_OTHER_P_VARS=${MANHATTAN_OTHER_P_VARS}" \
  --replace "MANHATTAN_GWAS_LABEL_NAMES=${MANHATTAN_GWAS_LABEL_NAMES}" \
  --replace "MANHATTAN_FIG_WIDTH=${MANHATTAN_FIG_WIDTH}" \
  --replace "MANHATTAN_FIG_HEIGHT=${MANHATTAN_FIG_HEIGHT}" \
  --replace "MANHATTAN_FONTSIZE=${MANHATTAN_FONTSIZE}" \
  --replace "MANHATTAN_Y_AXIS_LABEL_SIZE=${MANHATTAN_Y_AXIS_LABEL_SIZE}" \
  --replace "MANHATTAN_Y_AXIS_VALUE_SIZE=${MANHATTAN_Y_AXIS_VALUE_SIZE}" \
  --replace "MANHATTAN_GWAS_LABEL_SIZE=${MANHATTAN_GWAS_LABEL_SIZE}" \
  --replace "MANHATTAN_GWAS_LABEL_HALO_SIZE=${MANHATTAN_GWAS_LABEL_HALO_SIZE}" \
  --replace "OUTPUT_PREFIX=${OUTPUT_PREFIX}" \
  --replace "HTML_TITLE=${HTML_TITLE}" \
  --replace-file "WIDE_IMPORT_BLOCK=${IMPORT_BLOCK_RENDERED}"

if [[ "${EMIT_LOCAL_SAS_DEBUG}" == "1" || "${LOCAL_SAS_DEBUG_ONLY}" == "1" ]]; then
  perl "${LOCAL_SAS_DEBUG_EMITTER}" \
    --mode manhattan \
    --input "${RUN_SAS_RENDERED}" \
    --output "${LOCAL_DEBUG_SAS_RENDERED}" \
    --workdir "${WORKDIR}" \
    --deps-dir "${DEPS_DIR}" \
    --data-gz "${DATA_GZ}" \
    --manhattan-macro "${MACRO_SAS}" \
    --output-html-basename "$(basename "${HTML_OUT}")" >/dev/null
  echo "[local-sas] Emitted local-SAS debug script: ${LOCAL_DEBUG_SAS_RENDERED}"
fi

if [[ "${LOCAL_SAS_DEBUG_ONLY}" == "1" ]]; then
  echo "[local-sas] LOCAL_SAS_DEBUG_ONLY=1, skipping SAS ODA submit."
  exit 0
fi

rm -f "${PNG_OUT}" "${HTML_OUT}"

echo "[1/5] Uploading PNG Manhattan macro to SAS ODA..."
# File uploads land in the shared SAS ODA home directory, so they do not need
# to flow through the persistent session server.
"${ODA_PERL[@]}" \
  --upload-file "${MACRO_SAS}" \
  --output-prefix "upload_manhattan_png_macro_${stamp}"

if [[ "${SKIP_DATA_UPLOAD}" == "1" ]]; then
  echo "[2/5] Reusing already-uploaded gzipped Manhattan subset in SAS ODA: ${REMOTE_DATA_BASENAME}"
elif [[ "${KEEP_REMOTE_PLOT_DATA}" == "1" ]] && remote_data_exists && remote_data_matches_local_size; then
  echo "[2/5] Keeping and reusing existing remote Manhattan subset in SAS ODA: ${REMOTE_DATA_BASENAME}"
else
  echo "[2/5] Uploading gzipped Manhattan subset to SAS ODA..."
  if remote_data_exists && ! remote_data_matches_local_size; then
    echo "[repair] Existing remote file has the wrong size and will be replaced."
    delete_partial_remote_data
  fi
  # Large GWAS subset uploads are more reliable through a one-shot ODA
  # connection than through the persistent session server.
  upload_data_with_integrity_check
fi

echo "[3/5] Running SAS Manhattan plot..."
"${ODA_PERL[@]}" \
  --file "${RUN_SAS_RENDERED}" \
  --output-prefix "run_manhattan_png_${stamp}"

echo "[4/5] Downloading PNG and small HTML wrapper..."
oda_download_many \
  "download_manhattan_outputs_${stamp}" \
  --download-file "~/${OUTPUT_PREFIX}.png" \
  --download-local-path "${PNG_OUT}" \
  --download-file "~/${OUTPUT_PREFIX}_png.html" \
  --download-local-path "${HTML_OUT}"

if [[ ! -s "${PNG_OUT}" ]]; then
  echo "ERROR: Expected downloaded PNG was not created or is empty: ${PNG_OUT}" >&2
  exit 1
fi

recenter_png_horizontally "${PNG_OUT}"

if [[ ! -s "${HTML_OUT}" ]]; then
  echo "ERROR: Expected downloaded HTML was not created or is empty: ${HTML_OUT}" >&2
  exit 1
fi

echo "Verified PNG:  ${PNG_OUT} ($(wc -c < "${PNG_OUT}") bytes)"
echo "Verified HTML: ${HTML_OUT} ($(wc -c < "${HTML_OUT}") bytes)"

if [[ "${KEEP_REMOTE_PLOT_DATA}" == "1" ]]; then
  echo "[5/5] Keeping uploaded gz input in SAS ODA because KEEP_REMOTE_PLOT_DATA=${KEEP_REMOTE_PLOT_DATA}."
elif [[ "${CLEAN_ODA_INPUT}" == "1" ]]; then
  echo "[5/5] Removing uploaded gz input from SAS ODA to save space..."
  oda_delete_many \
    "cleanup_manhattan_input_${stamp}" \
    --delete-file "${REMOTE_DATA_BASENAME}"
else
  echo "[5/5] Keeping uploaded gz input in SAS ODA because CLEAN_ODA_INPUT=${CLEAN_ODA_INPUT}."
fi

echo "Done."
echo "Downloaded PNG:  ${PNG_OUT}"
echo "Downloaded HTML: ${HTML_OUT}"

if [[ "${OPEN_RESULT}" == "1" ]]; then
  open_html_result "${HTML_OUT}" || true
else
  echo "Not opening result because OPEN_RESULT=${OPEN_RESULT}."
fi
