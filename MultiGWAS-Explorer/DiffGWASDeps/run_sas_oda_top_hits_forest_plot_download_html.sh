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
LOCAL_TOP_HITS_CSV_HELPER="${LOCAL_TOP_HITS_CSV_HELPER:-${DEPS_DIR}/generate_requested_top_hits_csv.pl}"
RENDER_SAS_HELPER="${RENDER_SAS_HELPER:-${DEPS_DIR}/render_sas_template.pl}"
RUN_SAS_TEMPLATE="${RUN_SAS_TEMPLATE:-${DEPS_DIR}/run_sas_oda_top_hits_forest_plot.sas}"
SESSION_ID="${SESSION_ID:-mysession}"
USE_PERSISTENT_SESSION="${USE_PERSISTENT_SESSION:-0}"
CLEAN_ODA_INPUT="${CLEAN_ODA_INPUT:-1}"
SKIP_DATA_UPLOAD="${SKIP_DATA_UPLOAD:-0}"
KEEP_REMOTE_PLOT_DATA="${KEEP_REMOTE_PLOT_DATA:-0}"
OPEN_RESULT="${OPEN_RESULT:-1}"
EMIT_LOCAL_SAS_DEBUG="${EMIT_LOCAL_SAS_DEBUG:-0}"
LOCAL_SAS_DEBUG_ONLY="${LOCAL_SAS_DEBUG_ONLY:-0}"
REUSE_FOREST_TOP_HITS_CSV="${REUSE_FOREST_TOP_HITS_CSV:-1}"

TARGET_SNP_LIST="${TARGET_SNP_LIST:-}"
TARGET_SNP_GENES="${TARGET_SNP_GENES:-}"
TOP_HIT_MODE="${TOP_HIT_MODE:-differential}"
TOP_HIT_FOCUS_PVAR="${TOP_HIT_FOCUS_PVAR:-ALL_STD_P}"
TOP_HIT_SIGNAL_THRSHD="${TOP_HIT_SIGNAL_THRSHD:-1e-6}"
TOP_HIT_SIGNAL_THRSHDS="${TOP_HIT_SIGNAL_THRSHDS:-${TOP_HIT_SIGNAL_THRSHD}}"
TOP_HIT_DIST_BP="${TOP_HIT_DIST_BP:-1e6}"
TOP_HIT_MAF_THRESHOLD="${TOP_HIT_MAF_THRESHOLD:-0.01}"
TOP_HIT_MAX_LOCI="${TOP_HIT_MAX_LOCI:-0}"
TOP_HIT_GNOMAD_FREQ_FILE="${TOP_HIT_GNOMAD_FREQ_FILE:-}"
TOP_HIT_GNOMAD_POP_MAP="${TOP_HIT_GNOMAD_POP_MAP:-}"

FOREST_OUTPUT_PREFIX="${FOREST_OUTPUT_PREFIX:-${PROJECT_TAG}_SAS_top_hits_forest}"
FOREST_HTML_TITLE="${FOREST_HTML_TITLE:-${PROJECT_TAG} top-hit forest plots}"
FOREST_TOP_HITS_CSV_BASENAME="${FOREST_TOP_HITS_CSV_BASENAME:-${FOREST_OUTPUT_PREFIX}_top_hits.csv}"
FOREST_OUTPUT_HTML_BASENAME="${FOREST_OUTPUT_HTML_BASENAME:-${FOREST_OUTPUT_PREFIX}.html}"
FOREST_OUTPUT_MANIFEST_BASENAME="${FOREST_OUTPUT_MANIFEST_BASENAME:-${FOREST_OUTPUT_PREFIX}.manifest.tsv}"
FOREST_TRACK_IDS="${FOREST_TRACK_IDS:-}"
FOREST_TRACK_LABELS="${FOREST_TRACK_LABELS:-}"
FOREST_TRACK_BETA_VARS="${FOREST_TRACK_BETA_VARS:-}"
FOREST_TRACK_SE_VARS="${FOREST_TRACK_SE_VARS:-}"
FOREST_TRACK_P_VARS="${FOREST_TRACK_P_VARS:-}"
FOREST_TRACK_COUNT="${FOREST_TRACK_COUNT:-0}"
FOREST_DEFAULT_HIT_CLASS="${FOREST_DEFAULT_HIT_CLASS:-DIFFERENTIAL}"
FOREST_FIG_WIDTH="${FOREST_FIG_WIDTH:-900}"
FOREST_FIG_HEIGHT="${FOREST_FIG_HEIGHT:-}"
FOREST_DOTSIZE="${FOREST_DOTSIZE:-8}"
FOREST_Y_FONT_SIZE="${FOREST_Y_FONT_SIZE:-12}"
FOREST_MIN_AXIS="${FOREST_MIN_AXIS:-0.4}"
FOREST_MAX_AXIS="${FOREST_MAX_AXIS:-1.6}"
FOREST_XAXIS_VALUE_RANGE="${FOREST_XAXIS_VALUE_RANGE:-0.4 to 1.6 by 0.2}"

RAND_MACRO_SAS="${DEPS_DIR}/RandBetween.sas"
MKFMT_MACRO_SAS="${DEPS_DIR}/mkfmt4grps_by_var.sas"
FOREST_MACRO_SAS="${DEPS_DIR}/beta2OR_forest_plot.sas"

cd "${WORKDIR}"

ODA_HELPER_SCRIPT="${ODA_HELPER_SCRIPT:-${DEPS_DIR}/run_sas_codes_or_script_in_ODA.pl}"
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

open_html_result() {
  local target="${1:-}"
  [[ -n "${target}" && -e "${target}" ]] || return 1
  if command -v cygstart >/dev/null 2>&1; then
    echo "Opening forest HTML result with cygstart..."
    cygstart "${target}" >/dev/null 2>&1 &
    return 0
  fi
  if [[ "$(uname -s 2>/dev/null || printf '')" == "Darwin" ]] && command -v open >/dev/null 2>&1; then
    echo "Opening forest HTML result with open..."
    open "${target}" >/dev/null 2>&1 &
    return 0
  fi
  if command -v xdg-open >/dev/null 2>&1; then
    echo "Opening forest HTML result with xdg-open..."
    xdg-open "${target}" >/dev/null 2>&1 &
    return 0
  fi
  if command -v python3 >/dev/null 2>&1; then
    echo "Opening forest HTML result with python3 webbrowser..."
    python3 -c 'import pathlib, sys, webbrowser; webbrowser.open(pathlib.Path(sys.argv[1]).resolve().as_uri())' "${target}" >/dev/null 2>&1 &
    return 0
  fi
  if command -v python >/dev/null 2>&1; then
    echo "Opening forest HTML result with python webbrowser..."
    python -c 'import pathlib, sys, webbrowser; webbrowser.open(pathlib.Path(sys.argv[1]).resolve().as_uri())' "${target}" >/dev/null 2>&1 &
    return 0
  fi
  echo "WARN: Could not auto-open ${target}. Install xdg-open, use macOS open, or set OPEN_RESULT=0." >&2
  return 1
}

stamp="$(date +%Y%m%d_%H%M%S)"

stable_hash_text() {
  perl -MDigest::MD5=md5_hex -e 'print md5_hex(join("\0", @ARGV));' "$@"
}

if [[ -n "${TARGET_SNP_LIST}" ]]; then
  IFS=',' read -r -a _forest_target_snp_names <<< "${TARGET_SNP_LIST}"
  if [[ "${#_forest_target_snp_names[@]}" -eq 1 ]]; then
    _forest_target_tag="$(printf '%s' "${_forest_target_snp_names[0]}" | tr -c 'A-Za-z0-9._-' '_')"
  else
    _forest_target_tag="targets_$(stable_hash_text "${TARGET_SNP_LIST}" | cut -c1-12)"
  fi
  FOREST_TOP_HITS_CSV_BASENAME="${FOREST_TOP_HITS_CSV_BASENAME%.csv}_${_forest_target_tag}.csv"
fi

HTML_OUT="${WORKDIR}/${FOREST_OUTPUT_HTML_BASENAME}"
CSV_OUT="${WORKDIR}/${FOREST_TOP_HITS_CSV_BASENAME}"
MANIFEST_OUT="${WORKDIR}/${FOREST_OUTPUT_MANIFEST_BASENAME}"
RUN_SAS_RENDERED="${WORKDIR}/run_sas_oda_top_hits_forest_plot.${stamp}.sas"
LOCAL_DEBUG_SAS_RENDERED="${WORKDIR}/run_sas_local_debug_top_hits_forest_plot.${stamp}.sas"
REMOTE_TOP_HITS_BASENAME="$(basename "${CSV_OUT}")"
REMOTE_HTML_BASENAME="$(basename "${HTML_OUT}")"
REMOTE_MANIFEST_BASENAME="$(basename "${MANIFEST_OUT}")"

trap 'rm -f "${RUN_SAS_RENDERED}"' EXIT

run_oda_helper() {
  "${ODA_PERL[@]}" "$@"
}

oda_download_many() {
  local output_prefix="$1"
  shift
  run_oda_helper "$@" --output-prefix "${output_prefix}"
}

oda_delete_many() {
  local output_prefix="$1"
  shift
  run_oda_helper "$@" --output-prefix "${output_prefix}"
}

to_local_sas_path() {
  local p="$1"
  if command -v cygpath >/dev/null 2>&1; then
    cygpath -w "$p"
    return 0
  fi
  printf '%s\n' "$p"
}

remote_file_exists() {
  local remote_name="$1"
  local check_output
  check_output="$(
    run_oda_helper \
      --dir4listing '~' \
      --output-prefix "check_remote_forest_file_${stamp}" 2>&1 || true
  )"
  printf '%s\n' "${check_output}" | grep -Fxq "${remote_name}"
}

remote_file_size_bytes() {
  local remote_name="$1"
  local info_output
  info_output="$(
    run_oda_helper \
      --file-info "~/${remote_name}" \
      --output-prefix "check_remote_forest_size_${stamp}" 2>&1 || true
  )"
  printf '%s\n' "${info_output}" | awk -F '\t' '$1=="SIZE"{print $2}' | tail -n 1
}

remote_file_matches_local_size() {
  local local_path="$1"
  local remote_name="$2"
  [[ -s "${local_path}" ]] || return 1
  local local_size remote_size
  local_size="$(wc -c < "${local_path}" | tr -d '[:space:]')"
  remote_size="$(remote_file_size_bytes "${remote_name}" | tr -d '[:space:]')"
  [[ -n "${local_size}" && -n "${remote_size}" && "${local_size}" == "${remote_size}" ]]
}

delete_remote_file_quiet() {
  local remote_name="$1"
  oda_delete_many "cleanup_remote_forest_file_${stamp}" --delete-file "${remote_name}" >/dev/null 2>&1 || true
}

generate_requested_top_hits_csv_locally() {
  [[ -x "${LOCAL_TOP_HITS_CSV_HELPER}" || -f "${LOCAL_TOP_HITS_CSV_HELPER}" ]] || return 1
  echo "[prep] Generating requested forest-plot top-hit CSV locally..."
  local -a cmd=(
    perl "${LOCAL_TOP_HITS_CSV_HELPER}"
    --input "${DATA_GZ}"
    --output "${CSV_OUT}"
    --top-hit-mode "${TOP_HIT_MODE}"
    --top-hit-focus-pvar "${TOP_HIT_FOCUS_PVAR}"
    --top-hit-signal-thrshd "${TOP_HIT_SIGNAL_THRSHD}"
    --top-hit-signal-thrshds "${TOP_HIT_SIGNAL_THRSHDS}"
    --top-hit-dist-bp "${TOP_HIT_DIST_BP}"
    --maf-threshold "${TOP_HIT_MAF_THRESHOLD}"
    --max-hits "${TOP_HIT_MAX_LOCI}"
  )
  if [[ -n "${RUNNER_CONFIG_JSON}" ]]; then
    cmd+=(--runner-config "${RUNNER_CONFIG_JSON}")
  fi
  if [[ -n "${TARGET_SNP_LIST}" ]]; then
    cmd+=(--target-snps "${TARGET_SNP_LIST}")
  fi
  if [[ -n "${TARGET_SNP_GENES}" ]]; then
    cmd+=(--target-snp-genes "${TARGET_SNP_GENES}")
  fi
  if [[ -n "${TOP_HIT_GNOMAD_FREQ_FILE}" ]]; then
    cmd+=(--gnomad-freq-file "${TOP_HIT_GNOMAD_FREQ_FILE}")
  fi
  if [[ -n "${TOP_HIT_GNOMAD_POP_MAP}" ]]; then
    cmd+=(--gnomad-pop-map "${TOP_HIT_GNOMAD_POP_MAP}")
  fi
  "${cmd[@]}"
}

render_forest_template() {
  local output_sas="$1"
  local rand_include="$2"
  local mkfmt_include="$3"
  local forest_include="$4"
  local input_csv_path="$5"
  local output_image_prefix_path="$6"
  local output_image_basename_prefix="$7"
  local output_html_path="$8"
  local output_manifest_path="$9"

  perl "${RENDER_SAS_HELPER}" \
    --template "${RUN_SAS_TEMPLATE}" \
    --output "${output_sas}" \
    --replace "RAND_INCLUDE_PATH=${rand_include}" \
    --replace "MKFMT_INCLUDE_PATH=${mkfmt_include}" \
    --replace "FOREST_MACRO_INCLUDE_PATH=${forest_include}" \
    --replace "INPUT_CSV_PATH=${input_csv_path}" \
    --replace "OUTPUT_IMAGE_PREFIX_PATH=${output_image_prefix_path}" \
    --replace "OUTPUT_IMAGE_BASENAME_PREFIX=${output_image_basename_prefix}" \
    --replace "OUTPUT_HTML_PATH=${output_html_path}" \
    --replace "OUTPUT_MANIFEST_PATH=${output_manifest_path}" \
    --replace "FOREST_HTML_TITLE=${FOREST_HTML_TITLE}" \
    --replace "FOREST_TRACK_IDS=${FOREST_TRACK_IDS}" \
    --replace "FOREST_TRACK_LABELS=${FOREST_TRACK_LABELS}" \
    --replace "FOREST_TRACK_BETA_VARS=${FOREST_TRACK_BETA_VARS}" \
    --replace "FOREST_TRACK_SE_VARS=${FOREST_TRACK_SE_VARS}" \
    --replace "FOREST_TRACK_P_VARS=${FOREST_TRACK_P_VARS}" \
    --replace "FOREST_TRACK_COUNT=${FOREST_TRACK_COUNT}" \
    --replace "FOREST_FIG_WIDTH=${FOREST_FIG_WIDTH}" \
    --replace "FOREST_FIG_HEIGHT=${FOREST_FIG_HEIGHT}" \
    --replace "FOREST_DOTSIZE=${FOREST_DOTSIZE}" \
    --replace "FOREST_Y_FONT_SIZE=${FOREST_Y_FONT_SIZE}" \
    --replace "FOREST_MIN_AXIS=${FOREST_MIN_AXIS}" \
    --replace "FOREST_MAX_AXIS=${FOREST_MAX_AXIS}" \
    --replace "FOREST_XAXIS_VALUE_RANGE=${FOREST_XAXIS_VALUE_RANGE}" \
    --replace "FOREST_DEFAULT_HIT_CLASS=${FOREST_DEFAULT_HIT_CLASS}"
}

cleanup_remote_generated_outputs() {
  delete_remote_file_quiet "${REMOTE_HTML_BASENAME}"
  delete_remote_file_quiet "${REMOTE_MANIFEST_BASENAME}"
  local remote_pngs
  remote_pngs="$(
    run_oda_helper \
      --dir4listing '~' \
      --output-prefix "list_remote_forest_cleanup_${stamp}" 2>&1 || true
  )"
  while IFS= read -r remote_png; do
    [[ -z "${remote_png}" ]] && continue
    [[ "${remote_png}" =~ ^${FOREST_OUTPUT_PREFIX}_.+\.png$ ]] || continue
    delete_remote_file_quiet "${remote_png}"
  done <<< "${remote_pngs}"
}

[[ -n "${FOREST_TRACK_IDS}" ]] || {
  echo "ERROR: FOREST_TRACK_IDS is empty. The runner config did not resolve any single-GWAS forest panels." >&2
  exit 1
}
[[ "${FOREST_TRACK_COUNT}" =~ ^[0-9]+$ && "${FOREST_TRACK_COUNT}" -gt 0 ]] || {
  echo "ERROR: FOREST_TRACK_COUNT is invalid: ${FOREST_TRACK_COUNT}" >&2
  exit 1
}

if [[ "${REUSE_FOREST_TOP_HITS_CSV}" == "1" && -s "${CSV_OUT}" ]]; then
  echo "[prep] Reusing existing local forest top-hit CSV: ${CSV_OUT}"
else
  generate_requested_top_hits_csv_locally
fi

if [[ ! -s "${CSV_OUT}" ]]; then
  echo "ERROR: Forest top-hit CSV was not created: ${CSV_OUT}" >&2
  exit 1
fi

forest_hit_count="$(($(wc -l < "${CSV_OUT}") - 1))"
if (( forest_hit_count < 1 )); then
  echo "ERROR: Forest top-hit CSV does not contain any SNP rows: ${CSV_OUT}" >&2
  exit 1
fi

if [[ -z "${FOREST_FIG_HEIGHT}" ]]; then
  FOREST_FIG_HEIGHT=$((420 + (55 * forest_hit_count)))
  if (( FOREST_FIG_HEIGHT < 800 )); then
    FOREST_FIG_HEIGHT=800
  fi
  if (( FOREST_FIG_HEIGHT > 2200 )); then
    FOREST_FIG_HEIGHT=2200
  fi
fi

if [[ -z "${FOREST_Y_FONT_SIZE}" ]]; then
  if (( forest_hit_count <= 8 )); then
    FOREST_Y_FONT_SIZE=12
  elif (( forest_hit_count <= 16 )); then
    FOREST_Y_FONT_SIZE=11
  else
    FOREST_Y_FONT_SIZE=10
  fi
fi

render_forest_template \
  "${RUN_SAS_RENDERED}" \
  "~/RandBetween.sas" \
  "~/mkfmt4grps_by_var.sas" \
  "~/beta2OR_forest_plot.sas" \
  "~/${REMOTE_TOP_HITS_BASENAME}" \
  "~/${FOREST_OUTPUT_PREFIX}" \
  "${FOREST_OUTPUT_PREFIX}" \
  "~/${REMOTE_HTML_BASENAME}" \
  "~/${REMOTE_MANIFEST_BASENAME}"

if [[ "${EMIT_LOCAL_SAS_DEBUG}" == "1" || "${LOCAL_SAS_DEBUG_ONLY}" == "1" ]]; then
  local_debug_html_path="$(to_local_sas_path "${HTML_OUT}")"
  local_debug_manifest_path="$(to_local_sas_path "${MANIFEST_OUT}")"
  local_debug_csv_path="$(to_local_sas_path "${CSV_OUT}")"
  local_debug_image_prefix_path="$(to_local_sas_path "${WORKDIR}/${FOREST_OUTPUT_PREFIX}")"
  render_forest_template \
    "${LOCAL_DEBUG_SAS_RENDERED}" \
    "$(to_local_sas_path "${RAND_MACRO_SAS}")" \
    "$(to_local_sas_path "${MKFMT_MACRO_SAS}")" \
    "$(to_local_sas_path "${FOREST_MACRO_SAS}")" \
    "${local_debug_csv_path}" \
    "${local_debug_image_prefix_path}" \
    "${FOREST_OUTPUT_PREFIX}" \
    "${local_debug_html_path}" \
    "${local_debug_manifest_path}"
  echo "[local-sas] Emitted local-SAS debug script: ${LOCAL_DEBUG_SAS_RENDERED}"
fi

if [[ "${LOCAL_SAS_DEBUG_ONLY}" == "1" ]]; then
  echo "[local-sas] LOCAL_SAS_DEBUG_ONLY=1, skipping SAS ODA submit."
  exit 0
fi

rm -f "${HTML_OUT}" "${MANIFEST_OUT}"
for local_png_path in "${WORKDIR}/${FOREST_OUTPUT_PREFIX}"_*.png; do
  [[ -e "${local_png_path}" ]] || continue
  rm -f "${local_png_path}"
done

echo "[1/5] Uploading forest-plot SAS macros to SAS ODA..."
run_oda_helper \
  --upload-file "${RAND_MACRO_SAS}" \
  --upload-file "${MKFMT_MACRO_SAS}" \
  --upload-file "${FOREST_MACRO_SAS}" \
  --output-prefix "upload_forest_macros_${stamp}"

if [[ "${SKIP_DATA_UPLOAD}" == "1" ]]; then
  echo "[2/5] Reusing already-uploaded forest top-hit CSV in SAS ODA: ${REMOTE_TOP_HITS_BASENAME}"
elif [[ "${KEEP_REMOTE_PLOT_DATA}" == "1" ]] && remote_file_exists "${REMOTE_TOP_HITS_BASENAME}" && remote_file_matches_local_size "${CSV_OUT}" "${REMOTE_TOP_HITS_BASENAME}"; then
  echo "[2/5] Keeping and reusing existing remote forest top-hit CSV in SAS ODA: ${REMOTE_TOP_HITS_BASENAME}"
else
  if remote_file_exists "${REMOTE_TOP_HITS_BASENAME}" && ! remote_file_matches_local_size "${CSV_OUT}" "${REMOTE_TOP_HITS_BASENAME}"; then
    echo "[repair] Existing remote forest top-hit CSV has the wrong size and will be replaced."
    delete_remote_file_quiet "${REMOTE_TOP_HITS_BASENAME}"
  fi
  echo "[2/5] Uploading forest top-hit CSV to SAS ODA..."
  run_oda_helper \
    --upload-file "${CSV_OUT}" \
    --output-prefix "upload_forest_top_hits_csv_${stamp}"
fi

echo "[3/5] Running SAS forest plot..."
run_oda_helper \
  --file "${RUN_SAS_RENDERED}" \
  --output-prefix "run_top_hits_forest_plot_${stamp}"

echo "[4/5] Downloading forest HTML, manifest, and panel PNGs..."
oda_download_many \
  "download_top_hits_forest_support_${stamp}" \
  --download-file "~/${REMOTE_HTML_BASENAME}" \
  --download-local-path "${HTML_OUT}" \
  --download-file "~/${REMOTE_MANIFEST_BASENAME}" \
  --download-local-path "${MANIFEST_OUT}" || true

remote_pngs="$(
  run_oda_helper \
    --dir4listing '~' \
    --output-prefix "list_top_hits_forest_png_${stamp}" 2>&1 || true
)"

png_download_args=()
while IFS= read -r remote_png; do
  [[ -z "${remote_png}" ]] && continue
  [[ "${remote_png}" =~ ^${FOREST_OUTPUT_PREFIX}_.+\.png$ ]] || continue
  png_download_args+=(--download-file "~/${remote_png}" --download-local-path "${WORKDIR}/${remote_png}")
done <<< "${remote_pngs}"

if [[ ${#png_download_args[@]} -gt 0 ]]; then
  oda_download_many \
    "download_top_hits_forest_pngs_${stamp}" \
    "${png_download_args[@]}"
fi

if [[ ! -s "${MANIFEST_OUT}" ]]; then
  echo "ERROR: Expected forest manifest was not created or is empty: ${MANIFEST_OUT}" >&2
  exit 1
fi

tmp_manifest="${MANIFEST_OUT}.tmp"
awk -F '\t' '
  NR==1 { print; next }
  $1 ~ /^[[:space:]]*[0-9]+[[:space:]]*$/ && $4 ~ /\.(png|svg)$/ { print }
' "${MANIFEST_OUT}" > "${tmp_manifest}"
mv -f "${tmp_manifest}" "${MANIFEST_OUT}"

forest_png_count=0
while IFS= read -r manifest_png; do
  [[ -z "${manifest_png}" ]] && continue
  if [[ -s "${WORKDIR}/${manifest_png}" ]]; then
    forest_png_count=$((forest_png_count + 1))
  fi
done < <(awk -F '\t' 'NR>1 {print $4}' "${MANIFEST_OUT}")

if (( forest_png_count == 0 )); then
  echo "ERROR: No forest panel PNGs were downloaded for prefix ${FOREST_OUTPUT_PREFIX}." >&2
  exit 1
fi

{
  echo '<!doctype html>'
  echo '<html><head><meta charset="utf-8">'
  printf '<title>%s</title>\n' "${FOREST_HTML_TITLE}"
  echo '<style>'
  echo 'body{margin:0;padding:20px;font-family:Arial,sans-serif;background:#fff;color:#1f2937;}'
  echo '.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:24px;align-items:start;}'
  echo '.panel{margin:0;padding:14px;border:1px solid #d9dee7;border-radius:10px;background:#fff;box-shadow:0 1px 4px rgba(0,0,0,0.06);}'
  echo '.panel figcaption{font-size:18px;font-weight:700;margin-bottom:10px;text-align:center;}'
  echo '.panel img{width:100%;height:auto;display:block;}'
  echo '</style></head><body><div class="grid">'
  awk -F '\t' 'NR>1 {printf "<figure class=\"panel\"><figcaption>%s</figcaption><img src=\"%s\" alt=\"%s forest plot\"></figure>\n", $3, $4, $3}' "${MANIFEST_OUT}"
  echo '</div></body></html>'
} > "${HTML_OUT}"

if [[ ! -s "${HTML_OUT}" ]]; then
  echo "ERROR: Expected forest HTML was not created or is empty: ${HTML_OUT}" >&2
  exit 1
fi

echo "Verified forest HTML:     ${HTML_OUT} ($(wc -c < "${HTML_OUT}") bytes)"
echo "Verified forest manifest: ${MANIFEST_OUT} ($(wc -c < "${MANIFEST_OUT}") bytes)"
echo "Verified forest CSV:      ${CSV_OUT} ($(wc -c < "${CSV_OUT}") bytes)"
echo "Verified forest panels:   ${forest_png_count}"

if [[ "${KEEP_REMOTE_PLOT_DATA}" == "1" ]]; then
  echo "[5/5] Keeping remote forest outputs and uploaded inputs in SAS ODA because KEEP_REMOTE_PLOT_DATA=${KEEP_REMOTE_PLOT_DATA}."
elif [[ "${CLEAN_ODA_INPUT}" == "1" ]]; then
  echo "[5/5] Removing remote forest outputs and uploaded inputs from SAS ODA to save space..."
  cleanup_remote_generated_outputs
  oda_delete_many \
    "cleanup_top_hits_forest_inputs_${stamp}" \
    --delete-file "${REMOTE_TOP_HITS_BASENAME}" \
    --delete-file "RandBetween.sas" \
    --delete-file "mkfmt4grps_by_var.sas" \
    --delete-file "beta2OR_forest_plot.sas" || true
else
  echo "[5/5] Keeping uploaded forest inputs in SAS ODA because CLEAN_ODA_INPUT=${CLEAN_ODA_INPUT}."
  cleanup_remote_generated_outputs
fi

echo "Done."
echo "Downloaded forest HTML: ${HTML_OUT}"
echo "Downloaded forest manifest: ${MANIFEST_OUT}"
echo "Downloaded forest CSV: ${CSV_OUT}"

if [[ "${OPEN_RESULT}" == "1" ]]; then
  open_html_result "${HTML_OUT}" || echo "Forest HTML result is ready: ${HTML_OUT}"
else
  echo "Not opening forest result because OPEN_RESULT=${OPEN_RESULT}."
fi
