#!/usr/bin/env bash
set -euo pipefail

DEPS_DIR="${DEPS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)}"
WORKDIR="${WORKDIR:-$(cd "${DEPS_DIR}/.." && pwd -P)}"
ODA_HELPER_SCRIPT="${ODA_HELPER_SCRIPT:-${DEPS_DIR}/run_sas_codes_or_script_in_ODA.pl}"
if [[ ! -f "${ODA_HELPER_SCRIPT}" ]]; then
  ODA_HELPER_SCRIPT="${WORKDIR}/run_sas_codes_or_script_in_ODA.pl"
fi
if [[ -f "${ODA_HELPER_SCRIPT}" ]]; then
  ODA_PERL_BASE=(perl "${ODA_HELPER_SCRIPT}")
else
  ODA_PERL_BASE=(perl -S run_sas_codes_or_script_in_ODA.pl)
fi
RUNNER_CONFIG_JSON="${RUNNER_CONFIG_JSON:-}"
CALLER_LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES="${LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES-__UNSET__}"
if [[ -n "${RUNNER_CONFIG_JSON}" ]]; then
  eval "$("perl" "${DEPS_DIR}/emit_diff_gwas_runner_env.pl" --config "${RUNNER_CONFIG_JSON}")"
fi
if [[ "${CALLER_LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES}" != "__UNSET__" ]]; then
  LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES="${CALLER_LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES}"
fi
PROJECT_TAG="${PROJECT_TAG:-PGC_SCZ}"
DEFAULT_DATA_GZ="/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects.stdized.wide_beta_se_p_p_lt_0p05.final.tsv.gz"
DATA_GZ="${DATA_GZ:-${DEFAULT_DATA_GZ}}"
DEFAULT_REMOTE_DATA_BASENAME="PGC_SCZ_female_vs_male_diff_effects.stdized.wide_beta_se_p_p_lt_0p05.final.tsv.gz"
REMOTE_DATA_BASENAME="${REMOTE_DATA_BASENAME:-$(basename "${DATA_GZ}")}"
TOP_HIT_FOCUS_PVAR="${TOP_HIT_FOCUS_PVAR:-ASN_STD_P}"
TOP_HIT_FILTER_EXPR="${TOP_HIT_FILTER_EXPR:-((ASN_STD_P>0) and (ASN_STD_P<1E-6)) or ((EUR_STD_P>0) and (EUR_STD_P<1E-6))}"
TOP_HIT_SIGNAL_THRSHD="${TOP_HIT_SIGNAL_THRSHD:-1e-6}"
TOP_HIT_DIST_BP="${TOP_HIT_DIST_BP:-1e6}"
LOCAL_WINDOW_BP="${LOCAL_WINDOW_BP:-1e7}"
TARGET_SNP_LIST="${TARGET_SNP_LIST:-}"
TARGET_SNP_GENES="${TARGET_SNP_GENES:-}"

normalize_reference_build_shell() {
  local value="${1:-hg38}"
  value="$(printf '%s' "${value}" | tr '[:upper:]' '[:lower:]')"
  case "${value}" in
    hg19|grch37|b37|build37|lift37) printf '%s' 'hg19' ;;
    t2t|hs1|chm13|chm13v2|chm13v2.0|t2t-chm13) printf '%s' 't2t' ;;
    *) printf '%s' 'hg38' ;;
  esac
}
REFERENCE_BUILD="$(normalize_reference_build_shell "${REFERENCE_BUILD:-hg38}")"
case "${REFERENCE_BUILD}" in
  hg19)
    DEFAULT_GTF_DSD="FM.GTF_HG19"
    DEFAULT_GTF_LOCAL_DSD="gtf_hg19"
    DEFAULT_GTF_GZ_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh37_mapping/gencode.v49lift37.annotation.gtf.gz"
    ;;
  t2t)
    DEFAULT_GTF_DSD="FM.GTF_T2T"
    DEFAULT_GTF_LOCAL_DSD="gtf_t2t"
    DEFAULT_GTF_GZ_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/genes/hs1.ncbiRefSeq.gtf.gz"
    ;;
  *)
    DEFAULT_GTF_DSD="FM.GTF_HG38"
    DEFAULT_GTF_LOCAL_DSD="gtf_local_hits_manhattan"
    DEFAULT_GTF_GZ_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz"
    ;;
esac

GTF_DSD="${GTF_DSD:-${DEFAULT_GTF_DSD}}"
FM_LIBPATH="${FM_LIBPATH:-/home/cheng.zhong.shan/my_shared_file_links/cheng.zhong.shan/F_vs_M_Covid19_Hosp}"
GTF_LOCAL_DSD="${GTF_LOCAL_DSD:-${DEFAULT_GTF_LOCAL_DSD}}"
GTF_GZ_URL="${GTF_GZ_URL:-${DEFAULT_GTF_GZ_URL}}"
GET_GTF_MACRO_SAS="${GET_GTF_MACRO_SAS:-${DEPS_DIR}/get_genecode_gtf_data.sas}"
LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES="${LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES:-1}"
MACRO_SAS="${DEPS_DIR}/Manhattan4DiffGWASs_png.sas"
TOP_HIT_DIST_MACRO_SAS="${DEPS_DIR}/get_top_signal_within_dist.sas"
RUN_SAS_TEMPLATE="${DEPS_DIR}/run_sas_oda_local_top_hits_manhattan.sas"
SCHEMA_CONFIG_JSON="${SCHEMA_CONFIG_JSON:-${EXTRACTOR_CONFIG_JSON:-${WORKDIR}/configs/preset_pgc_scz_sex_diff.json}}"
SCHEMA_INCLUDE_HELPER="${SCHEMA_INCLUDE_HELPER:-${DEPS_DIR}/generate_sas_wide_import_include.pl}"
RENDER_SAS_HELPER="${RENDER_SAS_HELPER:-${DEPS_DIR}/render_sas_template.pl}"
LOCAL_SAS_DEBUG_EMITTER="${LOCAL_SAS_DEBUG_EMITTER:-${DEPS_DIR}/emit_local_sas_debug_script.pl}"
LOCAL_TOP_HITS_CSV_HELPER="${LOCAL_TOP_HITS_CSV_HELPER:-${DEPS_DIR}/generate_requested_top_hits_csv.pl}"
SINGLE_SNP_WIDE_HELPER="${SINGLE_SNP_WIDE_HELPER:-${DEPS_DIR}/extract_single_snp_wide_diff_gwas.pl}"
SESSION_ID="${SESSION_ID:-mysession}"
USE_PERSISTENT_SESSION="${USE_PERSISTENT_SESSION:-0}"
CLEAN_ODA_INPUT="${CLEAN_ODA_INPUT:-1}"
CLEAN_ODA_OUTPUT="${CLEAN_ODA_OUTPUT:-1}"
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
  DEFAULT_MANHATTAN_FIG_HEIGHT="760"
else
  DEFAULT_MANHATTAN_P_VAR="${DEFAULT_MULTI_P_VAR}"
  DEFAULT_MANHATTAN_OTHER_P_VARS="${DEFAULT_MULTI_OTHER_P_VARS}"
  DEFAULT_MANHATTAN_GWAS_LABEL_NAMES="${DEFAULT_MULTI_GWAS_LABEL_NAMES}"
  DEFAULT_MANHATTAN_FIG_HEIGHT="1000"
fi
MANHATTAN_P_VAR="${MANHATTAN_P_VAR:-${DEFAULT_MANHATTAN_P_VAR}}"
MANHATTAN_OTHER_P_VARS="${MANHATTAN_OTHER_P_VARS:-${DEFAULT_MANHATTAN_OTHER_P_VARS}}"
MANHATTAN_GWAS_LABEL_NAMES="${MANHATTAN_GWAS_LABEL_NAMES:-${DEFAULT_MANHATTAN_GWAS_LABEL_NAMES}}"
DEFAULT_MANHATTAN_FIG_WIDTH="${DEFAULT_MANHATTAN_FIG_WIDTH:-1800}"
MANHATTAN_FIG_WIDTH="${LOCAL_MANHATTAN_FIG_WIDTH:-${MANHATTAN_FIG_WIDTH:-${DEFAULT_MANHATTAN_FIG_WIDTH}}}"
MANHATTAN_FIG_HEIGHT="${LOCAL_MANHATTAN_FIG_HEIGHT:-${MANHATTAN_FIG_HEIGHT:-${DEFAULT_MANHATTAN_FIG_HEIGHT}}}"
LOCAL_MAX_HITS_PER_FIG="${LOCAL_MAX_HITS_PER_FIG:-6}"
LOCAL_MANHATTAN_ANGLE4XAXIS_LABEL="${LOCAL_MANHATTAN_ANGLE4XAXIS_LABEL:-}"
LOCAL_MANHATTAN_XGRP_Y_POS="${LOCAL_MANHATTAN_XGRP_Y_POS:-}"
LOCAL_MANHATTAN_YOFFSET_TOP="${LOCAL_MANHATTAN_YOFFSET_TOP:-}"
LOCAL_MANHATTAN_YOFFSET_BOTTOM="${LOCAL_MANHATTAN_YOFFSET_BOTTOM:-}"
LOCAL_MANHATTAN_FONTSIZE="${LOCAL_MANHATTAN_FONTSIZE:-}"
LOCAL_MANHATTAN_Y_AXIS_LABEL_SIZE="${LOCAL_MANHATTAN_Y_AXIS_LABEL_SIZE:-}"
LOCAL_MANHATTAN_Y_AXIS_VALUE_SIZE="${LOCAL_MANHATTAN_Y_AXIS_VALUE_SIZE:-}"
TOP_HIT_MAF_THRESHOLD="${TOP_HIT_MAF_THRESHOLD:-0.01}"
TOP_HIT_MAX_LOCI="${TOP_HIT_MAX_LOCI:-0}"

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
TOP_HIT_GNOMAD_FREQ_FILE="${TOP_HIT_GNOMAD_FREQ_FILE:-}"
TOP_HIT_GNOMAD_POP_MAP="${TOP_HIT_GNOMAD_POP_MAP:-}"
LOCAL_OUTPUT_PREFIX="${LOCAL_OUTPUT_PREFIX:-${OUTPUT_PREFIX:-${PROJECT_TAG}_SAS_local_top_hits_manhattan}}"
LOCAL_HTML_TITLE="${LOCAL_HTML_TITLE:-${HTML_TITLE:-${PROJECT_TAG} Local Top Hits Manhattan Plot}}"
LOCAL_TOP_HITS_CSV_BASENAME="${LOCAL_TOP_HITS_CSV_BASENAME:-${LOCAL_OUTPUT_PREFIX}_top_hits.csv}"
if [[ -n "${TARGET_SNP_LIST}" ]]; then
  IFS=',' read -r -a _target_snp_csv_names <<< "${TARGET_SNP_LIST}"
  if [[ "${#_target_snp_csv_names[@]}" -eq 1 ]]; then
    _single_target_csv_tag="$(printf '%s' "${_target_snp_csv_names[0]}" | tr -c 'A-Za-z0-9._-' '_')"
  else
    _single_target_csv_tag="targets_$(printf '%s' "${TARGET_SNP_LIST}" | perl -MDigest::MD5=md5_hex -ne 'print substr(md5_hex($_),0,12)')"
  fi
  LOCAL_TOP_HITS_CSV_BASENAME="${LOCAL_TOP_HITS_CSV_BASENAME%.csv}_${_single_target_csv_tag}.csv"
  echo "[prep] TARGET_SNP_LIST is set, so a target-specific local-top-hit CSV will be used: ${LOCAL_TOP_HITS_CSV_BASENAME}"
fi
OPEN_RESULT="${OPEN_RESULT:-1}"
platform_is_linux=0
if [[ "$(uname -s)" == "Linux" ]]; then
  platform_is_linux=1
fi
if [[ -z "${LOCAL_MH_SUBMIT_TIMEOUT_SECONDS:-}" ]]; then
  if [[ "${platform_is_linux}" == "1" ]]; then
    LOCAL_MH_SUBMIT_TIMEOUT_SECONDS=3600
  else
    LOCAL_MH_SUBMIT_TIMEOUT_SECONDS=1200
  fi
fi
if [[ -z "${LOCAL_MH_SUBMIT_TIMEOUT_GRACE_SECONDS:-}" ]]; then
  LOCAL_MH_SUBMIT_TIMEOUT_GRACE_SECONDS=30
fi
if [[ -z "${ODA_HELPER_TIMEOUT_SECONDS:-}" ]]; then
  if [[ "${platform_is_linux}" == "1" ]]; then
    ODA_HELPER_TIMEOUT_SECONDS=1200
  else
    ODA_HELPER_TIMEOUT_SECONDS=300
  fi
fi
if [[ -z "${ODA_HELPER_TIMEOUT_GRACE_SECONDS:-}" ]]; then
  ODA_HELPER_TIMEOUT_GRACE_SECONDS=20
fi
if [[ -z "${ODA_DELETE_TIMEOUT_SECONDS:-}" ]]; then
  if [[ "${platform_is_linux}" == "1" ]]; then
    ODA_DELETE_TIMEOUT_SECONDS=180
  else
    ODA_DELETE_TIMEOUT_SECONDS="${ODA_HELPER_TIMEOUT_SECONDS}"
  fi
fi
if [[ -z "${ODA_DELETE_TIMEOUT_GRACE_SECONDS:-}" ]]; then
  ODA_DELETE_TIMEOUT_GRACE_SECONDS=20
fi
INCLUDE_PREFLIGHT_STANDALONE_DEBUG="${INCLUDE_PREFLIGHT_STANDALONE_DEBUG:-0}"
export SAS_ODA_RUN_TIMEOUT_SECONDS="${SAS_ODA_RUN_TIMEOUT_SECONDS:-${LOCAL_MH_SUBMIT_TIMEOUT_SECONDS}}"
export SAS_ODA_RUN_TIMEOUT_GRACE_SECONDS="${SAS_ODA_RUN_TIMEOUT_GRACE_SECONDS:-${LOCAL_MH_SUBMIT_TIMEOUT_GRACE_SECONDS}}"
export STANDALONE_INCLUDE_TARGET_DEBUG="${INCLUDE_PREFLIGHT_STANDALONE_DEBUG}"
LOCAL_MH_REUSE_CACHE_DIR="${LOCAL_MH_REUSE_CACHE_DIR:-${WORKDIR}/cache/local_manhattan_reuse}"

cd "${WORKDIR}"

if [[ "${USE_PERSISTENT_SESSION}" == "1" ]]; then
  ODA_PERL=("${ODA_PERL_BASE[@]}" --persistent --session-id "${SESSION_ID}")
else
  ODA_PERL=("${ODA_PERL_BASE[@]}")
fi

stamp="$(date +%Y%m%d_%H%M%S)"
PNG_OUT="${WORKDIR}/${LOCAL_OUTPUT_PREFIX}.png"
HTML_OUT="${WORKDIR}/${LOCAL_OUTPUT_PREFIX}.html"
CSV_OUT="${WORKDIR}/${LOCAL_TOP_HITS_CSV_BASENAME}"
RUN_SAS_RENDERED="${WORKDIR}/run_sas_oda_local_top_hits_manhattan.${stamp}.sas"
LOCAL_DEBUG_SAS_RENDERED="${WORKDIR}/run_sas_local_debug_local_top_hits_manhattan.${stamp}.sas"
IMPORT_BLOCK_RENDERED="${WORKDIR}/auto_wide_import_local_hits.${stamp}.sas"
GET_GTF_MACRO_UPLOAD_DIR="${WORKDIR}/.autogen_get_gtf_macro_local_mh_${stamp}"
GET_GTF_MACRO_BASENAME="get_genecode_gtf_data_local_mh_${stamp}.sas"
GET_GTF_MACRO_NAME="gtfdlmh_${stamp}"
GET_GTF_MACRO_UPLOAD="${GET_GTF_MACRO_UPLOAD_DIR}/${GET_GTF_MACRO_BASENAME}"
RUN_PREFIX="run_local_hits_manhattan_png_${stamp}"
RUN_LOG_DIR="${WORKDIR}/${RUN_PREFIX}"
RUN_LOG_FILE="${RUN_LOG_DIR}/output.html.info.txt"
TARGET_SNP_AUG_DIR="${WORKDIR}/.target_snp_wide_aug_local_mh_${stamp}"
TARGET_SNP_AUG_GZ="${WORKDIR}/target_snp_augmented_local_mh_${stamp}.tsv.gz"
TARGET_SNP_AUG_CACHE_MANAGED=0

if [[ ! -f "${ODA_HELPER_SCRIPT}" ]]; then
  echo "ERROR: SAS ODA helper script was not found: ${ODA_HELPER_SCRIPT}" >&2
  exit 1
fi

echo "[helper] Using SAS ODA helper: ${ODA_HELPER_SCRIPT}"
echo "[helper] Include preflight standalone target debug: ${STANDALONE_INCLUDE_TARGET_DEBUG}"

cleanup_generated_artifacts() {
  rm -f "${RUN_SAS_RENDERED}" "${IMPORT_BLOCK_RENDERED}"
  if [[ "${TARGET_SNP_AUG_CACHE_MANAGED}" != "1" ]]; then
    rm -f "${TARGET_SNP_AUG_GZ}"
  fi
  rm -rf "${GET_GTF_MACRO_UPLOAD_DIR}" "${TARGET_SNP_AUG_DIR}"
}

trap cleanup_generated_artifacts EXIT

run_oda_helper() {
  run_oda_helper_with_timeout "${ODA_HELPER_TIMEOUT_SECONDS}" "${ODA_HELPER_TIMEOUT_GRACE_SECONDS}" "$@"
}

run_oda_helper_with_timeout() {
  local timeout_seconds="$1"
  local grace_seconds="$2"
  shift 2
  if [[ -x /usr/bin/timeout && "${timeout_seconds}" =~ ^[0-9]+$ && "${timeout_seconds}" -gt 0 ]]; then
    /usr/bin/timeout --kill-after="${grace_seconds}s" "${timeout_seconds}s" \
      "${ODA_PERL[@]}" "$@"
    return $?
  fi
  "${ODA_PERL[@]}" "$@"
}

oda_upload_many() {
  local output_prefix="$1"
  shift
  run_oda_helper "$@" --output-prefix "${output_prefix}"
}

oda_download_many() {
  local output_prefix="$1"
  shift
  run_oda_helper "$@" --output-prefix "${output_prefix}"
}

oda_delete_many() {
  local output_prefix="$1"
  shift
  run_oda_helper_with_timeout "${ODA_DELETE_TIMEOUT_SECONDS}" "${ODA_DELETE_TIMEOUT_GRACE_SECONDS}" \
    "$@" --output-prefix "${output_prefix}"
}

manhattan_submit_needs_retry() {
  [[ ! -s "${RUN_LOG_FILE}" ]] && return 0
  if grep -Eiq 'We failed in getConnection|The application could not log on to the server|server configuration is invalid|SAS process has terminated unexpectedly|SAS submit timed out' "${RUN_LOG_FILE}"; then
    return 0
  fi
  if grep -Eq 'HTML output saved to:|SAS job is completed!' "${RUN_LOG_FILE}"; then
    return 1
  fi
  local bytes
  bytes="$(wc -c < "${RUN_LOG_FILE}")"
  [[ "${bytes}" -lt 500 ]] && return 0
  return 1
}

run_manhattan_submit() {
  if [[ -x /usr/bin/timeout ]]; then
    /usr/bin/timeout --kill-after="${LOCAL_MH_SUBMIT_TIMEOUT_GRACE_SECONDS}s" "${LOCAL_MH_SUBMIT_TIMEOUT_SECONDS}s" \
      "${ODA_PERL[@]}" \
      --file "${RUN_SAS_RENDERED}" \
      --output-prefix "${RUN_PREFIX}"
    return $?
  fi

  "${ODA_PERL[@]}" \
    --file "${RUN_SAS_RENDERED}" \
    --output-prefix "${RUN_PREFIX}"
}

cleanup_remote_generated_outputs() {
  if [[ "${CLEAN_ODA_OUTPUT}" != "1" ]]; then
    echo "[cleanup] Keeping generated remote Manhattan outputs because CLEAN_ODA_OUTPUT=${CLEAN_ODA_OUTPUT}."
    return 0
  fi

  local remote_png seen
  seen=""
  echo "[cleanup] Removing generated remote Manhattan outputs from SAS ODA..."
  oda_delete_many "cleanup_local_hits_manhattan_output_html_${stamp}" --delete-file "${LOCAL_OUTPUT_PREFIX}.html" || true
  oda_delete_many "cleanup_local_hits_manhattan_output_csv_${stamp}" --delete-file "${LOCAL_TOP_HITS_CSV_BASENAME}" || true

  while IFS= read -r remote_png; do
    [[ -z "${remote_png}" ]] && continue
    [[ "${remote_png}" =~ ^${LOCAL_OUTPUT_PREFIX}(_part[0-9]+)?\.png$ ]] || continue
    case "|${seen}|" in
      *"|${remote_png}|"*) continue ;;
    esac
    seen="${seen}|${remote_png}"
    oda_delete_many "cleanup_local_hits_manhattan_output_png_${stamp}" --delete-file "${remote_png}" || true
  done <<EOF
${remote_pngs:-}
EOF
}

remote_data_exists() {
  local check_output
  check_output="$(
    run_oda_helper \
      --dir4listing '~' \
      --output-prefix "check_remote_local_manhattan_data_${stamp}" 2>&1 || true
  )"
  printf '%s\n' "${check_output}" | grep -Fxq "${REMOTE_DATA_BASENAME}"
}

remote_data_size_bytes() {
  local info_output
  info_output="$(
    run_oda_helper \
      --file-info "~/${REMOTE_DATA_BASENAME}" \
      --output-prefix "check_remote_local_manhattan_size_${stamp}" 2>&1 || true
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

remote_data_known_to_oda() {
  local remote_size
  remote_size="$(remote_data_size_bytes | tr -d '[:space:]')"
  [[ -n "${remote_size}" ]] && return 0
  remote_data_exists
}

remote_home_file_size_bytes() {
  local remote_basename="$1"
  local info_output
  info_output="$(
    run_oda_helper \
      --file-info "~/${remote_basename}" \
      --output-prefix "check_remote_home_file_${stamp}" 2>&1 || true
  )"
  printf '%s\n' "${info_output}" | awk -F '\t' '$1=="SIZE"{print $2}' | tail -n 1
}

remote_home_file_matches_local_size() {
  local local_path="$1"
  local remote_basename="$2"
  [[ -s "${local_path}" ]] || return 1
  local local_size remote_size
  local_size="$(wc -c < "${local_path}" | tr -d '[:space:]')"
  remote_size="$(remote_home_file_size_bytes "${remote_basename}" | tr -d '[:space:]')"
  [[ -n "${local_size}" && -n "${remote_size}" && "${local_size}" == "${remote_size}" ]]
}

upload_home_file_if_needed() {
  local step_label="$1"
  local local_path="$2"
  local remote_basename="$3"
  local output_prefix="$4"
  local upload_path="${local_path}"
  [[ -s "${local_path}" ]] || return 1
  if remote_home_file_matches_local_size "${local_path}" "${remote_basename}"; then
    echo "${step_label} Reusing existing remote file in SAS ODA home: ${remote_basename}"
    return 0
  fi
  echo "${step_label} Uploading $(basename "${local_path}") to SAS ODA home as ${remote_basename}..."
  if [[ "$(basename "${local_path}")" != "${remote_basename}" ]]; then
    mkdir -p "${WORKDIR}/.oda_upload_aliases"
    upload_path="${WORKDIR}/.oda_upload_aliases/${remote_basename}"
    cp -f "${local_path}" "${upload_path}"
  fi
  run_oda_helper \
    --upload-file "${upload_path}" \
    --output-prefix "${output_prefix}"
  if [[ "${upload_path}" != "${local_path}" ]]; then
    rm -f "${upload_path}"
  fi
}

delete_partial_remote_data() {
  echo "[repair] Removing partial or stale remote data file: ${REMOTE_DATA_BASENAME}"
  run_oda_helper \
    --delete-file "${REMOTE_DATA_BASENAME}" \
    --output-prefix "delete_partial_local_hits_subset_${stamp}" >/dev/null 2>&1 || true
}

upload_data_with_integrity_check() {
  local local_size remote_size verify_attempt upload_attempt max_upload_attempts max_verify_attempts
  local_size="$(wc -c < "${DATA_GZ}" | tr -d '[:space:]')"
  max_upload_attempts="${ODA_UPLOAD_MAX_ATTEMPTS:-2}"
  max_verify_attempts="${ODA_UPLOAD_VERIFY_ATTEMPTS:-4}"

  upload_attempt=1
  while [[ "${upload_attempt}" -le "${max_upload_attempts}" ]]; do
    run_oda_helper --upload-file "${DATA_GZ}" --output-prefix "upload_local_hits_subset_${stamp}_try${upload_attempt}"

    verify_attempt=1
    while [[ "${verify_attempt}" -le "${max_verify_attempts}" ]]; do
      remote_size="$(remote_data_size_bytes | tr -d '[:space:]')"
      if [[ -n "${remote_size}" && "${remote_size}" == "${local_size}" ]]; then
        return 0
      fi
      if [[ "${verify_attempt}" -lt "${max_verify_attempts}" ]]; then
        echo "[retry] Remote size check for ${REMOTE_DATA_BASENAME} did not stabilize yet (attempt ${verify_attempt}/${max_verify_attempts}, remote=${remote_size:-missing}). Reconnecting and retrying..."
        sleep "${ODA_UPLOAD_VERIFY_RETRY_SLEEP_SECONDS:-3}"
      fi
      verify_attempt=$((verify_attempt+1))
    done

    if [[ "${upload_attempt}" -lt "${max_upload_attempts}" ]]; then
      echo "[retry] Remote upload integrity check failed for ${REMOTE_DATA_BASENAME} after upload attempt ${upload_attempt}/${max_upload_attempts}. Re-uploading..."
      delete_partial_remote_data
      sleep "${ODA_UPLOAD_RETRY_SLEEP_SECONDS:-5}"
    fi
    upload_attempt=$((upload_attempt+1))
  done

  echo "ERROR: Remote upload size mismatch for ${REMOTE_DATA_BASENAME}. local=${local_size} remote=${remote_size:-missing}" >&2
  delete_partial_remote_data
  return 1
}

generate_requested_top_hits_csv_locally() {
  [[ -x "${LOCAL_TOP_HITS_CSV_HELPER}" || -f "${LOCAL_TOP_HITS_CSV_HELPER}" ]] || return 1
  echo "[prep] Generating MAF-filtered requested local-top-hit CSV locally..."
  local -a cmd=(
    perl "${LOCAL_TOP_HITS_CSV_HELPER}"
    --input "${DATA_GZ}"
    --output "${CSV_OUT}"
    --top-hit-mode "${TOP_HIT_MODE:-differential}"
    --top-hit-focus-pvar "${TOP_HIT_FOCUS_PVAR}"
    --top-hit-signal-thrshd "${TOP_HIT_SIGNAL_THRSHD}"
    --top-hit-signal-thrshds "${TOP_HIT_SIGNAL_THRSHDS:-${TOP_HIT_SIGNAL_THRSHD}}"
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

stable_hash_text() {
  perl -MDigest::MD5=md5_hex -e 'print md5_hex(join("\0", @ARGV));' "$@"
}

augment_data_gz_with_target_snp_windows() {
  [[ -n "${TARGET_SNP_LIST}" ]] || return 0

  if [[ -z "${SOURCE_LONG_GZ:-}" || ! -s "${SOURCE_LONG_GZ}" ]]; then
    echo "WARNING: TARGET_SNP_LIST was provided, but SOURCE_LONG_GZ is unavailable for building a compact local Manhattan subset." >&2
    return 0
  fi
  if [[ -z "${SCHEMA_CONFIG_JSON:-}" || ! -s "${SCHEMA_CONFIG_JSON}" ]]; then
    echo "WARNING: TARGET_SNP_LIST was provided, but SCHEMA_CONFIG_JSON is unavailable for building a compact local Manhattan subset." >&2
    return 0
  fi

  mkdir -p "${LOCAL_MH_REUSE_CACHE_DIR}"
  local target_cache_key target_cache_base
  target_cache_key="$(
    stable_hash_text \
      "${PROJECT_TAG}" \
      "${TARGET_SNP_LIST}" \
      "${LOCAL_WINDOW_BP}" \
      "${SOURCE_LONG_GZ}" \
      "${SCHEMA_CONFIG_JSON}"
  )"
  target_cache_base="${LOCAL_MH_REUSE_CACHE_DIR}/target_snp_augmented_${SAFE_PROJECT_TAG}_${target_cache_key}"
  if [[ -s "${target_cache_base}.tsv.gz" ]]; then
    TARGET_SNP_AUG_GZ="${target_cache_base}.tsv.gz"
    DATA_GZ="${TARGET_SNP_AUG_GZ}"
    REMOTE_DATA_BASENAME="$(basename "${TARGET_SNP_AUG_GZ}")"
    TARGET_SNP_AUG_CACHE_MANAGED=1
    echo "[prep] Reusing cached compact target-SNP local Manhattan subset: ${DATA_GZ}"
    return 0
  fi

  mkdir -p "${TARGET_SNP_AUG_DIR}"
  local snp safe_snp single_out single_manifest
  local -a target_snp_array=()
  local -a extra_args=()
  echo "[prep] Building a compact local Manhattan subset directly from SOURCE_LONG_GZ for the requested target SNP(s)..."
  IFS=',' read -r -a target_snp_array <<< "${TARGET_SNP_LIST}"
  for snp in "${target_snp_array[@]}"; do
    snp="$(printf '%s' "${snp}" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')"
    [[ -n "${snp}" ]] || continue
    safe_snp="$(printf '%s' "${snp}" | tr -c 'A-Za-z0-9._-' '_')"
    single_out="${TARGET_SNP_AUG_DIR}/${safe_snp}.wide.tsv.gz"
    single_manifest="${TARGET_SNP_AUG_DIR}/${safe_snp}.manifest.tsv"
    perl "${SINGLE_SNP_WIDE_HELPER}" \
      --config "${SCHEMA_CONFIG_JSON}" \
      --input "${SOURCE_LONG_GZ}" \
      --target-snp "${snp}" \
      --window-bp "${LOCAL_WINDOW_BP}" \
      --output "${single_out}" \
      --manifest "${single_manifest}"
    if [[ -s "${single_out}" ]]; then
      extra_args+=("${single_out}")
    else
      echo "WARNING: Failed to build a supplemental wide subset for target SNP ${snp}." >&2
    fi
  done

  [[ ${#extra_args[@]} -gt 0 ]] || return 0

  perl -MIO::Uncompress::Gunzip=gunzip,\$GunzipError -MIO::Compress::Gzip=gzip,\$GzipError -e '
    use strict;
    use warnings;
    my ($out_path, @inputs) = @ARGV;
    die "Need at least one input\n" unless @inputs;
    my $out = IO::Compress::Gzip->new($out_path) or die "Cannot write $out_path: $GzipError\n";
    my %seen;
    my @master_header;
    my %master_idx;
    my @datasets;
    for my $path (@inputs) {
      my $fh = IO::Uncompress::Gunzip->new($path) or die "Cannot open $path: $GunzipError\n";
      my $hdr = <$fh>;
      die "Missing header in $path\n" unless defined $hdr;
      chomp $hdr;
      my @h = split /\t/, $hdr, -1;
      for my $col (@h) {
        next if exists $master_idx{$col};
        $master_idx{$col} = scalar @master_header;
        push @master_header, $col;
      }
      push @datasets, [$path, \@h];
    }
    print {$out} join("\t", @master_header), "\n";
    for my $ds (@datasets) {
      my ($path, $header_ref) = @{$ds};
      my $fh = IO::Uncompress::Gunzip->new($path) or die "Cannot reopen $path: $GunzipError\n";
      <$fh>;
      my %idx;
      @idx{@{$header_ref}} = (0 .. $#{$header_ref});
      die "Required base columns missing in $path\n" unless exists $idx{CHR} && exists $idx{BP} && exists $idx{A1} && exists $idx{A2} && exists $idx{SNP};
      while (my $line = <$fh>) {
        chomp $line;
        next unless length $line;
        my @f = split /\t/, $line, -1;
        my $key = join("\t", @f[@idx{qw(CHR BP A1 A2 SNP)}]);
        next if $seen{$key}++;
        my %row;
        @row{@{$header_ref}} = @f;
        print {$out} join("\t", map { defined $row{$_} ? $row{$_} : q{} } @master_header), "\n";
      }
    }
    close $out or die "Failed closing $out_path: $!\n";
  ' "${TARGET_SNP_AUG_GZ}" "${extra_args[@]}"

  if [[ -s "${TARGET_SNP_AUG_GZ}" ]]; then
    cp -f "${TARGET_SNP_AUG_GZ}" "${target_cache_base}.tsv.gz"
    TARGET_SNP_AUG_GZ="${target_cache_base}.tsv.gz"
    TARGET_SNP_AUG_CACHE_MANAGED=1
    DATA_GZ="${TARGET_SNP_AUG_GZ}"
    REMOTE_DATA_BASENAME="$(basename "${TARGET_SNP_AUG_GZ}")"
    echo "[prep] Using compact target-SNP local Manhattan subset: ${DATA_GZ}"
  fi
}

SAFE_PROJECT_TAG="$(printf '%s' "${PROJECT_TAG}" | tr -c 'A-Za-z0-9._-' '_')"
augment_data_gz_with_target_snp_windows

perl "${SCHEMA_INCLUDE_HELPER}" \
  --config "${SCHEMA_CONFIG_JSON}" \
  --dataset scz_mh \
  --source-type gzip \
  --remote-basename "${REMOTE_DATA_BASENAME}" > "${IMPORT_BLOCK_RENDERED}"

perl "${RENDER_SAS_HELPER}" \
  --template "${RUN_SAS_TEMPLATE}" \
  --output "${RUN_SAS_RENDERED}" \
  --replace "TOP_HIT_FOCUS_PVAR=${TOP_HIT_FOCUS_PVAR}" \
  --replace "TOP_HIT_MODE=${TOP_HIT_MODE:-differential}" \
  --replace "TOP_HIT_FILTER_EXPR=${TOP_HIT_FILTER_EXPR}" \
  --replace "TOP_HIT_SIGNAL_THRSHD=${TOP_HIT_SIGNAL_THRSHD}" \
  --replace "TOP_HIT_SIGNAL_THRSHDS=${TOP_HIT_SIGNAL_THRSHDS:-${TOP_HIT_SIGNAL_THRSHD}}" \
  --replace "TOP_HIT_DIST_BP=${TOP_HIT_DIST_BP}" \
  --replace "TARGET_SNP_LIST=${TARGET_SNP_LIST}" \
  --replace "TARGET_SNP_GENES=${TARGET_SNP_GENES}" \
  --replace "COMMON_ASSOC_P_VARS=${COMMON_ASSOC_P_VARS:-}" \
  --replace "LOCAL_WINDOW_BP=${LOCAL_WINDOW_BP}" \
  --replace "GTF_DSD=${GTF_DSD}" \
  --replace "FM_LIBPATH=${FM_LIBPATH}" \
  --replace "GTF_LOCAL_DSD=${GTF_LOCAL_DSD}" \
  --replace "GTF_GZ_URL=${GTF_GZ_URL}" \
  --replace "GTF_INCLUDE_NON_PROTEIN_CODING=${LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES}" \
  --replace "GET_GTF_MACRO_BASENAME=${GET_GTF_MACRO_BASENAME}" \
  --replace "GET_GTF_MACRO_NAME=${GET_GTF_MACRO_NAME}" \
  --replace "MANHATTAN_P_VAR=${MANHATTAN_P_VAR}" \
  --replace "MANHATTAN_OTHER_P_VARS=${MANHATTAN_OTHER_P_VARS}" \
  --replace "MANHATTAN_GWAS_LABEL_NAMES=${MANHATTAN_GWAS_LABEL_NAMES}" \
  --replace "MANHATTAN_FIG_WIDTH=${MANHATTAN_FIG_WIDTH}" \
  --replace "MANHATTAN_FIG_HEIGHT=${MANHATTAN_FIG_HEIGHT}" \
  --replace "LOCAL_MAX_HITS_PER_FIG=${LOCAL_MAX_HITS_PER_FIG}" \
  --replace "LOCAL_MANHATTAN_ANGLE4XAXIS_LABEL=${LOCAL_MANHATTAN_ANGLE4XAXIS_LABEL}" \
  --replace "LOCAL_MANHATTAN_XGRP_Y_POS=${LOCAL_MANHATTAN_XGRP_Y_POS}" \
  --replace "LOCAL_MANHATTAN_YOFFSET_TOP=${LOCAL_MANHATTAN_YOFFSET_TOP}" \
  --replace "LOCAL_MANHATTAN_YOFFSET_BOTTOM=${LOCAL_MANHATTAN_YOFFSET_BOTTOM}" \
  --replace "LOCAL_MANHATTAN_FONTSIZE=${LOCAL_MANHATTAN_FONTSIZE}" \
  --replace "LOCAL_MANHATTAN_Y_AXIS_LABEL_SIZE=${LOCAL_MANHATTAN_Y_AXIS_LABEL_SIZE}" \
  --replace "LOCAL_MANHATTAN_Y_AXIS_VALUE_SIZE=${LOCAL_MANHATTAN_Y_AXIS_VALUE_SIZE}" \
  --replace "LOCAL_TOP_HITS_CSV_BASENAME=${LOCAL_TOP_HITS_CSV_BASENAME}" \
  --replace "OUTPUT_PREFIX=${LOCAL_OUTPUT_PREFIX}" \
  --replace "HTML_TITLE=${LOCAL_HTML_TITLE}" \
  --replace-file "WIDE_IMPORT_BLOCK=${IMPORT_BLOCK_RENDERED}"

rm -f "${PNG_OUT}" "${HTML_OUT}" "${CSV_OUT}"
mkdir -p "${GET_GTF_MACRO_UPLOAD_DIR}"

cat > "${GET_GTF_MACRO_UPLOAD}" <<'EOF'
%macro __GET_GTF_MACRO_NAME__(gtf_gz_url=,outdsd=,region_chrs=,region_starts=,region_ends=);
  %local _basename _home _cache_path _download_rc _n_regions _ri;
  %if %sysevalf(%superq(gtf_gz_url)=,boolean) %then %do;
    %put ERROR: get_genecode_gtf_data requires gtf_gz_url=;
    %abort 255;
  %end;
  %if %sysevalf(%superq(outdsd)=,boolean) %then %let outdsd=gtf_hg38;
  %let _n_regions=%sysfunc(countw(%superq(region_chrs),|));

  %let _basename=%scan(%superq(gtf_gz_url),-1,/);
  %if %sysevalf(%superq(_basename)=,boolean) %then %let _basename=&outdsd..gtf.gz;
  %let _home=%sysfunc(coalescec(%sysget(HOME),%sysfunc(pathname(work))));
  %let _cache_path=&_home/%superq(_basename);

  %if %sysfunc(fileexist("&_cache_path")) %then %do;
    %put NOTE: Reusing cached GTF.gz from SAS ODA home: &_cache_path;
  %end;
  %else %do;
    %put NOTE: No cached GTF.gz found in SAS ODA home. Downloading &_basename from &gtf_gz_url;
    filename _gtfraw "&_cache_path" recfm=n;
    proc http url="&gtf_gz_url" method="GET" out=_gtfraw;
    run;
    %let _download_rc=&SYS_PROCHTTP_STATUS_CODE;
    filename _gtfraw clear;
    %if %sysevalf(%superq(_download_rc)=,boolean) or &_download_rc >= 400 %then %do;
      %put ERROR: Failed to download &gtf_gz_url into &_cache_path (HTTP=&_download_rc).;
      %abort 255;
    %end;
  %end;

  filename _gtfgz zip "&_cache_path" gzip;
  data &outdsd;
    length seqname chr chr_raw source $64 feature type $32 score $16 strand $4 frame $8 attribute $1;
    length gene_id gene_name transcript_id transcript_name gene_type transcript_type gene genesymbol ensembl $256;
    length start end st en bp1 bp2 txStart txEnd protein_coding 8 attribute_raw $32767;
    infile _gtfgz dlm='09'x dsd truncover lrecl=1048576 firstobs=1;
    input seqname :$64.
          source :$64.
          feature :$32.
          start
          end
          score :$16.
          strand :$4.
          frame :$8.
          attribute_raw :$32767.;
    if missing(seqname) then delete;
    if substr(seqname,1,1)='#' then delete;
    if lowcase(feature) ne 'gene' then delete;

    chr_raw=seqname;
    chr=prxchange('s/^chr//i',1,strip(seqname));
    st=start;
    en=end;
    bp1=start;
    bp2=end;
    txStart=start;
    txEnd=end;

    %if &_n_regions > 0 %then %do;
    _keep_region=0;
    %do _ri=1 %to &_n_regions;
      if upcase(compress(strip(chr),'CHR')) = "%upcase(%scan(%superq(region_chrs),&_ri,|))"
         and end >= %scan(%superq(region_starts),&_ri,|)
         and start <= %scan(%superq(region_ends),&_ri,|) then _keep_region=1;
    %end;
    if _keep_region=0 then delete;
    %end;

    gene_id=prxchange('s/.*gene_id "([^"]+)".*/$1/i',1,attribute_raw);
    if gene_id=attribute_raw then gene_id='';
    gene_name=prxchange('s/.*gene_name "([^"]+)".*/$1/i',1,attribute_raw);
    if gene_name=attribute_raw then gene_name='';
    transcript_id='';
    transcript_name='';
    gene_type=prxchange('s/.*gene_type "([^"]+)".*/$1/i',1,attribute_raw);
    if gene_type=attribute_raw then gene_type='';
    transcript_type='';
    gene=coalescec(gene_name,gene_id,transcript_name,transcript_id,feature);
    genesymbol=coalescec(gene_name,gene,transcript_name,gene_id,transcript_id,feature);
    ensembl=coalescec(source,'gencode');
    type='gene';
    attribute='';
    protein_coding=(index(lowcase(coalescec(gene_type,attribute_raw,'')),'protein_coding')>0);
  run;
  filename _gtfgz clear;
%mend;
EOF
perl -0pi -e 's/__GET_GTF_MACRO_NAME__/\Q'"${GET_GTF_MACRO_NAME}"'\E/g' "${GET_GTF_MACRO_UPLOAD}"

if ! generate_requested_top_hits_csv_locally; then
  echo "WARNING: Local MAF-aware top-hit CSV generation did not succeed. The SAS script will fall back to its internal top-hit selection." >&2
fi

if [[ "${EMIT_LOCAL_SAS_DEBUG}" == "1" || "${LOCAL_SAS_DEBUG_ONLY}" == "1" ]]; then
  perl "${LOCAL_SAS_DEBUG_EMITTER}" \
    --mode local_manhattan \
    --input "${RUN_SAS_RENDERED}" \
    --output "${LOCAL_DEBUG_SAS_RENDERED}" \
    --workdir "${WORKDIR}" \
    --deps-dir "${DEPS_DIR}" \
    --data-gz "${DATA_GZ}" \
    --top-hits-csv "${CSV_OUT}" \
    --gtf-macro-upload "${GET_GTF_MACRO_UPLOAD}" \
    --gtf-local-dataset "${GTF_LOCAL_DSD}" \
    --manhattan-macro "${MACRO_SAS}" \
    --output-html-basename "$(basename "${HTML_OUT}")" >/dev/null
  echo "[local-sas] Emitted local-SAS debug script: ${LOCAL_DEBUG_SAS_RENDERED}"
fi

if [[ "${LOCAL_SAS_DEBUG_ONLY}" == "1" ]]; then
  echo "[local-sas] LOCAL_SAS_DEBUG_ONLY=1, skipping SAS ODA submit."
  exit 0
fi

echo "[1/5] Uploading PNG Manhattan macro to SAS ODA..."
# File uploads land in the shared SAS ODA home directory, so they do not need
# to flow through the persistent session server.
oda_upload_many \
  "upload_local_hits_support_${stamp}" \
  --upload-file "${MACRO_SAS}" \
  --upload-file "${TOP_HIT_DIST_MACRO_SAS}" \
  --upload-file "${GET_GTF_MACRO_UPLOAD}"

if [[ -s "${CSV_OUT}" ]]; then
  upload_home_file_if_needed \
    "[1b/5]" \
    "${CSV_OUT}" \
    "${LOCAL_TOP_HITS_CSV_BASENAME}" \
    "upload_local_hits_requested_csv_${stamp}"
fi

if [[ "${SKIP_DATA_UPLOAD}" == "1" ]]; then
  echo "[2/5] Reusing already-uploaded gzipped local-top-hit Manhattan subset in SAS ODA: ${REMOTE_DATA_BASENAME}"
elif [[ "${KEEP_REMOTE_PLOT_DATA}" == "1" ]] && remote_data_matches_local_size; then
  echo "[2/5] Keeping and reusing existing remote local-top-hit Manhattan subset in SAS ODA: ${REMOTE_DATA_BASENAME}"
else
  echo "[2/5] Uploading gzipped local-top-hit Manhattan subset to SAS ODA..."
  if remote_data_known_to_oda && ! remote_data_matches_local_size; then
    echo "[repair] Existing remote file has the wrong size and will be replaced."
    delete_partial_remote_data
  fi
  # Large GWAS subset uploads are more reliable through a one-shot ODA
  # connection than through the persistent session server.
  upload_data_with_integrity_check
fi

echo "[3/5] Running SAS local top-hits Manhattan plot..."
LOCAL_MH_SUBMIT_MAX_ATTEMPTS="${LOCAL_MH_SUBMIT_MAX_ATTEMPTS:-2}"
LOCAL_MH_SUBMIT_RETRY_SLEEP_SECONDS="${LOCAL_MH_SUBMIT_RETRY_SLEEP_SECONDS:-10}"
local_mh_submit_attempt=1
local_mh_submit_log_incomplete_after_success=0
while :; do
  rm -f "${RUN_LOG_FILE}"
  local_mh_submit_rc=0
  if run_manhattan_submit; then
    local_mh_submit_rc=0
  else
    local_mh_submit_rc=$?
    echo "[3b/5] Local Manhattan SAS submit attempt ${local_mh_submit_attempt} exited with status ${local_mh_submit_rc}."
  fi
  if [[ "${local_mh_submit_rc}" -eq 0 ]] && ! manhattan_submit_needs_retry; then
    break
  fi
  if [[ "${local_mh_submit_attempt}" -ge "${LOCAL_MH_SUBMIT_MAX_ATTEMPTS}" ]]; then
    if [[ "${local_mh_submit_rc}" -ne 0 ]]; then
      echo "ERROR: Local Manhattan SAS submit failed with status ${local_mh_submit_rc} after ${local_mh_submit_attempt} attempt(s)." >&2
      exit 1
    else
      echo "WARNING: Local Manhattan SAS submit log still looks incomplete after ${local_mh_submit_attempt} attempt(s): ${RUN_LOG_FILE}" >&2
      echo "WARNING: Proceeding to remote artifact download/recovery because the SAS submit itself exited with status 0." >&2
      local_mh_submit_log_incomplete_after_success=1
      break
    fi
  fi
  next_attempt=$((local_mh_submit_attempt + 1))
  if [[ "${local_mh_submit_rc}" -ne 0 ]]; then
    echo "[3b/5] Local Manhattan SAS submit attempt ${local_mh_submit_attempt} failed or timed out; retrying attempt ${next_attempt}/${LOCAL_MH_SUBMIT_MAX_ATTEMPTS} after ${LOCAL_MH_SUBMIT_RETRY_SLEEP_SECONDS}s..."
  else
    echo "[3b/5] Local Manhattan SAS submit attempt ${local_mh_submit_attempt} looked incomplete; retrying attempt ${next_attempt}/${LOCAL_MH_SUBMIT_MAX_ATTEMPTS} after ${LOCAL_MH_SUBMIT_RETRY_SLEEP_SECONDS}s..."
  fi
  sleep "${LOCAL_MH_SUBMIT_RETRY_SLEEP_SECONDS}"
  local_mh_submit_attempt="${next_attempt}"
done

echo "[4/5] Downloading PNG and small HTML wrapper..."
oda_download_many \
  "download_local_hits_manhattan_support_${stamp}" \
  --download-file "~/${LOCAL_OUTPUT_PREFIX}.html" \
  --download-local-path "${HTML_OUT}" \
  --download-file "~/${LOCAL_TOP_HITS_CSV_BASENAME}" \
  --download-local-path "${CSV_OUT}" || true

remote_pngs="$(
  run_oda_helper \
    --dir4listing '~' \
    --output-prefix "list_local_hits_manhattan_png_${stamp}" 2>&1 || true
)"

downloaded_pngs=()
png_download_args=()
while IFS= read -r remote_png; do
  [[ -z "${remote_png}" ]] && continue
  [[ "${remote_png}" =~ ^${LOCAL_OUTPUT_PREFIX}(_part[0-9]+)?\.png$ ]] || continue
  png_download_args+=(--download-file "~/${remote_png}" --download-local-path "${WORKDIR}/${remote_png}")
done <<< "${remote_pngs}"

if [[ ${#png_download_args[@]} -gt 0 ]]; then
  oda_download_many \
    "download_local_hits_manhattan_pngs_${stamp}" \
    "${png_download_args[@]}"
fi

downloaded_pngs=()
while IFS= read -r remote_png; do
  [[ -z "${remote_png}" ]] && continue
  [[ "${remote_png}" =~ ^${LOCAL_OUTPUT_PREFIX}(_part[0-9]+)?\.png$ ]] || continue
  if [[ -s "${WORKDIR}/${remote_png}" ]]; then
    downloaded_pngs+=("${remote_png}")
  fi
done <<< "${remote_pngs}"

if [[ ! -s "${PNG_OUT}" ]]; then
  if [[ "${local_mh_submit_log_incomplete_after_success}" == "1" ]]; then
    echo "ERROR: The SAS submit exited successfully, but no PNG artifact could be downloaded after log-based recovery fallback." >&2
    echo "ERROR: Check the saved SAS submit log: ${RUN_LOG_FILE}" >&2
  fi
  echo "ERROR: Expected downloaded PNG was not created or is empty: ${PNG_OUT}" >&2
  exit 1
fi

if [[ "${#downloaded_pngs[@]}" -gt 0 ]]; then
  ordered_pngs=()
  if [[ -s "${PNG_OUT}" ]]; then
    ordered_pngs+=("$(basename "${PNG_OUT}")")
  fi
  while IFS= read -r local_png; do
    [[ -z "${local_png}" ]] && continue
    [[ "${local_png}" == "$(basename "${PNG_OUT}")" ]] && continue
    ordered_pngs+=("${local_png}")
  done < <(printf '%s\n' "${downloaded_pngs[@]}" | grep '_part[0-9]\+\.png$' | sort -V)

  {
    echo '<!doctype html>'
    echo '<html><head><meta charset="utf-8">'
    echo "<title>${LOCAL_HTML_TITLE}</title>"
    echo '<style>body{margin:0;padding:16px;font-family:Arial,sans-serif;background:#fff;} img{max-width:100%;height:auto;display:block;}</style>'
    echo '</head><body>'
    part_idx=1
    for local_png in "${ordered_pngs[@]}"; do
      echo "<section style=\"margin-bottom:24px\"><img src=\"${local_png}\" alt=\"${LOCAL_HTML_TITLE} part ${part_idx} of ${#ordered_pngs[@]}\"></section>"
      part_idx=$((part_idx+1))
    done
    echo '</body></html>'
  } > "${HTML_OUT}"
fi

if [[ ! -s "${HTML_OUT}" ]]; then
  if [[ "${local_mh_submit_log_incomplete_after_success}" == "1" ]]; then
    echo "ERROR: The SAS submit exited successfully, but no HTML wrapper could be downloaded or rebuilt after log-based recovery fallback." >&2
    echo "ERROR: Check the saved SAS submit log: ${RUN_LOG_FILE}" >&2
  fi
  echo "ERROR: Expected downloaded or rebuilt HTML was not created or is empty: ${HTML_OUT}" >&2
  exit 1
fi

if [[ ! -s "${CSV_OUT}" ]]; then
  if [[ "${local_mh_submit_log_incomplete_after_success}" == "1" ]]; then
    echo "ERROR: The SAS submit exited successfully, but no top-hit CSV could be downloaded after log-based recovery fallback." >&2
    echo "ERROR: Check the saved SAS submit log: ${RUN_LOG_FILE}" >&2
  fi
  echo "ERROR: Expected downloaded top-hit CSV was not created or is empty: ${CSV_OUT}" >&2
  exit 1
fi

echo "Verified PNG:  ${PNG_OUT} ($(wc -c < "${PNG_OUT}") bytes)"
echo "Verified HTML: ${HTML_OUT} ($(wc -c < "${HTML_OUT}") bytes)"
echo "Verified CSV:  ${CSV_OUT} ($(wc -c < "${CSV_OUT}") bytes)"

cleanup_remote_generated_outputs

if [[ "${KEEP_REMOTE_PLOT_DATA}" == "1" ]]; then
  echo "[5/5] Keeping uploaded gz input in SAS ODA because KEEP_REMOTE_PLOT_DATA=${KEEP_REMOTE_PLOT_DATA}."
elif [[ "${CLEAN_ODA_INPUT}" == "1" ]]; then
  echo "[5/5] Removing uploaded gz input from SAS ODA to save space..."
  oda_delete_many \
    "cleanup_local_hits_manhattan_input_${stamp}" \
    --delete-file "${REMOTE_DATA_BASENAME}"
else
  echo "[5/5] Keeping uploaded gz input in SAS ODA because CLEAN_ODA_INPUT=${CLEAN_ODA_INPUT}."
fi

echo "Done."
echo "Downloaded PNG:  ${PNG_OUT}"
echo "Downloaded HTML: ${HTML_OUT}"
echo "Downloaded CSV:  ${CSV_OUT}"

if [[ "${OPEN_RESULT}" == "1" ]]; then
  open_html_result "${HTML_OUT}" || true
else
  echo "Not opening result because OPEN_RESULT=${OPEN_RESULT}."
fi
