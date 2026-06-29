#!/usr/bin/env bash
set -euo pipefail

DEPS_DIR="${DEPS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)}"
WORKDIR="${WORKDIR:-$(cd "${DEPS_DIR}/.." && pwd -P)}"
RUNNER_CONFIG_JSON="${RUNNER_CONFIG_JSON:-}"
CALLER_LOCAL_GTF_WINDOW_BP="${LOCAL_GTF_WINDOW_BP-__UNSET__}"
CALLER_LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES="${LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES-__UNSET__}"
if [[ -n "${RUNNER_CONFIG_JSON}" ]]; then
  eval "$("perl" "${DEPS_DIR}/emit_diff_gwas_runner_env.pl" --config "${RUNNER_CONFIG_JSON}")"
fi
if [[ "${CALLER_LOCAL_GTF_WINDOW_BP}" != "__UNSET__" ]]; then
  LOCAL_GTF_WINDOW_BP="${CALLER_LOCAL_GTF_WINDOW_BP}"
fi
if [[ "${CALLER_LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES}" != "__UNSET__" ]]; then
  LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES="${CALLER_LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES}"
fi
PROJECT_TAG="${PROJECT_TAG:-PGC_SCZ}"
DEFAULT_DATA_GZ="/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects.stdized.wide_beta_se_p_p_lt_0p05.final.tsv.gz"
DATA_GZ="${DATA_GZ:-${DEFAULT_DATA_GZ}}"
DEFAULT_REMOTE_DATA_BASENAME="PGC_SCZ_female_vs_male_diff_effects.stdized.wide_beta_se_p_p_lt_0p05.final.tsv.gz"
REMOTE_DATA_BASENAME="${REMOTE_DATA_BASENAME:-$(basename "${DATA_GZ}")}"
SCHEMA_CONFIG_JSON="${SCHEMA_CONFIG_JSON:-${EXTRACTOR_CONFIG_JSON:-${WORKDIR}/configs/preset_pgc_scz_sex_diff.json}}"
SCHEMA_INCLUDE_HELPER="${SCHEMA_INCLUDE_HELPER:-${DEPS_DIR}/generate_sas_wide_import_include.pl}"
GTF_IMPORT_INCLUDE_HELPER="${GTF_IMPORT_INCLUDE_HELPER:-${DEPS_DIR}/generate_sas_gtf_import_include.pl}"
GTF_SUBSET_HELPER="${GTF_SUBSET_HELPER:-${DEPS_DIR}/extract_gencode_gtf_subset.pl}"
SINGLE_SNP_WIDE_HELPER="${SINGLE_SNP_WIDE_HELPER:-${DEPS_DIR}/extract_single_snp_wide_diff_gwas.pl}"
RENDER_SAS_HELPER="${RENDER_SAS_HELPER:-${DEPS_DIR}/render_sas_template.pl}"
LOCAL_SAS_DEBUG_EMITTER="${LOCAL_SAS_DEBUG_EMITTER:-${DEPS_DIR}/emit_local_sas_debug_script.pl}"
LOCAL_TOP_HITS_CSV_HELPER="${LOCAL_TOP_HITS_CSV_HELPER:-${DEPS_DIR}/generate_requested_top_hits_csv.pl}"

TOP_HIT_FOCUS_PVAR="${TOP_HIT_FOCUS_PVAR:-ASN_STD_P}"
TOP_HIT_FILTER_EXPR="${TOP_HIT_FILTER_EXPR:-((ASN_STD_P>0) and (ASN_STD_P<1E-6)) or ((EUR_STD_P>0) and (EUR_STD_P<1E-6))}"
TOP_HIT_SIGNAL_THRSHD="${TOP_HIT_SIGNAL_THRSHD:-1e-6}"
TOP_HIT_DIST_BP="${TOP_HIT_DIST_BP:-1e6}"
TOP_HIT_MAF_THRESHOLD="${TOP_HIT_MAF_THRESHOLD:-0.01}"
TOP_HIT_GNOMAD_FREQ_FILE="${TOP_HIT_GNOMAD_FREQ_FILE:-}"
TOP_HIT_GNOMAD_POP_MAP="${TOP_HIT_GNOMAD_POP_MAP:-}"
LOCAL_WINDOW_BP="${LOCAL_WINDOW_BP:-1e7}"
LOCAL_GTF_WINDOW_BP="${LOCAL_GTF_WINDOW_BP:-${LOCAL_WINDOW_BP}}"
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
    DEFAULT_GTF_LOCAL_DSD="gtf_hg38"
    DEFAULT_GTF_GZ_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz"
    ;;
esac

GTF_DSD="${GTF_DSD:-${DEFAULT_GTF_DSD}}"
FM_LIBPATH="${FM_LIBPATH:-/home/cheng.zhong.shan/my_shared_file_links/cheng.zhong.shan/F_vs_M_Covid19_Hosp}"
GTF_LOCAL_DSD="${GTF_LOCAL_DSD:-${DEFAULT_GTF_LOCAL_DSD}}"
GTF_GZ_URL="${GTF_GZ_URL:-${DEFAULT_GTF_GZ_URL}}"
GET_GTF_MACRO_SAS="${GET_GTF_MACRO_SAS:-${DEPS_DIR}/get_genecode_gtf_data.sas}"
SNP_LOCAL_MACRO_SAS="${SNP_LOCAL_MACRO_SAS:-${DEPS_DIR}/SNP_Local_Manhattan_With_GTF.sas}"
PATCHED_LATTICE_MACRO_SAS="${PATCHED_LATTICE_MACRO_SAS:-${DEPS_DIR}/Lattice_gscatter_over_bed_track.sas}"
MAP_GRP_ASSOC_MACRO_SAS="${MAP_GRP_ASSOC_MACRO_SAS:-${DEPS_DIR}/map_grp_assoc2gene4covidsexgwas.sas}"
MULT_GSCATTER_GENE_MACRO_SAS="${MULT_GSCATTER_GENE_MACRO_SAS:-${DEPS_DIR}/Multgscatter_with_gene_exons.sas}"
ADJ_CLOSE_GENE_GRP_MACRO_SAS="${ADJ_CLOSE_GENE_GRP_MACRO_SAS:-${DEPS_DIR}/adj_grpnum4close_gene_bed_regs.sas}"
TOP_HIT_DIST_MACRO_SAS="${TOP_HIT_DIST_MACRO_SAS:-${DEPS_DIR}/get_top_signal_within_dist.sas}"
GTF_CACHE_DIR="${GTF_CACHE_DIR:-${WORKDIR}/cache/gtf}"
LOCAL_GTF_REUSE_CACHE_DIR="${LOCAL_GTF_REUSE_CACHE_DIR:-${WORKDIR}/cache/local_gtf_reuse}"
if [[ -z "${GTF_ASSOC_PVARS:-}" || -z "${GTF_ZSCORE_VARS:-}" || -z "${GTF_LABELS:-}" ]]; then
  case "${TOP_HIT_FOCUS_PVAR}" in
    ASN_*)
      DEFAULT_GTF_ASSOC_PVARS="ASN_FEMALE_P ASN_MALE_P ASN_DIFF_P"
      DEFAULT_GTF_ZSCORE_VARS="ASN_FEMALE_Z ASN_MALE_Z ASN_DIFF_Z"
      DEFAULT_GTF_LABELS="Asian_Female Asian_Male Asian_Diff"
      ;;
    EUR_*)
      DEFAULT_GTF_ASSOC_PVARS="EUR_FEMALE_P EUR_MALE_P EUR_DIFF_P"
      DEFAULT_GTF_ZSCORE_VARS="EUR_FEMALE_Z EUR_MALE_Z EUR_DIFF_Z"
      DEFAULT_GTF_LABELS="European_Female European_Male European_Diff"
      ;;
    *)
      DEFAULT_GTF_ASSOC_PVARS="ALL_FEMALE_P ALL_MALE_P ALL_DIFF_P"
      DEFAULT_GTF_ZSCORE_VARS="ALL_FEMALE_Z ALL_MALE_Z ALL_DIFF_Z"
      DEFAULT_GTF_LABELS="All_Female All_Male All_Diff"
      ;;
  esac
fi
GTF_ASSOC_PVARS="${GTF_ASSOC_PVARS:-${DEFAULT_GTF_ASSOC_PVARS}}"
GTF_ZSCORE_VARS="${GTF_ZSCORE_VARS:-${DEFAULT_GTF_ZSCORE_VARS}}"
GTF_LABELS="${GTF_LABELS:-${DEFAULT_GTF_LABELS}}"

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
infer_effect_metric_label_from_shell_words() {
  local text="${1:-}"
  local token u seen=0 all_z=1 all_beta=1 all_or=1
  for token in ${text}; do
    seen=1
    u="$(printf '%s' "${token}" | tr '[:lower:]' '[:upper:]')"
    [[ "${u}" =~ (^|_)Z(SCORE)?($|_) ]] || all_z=0
    [[ "${u}" == *BETA* ]] || all_beta=0
    [[ "${u}" =~ (^|_)OR($|_) || "${u}" == *ODDSRATIO* ]] || all_or=0
  done
  if [[ "${seen}" -eq 0 ]]; then
    printf '%s' 'Effect metric'
  elif [[ "${all_z}" -eq 1 ]]; then
    printf '%s' 'Z score'
  elif [[ "${all_beta}" -eq 1 ]]; then
    printf '%s' 'Beta'
  elif [[ "${all_or}" -eq 1 ]]; then
    printf '%s' 'Odds ratio'
  else
    printf '%s' 'Effect metric'
  fi
}
GTF_YAXIS_LABEL="${GTF_YAXIS_LABEL:--log10(P)}"
GTF_COLORBAR_LABEL="${GTF_COLORBAR_LABEL:-$(infer_effect_metric_label_from_shell_words "${GTF_ZSCORE_VARS}")}"
VISIBLE_YLABEL_ENABLED="${VISIBLE_YLABEL_ENABLED:-0}"
VISIBLE_YLABEL_LEFT_PAD_PX="${VISIBLE_YLABEL_LEFT_PAD_PX:-110}"
VISIBLE_YLABEL_POINTSIZE="${VISIBLE_YLABEL_POINTSIZE:-30}"
VISIBLE_YLABEL_FILL="${VISIBLE_YLABEL_FILL:-#20304f}"
VISIBLE_COLORBAR_LABEL_ENABLED="${VISIBLE_COLORBAR_LABEL_ENABLED:-0}"
VISIBLE_COLORBAR_LABEL_RIGHT_PAD_PX="${VISIBLE_COLORBAR_LABEL_RIGHT_PAD_PX:-90}"
VISIBLE_COLORBAR_LABEL_POINTSIZE="${VISIBLE_COLORBAR_LABEL_POINTSIZE:-24}"
VISIBLE_COLORBAR_LABEL_FILL="${VISIBLE_COLORBAR_LABEL_FILL:-#20304f}"
MAGICK_BIN="${MAGICK_BIN:-}"
# For local GTF plots, the requested local GTF window should define the actual
# displayed genomic half-window in the final figure.
GTF_DIST2SNP="${LOCAL_GTF_WINDOW_BP}"
GTF_DESIGN_WIDTH="${GTF_DESIGN_WIDTH:-950}"
count_shell_words() {
  local text="${1:-}"
  local count=0
  local token
  for token in ${text}; do
    count=$((count + 1))
  done
  printf '%s' "${count}"
}
GTF_SCATTER_TRACK_COUNT="$(count_shell_words "${GTF_ASSOC_PVARS}")"
if [[ -z "${GTF_SCATTER_TRACK_COUNT}" || "${GTF_SCATTER_TRACK_COUNT}" -le 0 ]]; then
  GTF_SCATTER_TRACK_COUNT=3
fi
if [[ -n "${GTF_DESIGN_HEIGHT:-}" ]]; then
  GTF_DESIGN_HEIGHT="${GTF_DESIGN_HEIGHT}"
else
  GTF_DESIGN_HEIGHT="$((700 + (45 * GTF_SCATTER_TRACK_COUNT)))"
  if (( GTF_DESIGN_HEIGHT < 1000 )); then
    GTF_DESIGN_HEIGHT=1000
  fi
  if (( GTF_DESIGN_HEIGHT > 1400 )); then
    GTF_DESIGN_HEIGHT=1400
  fi
fi
GTF_DIST2SEP_GENES="${GTF_DIST2SEP_GENES:-100000}"
GTF_SHIFT_TEXT_YVAL="${GTF_SHIFT_TEXT_YVAL:-0.2}"
# Give the bottom gene track a bit more vertical share by default so the
# manuscript local-GTF panels keep nearby genes readable without requiring
# per-run tuning.
GTF_PCT4NEG_Y="${GTF_PCT4NEG_Y:-1.4}"
GTF_ADJVAL4HEADER="${GTF_ADJVAL4HEADER:--0.6}"
GTF_YAXIS_OFFSET4MAX="${GTF_YAXIS_OFFSET4MAX:-}"
GTF_YOFFSET4TEXTLABELS="${GTF_YOFFSET4TEXTLABELS:-2.5}"
GTF_YOFFSET4MAX_DRAWMARKERSONTOP="${GTF_YOFFSET4MAX_DRAWMARKERSONTOP:-0.25}"
GTF_LABEL_SNPS="${GTF_LABEL_SNPS:-}"
GTF_LABEL_LAYOUT="${GTF_LABEL_LAYOUT:-auto}"
LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES="${LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES:-0}"
LOCAL_GTF_MAX_HITS_PER_FIG="${LOCAL_GTF_MAX_HITS_PER_FIG:-${LOCAL_MAX_HITS_PER_FIG:-1}}"
LOCAL_MAX_HITS_PER_FIG="${LOCAL_GTF_MAX_HITS_PER_FIG}"
TOP_HIT_MAX_LOCI="${TOP_HIT_MAX_LOCI:-0}"
LOCAL_TOP_HITS_CSV_BASENAME="${LOCAL_TOP_HITS_CSV_BASENAME:-${PROJECT_TAG}_SAS_local_top_hits_manhattan_top_hits.csv}"
LOCAL_TOP_HITS_INPUT_CSV_BASENAME="${LOCAL_TOP_HITS_INPUT_CSV_BASENAME:-}"
if [[ -n "${TARGET_SNP_LIST}" && -z "${LOCAL_TOP_HITS_INPUT_CSV_BASENAME}" ]]; then
  IFS=',' read -r -a _target_snp_csv_names <<< "${TARGET_SNP_LIST}"
  if [[ "${#_target_snp_csv_names[@]}" -eq 1 ]]; then
    _single_target_csv_tag="$(printf '%s' "${_target_snp_csv_names[0]}" | tr -c 'A-Za-z0-9._-' '_')"
  else
    _single_target_csv_tag="targets_$(printf '%s' "${TARGET_SNP_LIST}" | perl -MDigest::MD5=md5_hex -ne 'print substr(md5_hex($_),0,12)')"
  fi
  LOCAL_TOP_HITS_CSV_BASENAME="${LOCAL_TOP_HITS_CSV_BASENAME%.csv}_${_single_target_csv_tag}.csv"
  LOCAL_TOP_HITS_INPUT_CSV_BASENAME="${LOCAL_TOP_HITS_CSV_BASENAME}"
  echo "[prep] TARGET_SNP_LIST is set, so a target-specific local-top-hit CSV will be used: ${LOCAL_TOP_HITS_CSV_BASENAME}"
fi

RUN_SAS_TEMPLATE="${DEPS_DIR}/run_sas_oda_local_top_hits_with_gtf.sas"
SESSION_ID="${SESSION_ID:-mysession}"
USE_PERSISTENT_SESSION="${USE_PERSISTENT_SESSION:-0}"
CLEAN_ODA_INPUT="${CLEAN_ODA_INPUT:-1}"
CLEAN_ODA_OUTPUT="${CLEAN_ODA_OUTPUT:-1}"
SKIP_DATA_UPLOAD="${SKIP_DATA_UPLOAD:-0}"
KEEP_REMOTE_PLOT_DATA="${KEEP_REMOTE_PLOT_DATA:-0}"
OPEN_RESULT="${OPEN_RESULT:-1}"
OUTPUT_HTML_BASENAME="${OUTPUT_HTML_BASENAME:-${PROJECT_TAG}_SAS_local_top_hits_with_gtf.html}"
EMIT_LOCAL_SAS_DEBUG="${EMIT_LOCAL_SAS_DEBUG:-0}"
LOCAL_SAS_DEBUG_ONLY="${LOCAL_SAS_DEBUG_ONLY:-0}"
ASSUME_REMOTE_GTF_SUPPORT_READY="${ASSUME_REMOTE_GTF_SUPPORT_READY:-0}"
ASSUME_REMOTE_GTF_DYNAMIC_INPUTS_READY="${ASSUME_REMOTE_GTF_DYNAMIC_INPUTS_READY:-0}"

cd "${WORKDIR}"

ODA_HELPER_SCRIPT="${ODA_HELPER_SCRIPT:-${DEPS_DIR}/run_sas_codes_or_script_in_ODA.pl}"
if [[ ! -f "${ODA_HELPER_SCRIPT}" ]]; then
  ODA_HELPER_SCRIPT="${WORKDIR}/run_sas_codes_or_script_in_ODA.pl"
fi
if [[ -f "${ODA_HELPER_SCRIPT}" ]]; then
  ODA_PERL_BASE=(perl "${ODA_HELPER_SCRIPT}")
else
  ODA_PERL_BASE=(perl -S run_sas_codes_or_script_in_ODA.pl)
fi
# Default to one-shot SAS ODA submits/uploads because that path has proven
# much faster and more reliable than the local persistent-session relay for
# this pipeline. Users can still opt back into session reuse explicitly.
if [[ "${USE_PERSISTENT_SESSION}" == "1" ]]; then
  ODA_PERL=("${ODA_PERL_BASE[@]}" --persistent --session-id "${SESSION_ID}")
else
  ODA_PERL=("${ODA_PERL_BASE[@]}")
fi

stamp="$(date +%Y%m%d_%H%M%S)"
HTML_OUT="${WORKDIR}/${OUTPUT_HTML_BASENAME}"
CSV_OUT="${WORKDIR}/${LOCAL_TOP_HITS_CSV_BASENAME}"
RAW_HTML_OUT="${HTML_OUT%.html}.sasraw.html"
PNG_OUT="${HTML_OUT%.html}.png"
VERIFY_TOP_HITS_TSV=""
if [[ "${TOP_HIT_MODE:-differential}" == "common_association" ]]; then
  data_base="$(basename "${DATA_GZ}")"
  data_prefix="${data_base%.stdized.wide_beta_se_p_p_lt_0p05.final.tsv.gz}"
  VERIFY_TOP_HITS_TSV="$(dirname "${DATA_GZ}")/${data_prefix}.common_assoc_verify.tsv"
fi
RUN_SAS_RENDERED="${WORKDIR}/run_sas_oda_local_top_hits_with_gtf.${stamp}.sas"
LOCAL_DEBUG_SAS_RENDERED="${WORKDIR}/run_sas_local_debug_local_top_hits_with_gtf.${stamp}.sas"
IMPORT_BLOCK_RENDERED="${WORKDIR}/auto_wide_import_local_hits_with_gtf.${stamp}.sas"
GTF_IMPORT_BLOCK_RENDERED="${WORKDIR}/auto_gtf_import_local_hits_with_gtf.${stamp}.sas"
SAFE_PROJECT_TAG="$(printf '%s' "${PROJECT_TAG}" | tr -c 'A-Za-z0-9._-' '_')"
SAFE_LOCAL_GTF_WINDOW_BP="$(printf '%s' "${LOCAL_GTF_WINDOW_BP}" | tr -c 'A-Za-z0-9._-' '_')"
LOCAL_GTF_SUBSET="${WORKDIR}/local_gtf_subset_local_hits_${stamp}.tsv"
LOCAL_GTF_SUBSET_GZ="${LOCAL_GTF_SUBSET}.gz"
REMOTE_GTF_BASENAME="local_gtf_subset_${SAFE_PROJECT_TAG}_window_${SAFE_LOCAL_GTF_WINDOW_BP}_npc${LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES}.tsv.gz"
GET_GTF_MACRO_UPLOAD_DIR="${WORKDIR}/.autogen_get_gtf_macro_${stamp}"
GET_GTF_MACRO_BASENAME="get_genecode_gtf_data_local_top_hits.sas"
GET_GTF_MACRO_NAME="get_genecode_gtf_data"
GET_GTF_MACRO_UPLOAD="${GET_GTF_MACRO_UPLOAD_DIR}/${GET_GTF_MACRO_BASENAME}"
RUN_PREFIX="run_local_hits_with_gtf_${stamp}"
RUN_LOG_DIR="${WORKDIR}/${RUN_PREFIX}"
RUN_LOG_FILE="${RUN_LOG_DIR}/output.html.info.txt"
TARGET_SNP_AUG_DIR="${WORKDIR}/.target_snp_wide_aug_${stamp}"
TARGET_SNP_AUG_GZ="${WORKDIR}/target_snp_augmented_local_gtf_${stamp}.tsv.gz"
platform_is_linux=0
if [[ "$(uname -s)" == "Linux" ]]; then
  platform_is_linux=1
fi
if [[ -z "${GTF_SUBMIT_TIMEOUT_SECONDS:-}" ]]; then
  if [[ "${platform_is_linux}" == "1" ]]; then
    GTF_SUBMIT_TIMEOUT_SECONDS=3600
  else
    GTF_SUBMIT_TIMEOUT_SECONDS=1200
  fi
fi
if [[ -z "${GTF_SUBMIT_TIMEOUT_GRACE_SECONDS:-}" ]]; then
  GTF_SUBMIT_TIMEOUT_GRACE_SECONDS=30
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
if [[ -z "${ODA_PNG_DOWNLOAD_TIMEOUT_SECONDS:-}" ]]; then
  if [[ "${platform_is_linux}" == "1" ]]; then
    ODA_PNG_DOWNLOAD_TIMEOUT_SECONDS=180
  else
    ODA_PNG_DOWNLOAD_TIMEOUT_SECONDS="${ODA_HELPER_TIMEOUT_SECONDS}"
  fi
fi
if [[ -z "${ODA_PNG_DOWNLOAD_TIMEOUT_GRACE_SECONDS:-}" ]]; then
  ODA_PNG_DOWNLOAD_TIMEOUT_GRACE_SECONDS=20
fi
if [[ -z "${ODA_RESULT_DOWNLOAD_TIMEOUT_SECONDS:-}" ]]; then
  if [[ "${platform_is_linux}" == "1" ]]; then
    ODA_RESULT_DOWNLOAD_TIMEOUT_SECONDS=180
  else
    ODA_RESULT_DOWNLOAD_TIMEOUT_SECONDS="${ODA_HELPER_TIMEOUT_SECONDS}"
  fi
fi
if [[ -z "${ODA_RESULT_DOWNLOAD_TIMEOUT_GRACE_SECONDS:-}" ]]; then
  ODA_RESULT_DOWNLOAD_TIMEOUT_GRACE_SECONDS=20
fi
if [[ -z "${ODA_REMOTE_INFO_TIMEOUT_SECONDS:-}" ]]; then
  if [[ "${platform_is_linux}" == "1" ]]; then
    ODA_REMOTE_INFO_TIMEOUT_SECONDS=90
  else
    ODA_REMOTE_INFO_TIMEOUT_SECONDS="${ODA_HELPER_TIMEOUT_SECONDS}"
  fi
fi
if [[ -z "${ODA_REMOTE_INFO_TIMEOUT_GRACE_SECONDS:-}" ]]; then
  ODA_REMOTE_INFO_TIMEOUT_GRACE_SECONDS=15
fi
if [[ -z "${FORCE_DYNAMIC_GTF_SUPPORT_UPLOADS_ON_LINUX:-}" ]]; then
  if [[ "${platform_is_linux}" == "1" ]]; then
    FORCE_DYNAMIC_GTF_SUPPORT_UPLOADS_ON_LINUX=1
  else
    FORCE_DYNAMIC_GTF_SUPPORT_UPLOADS_ON_LINUX=0
  fi
fi
INCLUDE_PREFLIGHT_STANDALONE_DEBUG="${INCLUDE_PREFLIGHT_STANDALONE_DEBUG:-0}"
INCLUDE_PREFLIGHT_REFRESH_REMOTE="${INCLUDE_PREFLIGHT_REFRESH_REMOTE:-0}"
BATCH_INCLUDE_PREFLIGHT_REFRESH_REMOTE="${BATCH_INCLUDE_PREFLIGHT_REFRESH_REMOTE:-0}"
BATCH1_INCLUDE_PREFLIGHT_REFRESH_REMOTE="${BATCH1_INCLUDE_PREFLIGHT_REFRESH_REMOTE:-1}"
BATCH_INCLUDE_PREFLIGHT_ENABLED="${BATCH_INCLUDE_PREFLIGHT_ENABLED:-0}"
BATCH1_INCLUDE_PREFLIGHT_ENABLED="${BATCH1_INCLUDE_PREFLIGHT_ENABLED:-1}"
KEEP_RENDERED_DEBUG_FILES="${KEEP_RENDERED_DEBUG_FILES:-0}"
export SAS_ODA_RUN_TIMEOUT_SECONDS="${SAS_ODA_RUN_TIMEOUT_SECONDS:-${GTF_SUBMIT_TIMEOUT_SECONDS}}"
export SAS_ODA_RUN_TIMEOUT_GRACE_SECONDS="${SAS_ODA_RUN_TIMEOUT_GRACE_SECONDS:-${GTF_SUBMIT_TIMEOUT_GRACE_SECONDS}}"
export STANDALONE_INCLUDE_TARGET_DEBUG="${INCLUDE_PREFLIGHT_STANDALONE_DEBUG}"
export INCLUDE_PREFLIGHT_REFRESH_REMOTE="${INCLUDE_PREFLIGHT_REFRESH_REMOTE}"

mkdir -p "${LOCAL_GTF_REUSE_CACHE_DIR}"

LOCAL_GTF_SUBSET_CACHE_MANAGED=0
TARGET_SNP_AUG_CACHE_MANAGED=0

cleanup_generated_artifacts() {
  if [[ "${KEEP_RENDERED_DEBUG_FILES}" != "1" ]]; then
    rm -f "${RUN_SAS_RENDERED}" "${IMPORT_BLOCK_RENDERED}" "${GTF_IMPORT_BLOCK_RENDERED}"
    if [[ "${LOCAL_GTF_SUBSET_CACHE_MANAGED}" != "1" ]]; then
      rm -f "${LOCAL_GTF_SUBSET}" "${LOCAL_GTF_SUBSET_GZ}"
    fi
    if [[ "${TARGET_SNP_AUG_CACHE_MANAGED}" != "1" ]]; then
      rm -f "${TARGET_SNP_AUG_GZ}"
    fi
    rm -rf "${GET_GTF_MACRO_UPLOAD_DIR}" "${WORKDIR}/.oda_upload_aliases" "${TARGET_SNP_AUG_DIR}"
  else
    rm -rf "${WORKDIR}/.oda_upload_aliases"
  fi
}

trap cleanup_generated_artifacts EXIT

echo "[helper] Include preflight standalone target debug: ${STANDALONE_INCLUDE_TARGET_DEBUG}"
echo "[helper] Include preflight remote refresh: ${INCLUDE_PREFLIGHT_REFRESH_REMOTE}"
echo "[helper] Batch-1 include preflight enabled: ${BATCH1_INCLUDE_PREFLIGHT_ENABLED}"
echo "[helper] Batch include preflight enabled: ${BATCH_INCLUDE_PREFLIGHT_ENABLED}"
echo "[helper] Batch-1 include preflight remote refresh: ${BATCH1_INCLUDE_PREFLIGHT_REFRESH_REMOTE}"
echo "[helper] Batch include preflight remote refresh: ${BATCH_INCLUDE_PREFLIGHT_REFRESH_REMOTE}"

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

oda_download_many_with_timeout() {
  local timeout_seconds="$1"
  local grace_seconds="$2"
  local output_prefix="$3"
  shift 3
  run_oda_helper_with_timeout "${timeout_seconds}" "${grace_seconds}" "$@" --output-prefix "${output_prefix}"
}

oda_delete_many() {
  local output_prefix="$1"
  shift
  run_oda_helper_with_timeout "${ODA_DELETE_TIMEOUT_SECONDS}" "${ODA_DELETE_TIMEOUT_GRACE_SECONDS}" \
    "$@" --output-prefix "${output_prefix}"
}

cleanup_remote_generated_outputs() {
  if [[ "${CLEAN_ODA_OUTPUT}" != "1" ]]; then
    echo "[cleanup] Keeping generated remote local-GTF outputs because CLEAN_ODA_OUTPUT=${CLEAN_ODA_OUTPUT}."
    return 0
  fi

  local remote_listing remote_name remote_png_base
  remote_listing="$(
    run_oda_helper \
      --dir4listing '~' \
      --output-prefix "list_local_hits_with_gtf_outputs_${stamp}" 2>&1 || true
  )"

  echo "[cleanup] Removing generated remote local-GTF outputs from SAS ODA..."
  oda_delete_many "cleanup_local_hits_with_gtf_output_main_html_${stamp}" --delete-file "${OUTPUT_HTML_BASENAME}" || true
  if [[ -n "${LOCAL_TOP_HITS_INPUT_CSV_BASENAME:-}" && "${LOCAL_TOP_HITS_INPUT_CSV_BASENAME}" == "${LOCAL_TOP_HITS_CSV_BASENAME}" ]]; then
    echo "[cleanup] Preserving remote requested top-hit CSV because it is also the active local-GTF input: ${LOCAL_TOP_HITS_CSV_BASENAME}"
  else
    oda_delete_many "cleanup_local_hits_with_gtf_output_main_csv_${stamp}" --delete-file "${LOCAL_TOP_HITS_CSV_BASENAME}" || true
  fi
  oda_delete_many "cleanup_local_hits_with_gtf_output_prep_html_${stamp}" --delete-file "${OUTPUT_HTML_BASENAME%.html}.prep.html" || true

  while IFS= read -r remote_name; do
    [[ -z "${remote_name}" ]] && continue
    if [[ "${remote_name}" =~ ^${OUTPUT_HTML_BASENAME%.html}_part[0-9]+\.html$ ]]; then
      oda_delete_many "cleanup_local_hits_with_gtf_part_html_${stamp}" --delete-file "${remote_name}" || true
    fi
    if [[ "${remote_name}" =~ ^${LOCAL_TOP_HITS_CSV_BASENAME%.csv}_part[0-9]+\.csv$ ]]; then
      oda_delete_many "cleanup_local_hits_with_gtf_part_csv_${stamp}" --delete-file "${remote_name}" || true
    fi
  done <<EOF
${remote_listing}
EOF

  if [[ -n "${remote_png_path:-}" ]]; then
    remote_png_base="$(basename "${remote_png_path}")"
    if [[ -n "${remote_png_base}" ]]; then
      oda_delete_many "cleanup_local_hits_with_gtf_output_png_${stamp}" --delete-file "${remote_png_base}" || true
    fi
  fi
}

clear_stale_remote_expected_outputs() {
  if [[ "${SKIP_PRECLEAR_REMOTE_GTF_OUTPUTS:-0}" == "1" ]]; then
    echo "[cleanup] Skipping pre-submit remote local-GTF output cleanup because SKIP_PRECLEAR_REMOTE_GTF_OUTPUTS=${SKIP_PRECLEAR_REMOTE_GTF_OUTPUTS}."
    return 0
  fi
  echo "[cleanup] Clearing stale remote local-GTF outputs before submit..."
  oda_delete_many "preclear_local_hits_with_gtf_output_main_html_${stamp}" --delete-file "${OUTPUT_HTML_BASENAME}" || true
  if [[ -n "${LOCAL_TOP_HITS_INPUT_CSV_BASENAME:-}" && "${LOCAL_TOP_HITS_INPUT_CSV_BASENAME}" == "${LOCAL_TOP_HITS_CSV_BASENAME}" ]]; then
    echo "[cleanup] Preserving remote requested top-hit CSV before submit because it is also the active local-GTF input: ${LOCAL_TOP_HITS_CSV_BASENAME}"
  else
    oda_delete_many "preclear_local_hits_with_gtf_output_main_csv_${stamp}" --delete-file "${LOCAL_TOP_HITS_CSV_BASENAME}" || true
  fi
  oda_delete_many "preclear_local_hits_with_gtf_output_prep_html_${stamp}" --delete-file "${OUTPUT_HTML_BASENAME%.html}.prep.html" || true
}

cleanup_remote_local_gtf_subset() {
  if [[ ! -s "${LOCAL_GTF_SUBSET_GZ}" ]]; then
    return 0
  fi
  if [[ "${CLEAN_ODA_INPUT}" != "1" ]]; then
    echo "[cleanup] Keeping uploaded local GTF subset in SAS ODA because CLEAN_ODA_INPUT=${CLEAN_ODA_INPUT}."
    return 0
  fi
  echo "[cleanup] Removing uploaded local GTF subset from SAS ODA..."
  oda_delete_many \
    "cleanup_local_hits_with_gtf_subset_${stamp}" \
    --delete-file "${REMOTE_GTF_BASENAME}" || true
}

gtf_submit_needs_retry() {
  [[ ! -s "${RUN_LOG_FILE}" ]] && return 0
  if grep -Eiq 'We failed in getConnection|The application could not log on to the server|server configuration is invalid|SAS process has terminated unexpectedly' "${RUN_LOG_FILE}"; then
    return 0
  fi
  if grep -Eq 'HTML output saved to:|The final figure is put here:' "${RUN_LOG_FILE}"; then
    return 1
  fi
  local bytes
  bytes="$(wc -c < "${RUN_LOG_FILE}")"
  [[ "${bytes}" -lt 500 ]] && return 0
  return 1
}

gtf_submit_needs_retry_for_log() {
  local logfile="$1"
  [[ ! -s "${logfile}" ]] && return 0
  if grep -Eiq 'We failed in getConnection|The application could not log on to the server|server configuration is invalid|SAS process has terminated unexpectedly' "${logfile}"; then
    return 0
  fi
  if grep -Eq 'HTML output saved to:|The final figure is put here:' "${logfile}"; then
    return 1
  fi
  local bytes
  bytes="$(wc -c < "${logfile}")"
  [[ "${bytes}" -lt 500 ]] && return 0
  return 1
}

gtf_log_has_terminal_failure() {
  local logfile="$1"
  [[ -s "${logfile}" ]] || return 1
  if grep -Eiq 'ERROR: Insufficient space in file WORK\.|ERROR: File WORK\..* is damaged|ERROR: Sort initialization failure|We failed in getConnection|The application could not log on to the server|server configuration is invalid|No SAS process attached|SAS process has terminated unexpectedly|ERROR: The SAS job likely failed before producing the final figure' "${logfile}"; then
    return 0
  fi
  return 1
}

run_gtf_submit() {
  if [[ -x /usr/bin/timeout ]]; then
    /usr/bin/timeout --kill-after="${GTF_SUBMIT_TIMEOUT_GRACE_SECONDS}s" "${GTF_SUBMIT_TIMEOUT_SECONDS}s" \
      "${ODA_PERL[@]}" \
      --file "${RUN_SAS_RENDERED}" \
      --output-prefix "${RUN_PREFIX}"
    return $?
  fi

  "${ODA_PERL[@]}" \
    --file "${RUN_SAS_RENDERED}" \
    --output-prefix "${RUN_PREFIX}"
}

run_gtf_submit_for_file() {
  local sas_file="$1"
  local output_prefix="$2"
  local include_preflight_enabled="${3:-1}"
  local include_preflight_refresh_remote="${4:-${INCLUDE_PREFLIGHT_REFRESH_REMOTE}}"
  if [[ -x /usr/bin/timeout ]]; then
    env INCLUDE_PREFLIGHT_ENABLED="${include_preflight_enabled}" INCLUDE_PREFLIGHT_REFRESH_REMOTE="${include_preflight_refresh_remote}" \
      /usr/bin/timeout --kill-after="${GTF_SUBMIT_TIMEOUT_GRACE_SECONDS}s" "${GTF_SUBMIT_TIMEOUT_SECONDS}s" \
      "${ODA_PERL[@]}" \
      --file "${sas_file}" \
      --output-prefix "${output_prefix}"
    return $?
  fi

  env INCLUDE_PREFLIGHT_ENABLED="${include_preflight_enabled}" INCLUDE_PREFLIGHT_REFRESH_REMOTE="${include_preflight_refresh_remote}" \
    "${ODA_PERL[@]}" \
    --file "${sas_file}" \
    --output-prefix "${output_prefix}"
}

recover_html_from_submit_artifacts() {
  [[ -s "${HTML_OUT}" ]] && return 0
  local runner_html
  runner_html="$(find "${RUN_LOG_DIR}" -maxdepth 1 -type f -name 'sas_res_*.html' | head -n 1)"
  if [[ -n "${runner_html}" && -s "${runner_html}" ]]; then
    cp -f "${runner_html}" "${HTML_OUT}"
    echo "[recover] Reused the HTML artifact already saved by the submit helper: ${runner_html}"
    return 0
  fi
  return 1
}

recover_html_from_submit_artifacts_for_logdir() {
  local logdir="$1"
  local target_html="$2"
  [[ -s "${target_html}" ]] && return 0
  local runner_html
  runner_html="$(find "${logdir}" -maxdepth 1 -type f -name 'sas_res_*.html' | head -n 1)"
  if [[ -n "${runner_html}" && -s "${runner_html}" ]]; then
    cp -f "${runner_html}" "${target_html}"
    echo "[recover] Reused the HTML artifact already saved by the submit helper: ${runner_html}"
    return 0
  fi
  return 1
}

recover_csv_from_existing_local_copy() {
  [[ -s "${CSV_OUT}" ]] && return 0
  return 1
}

extract_remote_png_path_from_log_file() {
  local logfile="$1"
  [[ -s "${logfile}" ]] || return 1
  awk '
    /The final figure is put here:/ {
      for (i = 0; i < 40 && getline > 0; i++) {
        line = $0;
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", line);
        if (line ~ /^\/[^[:space:]]+\.png$/) {
          print line;
          exit;
        }
        if (match(line, /\/[^[:space:]]+\.png/)) {
          print substr(line, RSTART, RLENGTH);
          exit;
        }
      }
    }
  ' "${logfile}"
}

extract_remote_png_path_from_run_log() {
  extract_remote_png_path_from_log_file "${RUN_LOG_FILE}"
}

download_remote_png_to_path_if_reported() {
  local remote_png="$1"
  local png_path="$2"
  local output_prefix="$3"
  local remote_png_basename=""
  local remote_png_home=""
  [[ -n "${remote_png}" ]] || return 0
  rm -f "${png_path}"
  echo "[recover] Downloading remote PNG reported by SAS log: ${remote_png}"
  oda_download_many_with_timeout \
    "${ODA_PNG_DOWNLOAD_TIMEOUT_SECONDS}" \
    "${ODA_PNG_DOWNLOAD_TIMEOUT_GRACE_SECONDS}" \
    "${output_prefix}" \
    --download-file "${remote_png}" \
    --download-local-path "${png_path}" || true
  [[ -s "${png_path}" ]] && return 0

  remote_png_basename="$(basename "${remote_png}")"
  if [[ -n "${remote_png_basename}" && "${remote_png_basename}" != "${remote_png}" ]]; then
    echo "[recover] Direct remote PNG download did not succeed. Retrying with SAS ODA home-relative path: ~/${remote_png_basename}"
    oda_download_many_with_timeout \
      "${ODA_PNG_DOWNLOAD_TIMEOUT_SECONDS}" \
      "${ODA_PNG_DOWNLOAD_TIMEOUT_GRACE_SECONDS}" \
      "${output_prefix}_home" \
      --download-file "~/${remote_png_basename}" \
      --download-local-path "${png_path}" || true
    [[ -s "${png_path}" ]] && return 0
  fi

  if [[ -n "${remote_png_basename}" ]]; then
    echo "[recover] Direct remote PNG download still did not succeed. Retrying with basename only: ${remote_png_basename}"
    oda_download_many_with_timeout \
      "${ODA_PNG_DOWNLOAD_TIMEOUT_SECONDS}" \
      "${ODA_PNG_DOWNLOAD_TIMEOUT_GRACE_SECONDS}" \
      "${output_prefix}_base" \
      --download-file "${remote_png_basename}" \
      --download-local-path "${png_path}" || true
  fi
}

download_remote_png_if_reported() {
  local remote_png="$1"
  download_remote_png_to_path_if_reported "${remote_png}" "${PNG_OUT}" "download_local_hits_with_gtf_png_${stamp}"
}

extract_embedded_png_from_html_path_if_present() {
  local html_path="$1"
  local png_path="$2"
  [[ -s "${html_path}" ]] || return 1
  [[ -s "${png_path}" ]] && return 0
  perl -MMIME::Base64 -0777 -ne '
    if (m{data:image/png;base64,([^"'\'' ]+)}s) {
      open my $fh, ">:raw", $ARGV[1] or die "open $ARGV[1]: $!";
      print {$fh} MIME::Base64::decode_base64($1);
      close $fh;
      exit 0;
    }
    exit 1;
  ' "${html_path}" "${png_path}" 2>/dev/null || return 1
  [[ -s "${png_path}" ]]
}

extract_embedded_png_from_html_if_present() {
  extract_embedded_png_from_html_path_if_present "${HTML_OUT}" "${PNG_OUT}"
}

extract_embedded_png_data_uri_from_html_path() {
  local html_path="$1"
  [[ -s "${html_path}" ]] || return 1
  perl -0777 -ne '
    if (m{(data:image/png;base64,[^"'\'' ]+)}s) {
      print $1;
      exit 0;
    }
    exit 1;
  ' "${html_path}" 2>/dev/null
}

extract_embedded_png_data_uri_from_html() {
  extract_embedded_png_data_uri_from_html_path "${HTML_OUT}"
}

resolve_magick_bin() {
  local candidate resolved
  for candidate in \
    "${MAGICK_BIN:-}" \
    magick \
    magick.exe \
    /mnt/c/Users/cheng/Downloads/cygwin-portable-20210411/cygwin-portable/App/cygwin/bin/magick.exe \
    /cygdrive/c/Users/cheng/Downloads/cygwin-portable-20210411/cygwin-portable/App/cygwin/bin/magick.exe; do
    [[ -n "${candidate}" ]] || continue
    resolved="$(command -v "${candidate}" 2>/dev/null || true)"
    if [[ -n "${resolved}" && -x "${resolved}" ]]; then
      printf '%s\n' "${resolved}"
      return 0
    fi
    if [[ -x "${candidate}" ]]; then
      printf '%s\n' "${candidate}"
      return 0
    fi
  done
  return 1
}

annotate_png_with_visible_axis_labels_if_needed() {
  local png_path="$1"
  local ylabel="${2:-}"
  local colorbar_label="${3:-}"
  local magick_bin=""
  local canvas_png=""
  local label_png=""
  local final_png=""
  local y_enabled=0
  local cbar_enabled=0
  local left_pad=0
  local right_pad=0
  local y_label_png=""
  local cbar_label_png=""

  [[ -s "${png_path}" ]] || return 0

  if [[ "${VISIBLE_YLABEL_ENABLED}" == "1" && -n "${ylabel}" ]]; then
    y_enabled=1
    left_pad="${VISIBLE_YLABEL_LEFT_PAD_PX}"
  fi
  if [[ "${VISIBLE_COLORBAR_LABEL_ENABLED}" == "1" && -n "${colorbar_label}" ]]; then
    cbar_enabled=1
    right_pad="${VISIBLE_COLORBAR_LABEL_RIGHT_PAD_PX}"
  fi
  (( y_enabled == 1 || cbar_enabled == 1 )) || return 0

  magick_bin="$(resolve_magick_bin || true)"
  if [[ -z "${magick_bin}" ]]; then
    echo "[recover] Skipping visible axis annotation because ImageMagick is unavailable."
    return 0
  fi

  canvas_png="${png_path}.ylabel_canvas.png"
  final_png="${png_path}.ylabel_final.png"
  y_label_png="${png_path}.ylabel_label.png"
  cbar_label_png="${png_path}.cbar_label.png"

  if ! "${magick_bin}" "${png_path}" -background white -gravity center -splice "${left_pad}x0" -gravity east -splice "${right_pad}x0" "${canvas_png}" >/dev/null 2>&1; then
    rm -f "${canvas_png}" "${y_label_png}" "${cbar_label_png}" "${final_png}"
    echo "[recover] WARNING: Could not extend the PNG canvas for visible axis labels: ${png_path}" >&2
    return 0
  fi

  if (( y_enabled == 1 )); then
    if ! "${magick_bin}" -background none -fill "${VISIBLE_YLABEL_FILL}" -font Arial -pointsize "${VISIBLE_YLABEL_POINTSIZE}" label:"${ylabel}" -rotate 90 "${y_label_png}" >/dev/null 2>&1; then
      rm -f "${canvas_png}" "${y_label_png}" "${cbar_label_png}" "${final_png}"
      echo "[recover] WARNING: Could not render a visible y-axis label overlay: ${ylabel}" >&2
      return 0
    fi
    if ! "${magick_bin}" "${canvas_png}" "${y_label_png}" -gravity west -geometry +24+0 -composite "${final_png}" >/dev/null 2>&1; then
      rm -f "${canvas_png}" "${y_label_png}" "${cbar_label_png}" "${final_png}"
      echo "[recover] WARNING: Could not composite the visible y-axis label onto ${png_path}" >&2
      return 0
    fi
    mv -f "${final_png}" "${canvas_png}"
  fi

  if (( cbar_enabled == 1 )); then
    if ! "${magick_bin}" -background none -fill "${VISIBLE_COLORBAR_LABEL_FILL}" -font Arial -pointsize "${VISIBLE_COLORBAR_LABEL_POINTSIZE}" label:"${colorbar_label}" -rotate 270 "${cbar_label_png}" >/dev/null 2>&1; then
      rm -f "${canvas_png}" "${y_label_png}" "${cbar_label_png}" "${final_png}"
      echo "[recover] WARNING: Could not render a visible colorbar label overlay: ${colorbar_label}" >&2
      return 0
    fi
    if ! "${magick_bin}" "${canvas_png}" "${cbar_label_png}" -gravity east -geometry +18+0 -composite "${final_png}" >/dev/null 2>&1; then
      rm -f "${canvas_png}" "${y_label_png}" "${cbar_label_png}" "${final_png}"
      echo "[recover] WARNING: Could not composite the visible colorbar label onto ${png_path}" >&2
      return 0
    fi
    mv -f "${final_png}" "${canvas_png}"
  fi

  mv -f "${canvas_png}" "${png_path}"
  rm -f "${y_label_png}" "${cbar_label_png}" "${final_png}"
  if (( y_enabled == 1 )) && (( cbar_enabled == 1 )); then
    echo "[recover] Added visible y-axis label '${ylabel}' and colorbar label '${colorbar_label}' to local GTF PNG: ${png_path}"
  elif (( y_enabled == 1 )); then
    echo "[recover] Added visible y-axis label '${ylabel}' to local GTF PNG: ${png_path}"
  else
    echo "[recover] Added visible colorbar label '${colorbar_label}' to local GTF PNG: ${png_path}"
  fi
}

build_completed_html_from_png_assets_if_available() {
  local html_path="$1"
  local png_path="$2"
  local raw_html_path="$3"
  local csv_path="${4:-}"
  local figure_title="${5:-Local Manhattan and GTF plot for top hits}"
  local image_alt="${6:-Local top-hit GTF plot}"
  local image_src=""
  local use_iframe=0

  if [[ -s "${png_path}" ]]; then
    annotate_png_with_visible_axis_labels_if_needed "${png_path}" "${GTF_YAXIS_LABEL}" "${GTF_COLORBAR_LABEL}"
    image_src="$(basename "${png_path}")"
  else
    if extract_embedded_png_from_html_path_if_present "${html_path}" "${png_path}"; then
      [[ -s "${png_path}" ]] || return 0
      annotate_png_with_visible_axis_labels_if_needed "${png_path}" "${GTF_YAXIS_LABEL}" "${GTF_COLORBAR_LABEL}"
      image_src="$(basename "${png_path}")"
    else
      image_src="$(extract_embedded_png_data_uri_from_html_path "${html_path}" || true)"
      if [[ -z "${image_src}" ]]; then
        use_iframe=1
      fi
    fi
  fi
  [[ -s "${html_path}" ]] || return 0

  cp -f "${html_path}" "${raw_html_path}"
  {
    echo '<!DOCTYPE html>'
    echo "<html lang=\"en\"><head><meta charset=\"utf-8\"><title>${figure_title}</title></head><body style=\"font-family:Arial,Helvetica,sans-serif;margin:24px\">"
    echo "<h1 style=\"font-size:20px\">${figure_title}</h1>"
    echo "<div style=\"display:flex;align-items:center;gap:14px\">"
    if [[ "${VISIBLE_YLABEL_ENABLED}" == "1" ]]; then
      echo "<div style=\"writing-mode:vertical-rl;transform:rotate(180deg);font-size:20px;font-weight:600;color:${VISIBLE_YLABEL_FILL};line-height:1\">${GTF_YAXIS_LABEL}</div>"
    fi
    if [[ "${use_iframe}" -eq 1 ]]; then
      echo "<div style=\"flex:1;min-width:0\"><iframe src=\"$(basename "${raw_html_path}")\" title=\"${image_alt}\" style=\"width:100%;height:1700px;border:1px solid #ccc;background:#fff\"></iframe></div>"
    else
      echo "<div><img src=\"${image_src}\" alt=\"${image_alt}\" style=\"max-width:100%;height:auto;border:1px solid #ccc\"></div>"
    fi
    echo '</div>'
    if [[ -n "${csv_path}" && -s "${csv_path}" ]]; then
      echo "<p><a href=\"$(basename "${csv_path}")\">Top-hit CSV</a></p>"
    fi
    echo "<p><a href=\"$(basename "${raw_html_path}")\">Raw SAS HTML output</a></p>"
    echo '</body></html>'
  } > "${html_path}"
  if [[ "${use_iframe}" -eq 1 ]]; then
    if [[ "${VISIBLE_YLABEL_ENABLED}" == "1" ]]; then
      echo "[recover] Wrapped the raw SAS HTML with a visible y-axis label because no standalone PNG could be extracted: ${html_path}"
    else
      echo "[recover] Wrapped the raw SAS HTML because no standalone PNG could be extracted: ${html_path}"
    fi
  else
    echo "[recover] Replaced the opened HTML with a figure-first wrapper because a final PNG was generated: ${html_path}"
  fi
}

build_completed_html_from_png_if_available() {
  build_completed_html_from_png_assets_if_available \
    "${HTML_OUT}" \
    "${PNG_OUT}" \
    "${RAW_HTML_OUT}" \
    "${CSV_OUT}" \
    "Local Manhattan and GTF plot for top hits" \
    "Local top-hit GTF plot"
}

html_has_rendered_plot() {
  [[ -s "${HTML_OUT}" ]] || return 1
  perl -0777 -ne '
    exit 0 if /<img\b/i;
    exit 0 if /data:image\/png;base64,/i;
    exit 0 if /<svg\b/i;
    exit 0 if /<iframe\b/i;
    exit 1;
  ' "${HTML_OUT}" 2>/dev/null
}

delivered_gtf_artifact_ready() {
  [[ -s "${HTML_OUT}" ]] || return 1
  if html_has_rendered_plot; then
    return 0
  fi
  [[ -s "${RAW_HTML_OUT}" ]]
}

remote_data_exists() {
  local check_output
  check_output="$(
    run_oda_helper_with_timeout \
      "${ODA_REMOTE_INFO_TIMEOUT_SECONDS}" \
      "${ODA_REMOTE_INFO_TIMEOUT_GRACE_SECONDS}" \
      --dir4listing '~' \
      --output-prefix "check_remote_local_gtf_data_${stamp}" 2>&1 || true
  )"
  printf '%s\n' "${check_output}" | grep -Fxq "${REMOTE_DATA_BASENAME}"
}

remote_home_file_size_bytes() {
  local remote_basename="$1"
  local info_output
  info_output="$(
    run_oda_helper_with_timeout \
      "${ODA_REMOTE_INFO_TIMEOUT_SECONDS}" \
      "${ODA_REMOTE_INFO_TIMEOUT_GRACE_SECONDS}" \
      --file-info "~/${remote_basename}" \
      --output-prefix "check_remote_home_size_${stamp}" 2>&1 || true
  )"
  printf '%s\n' "${info_output}" | awk -F '\t' '$1=="SIZE"{print $2}' | tail -n 1
}

remote_data_size_bytes() {
  remote_home_file_size_bytes "${REMOTE_DATA_BASENAME}"
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

remote_data_matches_local_size() {
  remote_home_file_matches_local_size "${DATA_GZ}" "${REMOTE_DATA_BASENAME}"
}

remote_data_known_to_oda() {
  local remote_size
  remote_size="$(remote_data_size_bytes | tr -d '[:space:]')"
  [[ -n "${remote_size}" ]] && return 0
  remote_data_exists
}

delete_partial_remote_data() {
  run_oda_helper \
    --delete-file "${REMOTE_DATA_BASENAME}" \
    --output-prefix "delete_partial_local_gtf_subset_${stamp}" >/dev/null 2>&1 || true
}

upload_data_with_integrity_check() {
  local local_size remote_size verify_attempt upload_attempt max_upload_attempts max_verify_attempts
  local_size="$(wc -c < "${DATA_GZ}" | tr -d '[:space:]')"
  max_upload_attempts="${ODA_UPLOAD_MAX_ATTEMPTS:-2}"
  max_verify_attempts="${ODA_UPLOAD_VERIFY_ATTEMPTS:-4}"

  upload_attempt=1
  while [[ "${upload_attempt}" -le "${max_upload_attempts}" ]]; do
    run_oda_helper --upload-file "${DATA_GZ}" --output-prefix "upload_local_hits_with_gtf_subset_${stamp}_try${upload_attempt}"

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

upload_home_file_if_needed() {
  local step_label="$1"
  local local_path="$2"
  local remote_basename="$3"
  local output_prefix="$4"
  local upload_path="${local_path}"
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

stage_alias_upload_if_needed() {
  local local_path="$1"
  local remote_basename="$2"
  if remote_home_file_matches_local_size "${local_path}" "${remote_basename}"; then
    echo "[bulk] Reusing existing remote file in SAS ODA home: ${remote_basename}"
    return 1
  fi
  mkdir -p "${WORKDIR}/.oda_upload_aliases"
  local upload_path="${WORKDIR}/.oda_upload_aliases/${remote_basename}"
  cp -f "${local_path}" "${upload_path}"
  bulk_upload_args+=(--upload-file "${upload_path}")
  return 0
}

stable_hash_text() {
  perl -MDigest::MD5=md5_hex -e 'print md5_hex(join("\n", @ARGV))' "$@"
}

stable_hash_file() {
  local path="$1"
  perl -MDigest::MD5=md5_hex -e '
    use strict;
    use warnings;
    my ($path) = @ARGV;
    open my $fh, q{<}, $path or die "Cannot open $path: $!\n";
    binmode $fh;
    my $ctx = Digest::MD5->new;
    $ctx->addfile($fh);
    print $ctx->hexdigest;
  ' "${path}"
}

find_missing_target_snps_in_data_gz() {
  [[ -n "${TARGET_SNP_LIST}" ]] || return 0
  [[ -s "${DATA_GZ}" ]] || return 0
  perl -MIO::Uncompress::Gunzip=gunzip,\$GunzipError -e '
    use strict;
    use warnings;
    my ($path, $target_text) = @ARGV;
    my @targets = grep { length } map { s/^\s+|\s+$//gr } split /,/, ($target_text // q{});
    exit 0 unless @targets;
    my %want = map { $_ => 1 } @targets;
    my $fh = IO::Uncompress::Gunzip->new($path) or die "Cannot open $path: $GunzipError\n";
    my $hdr = <$fh>;
    exit 0 unless defined $hdr;
    chomp $hdr;
    my @h = split /\t/, $hdr, -1;
    my %idx;
    @idx{@h} = (0 .. $#h);
    die "SNP column not found in $path\n" unless exists $idx{SNP};
    while (my $line = <$fh>) {
      chomp $line;
      my @f = split /\t/, $line, -1;
      my $snp = $f[$idx{SNP}] // q{};
      delete $want{$snp} if exists $want{$snp};
      last unless %want;
    }
    print "$_\n" for @targets ? grep { exists $want{$_} } @targets : ();
  ' "${DATA_GZ}" "${TARGET_SNP_LIST}"
}

augment_data_gz_with_missing_target_snps() {
  [[ -n "${TARGET_SNP_LIST}" ]] || return 0

  if [[ -z "${SOURCE_LONG_GZ:-}" || ! -s "${SOURCE_LONG_GZ}" ]]; then
    echo "WARNING: TARGET_SNP_LIST was provided, but SOURCE_LONG_GZ is unavailable for building a compact target-only plotting subset." >&2
    return 0
  fi
  if [[ -z "${SCHEMA_CONFIG_JSON:-}" || ! -s "${SCHEMA_CONFIG_JSON}" ]]; then
    echo "WARNING: TARGET_SNP_LIST was provided, but SCHEMA_CONFIG_JSON is unavailable for building a compact target-only plotting subset." >&2
    return 0
  fi

  local target_cache_key target_cache_base
  target_cache_key="$(
    stable_hash_text \
      "${PROJECT_TAG}" \
      "${TARGET_SNP_LIST}" \
      "${LOCAL_GTF_WINDOW_BP}" \
      "${SOURCE_LONG_GZ}" \
      "${SCHEMA_CONFIG_JSON}"
  )"
  target_cache_base="${LOCAL_GTF_REUSE_CACHE_DIR}/target_snp_augmented_${SAFE_PROJECT_TAG}_${target_cache_key}"
  if [[ -s "${target_cache_base}.tsv.gz" ]]; then
    TARGET_SNP_AUG_GZ="${target_cache_base}.tsv.gz"
    DATA_GZ="${TARGET_SNP_AUG_GZ}"
    REMOTE_DATA_BASENAME="$(basename "${TARGET_SNP_AUG_GZ}")"
    TARGET_SNP_AUG_CACHE_MANAGED=1
    echo "[prep] Reusing cached compact target-SNP local GTF plotting subset: ${DATA_GZ}"
    return 0
  fi

  mkdir -p "${TARGET_SNP_AUG_DIR}"
  local snp safe_snp single_out single_manifest extra_args=()
  local -a target_snp_array=()
  echo "[prep] Building a compact local GTF plotting subset directly from SOURCE_LONG_GZ for the requested target SNP(s)..."
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
      --window-bp "${LOCAL_GTF_WINDOW_BP}" \
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
    echo "[prep] Using compact target-SNP local GTF plotting subset: ${DATA_GZ}"
  fi
}

augment_data_gz_with_missing_target_snps

perl "${SCHEMA_INCLUDE_HELPER}" \
  --config "${SCHEMA_CONFIG_JSON}" \
  --dataset scz_mh \
  --source-type gzip \
  --remote-basename "${REMOTE_DATA_BASENAME}" > "${IMPORT_BLOCK_RENDERED}"

: > "${GTF_IMPORT_BLOCK_RENDERED}"

gtf_region_source=""
gtf_region_args=()
if [[ -s "${CSV_OUT}" ]]; then
  gtf_region_source="${CSV_OUT}"
  mapfile -t gtf_region_args < <(
    perl -e '
      use strict;
      use warnings;
      my ($csv, $win) = @ARGV;
      open my $fh, q{<}, $csv or die "Cannot open $csv: $!\n";
      my $header = <$fh>;
      defined $header or exit 0;
      chomp $header;
      my @h = split /,/, $header;
      my %idx;
      for my $i (0 .. $#h) { $idx{$h[$i]} = $i; }
      die "CSV is missing CHR/BP columns\n" unless exists $idx{CHR} && exists $idx{BP};
      my %seen;
      while (my $line = <$fh>) {
        chomp $line;
        next unless length $line;
        my @f = split /,/, $line;
        next unless defined $f[$idx{CHR}] && defined $f[$idx{BP}];
        my $chr = $f[$idx{CHR}];
        my $bp  = $f[$idx{BP}];
        next unless defined $chr && defined $bp;
        next unless $chr =~ /^\d+$/ && $bp =~ /^\d+(?:\.\d+)?$/;
        my $start = int($bp - $win);
        $start = 1 if $start < 1;
        my $end = int($bp + $win);
        my $region = $chr . q{:} . $start . q{:} . $end;
        next if $seen{$region}++;
        print "--region\n", $region, "\n";
      }
      close $fh;
    ' "${CSV_OUT}" "${LOCAL_GTF_WINDOW_BP}"
  )
elif [[ -n "${VERIFY_TOP_HITS_TSV}" && -s "${VERIFY_TOP_HITS_TSV}" ]]; then
  gtf_region_source="${VERIFY_TOP_HITS_TSV}"
  mapfile -t gtf_region_args < <(
    perl -e '
      use strict;
      use warnings;
      my ($tsv, $win) = @ARGV;
      open my $fh, q{<}, $tsv or die "Cannot open $tsv: $!\n";
      my $header = <$fh>;
      defined $header or exit 0;
      chomp $header;
      my @h = split /\t/, $header;
      my %idx;
      for my $i (0 .. $#h) { $idx{$h[$i]} = $i; }
      die "TSV is missing CHR/BP columns\n" unless exists $idx{CHR} && exists $idx{BP};
      my %seen;
      while (my $line = <$fh>) {
        chomp $line;
        next unless length $line;
        my @f = split /\t/, $line;
        next unless defined $f[$idx{CHR}] && defined $f[$idx{BP}];
        my $chr = $f[$idx{CHR}];
        my $bp  = $f[$idx{BP}];
        next unless defined $chr && defined $bp;
        next unless $chr =~ /^\d+$/ && $bp =~ /^\d+(?:\.\d+)?$/;
        my $start = int($bp - $win);
        $start = 1 if $start < 1;
        my $end = int($bp + $win);
        my $region = $chr . q{:} . $start . q{:} . $end;
        next if $seen{$region}++;
        print "--region\n", $region, "\n";
      }
      close $fh;
    ' "${VERIFY_TOP_HITS_TSV}" "${LOCAL_GTF_WINDOW_BP}"
  )
fi

if [[ ${#gtf_region_args[@]} -gt 0 ]]; then
  local_gft_region_hash_input_file=""
  if [[ -n "${gtf_region_source}" && -s "${gtf_region_source}" ]]; then
    local_gft_region_hash_input_file="${gtf_region_source}"
  fi
  if [[ -n "${local_gft_region_hash_input_file}" ]]; then
    gtf_region_key="$(
      stable_hash_text \
        "${PROJECT_TAG}" \
        "${LOCAL_GTF_WINDOW_BP}" \
        "${LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES}" \
        "${GTF_GZ_URL}" \
        "$(stable_hash_file "${local_gft_region_hash_input_file}")"
    )"
  else
    gtf_region_key="$(
      stable_hash_text \
        "${PROJECT_TAG}" \
        "${LOCAL_GTF_WINDOW_BP}" \
        "${LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES}" \
        "${GTF_GZ_URL}" \
        "${gtf_region_args[*]}"
    )"
  fi
  gtf_subset_cache_base="${LOCAL_GTF_REUSE_CACHE_DIR}/local_gtf_subset_${SAFE_PROJECT_TAG}_${gtf_region_key}"
  if [[ -s "${gtf_subset_cache_base}.tsv.gz" ]]; then
    LOCAL_GTF_SUBSET="${gtf_subset_cache_base}.tsv"
    LOCAL_GTF_SUBSET_GZ="${gtf_subset_cache_base}.tsv.gz"
    LOCAL_GTF_SUBSET_CACHE_MANAGED=1
    echo "[prep] Reusing cached local GTF subset from ${gtf_region_source:-region list} using window ${LOCAL_GTF_WINDOW_BP}: ${LOCAL_GTF_SUBSET_GZ}"
  else
    echo "[prep] Building local GTF subset from ${gtf_region_source:-region list} using window ${LOCAL_GTF_WINDOW_BP}..."
    gtf_subset_cmd=(
      perl "${GTF_SUBSET_HELPER}"
      --gtf-url "${GTF_GZ_URL}"
      --cache-dir "${GTF_CACHE_DIR}"
      --output "${LOCAL_GTF_SUBSET}"
      "${gtf_region_args[@]}"
    )
    if [[ "${LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES}" == "1" ]]; then
      gtf_subset_cmd+=(--include-non-protein-coding)
    else
      gtf_subset_cmd+=(--no-include-non-protein-coding)
    fi
    "${gtf_subset_cmd[@]}"

    if [[ -s "${LOCAL_GTF_SUBSET}" ]]; then
      gzip -c "${LOCAL_GTF_SUBSET}" > "${LOCAL_GTF_SUBSET_GZ}"
      cp -f "${LOCAL_GTF_SUBSET}" "${gtf_subset_cache_base}.tsv"
      cp -f "${LOCAL_GTF_SUBSET_GZ}" "${gtf_subset_cache_base}.tsv.gz"
      LOCAL_GTF_SUBSET="${gtf_subset_cache_base}.tsv"
      LOCAL_GTF_SUBSET_GZ="${gtf_subset_cache_base}.tsv.gz"
      LOCAL_GTF_SUBSET_CACHE_MANAGED=1
    fi
  fi

  if [[ -s "${LOCAL_GTF_SUBSET_GZ}" ]]; then
    perl "${GTF_IMPORT_INCLUDE_HELPER}" \
      --dataset "${GTF_LOCAL_DSD}" \
      --remote-basename "${REMOTE_GTF_BASENAME}" > "${GTF_IMPORT_BLOCK_RENDERED}"
  fi
fi

render_gtf_runner() {
  local output_sas="$1"
  local output_html_basename="$2"
  local output_csv_basename="$3"
  local input_csv_basename="${4:-}"
  local prep_only="${5:-0}"
  local label_snps="${GTF_LABEL_SNPS}"
  local label_text_rotate_angle=""

  label_snps="${label_snps//,/ }"
  case "${GTF_LABEL_LAYOUT}" in
    vertical|VERTICAL)
      label_text_rotate_angle="90"
      ;;
    horizontal|HORIZONTAL)
      label_text_rotate_angle="0"
      ;;
    *)
      label_text_rotate_angle=""
      ;;
  esac

  perl "${RENDER_SAS_HELPER}" \
    --template "${RUN_SAS_TEMPLATE}" \
    --output "${output_sas}" \
    --replace "TOP_HIT_FOCUS_PVAR=${TOP_HIT_FOCUS_PVAR}" \
    --replace "TOP_HIT_MODE=${TOP_HIT_MODE:-differential}" \
    --replace "TOP_HIT_FILTER_EXPR=${TOP_HIT_FILTER_EXPR}" \
    --replace "TOP_HIT_SIGNAL_THRSHD=${TOP_HIT_SIGNAL_THRSHD}" \
    --replace "TOP_HIT_SIGNAL_THRSHDS=${TOP_HIT_SIGNAL_THRSHDS:-${TOP_HIT_SIGNAL_THRSHD}}" \
    --replace "TOP_HIT_DIST_BP=${TOP_HIT_DIST_BP}" \
    --replace "TARGET_SNP_LIST=${TARGET_SNP_LIST}" \
    --replace "TARGET_SNP_GENES=${TARGET_SNP_GENES}" \
    --replace "LOCAL_MAX_HITS_PER_FIG=${LOCAL_MAX_HITS_PER_FIG}" \
    --replace "LOCAL_TOP_HITS_CSV_BASENAME=${output_csv_basename}" \
    --replace "LOCAL_TOP_HITS_INPUT_CSV_BASENAME=${input_csv_basename}" \
    --replace "COMMON_ASSOC_P_VARS=${COMMON_ASSOC_P_VARS:-}" \
    --replace "PREP_ONLY=${prep_only}" \
    --replace "LOCAL_WINDOW_BP=${LOCAL_GTF_WINDOW_BP}" \
    --replace "OUTPUT_HTML=${output_html_basename}" \
    --replace "GTF_DSD=${GTF_DSD}" \
    --replace "FM_LIBPATH=${FM_LIBPATH}" \
    --replace "GTF_LOCAL_DSD=${GTF_LOCAL_DSD}" \
    --replace "GTF_GZ_URL=${GTF_GZ_URL}" \
    --replace "GTF_ASSOC_PVARS=${GTF_ASSOC_PVARS}" \
    --replace "GTF_ZSCORE_VARS=${GTF_ZSCORE_VARS}" \
    --replace "GTF_LABELS=${GTF_LABELS}" \
    --replace "GTF_DIST2SNP=${GTF_DIST2SNP}" \
    --replace "GTF_DESIGN_WIDTH=${GTF_DESIGN_WIDTH}" \
    --replace "GTF_DESIGN_HEIGHT=${GTF_DESIGN_HEIGHT}" \
    --replace "GTF_DIST2SEP_GENES=${GTF_DIST2SEP_GENES}" \
    --replace "GTF_SHIFT_TEXT_YVAL=${GTF_SHIFT_TEXT_YVAL}" \
    --replace "GTF_PCT4NEG_Y=${GTF_PCT4NEG_Y}" \
    --replace "GTF_ADJVAL4HEADER=${GTF_ADJVAL4HEADER}" \
    --replace "GTF_YAXIS_LABEL=${GTF_YAXIS_LABEL}" \
    --replace "GTF_COLORBAR_LABEL=${GTF_COLORBAR_LABEL}" \
    --replace "GTF_YAXIS_OFFSET4MAX=${GTF_YAXIS_OFFSET4MAX}" \
    --replace "GTF_YOFFSET4TEXTLABELS=${GTF_YOFFSET4TEXTLABELS}" \
    --replace "GTF_YOFFSET4MAX_DRAWMARKERSONTOP=${GTF_YOFFSET4MAX_DRAWMARKERSONTOP}" \
    --replace "GTF_LABEL_SNPS=${label_snps}" \
    --replace "GTF_LABEL_TEXT_ROTATE_ANGLE=${label_text_rotate_angle}" \
    --replace "GTF_INCLUDE_NON_PROTEIN_CODING=${LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES}" \
    --replace "GET_GTF_MACRO_BASENAME=${GET_GTF_MACRO_BASENAME}" \
    --replace "GET_GTF_MACRO_NAME=${GET_GTF_MACRO_NAME}" \
    --replace-file "WIDE_IMPORT_BLOCK=${IMPORT_BLOCK_RENDERED}" \
    --replace-file "GTF_IMPORT_BLOCK=${GTF_IMPORT_BLOCK_RENDERED}"
}

render_gtf_runner \
  "${RUN_SAS_RENDERED}" \
  "${OUTPUT_HTML_BASENAME}" \
  "${LOCAL_TOP_HITS_CSV_BASENAME}" \
  "${LOCAL_TOP_HITS_INPUT_CSV_BASENAME}" \
  "0"

if ! generate_requested_top_hits_csv_locally; then
  echo "WARNING: Local MAF-aware top-hit CSV generation did not succeed. The wrapper will fall back to the prep-only SAS export path if needed." >&2
fi

rm -f "${HTML_OUT}"
mkdir -p "${GET_GTF_MACRO_UPLOAD_DIR}"

if [[ -f "${GET_GTF_MACRO_SAS}" ]]; then
  cp "${GET_GTF_MACRO_SAS}" "${GET_GTF_MACRO_UPLOAD}"
else
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
    length seqname chr chr_raw source $64 feature $32 score $32 strand $4 frame $8 attribute $32767;
    length gene_id gene_name transcript_id transcript_name gene_type transcript_type exon_id exon_number level $256;
    length gene havana_gene havana_transcript transcript_support_level tag $256;
    length start end st en bp1 bp2 txStart txEnd 8;
    infile _gtfgz dlm='09'x dsd truncover lrecl=1048576 firstobs=1;
    input seqname :$64.
          source :$64.
          feature :$32.
          start
          end
          score :$32.
          strand :$4.
          frame :$8.
          attribute :$32767.;
    if missing(seqname) then delete;
    if substr(seqname,1,1)='#' then delete;

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

    gene_id=prxchange('s/.*gene_id "([^"]+)".*/$1/i',1,attribute);
    if gene_id=attribute then gene_id='';
    gene_name=prxchange('s/.*gene_name "([^"]+)".*/$1/i',1,attribute);
    if gene_name=attribute then gene_name='';
    transcript_id=prxchange('s/.*transcript_id "([^"]+)".*/$1/i',1,attribute);
    if transcript_id=attribute then transcript_id='';
    transcript_name=prxchange('s/.*transcript_name "([^"]+)".*/$1/i',1,attribute);
    if transcript_name=attribute then transcript_name='';
    gene_type=prxchange('s/.*gene_type "([^"]+)".*/$1/i',1,attribute);
    if gene_type=attribute then gene_type='';
    transcript_type=prxchange('s/.*transcript_type "([^"]+)".*/$1/i',1,attribute);
    if transcript_type=attribute then transcript_type='';
    exon_id=prxchange('s/.*exon_id "([^"]+)".*/$1/i',1,attribute);
    if exon_id=attribute then exon_id='';
    exon_number=prxchange('s/.*exon_number "([^"]+)".*/$1/i',1,attribute);
    if exon_number=attribute then exon_number='';
    level=prxchange('s/.*level "([^"]+)".*/$1/i',1,attribute);
    if level=attribute then level='';
    havana_gene=prxchange('s/.*havana_gene "([^"]+)".*/$1/i',1,attribute);
    if havana_gene=attribute then havana_gene='';
    havana_transcript=prxchange('s/.*havana_transcript "([^"]+)".*/$1/i',1,attribute);
    if havana_transcript=attribute then havana_transcript='';
    transcript_support_level=prxchange('s/.*transcript_support_level "([^"]+)".*/$1/i',1,attribute);
    if transcript_support_level=attribute then transcript_support_level='';
    tag=prxchange('s/.*tag "([^"]+)".*/$1/i',1,attribute);
    if tag=attribute then tag='';
    gene=coalescec(gene_name,gene_id,transcript_name,transcript_id,feature);
  run;
  filename _gtfgz clear;
%mend;
EOF
  perl -0pi -e 's/__GET_GTF_MACRO_NAME__/\Q'"${GET_GTF_MACRO_NAME}"'\E/g' "${GET_GTF_MACRO_UPLOAD}"
  echo "WARNING: ${GET_GTF_MACRO_SAS} not found; uploading an auto-generated downloader macro instead." >&2
fi

if [[ "${EMIT_LOCAL_SAS_DEBUG}" == "1" || "${LOCAL_SAS_DEBUG_ONLY}" == "1" ]]; then
  perl "${LOCAL_SAS_DEBUG_EMITTER}" \
    --mode local_gtf \
    --input "${RUN_SAS_RENDERED}" \
    --output "${LOCAL_DEBUG_SAS_RENDERED}" \
    --workdir "${WORKDIR}" \
    --deps-dir "${DEPS_DIR}" \
    --data-gz "${DATA_GZ}" \
    --gtf-subset-gz "${LOCAL_GTF_SUBSET_GZ}" \
    --top-hits-csv "${CSV_OUT}" \
    --gtf-macro-upload "${GET_GTF_MACRO_UPLOAD}" \
    --gtf-local-dataset "${GTF_LOCAL_DSD}" \
    --output-html-basename "${OUTPUT_HTML_BASENAME}" >/dev/null
  echo "[local-sas] Emitted local-SAS debug script: ${LOCAL_DEBUG_SAS_RENDERED}"
fi

if [[ "${LOCAL_SAS_DEBUG_ONLY}" == "1" ]]; then
  echo "[local-sas] LOCAL_SAS_DEBUG_ONLY=1, skipping SAS ODA submit."
  exit 0
fi

generate_top_hits_csv_for_batching() {
  local prep_run_sas_rendered="${WORKDIR}/run_sas_oda_local_top_hits_with_gtf.${stamp}.prep.sas"
  local prep_output_html_basename="${OUTPUT_HTML_BASENAME%.html}.prep.html"
  local prep_output_html="${WORKDIR}/${prep_output_html_basename}"

  echo "[prep] No local top-hits CSV is available yet. Running a prep-only SAS pass to export the hit list before batching..."
  render_gtf_runner \
    "${prep_run_sas_rendered}" \
    "${prep_output_html_basename}" \
    "${LOCAL_TOP_HITS_CSV_BASENAME}" \
    "" \
    "1"

  run_oda_helper \
    --file "${prep_run_sas_rendered}" \
    --output-prefix "${RUN_PREFIX}_prep"

  run_oda_helper \
    --download-file "~/${LOCAL_TOP_HITS_CSV_BASENAME}" \
    --download-local-path "${CSV_OUT}" \
    --output-prefix "download_local_hits_with_gtf_prep_csv_${stamp}" || true

  rm -f "${prep_run_sas_rendered}" "${prep_output_html}"
  [[ -s "${CSV_OUT}" ]]
}

write_batched_html_index() {
  local completed_parts="$1"
  local idx
  : > "${HTML_OUT}"
  echo '<!doctype html>' >> "${HTML_OUT}"
  echo '<html><head><meta charset="utf-8"><title>Local Top-Hit GTF Plot Batches</title></head>' >> "${HTML_OUT}"
  echo '<body style="font-family:Arial,Helvetica,sans-serif;margin:24px">' >> "${HTML_OUT}"
  echo '<h1 style="font-size:20px">Local Manhattan and GTF plot for top hits</h1>' >> "${HTML_OUT}"
  echo "<p>Top hits were split into ${completed_parts} batch runs to keep the SAS ODA job size manageable.</p>" >> "${HTML_OUT}"
  echo '<ul>' >> "${HTML_OUT}"
  idx=1
  while [[ "${idx}" -le "${completed_parts}" ]]; do
    echo "<li><a href=\"${OUTPUT_HTML_BASENAME%.html}_part${idx}.html\">Part ${idx}</a></li>" >> "${HTML_OUT}"
    idx=$((idx+1))
  done
  echo '</ul>' >> "${HTML_OUT}"
  if [[ -s "${CSV_OUT}" ]]; then
    echo "<p><a href=\"$(basename "${CSV_OUT}")\">Combined top-hit CSV</a></p>" >> "${HTML_OUT}"
  fi
  echo '</body></html>' >> "${HTML_OUT}"
}

support_uploads_performed=0
upload_support_file_if_needed() {
  local local_path="$1"
  local remote_basename="$2"
  local label="$3"
  local upload_path
  if [[ ! -f "${local_path}" && ! -s "${local_path}" ]]; then
    return 1
  fi
  if [[ "${label}" == "[dynamic]" && "${platform_is_linux}" == "1" && "${FORCE_DYNAMIC_GTF_SUPPORT_UPLOADS_ON_LINUX}" == "1" && "${remote_basename}" == "${REMOTE_GTF_BASENAME}" ]]; then
    if [[ "${support_uploads_performed}" -eq 0 ]]; then
      echo "[1-2c/5] Uploading local-top-hit GTF support files to SAS ODA..."
    fi
    echo "${label} Linux portable mode: uploading $(basename "${local_path}") to SAS ODA home as ${remote_basename} without remote size preflight..."
    upload_path="${local_path}"
    if [[ "$(basename "${local_path}")" != "${remote_basename}" ]]; then
      mkdir -p "${WORKDIR}/.oda_upload_aliases"
      upload_path="${WORKDIR}/.oda_upload_aliases/${remote_basename}"
      cp -f "${local_path}" "${upload_path}"
    fi
    run_oda_helper \
      --upload-file "${upload_path}" \
      --output-prefix "upload_local_hits_with_gtf_support_${stamp}"
    if [[ "${upload_path}" != "${local_path}" ]]; then
      rm -f "${upload_path}"
    fi
    support_uploads_performed=1
    return 0
  fi
  if remote_home_file_matches_local_size "${local_path}" "${remote_basename}"; then
    echo "[bulk] Reusing existing remote file in SAS ODA home: ${remote_basename}"
    return 1
  fi
  if [[ "${support_uploads_performed}" -eq 0 ]]; then
    echo "[1-2c/5] Uploading local-top-hit GTF support files to SAS ODA..."
  fi
  upload_home_file_if_needed "${label}" "${local_path}" "${remote_basename}" "upload_local_hits_with_gtf_support_${stamp}"
  support_uploads_performed=1
  return 0
}

upload_dynamic_local_gtf_inputs_if_needed() {
  local dynamic_uploads_performed=0
  local -a dynamic_bulk_upload_args=()
  local -a dynamic_alias_paths=()

  if [[ "${ASSUME_REMOTE_GTF_DYNAMIC_INPUTS_READY}" == "1" ]]; then
    echo "[1-2c/5] Assuming dynamic local-top-hit GTF inputs are already present in SAS ODA home."
    return 0
  fi

  if [[ "${platform_is_linux}" == "1" && "${FORCE_DYNAMIC_GTF_SUPPORT_UPLOADS_ON_LINUX}" == "1" ]]; then
    queue_dynamic_upload_for_linux() {
      local local_path="$1"
      local remote_basename="$2"
      local upload_path="${local_path}"
      [[ -s "${local_path}" ]] || return 0
      echo "[dynamic] Linux portable mode: queueing $(basename "${local_path}") for upload to SAS ODA home as ${remote_basename}..."
      if [[ "$(basename "${local_path}")" != "${remote_basename}" ]]; then
        mkdir -p "${WORKDIR}/.oda_upload_aliases"
        upload_path="${WORKDIR}/.oda_upload_aliases/${remote_basename}"
        cp -f "${local_path}" "${upload_path}"
        dynamic_alias_paths+=("${upload_path}")
      fi
      dynamic_bulk_upload_args+=(--upload-file "${upload_path}")
      dynamic_uploads_performed=1
    }

    if [[ -s "${CSV_OUT}" && -n "${LOCAL_TOP_HITS_INPUT_CSV_BASENAME:-}" ]]; then
      queue_dynamic_upload_for_linux "${CSV_OUT}" "${LOCAL_TOP_HITS_INPUT_CSV_BASENAME}"
    fi
    if [[ -s "${GET_GTF_MACRO_UPLOAD}" ]]; then
      queue_dynamic_upload_for_linux "${GET_GTF_MACRO_UPLOAD}" "${GET_GTF_MACRO_BASENAME}"
    fi
    if [[ -s "${LOCAL_GTF_SUBSET_GZ}" ]]; then
      queue_dynamic_upload_for_linux "${LOCAL_GTF_SUBSET_GZ}" "${REMOTE_GTF_BASENAME}"
    fi

    if [[ "${dynamic_uploads_performed}" -eq 1 ]]; then
      if [[ "${support_uploads_performed}" -eq 0 ]]; then
        echo "[1-2c/5] Uploading local-top-hit GTF support files to SAS ODA..."
      fi
      run_oda_helper "${dynamic_bulk_upload_args[@]}" --output-prefix "upload_local_hits_with_gtf_support_${stamp}"
      support_uploads_performed=1
      if [[ "${#dynamic_alias_paths[@]}" -gt 0 ]]; then
        rm -f "${dynamic_alias_paths[@]}"
      fi
    fi
    return 0
  fi

  if [[ -s "${CSV_OUT}" && -n "${LOCAL_TOP_HITS_INPUT_CSV_BASENAME:-}" ]]; then
    upload_support_file_if_needed "${CSV_OUT}" "${LOCAL_TOP_HITS_INPUT_CSV_BASENAME}" "[dynamic]" && dynamic_uploads_performed=1 || true
  fi

  if [[ -s "${GET_GTF_MACRO_UPLOAD}" ]]; then
    upload_support_file_if_needed "${GET_GTF_MACRO_UPLOAD}" "${GET_GTF_MACRO_BASENAME}" "[dynamic]" && dynamic_uploads_performed=1 || true
  fi

  if [[ -s "${LOCAL_GTF_SUBSET_GZ}" ]]; then
    upload_support_file_if_needed "${LOCAL_GTF_SUBSET_GZ}" "${REMOTE_GTF_BASENAME}" "[dynamic]" && dynamic_uploads_performed=1 || true
  fi

  return 0
}

if [[ "${ASSUME_REMOTE_GTF_SUPPORT_READY}" == "1" ]]; then
  echo "[1-2c/5] Assuming static local-top-hit GTF support macros are already present in SAS ODA home."
else
  if [[ -f "${SNP_LOCAL_MACRO_SAS}" ]]; then
    upload_support_file_if_needed "${SNP_LOCAL_MACRO_SAS}" "$(basename "${SNP_LOCAL_MACRO_SAS}")" "[bulk]" || true
  else
    echo "WARNING: Local SNP_Local_Manhattan_With_GTF macro not found: ${SNP_LOCAL_MACRO_SAS}" >&2
  fi
  if [[ -f "${MAP_GRP_ASSOC_MACRO_SAS}" ]]; then
    upload_support_file_if_needed "${MAP_GRP_ASSOC_MACRO_SAS}" "$(basename "${MAP_GRP_ASSOC_MACRO_SAS}")" "[bulk]" || true
  fi
  if [[ -f "${MULT_GSCATTER_GENE_MACRO_SAS}" ]]; then
    upload_support_file_if_needed "${MULT_GSCATTER_GENE_MACRO_SAS}" "$(basename "${MULT_GSCATTER_GENE_MACRO_SAS}")" "[bulk]" || true
  fi
  if [[ -f "${ADJ_CLOSE_GENE_GRP_MACRO_SAS}" ]]; then
    upload_support_file_if_needed "${ADJ_CLOSE_GENE_GRP_MACRO_SAS}" "$(basename "${ADJ_CLOSE_GENE_GRP_MACRO_SAS}")" "[bulk]" || true
  fi
  upload_support_file_if_needed "${PATCHED_LATTICE_MACRO_SAS}" "$(basename "${PATCHED_LATTICE_MACRO_SAS}")" "[bulk]" || true
  upload_support_file_if_needed "${TOP_HIT_DIST_MACRO_SAS}" "$(basename "${TOP_HIT_DIST_MACRO_SAS}")" "[bulk]" || true
  if [[ "${support_uploads_performed}" -eq 0 ]]; then
    echo "[1-2c/5] Reusing all local-top-hit GTF support files already present in SAS ODA home."
  fi
fi

upload_dynamic_local_gtf_inputs_if_needed

if [[ "${SKIP_DATA_UPLOAD}" == "1" ]]; then
  echo "[3/5] Reusing already-uploaded gzipped local-top-hit subset in SAS ODA: ${REMOTE_DATA_BASENAME}"
elif [[ "${KEEP_REMOTE_PLOT_DATA}" == "1" ]] && remote_data_matches_local_size; then
  echo "[3/5] Keeping and reusing existing remote local-top-hit subset in SAS ODA: ${REMOTE_DATA_BASENAME}"
else
  echo "[3/5] Uploading gzipped local-top-hit subset to SAS ODA..."
  if remote_data_known_to_oda && ! remote_data_matches_local_size; then
    echo "[repair] Existing remote file has the wrong size and will be replaced."
    delete_partial_remote_data
  fi
  # Large GWAS subset uploads are more reliable through a one-shot ODA
  # connection than through the persistent session server.
  upload_data_with_integrity_check
fi

BATCH_SIZE="${LOCAL_MAX_HITS_PER_FIG:-4}"
if [[ ! -s "${CSV_OUT}" ]]; then
  generate_top_hits_csv_for_batching || true
fi
if [[ -f "${CSV_OUT}" && -s "${CSV_OUT}" ]]; then
  total_hits="$(($(wc -l < "${CSV_OUT}") - 1))"
  if [[ "${total_hits}" -gt "${BATCH_SIZE}" ]]; then
    echo "[batch] Detected ${total_hits} local-top-hit rows in ${CSV_OUT}. Splitting into batches of ${BATCH_SIZE} hits..."
    header_line="$(head -n1 "${CSV_OUT}")"
    mkdir -p "${RUN_LOG_DIR}"
    rm -f "${RUN_LOG_DIR}"/hits_part_*
    tail -n +2 "${CSV_OUT}" | split -l "${BATCH_SIZE}" - "${RUN_LOG_DIR}/hits_part_"

    part=1
    for partfile in "${RUN_LOG_DIR}"/hits_part_*; do
      [[ -f "${partfile}" ]] || continue
      batch_remote_csv="${LOCAL_TOP_HITS_CSV_BASENAME%.csv}_part${part}.csv"
      batch_csv="${WORKDIR}/${batch_remote_csv}"
      batch_output_html_basename="${OUTPUT_HTML_BASENAME%.html}_part${part}.html"
      batch_output_html="${WORKDIR}/${batch_output_html_basename}"
      batch_raw_html="${batch_output_html%.html}.sasraw.html"
      batch_png="${batch_output_html%.html}.png"
      batch_run_sas_rendered="${WORKDIR}/run_sas_oda_local_top_hits_with_gtf.${stamp}.part${part}.sas"
      batch_run_prefix="${RUN_PREFIX}_part${part}"
      batch_run_log_dir="${WORKDIR}/${batch_run_prefix}"
      batch_run_log_file="${batch_run_log_dir}/output.html.info.txt"
      batch_submit_attempt=1
      batch_include_preflight_enabled="${BATCH_INCLUDE_PREFLIGHT_ENABLED}"
      batch_include_preflight_refresh_remote="${BATCH_INCLUDE_PREFLIGHT_REFRESH_REMOTE}"
      if [[ "${part}" -eq 1 ]]; then
        batch_include_preflight_enabled="${BATCH1_INCLUDE_PREFLIGHT_ENABLED}"
        batch_include_preflight_refresh_remote="${BATCH1_INCLUDE_PREFLIGHT_REFRESH_REMOTE}"
      fi

      printf '%s\n' "${header_line}" > "${batch_csv}"
      cat "${partfile}" >> "${batch_csv}"

      upload_home_file_if_needed \
        "[batch ${part}]" \
        "${batch_csv}" \
        "${batch_remote_csv}" \
        "upload_top_hits_batch_${stamp}_${part}"

      render_gtf_runner \
        "${batch_run_sas_rendered}" \
        "${batch_output_html_basename}" \
        "${batch_remote_csv}" \
        "${batch_remote_csv}" \
        "0"

      while :; do
        rm -f "${batch_run_log_file}"
        batch_submit_rc=0
        if run_gtf_submit_for_file "${batch_run_sas_rendered}" "${batch_run_prefix}" "${batch_include_preflight_enabled}" "${batch_include_preflight_refresh_remote}"; then
          batch_submit_rc=0
        else
          batch_submit_rc=$?
          echo "[batch ${part}] GTF SAS submit attempt ${batch_submit_attempt} exited with status ${batch_submit_rc}."
        fi
        if [[ "${batch_submit_rc}" -eq 0 ]] && ! gtf_submit_needs_retry_for_log "${batch_run_log_file}"; then
          break
        fi
        if [[ "${batch_submit_attempt}" -ge "${GTF_SUBMIT_MAX_ATTEMPTS}" ]]; then
          if [[ "${batch_submit_rc}" -ne 0 ]]; then
            echo "ERROR: Batch ${part} GTF SAS submit failed with status ${batch_submit_rc} after ${batch_submit_attempt} attempt(s)." >&2
          else
            echo "ERROR: Batch ${part} GTF SAS submit log still looks incomplete after ${batch_submit_attempt} attempt(s): ${batch_run_log_file}" >&2
          fi
          exit 1
        fi
        next_attempt=$((batch_submit_attempt + 1))
        if [[ "${batch_submit_rc}" -ne 0 ]]; then
          echo "[batch ${part}] GTF SAS submit attempt ${batch_submit_attempt} failed or timed out; retrying attempt ${next_attempt}/${GTF_SUBMIT_MAX_ATTEMPTS} after ${GTF_SUBMIT_RETRY_SLEEP_SECONDS}s..."
        else
          echo "[batch ${part}] GTF SAS submit attempt ${batch_submit_attempt} looked incomplete; retrying attempt ${next_attempt}/${GTF_SUBMIT_MAX_ATTEMPTS} after ${GTF_SUBMIT_RETRY_SLEEP_SECONDS}s..."
        fi
        sleep "${GTF_SUBMIT_RETRY_SLEEP_SECONDS}"
        batch_submit_attempt=$next_attempt
      done

      if gtf_log_has_terminal_failure "${batch_run_log_file}"; then
        echo "ERROR: Batch ${part} SAS log shows a terminal failure before HTML download: ${batch_run_log_file}" >&2
        exit 1
      fi

      recover_html_from_submit_artifacts_for_logdir "${batch_run_log_dir}" "${batch_output_html}" || true
      if [[ ! -s "${batch_output_html}" ]]; then
        oda_download_many_with_timeout \
          "${ODA_RESULT_DOWNLOAD_TIMEOUT_SECONDS}" \
          "${ODA_RESULT_DOWNLOAD_TIMEOUT_GRACE_SECONDS}" \
          "download_local_hits_with_gtf_html_${stamp}_part${part}" \
          --download-file "~/${batch_output_html_basename}" \
          --download-local-path "${batch_output_html}" || true
      fi

      recover_html_from_submit_artifacts_for_logdir "${batch_run_log_dir}" "${batch_output_html}" || true
      batch_remote_png_path="$(extract_remote_png_path_from_log_file "${batch_run_log_file}" || true)"
      download_remote_png_to_path_if_reported \
        "${batch_remote_png_path}" \
        "${batch_png}" \
        "download_local_hits_with_gtf_png_${stamp}_part${part}" || true
      build_completed_html_from_png_assets_if_available \
        "${batch_output_html}" \
        "${batch_png}" \
        "${batch_raw_html}" \
        "${batch_csv}" \
        "Local Manhattan and GTF plot for top hits (Part ${part})" \
        "Local top-hit GTF plot (Part ${part})" || true

      if [[ ! -s "${batch_output_html}" ]]; then
        echo "ERROR: Batch ${part} expected downloaded HTML was not created or is empty: ${batch_output_html}" >&2
        echo "ERROR: Check the saved batch log: ${batch_run_log_file}" >&2
        exit 1
      fi

      echo "[batch ${part}] Wrote batch HTML: ${batch_output_html}"

      if [[ "${CLEAN_ODA_INPUT}" == "1" ]]; then
        run_oda_helper \
          --delete-file "${batch_remote_csv}" \
          --output-prefix "cleanup_top_hits_batch_${stamp}_${part}" >/dev/null 2>&1 || true
      fi

      part=$((part+1))
    done

    completed_parts=$((part-1))
    write_batched_html_index "${completed_parts}"

    if [[ "${KEEP_REMOTE_PLOT_DATA}" == "1" ]]; then
      echo "[cleanup] Keeping uploaded gz input in SAS ODA because KEEP_REMOTE_PLOT_DATA=${KEEP_REMOTE_PLOT_DATA}."
    elif [[ "${CLEAN_ODA_INPUT}" == "1" ]]; then
      echo "[cleanup] Removing uploaded gz input from SAS ODA to save space..."
      oda_delete_many \
        "cleanup_local_hits_with_gtf_input_${stamp}" \
        --delete-file "${REMOTE_DATA_BASENAME}" || true
    else
      echo "[cleanup] Keeping uploaded gz input in SAS ODA because CLEAN_ODA_INPUT=${CLEAN_ODA_INPUT}."
    fi

    cleanup_remote_generated_outputs
    cleanup_remote_local_gtf_subset
    echo "Done."
    echo "Downloaded HTML index: ${HTML_OUT}"
    if [[ "${OPEN_RESULT}" == "1" ]]; then
      open_html_result "${HTML_OUT}" || true
    else
      echo "Not opening result because OPEN_RESULT=${OPEN_RESULT}."
    fi
    exit 0
  fi
fi

echo "[4/5] Running SAS local top-hits gene-track plot..."
clear_stale_remote_expected_outputs
GTF_SUBMIT_MAX_ATTEMPTS="${GTF_SUBMIT_MAX_ATTEMPTS:-5}"
GTF_SUBMIT_RETRY_SLEEP_SECONDS="${GTF_SUBMIT_RETRY_SLEEP_SECONDS:-10}"
gtf_submit_attempt=1
while :; do
  rm -f "${RUN_LOG_FILE}"
  gtf_submit_rc=0
  if run_gtf_submit; then
    gtf_submit_rc=0
  else
    gtf_submit_rc=$?
    echo "[4b/5] GTF SAS submit attempt ${gtf_submit_attempt} exited with status ${gtf_submit_rc}."
  fi
  if [[ "${gtf_submit_rc}" -eq 0 ]] && ! gtf_submit_needs_retry; then
    break
  fi
  if [[ "${gtf_submit_attempt}" -ge "${GTF_SUBMIT_MAX_ATTEMPTS}" ]]; then
    if [[ "${gtf_submit_rc}" -ne 0 ]]; then
      echo "ERROR: GTF SAS submit failed with status ${gtf_submit_rc} after ${gtf_submit_attempt} attempt(s)." >&2
    else
      echo "ERROR: GTF SAS submit log still looks incomplete after ${gtf_submit_attempt} attempts: ${RUN_LOG_FILE}" >&2
    fi
    exit 1
  fi
  next_attempt=$((gtf_submit_attempt + 1))
  if [[ "${gtf_submit_rc}" -ne 0 ]]; then
    echo "[4b/5] GTF SAS submit attempt ${gtf_submit_attempt} failed or timed out; retrying attempt ${next_attempt}/${GTF_SUBMIT_MAX_ATTEMPTS} after ${GTF_SUBMIT_RETRY_SLEEP_SECONDS}s..."
  else
    echo "[4b/5] GTF SAS submit attempt ${gtf_submit_attempt} looked incomplete; retrying attempt ${next_attempt}/${GTF_SUBMIT_MAX_ATTEMPTS} after ${GTF_SUBMIT_RETRY_SLEEP_SECONDS}s..."
  fi
  sleep "${GTF_SUBMIT_RETRY_SLEEP_SECONDS}"
  gtf_submit_attempt="${next_attempt}"
done

if gtf_log_has_terminal_failure "${RUN_LOG_FILE}"; then
  echo "ERROR: GTF SAS log shows a terminal failure before HTML download: ${RUN_LOG_FILE}" >&2
  exit 1
fi

echo "[5/5] Downloading self-contained HTML result..."
rm -f "${HTML_OUT}" "${RAW_HTML_OUT}" "${PNG_OUT}"
recover_html_from_submit_artifacts || true
if [[ -s "${HTML_OUT}" ]]; then
  echo "[recover] Using the local SAS HTML artifact saved by the submit helper; remote HTML download will be skipped unless another file is still missing."
  if [[ ! -s "${RAW_HTML_OUT}" ]]; then
    build_completed_html_from_png_if_available || true
  fi
fi
recover_csv_from_existing_local_copy || true
download_outputs_args=()
if [[ ! -s "${HTML_OUT}" ]]; then
  download_outputs_args+=(--download-file "~/${OUTPUT_HTML_BASENAME}" --download-local-path "${HTML_OUT}")
fi
if [[ ! -s "${CSV_OUT}" ]]; then
  download_outputs_args+=(--download-file "~/${LOCAL_TOP_HITS_CSV_BASENAME}" --download-local-path "${CSV_OUT}")
fi
if [[ ${#download_outputs_args[@]} -gt 0 ]]; then
  oda_download_many_with_timeout \
    "${ODA_RESULT_DOWNLOAD_TIMEOUT_SECONDS}" \
    "${ODA_RESULT_DOWNLOAD_TIMEOUT_GRACE_SECONDS}" \
    "download_local_hits_with_gtf_outputs_${stamp}" \
    "${download_outputs_args[@]}" || true
fi

recover_html_from_submit_artifacts || true
if [[ -s "${HTML_OUT}" && ! -s "${RAW_HTML_OUT}" ]]; then
  build_completed_html_from_png_if_available || true
fi
recover_csv_from_existing_local_copy || true
if delivered_gtf_artifact_ready; then
  echo "[recover] Final delivered local-GTF artifact is already available locally; skipping remote PNG recovery."
else
  remote_png_path="$(extract_remote_png_path_from_run_log || true)"
  download_remote_png_if_reported "${remote_png_path}" || true
  build_completed_html_from_png_if_available || true
fi

if [[ ! -s "${HTML_OUT}" ]]; then
  echo "ERROR: Expected downloaded HTML was not created or is empty: ${HTML_OUT}" >&2
  exit 1
fi

if ! delivered_gtf_artifact_ready; then
  echo "ERROR: Downloaded HTML exists but contains no rendered plot content: ${HTML_OUT}" >&2
  echo "ERROR: The SAS job likely failed before producing the final figure. Check the SAS log saved under ${RUN_LOG_DIR}." >&2
  exit 1
fi

echo "Verified HTML: ${HTML_OUT} ($(wc -c < "${HTML_OUT}") bytes)"

cleanup_remote_generated_outputs

if [[ "${KEEP_REMOTE_PLOT_DATA}" == "1" ]]; then
  echo "[cleanup] Keeping uploaded gz input in SAS ODA because KEEP_REMOTE_PLOT_DATA=${KEEP_REMOTE_PLOT_DATA}."
elif [[ "${CLEAN_ODA_INPUT}" == "1" ]]; then
  echo "[cleanup] Removing uploaded gz input from SAS ODA to save space..."
  oda_delete_many \
    "cleanup_local_hits_with_gtf_input_${stamp}" \
    --delete-file "${REMOTE_DATA_BASENAME}"
else
  echo "[cleanup] Keeping uploaded gz input in SAS ODA because CLEAN_ODA_INPUT=${CLEAN_ODA_INPUT}."
fi

cleanup_remote_local_gtf_subset

echo "Done."
echo "Downloaded HTML: ${HTML_OUT}"

if [[ "${OPEN_RESULT}" == "1" ]]; then
  open_html_result "${HTML_OUT}" || true
else
  echo "Not opening result because OPEN_RESULT=${OPEN_RESULT}."
fi
