#!/usr/bin/env bash
set -euo pipefail

DEPS_DIR="${DEPS_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)}"
WORKDIR="${WORKDIR:-$(cd "${DEPS_DIR}/.." && pwd -P)}"
RUNNER_CONFIG_JSON="${RUNNER_CONFIG_JSON:-}"
CALLER_TARGET_SNP="${TARGET_SNP-__UNSET__}"
CALLER_OUTPUT_HTML_BASENAME="${OUTPUT_HTML_BASENAME-__UNSET__}"
CALLER_LOCAL_WINDOW_BP="${LOCAL_WINDOW_BP-__UNSET__}"
CALLER_GTF_LABEL_SNPS="${GTF_LABEL_SNPS-__UNSET__}"
CALLER_LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES="${LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES-__UNSET__}"
CALLER_OPEN_RESULT="${OPEN_RESULT-__UNSET__}"
CALLER_DATA_GZ="${DATA_GZ-__UNSET__}"
CALLER_REMOTE_DATA_BASENAME="${REMOTE_DATA_BASENAME-__UNSET__}"
if [[ -n "${RUNNER_CONFIG_JSON}" ]]; then
  eval "$("perl" "${DEPS_DIR}/emit_diff_gwas_runner_env.pl" --config "${RUNNER_CONFIG_JSON}")"
fi
if [[ "${CALLER_TARGET_SNP}" != "__UNSET__" ]]; then
  TARGET_SNP="${CALLER_TARGET_SNP}"
fi
if [[ "${CALLER_OUTPUT_HTML_BASENAME}" != "__UNSET__" ]]; then
  OUTPUT_HTML_BASENAME="${CALLER_OUTPUT_HTML_BASENAME}"
fi
if [[ "${CALLER_LOCAL_WINDOW_BP}" != "__UNSET__" ]]; then
  LOCAL_WINDOW_BP="${CALLER_LOCAL_WINDOW_BP}"
fi
if [[ "${CALLER_GTF_LABEL_SNPS}" != "__UNSET__" ]]; then
  GTF_LABEL_SNPS="${CALLER_GTF_LABEL_SNPS}"
fi
if [[ "${CALLER_LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES}" != "__UNSET__" ]]; then
  LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES="${CALLER_LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES}"
fi
if [[ "${CALLER_OPEN_RESULT}" != "__UNSET__" ]]; then
  OPEN_RESULT="${CALLER_OPEN_RESULT}"
fi
if [[ "${CALLER_DATA_GZ}" != "__UNSET__" ]]; then
  DATA_GZ="${CALLER_DATA_GZ}"
fi
if [[ "${CALLER_REMOTE_DATA_BASENAME}" != "__UNSET__" ]]; then
  REMOTE_DATA_BASENAME="${CALLER_REMOTE_DATA_BASENAME}"
fi
PROJECT_TAG="${PROJECT_TAG:-PGC_SCZ}"
DEFAULT_SOURCE_LONG_GZ="/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects.stdized.tsv.gz"
SOURCE_LONG_GZ="${SOURCE_LONG_GZ:-${DEFAULT_SOURCE_LONG_GZ}}"
DATA_GZ="${DATA_GZ:-}"
REMOTE_DATA_BASENAME="${REMOTE_DATA_BASENAME:-}"
LOCAL_WIDE_HELPER="${LOCAL_WIDE_HELPER:-${DEPS_DIR}/extract_single_snp_wide_diff_gwas.pl}"
EXTRACTOR_CONFIG_JSON="${EXTRACTOR_CONFIG_JSON:-}"
SCHEMA_CONFIG_JSON="${SCHEMA_CONFIG_JSON:-${EXTRACTOR_CONFIG_JSON:-${WORKDIR}/configs/preset_pgc_scz_sex_diff.json}}"
SCHEMA_INCLUDE_HELPER="${SCHEMA_INCLUDE_HELPER:-${DEPS_DIR}/generate_sas_wide_import_include.pl}"
GTF_IMPORT_INCLUDE_HELPER="${GTF_IMPORT_INCLUDE_HELPER:-${DEPS_DIR}/generate_sas_gtf_import_include.pl}"
GTF_SUBSET_HELPER="${GTF_SUBSET_HELPER:-${DEPS_DIR}/extract_gencode_gtf_subset.pl}"
RENDER_SAS_HELPER="${RENDER_SAS_HELPER:-${DEPS_DIR}/render_sas_template.pl}"
LOCAL_WIDE_MANIFEST=""
LOCAL_WIDE_AUTOGEN=0

TARGET_SNP="${TARGET_SNP:-}"
if [[ -z "${TARGET_SNP}" ]]; then
  echo "ERROR: TARGET_SNP is required, for example TARGET_SNP=rs42067 ./run_sas_oda_single_snp_with_gtf_download_html.sh" >&2
  exit 2
fi
SAFE_TARGET_SNP="$(printf '%s' "${TARGET_SNP}" | tr -c 'A-Za-z0-9._-' '_')"
OUTPUT_HTML_BASENAME="${OUTPUT_HTML_BASENAME:-}"
LOCAL_WINDOW_BP="${LOCAL_WINDOW_BP:-1e7}"
GTF_LABEL_SNPS="${GTF_LABEL_SNPS:-${TARGET_SNP}}"
GTF_LABEL_SNPS="${GTF_LABEL_SNPS//,/ }"

# A single-SNP context run should not silently reuse the generic genome-wide
# wide subset or generic local-top-hits HTML basename from the multi-hit runner
# config. Fall back to SNP-specific helper output instead.
if [[ -n "${DATA_GZ}" && "$(basename "${DATA_GZ}")" == *.stdized.wide_beta_se_p_p_lt_0p05.final.tsv.gz ]]; then
  echo "[prep] Ignoring generic DATA_GZ from runner config for single-SNP mode: ${DATA_GZ}"
  DATA_GZ=""
fi
if [[ -n "${OUTPUT_HTML_BASENAME}" && "${OUTPUT_HTML_BASENAME}" == *_local_top_hits_with_gtf*.html ]]; then
  echo "[prep] Ignoring generic OUTPUT_HTML_BASENAME from runner config for single-SNP mode: ${OUTPUT_HTML_BASENAME}"
  OUTPUT_HTML_BASENAME=""
fi

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
    DEFAULT_GTF_LOCAL_DSD=""
    DEFAULT_GTF_GZ_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz"
    ;;
esac

GTF_DSD="${GTF_DSD:-${DEFAULT_GTF_DSD}}"
FM_LIBPATH="${FM_LIBPATH:-/home/cheng.zhong.shan/my_shared_file_links/cheng.zhong.shan/F_vs_M_Covid19_Hosp}"
GTF_LOCAL_DSD="${GTF_LOCAL_DSD:-${DEFAULT_GTF_LOCAL_DSD}}"
GTF_GZ_URL="${GTF_GZ_URL:-${DEFAULT_GTF_GZ_URL}}"
SNP_LOCAL_MACRO_SAS="${SNP_LOCAL_MACRO_SAS:-${DEPS_DIR}/SNP_Local_Manhattan_With_GTF.sas}"
PATCHED_LATTICE_MACRO_SAS="${PATCHED_LATTICE_MACRO_SAS:-${DEPS_DIR}/Lattice_gscatter_over_bed_track.sas}"
MAP_GRP_ASSOC_MACRO_SAS="${MAP_GRP_ASSOC_MACRO_SAS:-${DEPS_DIR}/map_grp_assoc2gene4covidsexgwas.sas}"
MULT_GSCATTER_GENE_MACRO_SAS="${MULT_GSCATTER_GENE_MACRO_SAS:-${DEPS_DIR}/Multgscatter_with_gene_exons.sas}"
ADJ_CLOSE_GENE_GRP_MACRO_SAS="${ADJ_CLOSE_GENE_GRP_MACRO_SAS:-${DEPS_DIR}/adj_grpnum4close_gene_bed_regs.sas}"
GTF_CACHE_DIR="${GTF_CACHE_DIR:-${WORKDIR}/cache/gtf}"

GTF_ASSOC_PVARS="${GTF_ASSOC_PVARS:-ALL_DIFF_P ASN_DIFF_P EUR_DIFF_P}"
GTF_ZSCORE_VARS="${GTF_ZSCORE_VARS:-ALL_DIFF_Z ASN_DIFF_Z EUR_DIFF_Z}"
GTF_LABELS="${GTF_LABELS:-All_Diff Asian_Diff European_Diff}"

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
GTF_DIST2SNP="${GTF_DIST2SNP:-500000}"
GTF_DESIGN_WIDTH="${GTF_DESIGN_WIDTH:-950}"
GTF_DESIGN_HEIGHT="${GTF_DESIGN_HEIGHT:-1000}"
GTF_DIST2SEP_GENES="${GTF_DIST2SEP_GENES:-100000}"
GTF_SHIFT_TEXT_YVAL="${GTF_SHIFT_TEXT_YVAL:-0.2}"
# Match the multi-hit local-GTF default so single-SNP manuscript reruns keep
# a larger, easier-to-read bottom gene track.
GTF_PCT4NEG_Y="${GTF_PCT4NEG_Y:-1.4}"
GTF_ADJVAL4HEADER="${GTF_ADJVAL4HEADER:--0.6}"
GTF_YOFFSET4TEXTLABELS="${GTF_YOFFSET4TEXTLABELS:-1.5}"
LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES="${LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES:-0}"

RUN_SAS_TEMPLATE="${DEPS_DIR}/run_sas_oda_single_snp_with_gtf.sas"
SESSION_ID="${SESSION_ID:-mysession}"
USE_PERSISTENT_SESSION="${USE_PERSISTENT_SESSION:-0}"
CLEAN_ODA_INPUT="${CLEAN_ODA_INPUT:-1}"
CLEAN_ODA_OUTPUT="${CLEAN_ODA_OUTPUT:-1}"
CLEAN_ODA_MACROS="${CLEAN_ODA_MACROS:-1}"
CLEAN_LOCAL_AUTOGEN="${CLEAN_LOCAL_AUTOGEN:-1}"
OPEN_RESULT="${OPEN_RESULT:-1}"
platform_is_linux=0
platform_is_cygwin=0
if [[ "$(uname -s)" == "Linux" ]]; then
  platform_is_linux=1
elif [[ "$(uname -s)" == CYGWIN* ]]; then
  platform_is_cygwin=1
fi
if [[ -z "${SINGLE_SNP_GTF_SUBMIT_TIMEOUT_SECONDS:-}" ]]; then
  if [[ "${platform_is_linux}" == "1" || "${platform_is_cygwin}" == "1" ]]; then
    SINGLE_SNP_GTF_SUBMIT_TIMEOUT_SECONDS=3600
  else
    SINGLE_SNP_GTF_SUBMIT_TIMEOUT_SECONDS=1200
  fi
fi
if [[ -z "${SINGLE_SNP_GTF_SUBMIT_TIMEOUT_GRACE_SECONDS:-}" ]]; then
  SINGLE_SNP_GTF_SUBMIT_TIMEOUT_GRACE_SECONDS=30
fi
if [[ -z "${ODA_HELPER_TIMEOUT_SECONDS:-}" ]]; then
  if [[ "${platform_is_linux}" == "1" || "${platform_is_cygwin}" == "1" ]]; then
    ODA_HELPER_TIMEOUT_SECONDS=1200
  else
    ODA_HELPER_TIMEOUT_SECONDS=300
  fi
fi
if [[ -z "${ODA_HELPER_TIMEOUT_GRACE_SECONDS:-}" ]]; then
  ODA_HELPER_TIMEOUT_GRACE_SECONDS=20
fi
if [[ -z "${ODA_DELETE_TIMEOUT_SECONDS:-}" ]]; then
  if [[ "${platform_is_linux}" == "1" || "${platform_is_cygwin}" == "1" ]]; then
    ODA_DELETE_TIMEOUT_SECONDS=180
  else
    ODA_DELETE_TIMEOUT_SECONDS="${ODA_HELPER_TIMEOUT_SECONDS}"
  fi
fi
if [[ -z "${ODA_DELETE_TIMEOUT_GRACE_SECONDS:-}" ]]; then
  ODA_DELETE_TIMEOUT_GRACE_SECONDS=20
fi
INCLUDE_PREFLIGHT_STANDALONE_DEBUG="${INCLUDE_PREFLIGHT_STANDALONE_DEBUG:-0}"
export SAS_ODA_RUN_TIMEOUT_SECONDS="${SAS_ODA_RUN_TIMEOUT_SECONDS:-${SINGLE_SNP_GTF_SUBMIT_TIMEOUT_SECONDS}}"
export SAS_ODA_RUN_TIMEOUT_GRACE_SECONDS="${SAS_ODA_RUN_TIMEOUT_GRACE_SECONDS:-${SINGLE_SNP_GTF_SUBMIT_TIMEOUT_GRACE_SECONDS}}"
export STANDALONE_INCLUDE_TARGET_DEBUG="${INCLUDE_PREFLIGHT_STANDALONE_DEBUG}"

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

echo "[helper] Include preflight standalone target debug: ${STANDALONE_INCLUDE_TARGET_DEBUG}"

stamp="$(date +%Y%m%d_%H%M%S)"
if [[ -z "${GTF_LOCAL_DSD}" ]]; then
  GTF_LOCAL_DSD="FM.GTFSP_${stamp//_/}"
fi
GWAS_DATASET="SCZMH_${stamp//_/}"
TARGET_HIT_DATASET="THIT_${stamp//_/}"
TARGET_LOCAL_DATASET="TLOC_${stamp//_/}"
RUN_SAS_RENDERED="${WORKDIR}/run_sas_oda_single_snp_with_gtf.${stamp}.sas"
RENDERED_SAS_BASENAME="$(basename "${RUN_SAS_RENDERED}")"
RUN_SAS_RENDERED_SUBMIT="./${RENDERED_SAS_BASENAME}"
IMPORT_BLOCK_RENDERED="${WORKDIR}/auto_wide_import_single_snp.${stamp}.sas"
GTF_IMPORT_BLOCK_RENDERED="${WORKDIR}/auto_gtf_import_single_snp.${stamp}.sas"
LOCAL_GTF_SUBSET="${WORKDIR}/local_gtf_subset_${SAFE_TARGET_SNP}_${stamp}.tsv"
LOCAL_GTF_SUBSET_GZ="${LOCAL_GTF_SUBSET}.gz"
REMOTE_GTF_BASENAME="$(basename "${LOCAL_GTF_SUBSET_GZ}")"
PNG_OUT=""
RAW_HTML_OUT=""
RUN_PREFIX="run_single_snp_with_gtf_${stamp}"
RUN_LOG_DIR="${WORKDIR}/${RUN_PREFIX}"
RUN_LOG_FILE="${RUN_LOG_DIR}/output.html.info.txt"

cleanup_local_artifacts() {
  rm -f "${RUN_SAS_RENDERED}" "${IMPORT_BLOCK_RENDERED}" "${GTF_IMPORT_BLOCK_RENDERED}"
  if [[ "${CLEAN_LOCAL_AUTOGEN}" == "1" ]]; then
    rm -f "${LOCAL_GTF_SUBSET}" "${LOCAL_GTF_SUBSET_GZ}"
  fi
  if [[ "${LOCAL_WIDE_AUTOGEN}" == "1" && "${CLEAN_LOCAL_AUTOGEN}" == "1" ]]; then
    [[ -n "${DATA_GZ}" ]] && rm -f "${DATA_GZ}"
    [[ -n "${LOCAL_WIDE_MANIFEST}" ]] && rm -f "${LOCAL_WIDE_MANIFEST}"
  fi
}

trap cleanup_local_artifacts EXIT

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

remote_home_file_size_bytes() {
  local remote_basename="$1"
  local info_output
  info_output="$(
    run_oda_helper \
      --file-info "~/${remote_basename}" \
      --output-prefix "check_remote_single_snp_home_size_${stamp}" 2>&1 || true
  )"
  printf '%s\n' "${info_output}" | awk -F '\t' '$1=="SIZE"{print $2}' | tail -n 1
}

remote_home_file_exists() {
  local remote_basename="$1"
  local info_output
  info_output="$(
    run_oda_helper \
      --file-info "~/${remote_basename}" \
      --output-prefix "check_remote_single_snp_home_exists_${stamp}" 2>&1 || true
  )"
  [[ "$(printf '%s\n' "${info_output}" | awk -F '\t' '$1=="EXISTS"{print $2}' | tail -n 1 | tr -d '[:space:]')" == "1" ]]
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

delete_partial_remote_home_file() {
  local remote_basename="$1"
  run_oda_helper \
    --delete-file "${remote_basename}" \
    --output-prefix "delete_partial_single_snp_home_${stamp}" >/dev/null 2>&1 || true
}

extract_target_row_from_local_subset() {
  local data_path="$1"
  local target_snp="$2"
  perl -MIO::Uncompress::Gunzip=gunzip,\$GunzipError -e '
    use strict;
    use warnings;
    my ($path, $target) = @ARGV;
    my $fh;
    if ($path =~ /\.gz$/i) {
      $fh = IO::Uncompress::Gunzip->new($path)
        or die "Cannot open $path: $GunzipError\n";
    } else {
      open $fh, q{<}, $path or die "Cannot open $path: $!\n";
    }
    my $header = <$fh>;
    while (my $line = <$fh>) {
      chomp $line;
      $line =~ s/\r$//;
      next if $line =~ /^\s*$/;
      my @f = split /\t/, $line, -1;
      next unless @f >= 5;
      next unless defined $f[4] && $f[4] eq $target;
      print $line;
      exit 0;
    }
    exit 1;
  ' "${data_path}" "${target_snp}" 2>/dev/null
}

manifest_metric_value() {
  local metric="$1"
  local manifest_path="$2"
  [[ -s "${manifest_path}" ]] || return 1
  awk -F '\t' -v key="${metric}" '$1==key{print $2; exit}' "${manifest_path}"
}

log_local_subset_target_summary_if_available() {
  [[ -n "${LOCAL_WIDE_MANIFEST:-}" && -s "${LOCAL_WIDE_MANIFEST}" ]] || return 0
  local found present missing count all_prefixes
  found="$(manifest_metric_value "target_row_found_in_window" "${LOCAL_WIDE_MANIFEST}" || true)"
  present="$(manifest_metric_value "target_row_groups_present" "${LOCAL_WIDE_MANIFEST}" || true)"
  missing="$(manifest_metric_value "target_row_groups_missing" "${LOCAL_WIDE_MANIFEST}" || true)"
  count="$(manifest_metric_value "target_row_group_count" "${LOCAL_WIDE_MANIFEST}" || true)"
  all_prefixes="$(manifest_metric_value "target_row_has_all_prefixes" "${LOCAL_WIDE_MANIFEST}" || true)"
  [[ -n "${found}" ]] || return 0
  echo "[prep] Local single-SNP subset manifest for ${TARGET_SNP}: target_row_found=${found} group_count=${count:-0} present_groups=${present:-none} missing_groups=${missing:-none} all_prefixes=${all_prefixes:-0}"
  if [[ "${all_prefixes}" == "0" && -n "${missing}" ]]; then
    echo "[prep] Note: ${TARGET_SNP} is still a valid centered target even though some subgroup prefix blocks are blank in the emitted wide row: missing=${missing}."
  fi
}

run_log_reports_missing_target_snp() {
  [[ -s "${RUN_LOG_FILE}" ]] || return 1
  grep -Fq 'Target SNP &target_snp was not found in the uploaded GWAS subset.' "${RUN_LOG_FILE}"
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

build_completed_html_from_single_snp_assets_if_available() {
  local raw_html_path="$1"
  local png_path="$2"
  local final_html_path="$3"
  local figure_title="${4:-Single-SNP Local GTF Plot}"
  local image_alt="${5:-Single-SNP local GTF plot}"
  local image_src=""
  local use_iframe=0

  if [[ -s "${png_path}" ]]; then
    annotate_png_with_visible_axis_labels_if_needed "${png_path}" "${GTF_YAXIS_LABEL}" "${GTF_COLORBAR_LABEL}"
    image_src="$(basename "${png_path}")"
  else
    if extract_embedded_png_from_html_path_if_present "${raw_html_path}" "${png_path}"; then
      [[ -s "${png_path}" ]] || return 0
      annotate_png_with_visible_axis_labels_if_needed "${png_path}" "${GTF_YAXIS_LABEL}" "${GTF_COLORBAR_LABEL}"
      image_src="$(basename "${png_path}")"
    else
      image_src="$(extract_embedded_png_data_uri_from_html_path "${raw_html_path}" || true)"
      if [[ -z "${image_src}" ]]; then
        use_iframe=1
      fi
    fi
  fi
  [[ -s "${raw_html_path}" ]] || return 0

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
    echo "<p><a href=\"$(basename "${raw_html_path}")\">Raw SAS HTML output</a></p>"
    echo '</body></html>'
  } > "${final_html_path}"
  if [[ "${use_iframe}" -eq 1 ]]; then
    echo "[recover] Wrapped the raw single-SNP SAS HTML with a visible y-axis label because no standalone PNG could be extracted: ${final_html_path}"
  else
    echo "[recover] Replaced the opened HTML with a figure-first wrapper because a final PNG was generated: ${final_html_path}"
  fi
}

single_snp_submit_needs_retry() {
  [[ ! -s "${RUN_LOG_FILE}" ]] && return 0
  if grep -Eiq 'We failed in getConnection|The application could not log on to the server|server configuration is invalid|SAS process has terminated unexpectedly|SAS submit timed out' "${RUN_LOG_FILE}"; then
    return 0
  fi
  if grep -Eq 'HTML output saved to:|The final figure is put here:|SAS job is completed!' "${RUN_LOG_FILE}"; then
    return 1
  fi
  local bytes
  bytes="$(wc -c < "${RUN_LOG_FILE}")"
  [[ "${bytes}" -lt 500 ]] && return 0
  return 1
}

run_single_snp_submit() {
  if [[ -x /usr/bin/timeout ]]; then
    /usr/bin/timeout --kill-after="${SINGLE_SNP_GTF_SUBMIT_TIMEOUT_GRACE_SECONDS}s" "${SINGLE_SNP_GTF_SUBMIT_TIMEOUT_SECONDS}s" \
      "${ODA_PERL[@]}" \
      --file "${RUN_SAS_RENDERED_SUBMIT}" \
      --output-prefix "${RUN_PREFIX}"
    return $?
  fi

  "${ODA_PERL[@]}" \
    --file "${RUN_SAS_RENDERED_SUBMIT}" \
    --output-prefix "${RUN_PREFIX}"
}

cleanup_remote_generated_outputs() {
  if [[ "${CLEAN_ODA_OUTPUT}" != "1" ]]; then
    echo "[cleanup] Keeping generated remote single-SNP GTF outputs because CLEAN_ODA_OUTPUT=${CLEAN_ODA_OUTPUT}."
    return 0
  fi

  local remote_png_base=""
  echo "[cleanup] Removing generated remote single-SNP GTF outputs from SAS ODA..."
  oda_delete_many "cleanup_single_snp_with_gtf_output_html_${stamp}" --delete-file "${OUTPUT_HTML_BASENAME}" || true
  if [[ -n "${remote_png_path:-}" ]]; then
    remote_png_base="$(basename "${remote_png_path}")"
    if [[ -n "${remote_png_base}" ]]; then
      oda_delete_many "cleanup_single_snp_with_gtf_output_png_${stamp}" --delete-file "${remote_png_base}" || true
    fi
  fi
}

TARGET_CHR=""
TARGET_BP=""

if [[ -z "${DATA_GZ}" ]]; then
  if [[ -z "${EXTRACTOR_CONFIG_JSON}" && ! -s "${SOURCE_LONG_GZ}" ]]; then
    echo "ERROR: SOURCE_LONG_GZ does not exist or is empty: ${SOURCE_LONG_GZ}" >&2
    exit 2
  fi
  if [[ ! -f "${LOCAL_WIDE_HELPER}" ]]; then
    echo "ERROR: Local helper script not found: ${LOCAL_WIDE_HELPER}" >&2
    exit 2
  fi

  echo "[prep] Building a single-SNP wide local subset from the long standardized GWAS..."
  helper_cmd=(perl "${LOCAL_WIDE_HELPER}" --target-snp "${TARGET_SNP}" --window-bp "${LOCAL_WINDOW_BP}" --output-dir "${WORKDIR}")
  if [[ -n "${EXTRACTOR_CONFIG_JSON}" ]]; then
    helper_cmd+=(--config "${EXTRACTOR_CONFIG_JSON}")
  else
    helper_cmd+=(--input "${SOURCE_LONG_GZ}")
  fi
  helper_out="$("${helper_cmd[@]}")"

  DATA_GZ="$(printf '%s\n' "${helper_out}" | awk -F '\t' '$1=="OUTPUT"{print $2}')"
  LOCAL_WIDE_MANIFEST="$(printf '%s\n' "${helper_out}" | awk -F '\t' '$1=="MANIFEST"{print $2}')"
  TARGET_CHR="$(printf '%s\n' "${helper_out}" | awk -F '\t' '$1=="TARGET_CHR"{print $2}')"
  TARGET_BP="$(printf '%s\n' "${helper_out}" | awk -F '\t' '$1=="TARGET_BP"{print $2}')"
  LOCAL_WIDE_AUTOGEN=1

  if [[ -z "${DATA_GZ}" || ! -s "${DATA_GZ}" ]]; then
    echo "ERROR: Helper did not produce a valid local wide subset." >&2
    printf '%s\n' "${helper_out}" >&2
    exit 1
  fi
  if [[ -z "${TARGET_CHR}" || -z "${TARGET_BP}" ]]; then
    echo "ERROR: Helper did not emit target chr/bp metadata." >&2
    printf '%s\n' "${helper_out}" >&2
    exit 1
  fi
else
  TARGET_CHR="${TARGET_CHR:-NA}"
  TARGET_BP="${TARGET_BP:-NA}"
fi

local_target_row="$(
  extract_target_row_from_local_subset "${DATA_GZ}" "${TARGET_SNP}" || true
)"
if [[ -z "${local_target_row}" ]]; then
  echo "ERROR: The local single-SNP wide subset does not contain ${TARGET_SNP}: ${DATA_GZ}" >&2
  if [[ -n "${LOCAL_WIDE_MANIFEST:-}" && -s "${LOCAL_WIDE_MANIFEST}" ]]; then
    echo "Manifest: ${LOCAL_WIDE_MANIFEST}" >&2
  fi
  exit 1
fi
echo "[prep] Verified that the local single-SNP wide subset contains ${TARGET_SNP}."
log_local_subset_target_summary_if_available

REMOTE_DATA_BASENAME="${REMOTE_DATA_BASENAME:-$(basename "${DATA_GZ}")}"

if [[ "${TARGET_CHR}" == "NA" || "${TARGET_BP}" == "NA" || -z "${TARGET_CHR}" || -z "${TARGET_BP}" ]]; then
  locus_out="$(
    perl -MIO::Uncompress::Gunzip=gunzip,\$GunzipError -e '
      use strict;
      use warnings;
      my ($path, $target) = @ARGV;
      my $fh = IO::Uncompress::Gunzip->new($path)
        or die "Cannot open $path: $GunzipError\n";
      my $header = <$fh>;
      while (my $line = <$fh>) {
        chomp $line;
        my @f = split /\t/, $line;
        next unless @f >= 5;
        next unless defined $f[4] && $f[4] eq $target;
        print $f[0], "\t", $f[1], "\n";
        exit 0;
      }
      exit 1;
    ' "${DATA_GZ}" "${TARGET_SNP}" 2>/dev/null || true
  )"
  if [[ -n "${locus_out}" ]]; then
    TARGET_CHR="$(printf '%s\n' "${locus_out}" | awk -F '\t' 'NR==1{print $1}')"
    TARGET_BP="$(printf '%s\n' "${locus_out}" | awk -F '\t' 'NR==1{print $2}')"
  fi
fi

if [[ -z "${OUTPUT_HTML_BASENAME}" ]]; then
  OUTPUT_HTML_BASENAME="${PROJECT_TAG}_SAS_single_snp_with_gtf_${SAFE_TARGET_SNP}_chr${TARGET_CHR}_bp${TARGET_BP}.html"
fi
HTML_OUT="${WORKDIR}/${OUTPUT_HTML_BASENAME}"
PNG_OUT="${WORKDIR}/${OUTPUT_HTML_BASENAME%.html}.png"
RAW_HTML_OUT="${WORKDIR}/${OUTPUT_HTML_BASENAME%.html}.sasraw.html"

if [[ -z "${TARGET_CHR}" || -z "${TARGET_BP}" || "${TARGET_CHR}" == "NA" || "${TARGET_BP}" == "NA" ]]; then
  echo "ERROR: Could not resolve target chr/bp for ${TARGET_SNP} from ${DATA_GZ}." >&2
  exit 1
fi

REGION_START="$(
  perl -e 'my ($bp, $win) = @ARGV; my $start = int($bp - $win); $start = 1 if $start < 1; print $start;' \
    "${TARGET_BP}" "${LOCAL_WINDOW_BP}"
)"
REGION_END="$(
  perl -e 'my ($bp, $win) = @ARGV; my $end = int($bp + $win); print $end;' \
    "${TARGET_BP}" "${LOCAL_WINDOW_BP}"
)"

echo "[prep] Building local GTF subset for ${TARGET_SNP} at chr${TARGET_CHR}:${REGION_START}-${REGION_END}..."
gtf_helper_out="$(
  perl "${GTF_SUBSET_HELPER}" \
    --gtf-url "${GTF_GZ_URL}" \
    --cache-dir "${GTF_CACHE_DIR}" \
    --output "${LOCAL_GTF_SUBSET}" \
    --region "${TARGET_CHR}:${REGION_START}:${REGION_END}" \
    $( [[ "${LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES}" == "1" ]] && printf '%s' "--include-non-protein-coding" || printf '%s' "--no-include-non-protein-coding" )
)"
if [[ ! -s "${LOCAL_GTF_SUBSET}" ]]; then
  echo "ERROR: Failed to build local GTF subset for ${TARGET_SNP}." >&2
  printf '%s\n' "${gtf_helper_out}" >&2
  exit 1
fi
gzip -c "${LOCAL_GTF_SUBSET}" > "${LOCAL_GTF_SUBSET_GZ}"
if [[ ! -s "${LOCAL_GTF_SUBSET_GZ}" ]]; then
  echo "ERROR: Failed to gzip local GTF subset for ${TARGET_SNP}." >&2
  exit 1
fi

perl "${SCHEMA_INCLUDE_HELPER}" \
  --config "${SCHEMA_CONFIG_JSON}" \
  --dataset "${GWAS_DATASET}" \
  --source-type gzip \
  --remote-basename "${REMOTE_DATA_BASENAME}" > "${IMPORT_BLOCK_RENDERED}"

perl "${GTF_IMPORT_INCLUDE_HELPER}" \
  --dataset "${GTF_LOCAL_DSD}" \
  --remote-basename "${REMOTE_GTF_BASENAME}" > "${GTF_IMPORT_BLOCK_RENDERED}"

perl "${RENDER_SAS_HELPER}" \
  --template "${RUN_SAS_TEMPLATE}" \
  --output "${RUN_SAS_RENDERED}" \
  --replace "TARGET_SNP=${TARGET_SNP}" \
  --replace "LOCAL_WINDOW_BP=${LOCAL_WINDOW_BP}" \
  --replace "GTF_LABEL_SNPS=${GTF_LABEL_SNPS}" \
  --replace "OUTPUT_HTML=${OUTPUT_HTML_BASENAME}" \
  --replace "GWAS_DATASET=${GWAS_DATASET}" \
  --replace "TARGET_HIT_DATASET=${TARGET_HIT_DATASET}" \
  --replace "TARGET_LOCAL_DATASET=${TARGET_LOCAL_DATASET}" \
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
  --replace "GTF_YOFFSET4TEXTLABELS=${GTF_YOFFSET4TEXTLABELS}" \
  --replace "GTF_INCLUDE_NON_PROTEIN_CODING=${LOCAL_GTF_INCLUDE_NON_PROTEIN_CODING_GENES}" \
  --replace-file "WIDE_IMPORT_BLOCK=${IMPORT_BLOCK_RENDERED}" \
  --replace-file "GTF_IMPORT_BLOCK=${GTF_IMPORT_BLOCK_RENDERED}"

rm -f "${HTML_OUT}"

upload_support_args=(--upload-file "${PATCHED_LATTICE_MACRO_SAS}" --upload-file "${LOCAL_GTF_SUBSET_GZ}" --upload-file "${DATA_GZ}")
if [[ -f "${MAP_GRP_ASSOC_MACRO_SAS}" ]]; then
  upload_support_args=(--upload-file "${MAP_GRP_ASSOC_MACRO_SAS}" "${upload_support_args[@]}")
fi
if [[ -f "${MULT_GSCATTER_GENE_MACRO_SAS}" ]]; then
  upload_support_args=(--upload-file "${MULT_GSCATTER_GENE_MACRO_SAS}" "${upload_support_args[@]}")
fi
if [[ -f "${ADJ_CLOSE_GENE_GRP_MACRO_SAS}" ]]; then
  upload_support_args=(--upload-file "${ADJ_CLOSE_GENE_GRP_MACRO_SAS}" "${upload_support_args[@]}")
fi
if [[ -f "${SNP_LOCAL_MACRO_SAS}" ]]; then
  echo "[1/6] Uploading SNP_Local_Manhattan_With_GTF macro and support files..."
  upload_support_args=(--upload-file "${SNP_LOCAL_MACRO_SAS}" "${upload_support_args[@]}")
else
  echo "[1/6] Local SNP_Local_Manhattan_With_GTF macro not found; relying on the SAS ODA built-in macro." >&2
  echo "[1/6] Uploading support files..."
fi
oda_upload_many \
  "upload_single_snp_with_gtf_support_${stamp}" \
  "${upload_support_args[@]}"

if ! remote_home_file_matches_local_size "${DATA_GZ}" "${REMOTE_DATA_BASENAME}"; then
  echo "[1b/6] Remote GWAS subset size check failed after bulk upload. Re-uploading ${REMOTE_DATA_BASENAME}..."
  delete_partial_remote_home_file "${REMOTE_DATA_BASENAME}"
  run_oda_helper --upload-file "${DATA_GZ}" --output-prefix "reupload_single_snp_gwas_subset_${stamp}"
  if ! remote_home_file_matches_local_size "${DATA_GZ}" "${REMOTE_DATA_BASENAME}"; then
    echo "ERROR: Remote GWAS subset size mismatch for ${REMOTE_DATA_BASENAME} after re-upload." >&2
    exit 1
  fi
fi
echo "[1b/6] Verified remote GWAS subset upload: ${REMOTE_DATA_BASENAME} ($(remote_home_file_size_bytes "${REMOTE_DATA_BASENAME}") bytes)."

if ! remote_home_file_matches_local_size "${LOCAL_GTF_SUBSET_GZ}" "${REMOTE_GTF_BASENAME}"; then
  echo "[1c/6] Remote local-GTF subset size check failed after bulk upload. Re-uploading ${REMOTE_GTF_BASENAME}..."
  delete_partial_remote_home_file "${REMOTE_GTF_BASENAME}"
  run_oda_helper --upload-file "${LOCAL_GTF_SUBSET_GZ}" --output-prefix "reupload_single_snp_local_gtf_subset_${stamp}"
  if ! remote_home_file_matches_local_size "${LOCAL_GTF_SUBSET_GZ}" "${REMOTE_GTF_BASENAME}"; then
    echo "ERROR: Remote local-GTF subset size mismatch for ${REMOTE_GTF_BASENAME} after re-upload." >&2
    exit 1
  fi
fi
echo "[1c/6] Verified remote local-GTF subset upload: ${REMOTE_GTF_BASENAME} ($(remote_home_file_size_bytes "${REMOTE_GTF_BASENAME}") bytes)."

echo "[5/6] Running SAS single-SNP local Manhattan gene-track plot..."
SINGLE_SNP_GTF_SUBMIT_MAX_ATTEMPTS="${SINGLE_SNP_GTF_SUBMIT_MAX_ATTEMPTS:-2}"
SINGLE_SNP_GTF_SUBMIT_RETRY_SLEEP_SECONDS="${SINGLE_SNP_GTF_SUBMIT_RETRY_SLEEP_SECONDS:-10}"
single_snp_submit_attempt=1
while :; do
  rm -f "${RUN_LOG_FILE}"
  single_snp_submit_rc=0
  if run_single_snp_submit; then
    single_snp_submit_rc=0
  else
    single_snp_submit_rc=$?
    echo "[5b/6] Single-SNP GTF SAS submit attempt ${single_snp_submit_attempt} exited with status ${single_snp_submit_rc}."
  fi
  if [[ "${single_snp_submit_rc}" -eq 0 ]] && ! single_snp_submit_needs_retry; then
    break
  fi
  if [[ "${single_snp_submit_rc}" -ne 0 ]] && remote_home_file_exists "${OUTPUT_HTML_BASENAME}"; then
    echo "[recover] The SAS helper exited with status ${single_snp_submit_rc}, but the remote HTML ${OUTPUT_HTML_BASENAME} already exists. Continuing with download and local recovery."
    break
  fi
  if [[ "${single_snp_submit_rc}" -ne 0 ]] && run_log_reports_missing_target_snp; then
    if [[ -n "${local_target_row}" ]]; then
      echo "[retry] SAS reported ${TARGET_SNP} missing from the uploaded GWAS subset, but the local subset still contains the target row. Re-uploading ${REMOTE_DATA_BASENAME} before the next retry..."
      delete_partial_remote_home_file "${REMOTE_DATA_BASENAME}"
      run_oda_helper --upload-file "${DATA_GZ}" --output-prefix "reupload_single_snp_missing_target_${stamp}_try${single_snp_submit_attempt}" || true
    else
      echo "ERROR: SAS reported ${TARGET_SNP} missing and the local subset also does not contain that target row: ${DATA_GZ}" >&2
      exit 1
    fi
  fi
  if [[ "${single_snp_submit_attempt}" -ge "${SINGLE_SNP_GTF_SUBMIT_MAX_ATTEMPTS}" ]]; then
    if [[ "${single_snp_submit_rc}" -ne 0 ]]; then
      echo "ERROR: Single-SNP GTF SAS submit failed with status ${single_snp_submit_rc} after ${single_snp_submit_attempt} attempt(s)." >&2
    else
      echo "ERROR: Single-SNP GTF SAS submit log still looks incomplete after ${single_snp_submit_attempt} attempts: ${RUN_LOG_FILE}" >&2
    fi
    exit 1
  fi
  next_attempt=$((single_snp_submit_attempt + 1))
  if [[ "${single_snp_submit_rc}" -ne 0 ]]; then
    echo "[5b/6] Single-SNP GTF SAS submit attempt ${single_snp_submit_attempt} failed or timed out; retrying attempt ${next_attempt}/${SINGLE_SNP_GTF_SUBMIT_MAX_ATTEMPTS} after ${SINGLE_SNP_GTF_SUBMIT_RETRY_SLEEP_SECONDS}s..."
  else
    echo "[5b/6] Single-SNP GTF SAS submit attempt ${single_snp_submit_attempt} looked incomplete; retrying attempt ${next_attempt}/${SINGLE_SNP_GTF_SUBMIT_MAX_ATTEMPTS} after ${SINGLE_SNP_GTF_SUBMIT_RETRY_SLEEP_SECONDS}s..."
  fi
  sleep "${SINGLE_SNP_GTF_SUBMIT_RETRY_SLEEP_SECONDS}"
  single_snp_submit_attempt="${next_attempt}"
done

echo "[6/6] Downloading HTML result..."
rm -f "${HTML_OUT}" "${RAW_HTML_OUT}" "${PNG_OUT}"
oda_download_many \
  "download_single_snp_with_gtf_html_${stamp}" \
  --download-file "~/${OUTPUT_HTML_BASENAME}" \
  --download-local-path "${RAW_HTML_OUT}" || true

if [[ ! -s "${RAW_HTML_OUT}" ]]; then
  fallback_raw_html="$(
    find "${RUN_LOG_DIR}" -maxdepth 1 -type f -name 'sas_res_*.html' 2>/dev/null | head -n 1
  )"
  if [[ -n "${fallback_raw_html}" && -s "${fallback_raw_html}" ]]; then
    cp "${fallback_raw_html}" "${RAW_HTML_OUT}"
    echo "[recover] Reused local raw SAS HTML sidecar because the remote HTML download was unavailable: ${RAW_HTML_OUT}"
  fi
fi

if [[ ! -s "${RAW_HTML_OUT}" ]]; then
  echo "ERROR: Expected downloaded HTML was not created or is empty: ${RAW_HTML_OUT}" >&2
  exit 1
fi

remote_png_path="$(
  perl -ne 'if(/The final figure is put here:/){$want=1; next} if($want && m{(/home/\S+\.(?:png|jpg|jpeg|svg))}){print $1; exit}' "${RUN_LOG_FILE}" 2>/dev/null || true
)"
if [[ -n "${remote_png_path}" ]]; then
  echo "[6b/6] Downloading generated local plot PNG..."
  rm -f "${PNG_OUT}"
  run_oda_helper \
    --download-file "${remote_png_path}" \
    --download-local-path "${PNG_OUT}" \
    --output-prefix "download_single_snp_with_gtf_png_${stamp}" || true
fi

build_completed_html_from_single_snp_assets_if_available \
  "${RAW_HTML_OUT}" \
  "${PNG_OUT}" \
  "${HTML_OUT}" \
  "Local Manhattan plot for ${TARGET_SNP}" \
  "Local GTF plot for ${TARGET_SNP}" || true

if [[ ! -s "${HTML_OUT}" ]]; then
  cp "${RAW_HTML_OUT}" "${HTML_OUT}"
fi

echo "Verified HTML: ${HTML_OUT} ($(wc -c < "${HTML_OUT}") bytes)"
if [[ -s "${PNG_OUT}" ]]; then
  echo "Verified PNG: ${PNG_OUT} ($(wc -c < "${PNG_OUT}") bytes)"
fi
if [[ -n "${LOCAL_WIDE_MANIFEST}" && -s "${LOCAL_WIDE_MANIFEST}" ]]; then
  echo "Local subset manifest: ${LOCAL_WIDE_MANIFEST}"
fi
if [[ "${LOCAL_WIDE_AUTOGEN}" == "1" ]]; then
  echo "Auto-generated local wide subset: ${DATA_GZ}"
fi

if [[ "${CLEAN_ODA_INPUT}" == "1" ]]; then
  echo "[cleanup] Removing uploaded GWAS and local GTF subsets from SAS ODA..."
  oda_delete_many \
    "cleanup_single_snp_with_gtf_inputs_${stamp}" \
    --delete-file "${REMOTE_DATA_BASENAME}" \
    --delete-file "${REMOTE_GTF_BASENAME}"
else
  echo "[cleanup] Keeping uploaded GWAS and local GTF subsets because CLEAN_ODA_INPUT=${CLEAN_ODA_INPUT}."
fi

cleanup_remote_generated_outputs

if [[ "${CLEAN_ODA_MACROS}" == "1" ]]; then
  echo "[cleanup] Removing uploaded SAS helper files from SAS ODA..."
  delete_macro_args=(--delete-file "$(basename "${PATCHED_LATTICE_MACRO_SAS}")" --delete-file "${RENDERED_SAS_BASENAME}")
  if [[ -f "${MAP_GRP_ASSOC_MACRO_SAS}" ]]; then
    delete_macro_args=(--delete-file "$(basename "${MAP_GRP_ASSOC_MACRO_SAS}")" "${delete_macro_args[@]}")
  fi
  if [[ -f "${MULT_GSCATTER_GENE_MACRO_SAS}" ]]; then
    delete_macro_args=(--delete-file "$(basename "${MULT_GSCATTER_GENE_MACRO_SAS}")" "${delete_macro_args[@]}")
  fi
  if [[ -f "${ADJ_CLOSE_GENE_GRP_MACRO_SAS}" ]]; then
    delete_macro_args=(--delete-file "$(basename "${ADJ_CLOSE_GENE_GRP_MACRO_SAS}")" "${delete_macro_args[@]}")
  fi
  if [[ -f "${SNP_LOCAL_MACRO_SAS}" ]]; then
    delete_macro_args=(--delete-file "$(basename "${SNP_LOCAL_MACRO_SAS}")" "${delete_macro_args[@]}")
  fi
  oda_delete_many \
    "cleanup_single_snp_macros_${stamp}" \
    "${delete_macro_args[@]}"
else
  echo "[cleanup] Keeping uploaded SAS helper files because CLEAN_ODA_MACROS=${CLEAN_ODA_MACROS}."
fi

echo "Done."
echo "Downloaded HTML: ${HTML_OUT}"

if [[ "${OPEN_RESULT}" == "1" ]]; then
  open_html_result "${HTML_OUT}" || true
else
  echo "Not opening result because OPEN_RESULT=${OPEN_RESULT}."
fi
