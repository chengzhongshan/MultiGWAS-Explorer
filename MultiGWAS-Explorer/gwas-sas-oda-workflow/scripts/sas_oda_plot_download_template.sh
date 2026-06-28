#!/usr/bin/env bash
set -euo pipefail

# Template for SAS ODA plot workflows.
# Copy this into a project and set the paths below before running.

WORKDIR="${WORKDIR:?Set WORKDIR to the project directory in Cygwin form}"
DATA_GZ="${DATA_GZ:?Set DATA_GZ to the gzipped plot subset}"
MACRO_SAS="${MACRO_SAS:?Set MACRO_SAS to the SAS plotting macro}"
RUN_SAS="${RUN_SAS:?Set RUN_SAS to the SAS driver script}"
REMOTE_PNG="${REMOTE_PNG:-PGC_SCZ_SAS_manhattan.png}"
REMOTE_HTML="${REMOTE_HTML:-PGC_SCZ_SAS_manhattan_png.html}"
SESSION_ID="${SESSION_ID:-mysession}"
CLEAN_ODA_INPUT="${CLEAN_ODA_INPUT:-1}"
OPEN_RESULT="${OPEN_RESULT:-1}"

cd "${WORKDIR}"

open_html_result() {
  local target="${1:-}"
  [[ -n "${target}" && -e "${target}" ]] || return 1
  if command -v cygstart >/dev/null 2>&1; then
    cygstart "${target}" >/dev/null 2>&1 &
    return 0
  fi
  if [[ "$(uname -s 2>/dev/null || printf '')" == "Darwin" ]] && command -v open >/dev/null 2>&1; then
    open "${target}" >/dev/null 2>&1 &
    return 0
  fi
  if command -v xdg-open >/dev/null 2>&1; then
    xdg-open "${target}" >/dev/null 2>&1 &
    return 0
  fi
  if command -v python3 >/dev/null 2>&1; then
    python3 -c 'import pathlib, sys, webbrowser; webbrowser.open(pathlib.Path(sys.argv[1]).resolve().as_uri())' "${target}" >/dev/null 2>&1 &
    return 0
  fi
  if command -v python >/dev/null 2>&1; then
    python -c 'import pathlib, sys, webbrowser; webbrowser.open(pathlib.Path(sys.argv[1]).resolve().as_uri())' "${target}" >/dev/null 2>&1 &
    return 0
  fi
  echo "WARN: Could not auto-open ${target}. Set OPEN_RESULT=0 to disable auto-open." >&2
  return 1
}

stamp="$(date +%Y%m%d_%H%M%S)"
PNG_OUT="${WORKDIR}/${REMOTE_PNG}"
HTML_OUT="${WORKDIR}/${REMOTE_HTML}"

rm -f "${PNG_OUT}" "${HTML_OUT}"

perl -S run_sas_codes_or_script_in_ODA.pl \
  --upload-file "${MACRO_SAS}" \
  --persistent \
  --session-id "${SESSION_ID}" \
  --output-prefix "upload_plot_macro_${stamp}"

perl -S run_sas_codes_or_script_in_ODA.pl \
  --upload-file "${DATA_GZ}" \
  --persistent \
  --session-id "${SESSION_ID}" \
  --output-prefix "upload_plot_data_${stamp}"

perl -S run_sas_codes_or_script_in_ODA.pl \
  --file "${RUN_SAS}" \
  --persistent \
  --session-id "${SESSION_ID}" \
  --output-prefix "run_plot_${stamp}"

perl -S run_sas_codes_or_script_in_ODA.pl \
  --download-file "~/${REMOTE_PNG}" \
  --persistent \
  --session-id "${SESSION_ID}" \
  --output-prefix "download_plot_png_${stamp}"

perl -S run_sas_codes_or_script_in_ODA.pl \
  --download-file "~/${REMOTE_HTML}" \
  --persistent \
  --session-id "${SESSION_ID}" \
  --output-prefix "download_plot_html_${stamp}"

[[ -s "${PNG_OUT}" ]] || { echo "ERROR: missing downloaded PNG: ${PNG_OUT}" >&2; exit 1; }
[[ -s "${HTML_OUT}" ]] || { echo "ERROR: missing downloaded HTML: ${HTML_OUT}" >&2; exit 1; }

if [[ "${CLEAN_ODA_INPUT}" == "1" ]]; then
  perl -S run_sas_codes_or_script_in_ODA.pl \
    --delete-file "$(basename "${DATA_GZ}")" \
    --persistent \
    --session-id "${SESSION_ID}" \
    --output-prefix "cleanup_plot_data_${stamp}"
fi

echo "Downloaded PNG:  ${PNG_OUT} ($(wc -c < "${PNG_OUT}") bytes)"
echo "Downloaded HTML: ${HTML_OUT} ($(wc -c < "${HTML_OUT}") bytes)"

if [[ "${OPEN_RESULT}" == "1" ]]; then
  open_html_result "${HTML_OUT}" || true
fi
