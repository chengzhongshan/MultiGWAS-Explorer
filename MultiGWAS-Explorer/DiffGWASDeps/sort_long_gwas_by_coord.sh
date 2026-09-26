#!/usr/bin/env bash
set -euo pipefail

INPUT_GZ="${INPUT_GZ:-}"
OUTPUT_GZ="${OUTPUT_GZ:-}"
EXCLUDED_GZ="${EXCLUDED_GZ:-}"
TMPDIR_SORT="${TMPDIR_SORT:-}"
HTSBIN="${HTSBIN:-}"

if [[ -z "${INPUT_GZ}" || -z "${OUTPUT_GZ}" || -z "${EXCLUDED_GZ}" || -z "${TMPDIR_SORT}" ]]; then
  echo "Required env vars: INPUT_GZ OUTPUT_GZ EXCLUDED_GZ TMPDIR_SORT" >&2
  exit 1
fi

mkdir -p "${TMPDIR_SORT}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
for candidate in "${HTSBIN}" "${SCRIPT_DIR}/../local/bin" \
                 "${SCRIPT_DIR}/../../local/bin" "${SCRIPT_DIR}/../../../local/bin"; do
  [[ -n "${candidate}" ]] || continue
  if { [[ -x "${candidate}/bgzip" ]] || [[ -x "${candidate}/bgzip.exe" ]]; } \
      && { [[ -x "${candidate}/tabix" ]] || [[ -x "${candidate}/tabix.exe" ]]; }; then
    export PATH="${candidate}:$PATH"
    break
  fi
done

HAS_BGZIP=0
HAS_TABIX=0
if command -v bgzip >/dev/null 2>&1; then
  HAS_BGZIP=1
fi
if command -v tabix >/dev/null 2>&1; then
  HAS_TABIX=1
fi

if [[ "${HAS_BGZIP}" -eq 0 ]]; then
  echo "bgzip not found on PATH; activate the pipeline environment or set HTSBIN" >&2
  exit 1
fi
if [[ "${HAS_TABIX}" -eq 0 ]]; then
  echo "tabix not found on PATH; activate the pipeline environment or set HTSBIN" >&2
  exit 1
fi

echo "Input:    ${INPUT_GZ}"
echo "Output:   ${OUTPUT_GZ}"
echo "Excluded: ${EXCLUDED_GZ}"
echo "Tmpdir:   ${TMPDIR_SORT}"
echo "Start:    $(date)"

COMPRESS_CMD=(bgzip -c)

{
  set +o pipefail
  gzip -dc "${INPUT_GZ}" | head -n 1 | sed 's/^/#/'
  set -o pipefail
  gzip -dc "${INPUT_GZ}" |
    tail -n +2 |
    awk -F $'\t' '$1 != "" && $2 ~ /^[0-9]+$/' |
    LC_ALL=C sort \
      -T "${TMPDIR_SORT}" \
      -S 50% \
      -t $'\t' \
      -k1,1V \
      -k2,2n
} | "${COMPRESS_CMD[@]}" > "${OUTPUT_GZ}"

gzip -dc "${INPUT_GZ}" |
  tail -n +2 |
  awk -F $'\t' '$1 == "" || $2 !~ /^[0-9]+$/' |
  gzip -c > "${EXCLUDED_GZ}"

# Windows htslib cannot open Cygwin /mnt or /cygdrive paths. Native Cygwin
# tabix accepts them, so convert only for PE binaries without cygwin1.dll.
tabix_input="${OUTPUT_GZ}"
if [[ "$(uname -s)" == CYGWIN* ]]; then
  tabix_bin="$(command -v tabix)"
  if ! cygcheck "$tabix_bin" 2>/dev/null | grep -qi 'cygwin1.dll'; then
    tabix_input="$(cygpath -w "${OUTPUT_GZ}")"
  fi
fi
tabix -f -s 1 -b 2 -e 2 -S 1 "$tabix_input"

echo "Done: $(date)"
ls -lh "${OUTPUT_GZ}" "${OUTPUT_GZ}.tbi" "${EXCLUDED_GZ}"
