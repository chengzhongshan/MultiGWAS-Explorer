#!/usr/bin/env bash
set -euo pipefail

INPUT_GZ="${INPUT_GZ:-}"
OUTPUT_GZ="${OUTPUT_GZ:-}"
EXCLUDED_GZ="${EXCLUDED_GZ:-}"
TMPDIR_SORT="${TMPDIR_SORT:-}"
HTSBIN="${HTSBIN:-/mnt/g/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/local/bin}"

if [[ -z "${INPUT_GZ}" || -z "${OUTPUT_GZ}" || -z "${EXCLUDED_GZ}" || -z "${TMPDIR_SORT}" ]]; then
  echo "Required env vars: INPUT_GZ OUTPUT_GZ EXCLUDED_GZ TMPDIR_SORT" >&2
  exit 1
fi

mkdir -p "${TMPDIR_SORT}"

if [[ -x "${HTSBIN}/bgzip" && -x "${HTSBIN}/tabix" ]]; then
  export PATH="${HTSBIN}:$PATH"
elif ! command -v bgzip >/dev/null 2>&1 || ! command -v tabix >/dev/null 2>&1; then
  export PATH="/mnt/e/plink_win64:$PATH"
fi

HAS_BGZIP=0
HAS_TABIX=0
if command -v bgzip >/dev/null 2>&1; then
  HAS_BGZIP=1
fi
if command -v tabix >/dev/null 2>&1; then
  HAS_TABIX=1
fi

if [[ "${HAS_BGZIP}" -eq 0 ]]; then
  echo "bgzip not found on PATH; falling back to gzip output without block indexing" >&2
fi
if [[ "${HAS_TABIX}" -eq 0 ]]; then
  echo "tabix not found on PATH; a placeholder ${OUTPUT_GZ}.tbi note will be written" >&2
fi

echo "Input:    ${INPUT_GZ}"
echo "Output:   ${OUTPUT_GZ}"
echo "Excluded: ${EXCLUDED_GZ}"
echo "Tmpdir:   ${TMPDIR_SORT}"
echo "Start:    $(date)"

if [[ "${HAS_BGZIP}" -eq 1 ]]; then
  COMPRESS_CMD=(bgzip -c)
else
  COMPRESS_CMD=(gzip -c)
fi

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

if [[ "${HAS_TABIX}" -eq 1 && "${HAS_BGZIP}" -eq 1 ]]; then
  tabix -f -s 1 -b 2 -e 2 -S 1 "${OUTPUT_GZ}"
else
  cat > "${OUTPUT_GZ}.tbi" <<EOF
placeholder_index
reason=$([[ "${HAS_BGZIP}" -eq 0 ]] && echo "bgzip_missing" || echo "tabix_missing")
file=${OUTPUT_GZ}
created=$(date)
EOF
fi

echo "Done: $(date)"
ls -lh "${OUTPUT_GZ}" "${OUTPUT_GZ}.tbi" "${EXCLUDED_GZ}"
