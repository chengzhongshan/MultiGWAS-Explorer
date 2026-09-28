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
if [[ "$(uname -s)" == CYGWIN* ]]; then
  PATH="/usr/local/bin:/usr/bin:${PATH}"
  export PATH
fi
resolve_hts_tool() {
  perl -I"${SCRIPT_DIR}" -MHTSToolResolver=resolve_hts_tool -e '
    my ($name, $explicit, $start) = @ARGV;
    my $path = resolve_hts_tool($name, explicit => $explicit, start_dir => $start);
    exit 1 unless defined($path) && length($path);
    print $path;
  ' "$1" "${HTSBIN}" "${SCRIPT_DIR}"
}
bgzip_bin="$(resolve_hts_tool bgzip)" || {
  echo "Native bgzip not found; run install/repair_and_test_cygwin.sh or set HTSBIN" >&2
  exit 1
}
tabix_bin="$(resolve_hts_tool tabix)" || {
  echo "Native tabix not found; run install/repair_and_test_cygwin.sh or set HTSBIN" >&2
  exit 1
}

echo "Input:    ${INPUT_GZ}"
echo "Output:   ${OUTPUT_GZ}"
echo "Excluded: ${EXCLUDED_GZ}"
echo "Tmpdir:   ${TMPDIR_SORT}"
echo "Start:    $(date)"

COMPRESS_CMD=("${bgzip_bin}" -c)

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

"${tabix_bin}" -f -s 1 -b 2 -e 2 -S 1 "${OUTPUT_GZ}"

echo "Done: $(date)"
ls -lh "${OUTPUT_GZ}" "${OUTPUT_GZ}.tbi" "${EXCLUDED_GZ}"
