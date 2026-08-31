#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXTRACTOR="${SCRIPT_DIR}/extract_haploreg_ld_for_candidates.pl"
COMBINER="${SCRIPT_DIR}/combine_haploreg_ld_caches.pl"
HAPLOREG_DATA_URL="${HAPLOREG_DATA_URL:-https://pubs.broadinstitute.org/mammals/haploreg/data}"
HAPLOREG_LD_ARCHIVE_DIR="${HAPLOREG_LD_ARCHIVE_DIR:-}"

CANDIDATES=""
OUTPUT=""
POPULATIONS="EUR ASN"
MIN_R2="0.2"
RETAIN_ALL_PROXIES=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --candidates) CANDIDATES="$2"; shift 2 ;;
    --output) OUTPUT="$2"; shift 2 ;;
    --populations) POPULATIONS="$2"; shift 2 ;;
    --min-r2) MIN_R2="$2"; shift 2 ;;
    --retain-all-proxies) RETAIN_ALL_PROXIES=1; shift ;;
    -h|--help)
      echo "Usage: $0 --candidates FILE --output FILE [--populations 'EUR ASN'] [--min-r2 0.2] [--retain-all-proxies]"
      echo "Set HAPLOREG_LD_ARCHIVE_DIR to reuse downloaded LD_POP.tsv.gz files."
      exit 0
      ;;
    *) echo "ERROR: Unknown argument: $1" >&2; exit 2 ;;
  esac
done

[[ -s "${CANDIDATES}" ]] || { echo "ERROR: Missing candidate file: ${CANDIDATES}" >&2; exit 2; }
[[ -n "${OUTPUT}" ]] || { echo "ERROR: --output is required" >&2; exit 2; }
perl -e 'exit(($ARGV[0]+0) < 0.2 || ($ARGV[0]+0) > 1 ? 1 : 0)' "${MIN_R2}" || {
  echo "ERROR: Downloadable HaploReg v4 LD archives support r2>=0.2 only." >&2
  exit 2
}

mkdir -p "$(dirname "${OUTPUT}")"
tmp_dir="$(mktemp -d "${TMPDIR:-/tmp}/haploreg_candidate_cache.XXXXXX")"
pids=()
cache_files=()
cleanup_cache_temp() {
  local path
  for path in "${cache_files[@]:-}"; do
    if [[ "${path}" == "${tmp_dir}/"* ]]; then
      rm -f -- "${path}"
    fi
  done
  rmdir -- "${tmp_dir}" 2>/dev/null || true
}
trap cleanup_cache_temp EXIT

for population in ${POPULATIONS}; do
  population="$(printf '%s' "${population}" | tr '[:lower:]' '[:upper:]')"
  case "${population}" in AFR|AMR|ASN|EUR) ;; *) echo "ERROR: Invalid population: ${population}" >&2; exit 2 ;; esac
  population_cache="${tmp_dir}/candidate_ld_${population}.tsv"
  population_summary="${OUTPUT}.${population}.summary.json"
  cache_files+=("${population_cache}")
  extractor_proxy_args=()
  if [[ "${RETAIN_ALL_PROXIES}" == "1" ]]; then
    extractor_proxy_args+=(--no-candidate-proxies-only)
  fi

  if [[ -n "${HAPLOREG_LD_ARCHIVE_DIR}" ]]; then
    archive="${HAPLOREG_LD_ARCHIVE_DIR}/LD_${population}.tsv.gz"
    [[ -s "${archive}" ]] || { echo "ERROR: Missing local archive: ${archive}" >&2; exit 2; }
    perl "${EXTRACTOR}" \
      --candidates "${CANDIDATES}" --population "${population}" \
      --input "${archive}" --output "${population_cache}" \
      --summary "${population_summary}" --min-r2 "${MIN_R2}" \
      "${extractor_proxy_args[@]}" &
  else
    (
      set -o pipefail
      curl -L --fail --retry 5 --silent --show-error \
        "${HAPLOREG_DATA_URL}/LD_${population}.tsv.gz" |
        gzip -dc |
        perl "${EXTRACTOR}" \
          --candidates "${CANDIDATES}" --population "${population}" \
          --input - --output "${population_cache}" \
          --summary "${population_summary}" --min-r2 "${MIN_R2}" \
          "${extractor_proxy_args[@]}"
    ) &
  fi
  pids+=("$!")
done

failed=0
for pid in "${pids[@]}"; do
  wait "${pid}" || failed=1
done
[[ "${failed}" -eq 0 ]] || { echo "ERROR: At least one HaploReg archive extraction failed." >&2; exit 1; }

combine_args=()
for cache in "${cache_files[@]}"; do
  combine_args+=(--input "${cache}")
done
perl "${COMBINER}" "${combine_args[@]}" --output "${OUTPUT}"

echo "Candidate-specific HaploReg cache: ${OUTPUT}"
echo "Archive minimum r2: ${MIN_R2}"
