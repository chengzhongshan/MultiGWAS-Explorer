#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf -- "${TMP_DIR}"' EXIT

perl "${PROJECT_ROOT}/DiffGWASDeps/clump_top_hits_with_haploreg_cache.pl" \
  --candidates "${SCRIPT_DIR}/fixtures/haploreg_cache_candidates.csv" \
  --cache "${SCRIPT_DIR}/fixtures/haploreg_cache_edges.tsv" \
  --sqlite "${TMP_DIR}/cache.sqlite" \
  --output-leads "${TMP_DIR}/leads.csv" \
  --output-audit "${TMP_DIR}/audit.tsv" \
  --signal-column P \
  --populations EUR \
  --r2-threshold 0.2 \
  --population-rule ANY

diff -u "${SCRIPT_DIR}/fixtures/haploreg_cache_clumped_leads.csv" "${TMP_DIR}/leads.csv"
diff -u "${SCRIPT_DIR}/fixtures/haploreg_cache_clumped_audit.tsv" "${TMP_DIR}/audit.tsv"
printf 'HaploReg cache LD-clumping fixture: PASS\n'
