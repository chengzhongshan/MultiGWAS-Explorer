#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP="$(mktemp -d)"
trap 'rm -rf -- "$TMP"' EXIT

perl "$ROOT/DiffGWASDeps/clump_top_hits_with_ld_cache.pl" \
  --candidates "$ROOT/benchmark/fixtures/plink2_1kg_ld_candidates.tsv" \
  --cache "$ROOT/benchmark/fixtures/plink2_1kg_ld_edges.tsv" \
  --sqlite "$TMP/ld.sqlite" \
  --output-leads "$TMP/leads.tsv" \
  --output-audit "$TMP/audit.tsv" \
  --signal-column P --signal-threshold 1e-6 \
  --populations "EUR EAS" --population-rule ANY \
  --r2-threshold 0.1 --cache-min-r2 0 \
  --cache-source PLINK2_1KG_DIRECT --fallback-distance-bp 0

[[ "$(wc -l < "$TMP/leads.tsv")" -eq 4 ]]
grep -q $'RS185665940\tRS185665940\tSELECTED_LEAD\tLD\tPLINK2_1KG_DIRECT\tPARTIAL_ESTIMABLE' "$TMP/audit.tsv"
grep -q $'RS7755143\tRS753634\tPRUNED_LD\tLD\tPLINK2_1KG_DIRECT' "$TMP/audit.tsv"
! grep -q 'PRUNED_DISTANCE_FALLBACK' "$TMP/audit.tsv"
echo "PLINK2 LD cache fixture: PASS"
