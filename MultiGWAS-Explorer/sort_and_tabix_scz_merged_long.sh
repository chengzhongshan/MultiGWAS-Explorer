#!/usr/bin/env bash
set -euo pipefail

workdir="/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs"
input="${workdir}/PGC_SCZ_sex_stratified_merged_long.tsv.gz"
output="${workdir}/PGC_SCZ_sex_stratified_merged_long.sorted.coord.tsv.gz"
excluded="${workdir}/PGC_SCZ_sex_stratified_merged_long.sorted.excluded_noncoord.tsv.gz"
tmpdir="${workdir}/sort_tmp"
htsbin="/mnt/g/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/local/bin"

mkdir -p "$tmpdir"
cd "$workdir"

if [[ ! -s "$input" ]]; then
  echo "Missing input: $input" >&2
  exit 1
fi

if [[ -x "${htsbin}/bgzip" && -x "${htsbin}/tabix" ]]; then
  export PATH="${htsbin}:$PATH"
elif ! command -v bgzip >/dev/null 2>&1 || ! command -v tabix >/dev/null 2>&1; then
  export PATH="/mnt/e/plink_win64:$PATH"
fi

if ! command -v bgzip >/dev/null 2>&1 || ! command -v tabix >/dev/null 2>&1; then
  echo "bgzip and/or tabix not found on PATH" >&2
  exit 1
fi

echo "Input:  $input"
echo "Output: $output"
echo "Excluded non-coordinate rows: $excluded"
echo "Tmpdir: $tmpdir"
echo "Start:  $(date)"
echo "bgzip:  $(command -v bgzip)"
echo "tabix:  $(command -v tabix)"

{
  set +o pipefail
  zcat "$input" | head -n 1 | sed 's/^/#/'
  set -o pipefail
  zcat "$input" |
    tail -n +2 |
    awk -F $'\t' '$1 != "" && $2 ~ /^[0-9]+$/' |
    LC_ALL=C sort \
      -T "$tmpdir" \
      -S 50% \
      -t $'\t' \
      -k1,1V \
      -k2,2n
} | bgzip -c > "$output"

zcat "$input" |
  tail -n +2 |
  awk -F $'\t' '$1 == "" || $2 !~ /^[0-9]+$/' |
  gzip -c > "$excluded"

tabix -f -s 1 -b 2 -e 2 -S 1 "$output"

echo "Done:   $(date)"
ls -lh "$output" "$output.tbi" "$excluded"
