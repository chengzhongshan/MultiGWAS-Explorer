#!/usr/bin/env bash
set -euo pipefail

workdir="/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs"
input="${workdir}/PGC_SCZ_sex_stratified_merged_long.tsv.gz"
output="${workdir}/PGC_SCZ_sex_stratified_merged_long.sorted.coord.tsv.gz"
excluded="${workdir}/PGC_SCZ_sex_stratified_merged_long.sorted.excluded_noncoord.tsv.gz"
tmpdir="${workdir}/sort_tmp"
htsbin="/mnt/g/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/local/bin"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

mkdir -p "$tmpdir"
cd "$workdir"

if [[ ! -s "$input" ]]; then
  echo "Missing input: $input" >&2
  exit 1
fi

if [[ "$(uname -s)" == CYGWIN* ]]; then
  PATH="/usr/local/bin:/usr/bin:${PATH}"
  export PATH
fi
resolve_hts_tool() {
  perl -I"${script_dir}/DiffGWASDeps" -MHTSToolResolver=resolve_hts_tool -e '
    my ($name, $explicit, $start) = @ARGV;
    my $path = resolve_hts_tool($name, explicit => $explicit, start_dir => $start);
    exit 1 unless defined($path) && length($path);
    print $path;
  ' "$1" "${HTSBIN:-${htsbin}}" "${script_dir}"
}
bgzip_bin="$(resolve_hts_tool bgzip)" || { echo "Native bgzip not found" >&2; exit 1; }
tabix_bin="$(resolve_hts_tool tabix)" || { echo "Native tabix not found" >&2; exit 1; }

echo "Input:  $input"
echo "Output: $output"
echo "Excluded non-coordinate rows: $excluded"
echo "Tmpdir: $tmpdir"
echo "Start:  $(date)"
echo "bgzip:  ${bgzip_bin}"
echo "tabix:  ${tabix_bin}"

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
} | "${bgzip_bin}" -c > "$output"

zcat "$input" |
  tail -n +2 |
  awk -F $'\t' '$1 == "" || $2 !~ /^[0-9]+$/' |
  gzip -c > "$excluded"

"${tabix_bin}" -f -s 1 -b 2 -e 2 -S 1 "$output"

echo "Done:   $(date)"
ls -lh "$output" "$output.tbi" "$excluded"
