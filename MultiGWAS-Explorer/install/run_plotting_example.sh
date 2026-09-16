#!/usr/bin/env bash
# A synthetic local-locus plot with no GWAS downloads or SAS credentials.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${SCRIPT_DIR}/common.sh"
activate_perl_env
activate_python_env
output_dir="${1:-${PIPELINE_ROOT}/example-output}"
mkdir -p "$output_dir"
output_dir="$(cd "$output_dir" && pwd)"
cat > "${output_dir}/example.tsv" <<'EOF'
CHR	BP	SNP	EUR_P	EUR_Z
1	100	rs100	1e-6	-4
1	120	rs200	1e-5	2
1	140	rs300	1e-4	4
EOF
perl "${PIPELINE_ROOT}/DiffGWASDeps/gnuplot/pdl_gunplot_local_locus.pl" \
  --data "${output_dir}/example.tsv" \
  --snp rs200 --window-bp 100 \
  --pcols EUR_P --zcols EUR_Z --labels 'Synthetic GWAS' \
  --out-prefix "${output_dir}/example"
[ -s "${output_dir}/example.png" ] || die "Example plot was not created"
log "Synthetic example plot: ${output_dir}/example.png"
