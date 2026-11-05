#!/usr/bin/env bash
set -euo pipefail
export PATH=/usr/bin:$PATH
project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$project_dir"
output_dir="$(mktemp -d "$project_dir/SAS_Output_gtex_cache_demo_XXXXXX")"
cat examples/sas_macros/Boxplots4GenesInGTExV8_Cached.sas \
    examples/sas/gtex_boxplot_demo.sas > "$output_dir/demo.sas"
printf 'SAS learning example (synthetic data):\n'
cat examples/sas/gtex_boxplot_demo.sas
printf '\nRunning: perl run_sas_codes_or_script_in_ODA.pl --file "%s" --output-prefix "%s" --run-timeout-seconds 180\n' \
  "$output_dir/demo.sas" "$output_dir/plot"
OPEN_RESULT=0 SAS_ODA_AUTOLOAD_MACROS=0 \
  perl run_sas_codes_or_script_in_ODA.pl \
    --file "$output_dir/demo.sas" --output-prefix "$output_dir/plot" \
    --run-timeout-seconds 180 2>&1 | tee "$output_dir/command.log"
printf '\nValidation:\n'
grep -E '^GTEX_DEMO_PASS:' "$output_dir/plot/output.html.info.txt"
printf '\nHTML results:\n'
find "$output_dir/plot" -type f -name '*.html' -print
printf '\nSAS log: %s\n' "$output_dir/plot/output.html.info.txt"
