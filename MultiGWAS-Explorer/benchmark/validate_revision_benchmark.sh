#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

die() { printf '[validate] ERROR: %s\n' "$*" >&2; exit 1; }
require_file() { [[ -s "$1" ]] || die "missing or empty file: $1"; }
require_text() {
  local pattern="$1" path="$2"
  grep -Eq "$pattern" "$path" || die "expected pattern '$pattern' in $path"
}

required_files=(
  "${PROJECT_ROOT}/auto_prepare_and_run_diff_gwas.pl"
  "${PROJECT_ROOT}/server.pl"
  "${PROJECT_ROOT}/DiffGWASDeps/diff_pairwise_gwas.pl"
  "${PROJECT_ROOT}/DiffGWASDeps/get_top_signal_with_ld.sas"
  "${PROJECT_ROOT}/DiffGWASDeps/QueryLD_SNPs_at_Haploreg4.sas"
  "${PROJECT_ROOT}/DiffGWASDeps/QueryMulti_LD_SNPs_at_Haploreg4.sas"
  "${PROJECT_ROOT}/DiffGWASDeps/build_haploreg_ld_candidate_cache.sh"
  "${PROJECT_ROOT}/DiffGWASDeps/clump_top_hits_with_haploreg_cache.pl"
  "${SCRIPT_DIR}/Dockerfile.comparators"
  "${SCRIPT_DIR}/run_agent_interface_benchmark.ps1"
  "${SCRIPT_DIR}/agent_interface/AI_Interface_Evaluation_Requests.jsonl"
  "${SCRIPT_DIR}/agent_interface/agent_interface_summary.tsv"
  "${SCRIPT_DIR}/agent_interface/failure_recovery_summary.tsv"
  "${SCRIPT_DIR}/fixture_all_100k/official_comparison.tsv"
  "${SCRIPT_DIR}/harmonization_audit/harmonization_audit.tsv"
  "${SCRIPT_DIR}/unfiltered_qc_harmonized/unfiltered_genomic_inflation.tsv"
  "${SCRIPT_DIR}/ld_top_hit_benchmark_summary.tsv"
)
for path in "${required_files[@]}"; do require_file "$path"; done

perl_scripts=(
  "${PROJECT_ROOT}/auto_prepare_and_run_diff_gwas.pl"
  "${PROJECT_ROOT}/server.pl"
  "${PROJECT_ROOT}/DiffGWASDeps/diff_pairwise_gwas.pl"
  "${PROJECT_ROOT}/DiffGWASDeps/audit_allele_harmonization.pl"
  "${PROJECT_ROOT}/DiffGWASDeps/benchmark_diff_methods.pl"
  "${PROJECT_ROOT}/DiffGWASDeps/clump_top_hits_with_haploreg_cache.pl"
  "${PROJECT_ROOT}/DiffGWASDeps/compare_official_benchmarks.pl"
  "${PROJECT_ROOT}/DiffGWASDeps/export_official_benchmark_fixture.pl"
  "${PROJECT_ROOT}/DiffGWASDeps/extract_haploreg_ld_for_candidates.pl"
  "${PROJECT_ROOT}/DiffGWASDeps/reviewer_qc_long_report.pl"
)
for path in "${perl_scripts[@]}"; do perl -c "$path" >/dev/null; done

require_text 'top_hit_selection_method.*ld' "${PROJECT_ROOT}/configs/spec_pgc_scz_sex_common_automation.json"
require_text 'exclude_strand_ambiguous.*1' "${PROJECT_ROOT}/configs/spec_pgc_scz_sex_common_automation.json"
require_text 'max_eaf_abs_diff.*0\.2' "${PROJECT_ROOT}/configs/spec_pgc_scz_sex_common_automation.json"
require_text 'get_top_signal_with_ld' "${PROJECT_ROOT}/DiffGWASDeps/get_top_hits4Manhattan.sas"
require_text 'PRUNED_DISTANCE_FALLBACK' "${PROJECT_ROOT}/DiffGWASDeps/get_top_signal_with_ld.sas"

perl -MJSON::PP -e '
  my $path=shift;
  open my $fh,"<",$path or die $!;
  my (%n,$total,$timed);
  while(<$fh>){ next unless /\S/; my $r=decode_json($_); $total++; $timed++ if $r->{timed}; $n{$r->{record_type}}++ }
  die "request count mismatch\n" unless $total==17 && $timed==16;
  die "request classes mismatch\n" unless $n{infrastructure_template}==1 && $n{successful_task}==10 && $n{intentionally_invalid_task}==3 && $n{corrected_retry}==3;
' "${SCRIPT_DIR}/agent_interface/AI_Interface_Evaluation_Requests.jsonl"

awk -F '\t' '
  NR==1 { for(i=1;i<=NF;i++) h[$i]=i; next }
  $h["PATH"]=="CLI" || $h["PATH"]=="MCP_AGENT" {
    if ($h["RUNS"] != 10 || $h["SUCCESS_RATE"] != 1 || $h["PARITY_RATE"] != 1 || $h["UNIQUE_ARTIFACT_HASHES"] != 1) exit 1;
    seen++
  }
  END { if (seen != 2) exit 1 }
' "${SCRIPT_DIR}/agent_interface/agent_interface_summary.tsv" || die 'agent summary validation failed'

awk -F '\t' '
  NR==2 { if ($1 != 3 || $2 != 1 || $3 != 1) exit 1; seen=1 }
  END { if (!seen) exit 1 }
' "${SCRIPT_DIR}/agent_interface/failure_recovery_summary.tsv" || die 'failure-recovery validation failed'

awk -F '\t' '
  NR==1 { for(i=1;i<=NF;i++) h[$i]=i; next }
  { if ($h["missing_rows"] != 0 || $h["discordant_p_lt_0p05"] != 0 || $h["discordant_p_lt_1e_5"] != 0 || $h["discordant_p_lt_5e_8"] != 0) exit 1; seen++ }
  END { if (seen != 2) exit 1 }
' "${SCRIPT_DIR}/fixture_all_100k/official_comparison.tsv" || die 'official comparator validation failed'

awk -F '\t' '
  NR==1 { for(i=1;i<=NF;i++) h[$i]=i; next }
  { direct += $h["DIRECT"]; ambiguous += $h["STRAND_AMBIGUOUS"]; eligible += $h["REVISED_HARMONIZATION_ELIGIBLE"] }
  END { if (direct != 22292345 || ambiguous != 3386012 || eligible != 18906333) exit 1 }
' "${SCRIPT_DIR}/harmonization_audit/harmonization_audit.tsv" || die 'harmonization totals validation failed'

"${BASH:-bash}" "${SCRIPT_DIR}/test_haploreg_cache_fixture.sh"

printf '[validate] PASS: syntax, exact request records, comparator agreement, harmonization totals, and LD-clumping assets validated.\n'
