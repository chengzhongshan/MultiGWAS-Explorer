#!/usr/bin/env bash
# Reproduce the technical measurements summarized in benchmark/README.md.
#
# The default output is benchmark/reproduction so versioned benchmark fixtures
# are not overwritten. Run this script from Cygwin,
# WSL, Linux, or another Bash environment with access to the source GWAS files.
# SAS OnDemand stages also require the repository SASPy/ODA setup and account.

set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DEPS_DIR="${PROJECT_ROOT}/DiffGWASDeps"
RUN_ROOT="${BENCHMARK_OUT_DIR:-${SCRIPT_DIR}/reproduction}"
FIXTURE_DIR="${RUN_ROOT}/fixture_all_100k"

# Keep compiled Perl dependencies platform-specific when the same checkout is
# used from Cygwin and Linux/Docker.
case "$(uname -s)" in
  CYGWIN*) PERL_LOCAL_TREE="${PROJECT_ROOT}/local/perl5-cygwin/lib/perl5" ;;
  Darwin*) PERL_LOCAL_TREE="${PROJECT_ROOT}/local/perl5-darwin/lib/perl5" ;;
  *)       PERL_LOCAL_TREE="${PROJECT_ROOT}/local/perl5-linux/lib/perl5" ;;
esac
if [[ -d "${PERL_LOCAL_TREE}" ]]; then
  export PERL5LIB="${PERL_LOCAL_TREE}${PERL5LIB:+:${PERL5LIB}}"
fi

INPUT_DIR_DEFAULT="/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs"
INPUT_DIR="${PGC_INPUT_DIR:-${INPUT_DIR_DEFAULT}}"
MERGED_LONG_GZ="${MERGED_LONG_GZ:-${INPUT_DIR}/PGC_SCZ_female_vs_male_diff_effects_merged_long.sorted.coord.tsv.gz}"
[[ -n "${QC_STANDARDIZED_LONG_GZ+x}" ]] && QC_INPUT_EXPLICIT=1 || QC_INPUT_EXPLICIT=0
QC_STANDARDIZED_LONG_GZ="${QC_STANDARDIZED_LONG_GZ:-${INPUT_DIR}/PGC_SCZ_female_vs_male_diff_effects.stdized.tsv.gz}"
WIDE_GZ="${WIDE_GZ:-${INPUT_DIR}/PGC_SCZ_female_vs_male_diff_effects.stdized.wide_beta_se_p_p_lt_0p05.final.tsv.gz}"

SPEC_JSON="${SPEC_JSON:-${PROJECT_ROOT}/configs/spec_pgc_scz_sex_common_automation.json}"
MERGE_CONFIG="${MERGE_CONFIG:-${PROJECT_ROOT}/configs/auto_PGC_SCZ_female_vs_male_diff_effects_merge.json}"
DIFF_CONFIG="${DIFF_CONFIG:-${PROJECT_ROOT}/configs/auto_PGC_SCZ_female_vs_male_diff_effects_diff.json}"
RUNNER_CONFIG="${RUNNER_CONFIG:-${PROJECT_ROOT}/configs/auto_PGC_SCZ_female_vs_male_diff_effects_runner.json}"

COMPARATOR_IMAGE="${COMPARATOR_IMAGE:-multigwas-comparators:20260819}"
GNU_TIME="${GNU_TIME:-/usr/bin/time}"
SAS_SESSION_ID="${SAS_SESSION_ID:-pipeline_revision_benchmark}"
SAS_RUNNER="${SAS_RUNNER:-${PROJECT_ROOT}/run_sas_codes_or_script_in_ODA.pl}"
QUERY_LD_MACRO="${QUERY_LD_MACRO:-${DEPS_DIR}/QueryLD_SNPs_at_Haploreg4.sas}"
# Repository copy of the multi-SNP macro used by this benchmark.
QUERY_MULTI_LD_MACRO="${QUERY_MULTI_LD_MACRO:-${DEPS_DIR}/QueryMulti_LD_SNPs_at_Haploreg4.sas}"
GET_TOP_SIGNAL_LD_MACRO="${GET_TOP_SIGNAL_LD_MACRO:-${DEPS_DIR}/get_top_signal_with_ld.sas}"
GET_TOP_SIGNAL_DISTANCE_MACRO="${GET_TOP_SIGNAL_DISTANCE_MACRO:-${DEPS_DIR}/get_top_signal_within_dist.sas}"

log() { printf '[benchmark] %s\n' "$*"; }
die() { printf '[benchmark] ERROR: %s\n' "$*" >&2; exit 1; }
require_command() { command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"; }
require_file() { [[ -s "$1" ]] || die "Required input is absent or empty: $1"; }

docker_path() {
  case "$(uname -s)" in
    CYGWIN*|MINGW*|MSYS*) require_command cygpath; cygpath -w "$1" ;;
    *) printf '%s\n' "$1" ;;
  esac
}

powershell_path() {
  case "$(uname -s)" in
    CYGWIN*|MINGW*|MSYS*) require_command cygpath; cygpath -w "$1" ;;
    *) printf '%s\n' "$1" ;;
  esac
}

run_powershell_file() {
  local script_path="$1"
  shift
  if command -v powershell.exe >/dev/null 2>&1; then
    powershell.exe -NoProfile -ExecutionPolicy Bypass -File "$(powershell_path "${script_path}")" "$@"
  elif command -v pwsh >/dev/null 2>&1; then
    pwsh -NoProfile -File "${script_path}" "$@"
  else
    die "PowerShell is required for this stage (powershell.exe or pwsh)."
  fi
}

elapsed_from_gnu_time() {
  perl -ne '
    if (/Elapsed \(wall clock\).*?:\s*([0-9:.]+)\s*$/) {
      @p=split /:/,$1;
      $s=@p==3 ? $p[0]*3600+$p[1]*60+$p[2] : $p[0]*60+$p[1];
      printf "%.3f",$s; exit
    }
  ' "$1"
}

rss_from_gnu_time() {
  awk -F': ' '/Maximum resident set size/ { print $NF; exit }' "$1"
}

data_rows() {
  awk 'END { print NR > 0 ? NR - 1 : 0 }' "$1"
}

show_inputs() {
  cat <<EOF
PROJECT_ROOT=${PROJECT_ROOT}
RUN_ROOT=${RUN_ROOT}
MERGED_LONG_GZ=${MERGED_LONG_GZ}
QC_STANDARDIZED_LONG_GZ=${QC_STANDARDIZED_LONG_GZ}
WIDE_GZ=${WIDE_GZ}
COMPARATOR_IMAGE=${COMPARATOR_IMAGE}
SAS_SESSION_ID=${SAS_SESSION_ID}
HAPLOREG_LD_ARCHIVE_DIR=${HAPLOREG_LD_ARCHIVE_DIR:-<stream from Broad>}
EOF
}

usage() {
  cat <<'EOF'
Usage: bash benchmark/run_differential_gwas_revision_benchmarks.sh TARGET

Targets:
  show-inputs          Print resolved source, output, image, and SAS paths.
  prepare-pipeline     Rebuild merged/differential/plot data from raw PGC files.
  prepare-qc-source    Build an unfiltered paired table retaining ambiguous SNPs
                       so the QC stage, rather than an upstream step, records the
                       3,386,012 conservative ambiguity exclusions.
  comparators          Run the shared 100,000-row MultiGWAS-Explorer, GWAMA 2.2.2,
                       and EasyStrata 8.6 benchmark and numerical comparisons.
  full-qc              Run full-source allele audit and both unfiltered QC scans.
  power                Recalculate observed-SE detectable effects and power.
  agent                Run 10 CLI/agent parity trials and 3 recovery trials.
  ld-archive           Generate common candidates, stream/reuse HaploReg EUR/ASN
                       archives, build SQLite, and perform three local clumps.
  ld-sas               Run nominated-SNP LD, pairwise LD, four-SNP regression,
                       and real differential clumping in SAS ODA.
  ld-sas-common        Run the 11,701-candidate SAS cache-plus-web hybrid stage.
  summarize-ld         Consolidate completed LD audit files into the report TSV.
  web-only-control     Opt-in web-only 11,701-candidate negative control; expected
                       to time out/non-complete, as documented in the report.
  validate-results     Validate scripts and the compact published result tables.
  local-report         Run comparators, full-qc, power, agent, and ld-archive.
  all                  Run local-report plus the successful SAS ODA stages.

Important environment overrides:
  PGC_INPUT_DIR, MERGED_LONG_GZ, QC_STANDARDIZED_LONG_GZ, WIDE_GZ
  BENCHMARK_OUT_DIR, COMPARATOR_IMAGE, GNU_TIME
  HAPLOREG_LD_ARCHIVE_DIR (directory containing LD_EUR.tsv.gz and LD_ASN.tsv.gz)
  SAS_SESSION_ID, SAS_RUNNER, QUERY_LD_MACRO, QUERY_MULTI_LD_MACRO

No target is run when TARGET is omitted. Outputs default to benchmark/reproduction.
The intentional web-only timeout is never included in "all"; invoke it explicitly.
EOF
}

prepare_pipeline() {
  require_command perl
  require_file "${SPEC_JSON}"
  log "Building the raw-source pipeline with plots disabled."
  (
    cd "${PROJECT_ROOT}"
    perl auto_prepare_and_run_diff_gwas.pl \
      --spec "${SPEC_JSON}" \
      --mode full \
      --skip-plots
  )
}

prepare_qc_source() {
  require_command perl
  require_file "${MERGED_LONG_GZ}"
  require_file "${DIFF_CONFIG}"
  mkdir -p "${RUN_ROOT}/qc_source"
  local qc_config="${RUN_ROOT}/qc_source/unfiltered_pairwise.config.json"
  local qc_diff="${RUN_ROOT}/qc_source/unfiltered_pairwise.tsv.gz"
  local qc_diff_manifest="${RUN_ROOT}/qc_source/unfiltered_pairwise.manifest.tsv"
  local qc_std="${RUN_ROOT}/qc_source/unfiltered_pairwise.stdized.tsv.gz"
  local qc_std_manifest="${RUN_ROOT}/qc_source/unfiltered_pairwise.stdized.manifest.tsv"

  BASE_DIFF_CONFIG="${DIFF_CONFIG}" \
  QC_INPUT="${MERGED_LONG_GZ}" \
  QC_OUTPUT="${qc_diff}" \
  QC_MANIFEST="${qc_diff_manifest}" \
  perl -MJSON::PP -0777 -e '
    open my $fh,"<",$ENV{BASE_DIFF_CONFIG} or die $!;
    my $cfg=decode_json(<$fh>); close $fh;
    $cfg->{input}=$ENV{QC_INPUT};
    $cfg->{output}=$ENV{QC_OUTPUT};
    $cfg->{manifest}=$ENV{QC_MANIFEST};
    $cfg->{exclude_strand_ambiguous}=JSON::PP::false;
    $cfg->{max_eaf_abs_diff}=0.20;
    print JSON::PP->new->canonical->pretty->encode($cfg);
  ' > "${qc_config}"

  perl "${DEPS_DIR}/diff_pairwise_gwas.pl" --config "${qc_config}"
  perl "${DEPS_DIR}/standardize_diff_gwas_zscore.pl" \
    --input "${qc_diff}" \
    --output "${qc_std}" \
    --manifest "${qc_std_manifest}" \
    --z-col DIFF_Z
  log "QC input retaining direct strand-ambiguous pairs: ${qc_std}"
  log "A later full-qc target will automatically prefer this generated file unless QC_STANDARDIZED_LONG_GZ is explicitly set."
}

run_comparators() {
  require_command perl
  require_command docker
  require_file "${MERGED_LONG_GZ}"
  require_file "${SCRIPT_DIR}/Dockerfile.comparators"
  mkdir -p "${FIXTURE_DIR}"

  log "Exporting one shared fixture with the first 100,000 valid ALL female/male pairs."
  perl "${DEPS_DIR}/export_official_benchmark_fixture.pl" \
    --input "${MERGED_LONG_GZ}" \
    --output-dir "${FIXTURE_DIR}" \
    --prefix ALL \
    --rows 100000
  cp "${SCRIPT_DIR}/fixture_all_100k/run_easystrata.R" "${FIXTURE_DIR}/run_easystrata.R"

  log "Building the pinned Ubuntu comparator image."
  docker build \
    --tag "${COMPARATOR_IMAGE}" \
    --file "$(docker_path "${SCRIPT_DIR}/Dockerfile.comparators")" \
    "$(docker_path "${SCRIPT_DIR}")"

  local docker_project docker_run
  docker_project="$(docker_path "${PROJECT_ROOT}")"
  docker_run="$(docker_path "${RUN_ROOT}")"
  log "Running the direct audit, GWAMA, and EasyStrata on the identical fixture."
  docker run --rm \
    --volume "${docker_project}:/project:ro" \
    --volume "${docker_run}:/results" \
    --workdir /results/fixture_all_100k \
    "${COMPARATOR_IMAGE}" bash -lc '
      set -euo pipefail
      /usr/bin/time -v -o multigwas.time.txt \
        perl /project/DiffGWASDeps/benchmark_diff_methods.pl \
          --input multigwas.tsv --prefix ALL --rho 0 --limit 100000 \
          --output multigwas_numerical.tsv \
          >multigwas.stdout.txt 2>multigwas.stderr.txt
      /usr/bin/time -v -o gwama.time.txt \
        GWAMA --filelist gwama.in --output gwama --quantitative --sex \
          >gwama.stdout.txt 2>gwama.stderr.txt
      /usr/bin/time -v -o easystrata.time.txt \
        Rscript run_easystrata.R \
          >easystrata.stdout.txt 2>easystrata.stderr.txt
      printf "GWAMA\t" > comparator_versions.tsv
      dpkg-query -W gwama | awk "{print \$2}" >> comparator_versions.tsv
      printf "\nEasyStrata\t" >> comparator_versions.tsv
      Rscript -e "cat(as.character(packageVersion(\"EasyStrata\")))" >> comparator_versions.tsv
      printf "\n" >> comparator_versions.tsv
    '

  perl "${DEPS_DIR}/compare_official_benchmarks.pl" \
    --expected "${FIXTURE_DIR}/expected.tsv" \
    --gwama "${FIXTURE_DIR}/gwama.out" \
    --easystrata "${FIXTURE_DIR}/EasyStrata.easystrata.results.txt" \
    --out "${FIXTURE_DIR}/official_comparison.tsv"
  perl "${DEPS_DIR}/compare_official_benchmarks_by_allele_class.pl" \
    --expected "${FIXTURE_DIR}/expected.tsv" \
    --gwama "${FIXTURE_DIR}/gwama.out" \
    --easystrata "${FIXTURE_DIR}/EasyStrata.easystrata.results.txt" \
    --out "${FIXTURE_DIR}/official_comparison_by_allele_class.tsv"

  local direct_elapsed direct_rss gwama_elapsed gwama_rss easy_elapsed easy_rss
  local direct_rows gwama_rows easy_rows
  direct_elapsed="$(elapsed_from_gnu_time "${FIXTURE_DIR}/multigwas.time.txt")"
  direct_rss="$(rss_from_gnu_time "${FIXTURE_DIR}/multigwas.time.txt")"
  gwama_elapsed="$(elapsed_from_gnu_time "${FIXTURE_DIR}/gwama.time.txt")"
  gwama_rss="$(rss_from_gnu_time "${FIXTURE_DIR}/gwama.time.txt")"
  easy_elapsed="$(elapsed_from_gnu_time "${FIXTURE_DIR}/easystrata.time.txt")"
  easy_rss="$(rss_from_gnu_time "${FIXTURE_DIR}/easystrata.time.txt")"
  direct_rows="$(awk -F '\t' '$1=="rows_compared" {print $2}' "${FIXTURE_DIR}/multigwas_numerical.tsv")"
  gwama_rows="$(data_rows "${FIXTURE_DIR}/gwama.out")"
  easy_rows="$(data_rows "${FIXTURE_DIR}/EasyStrata.easystrata.results.txt")"
  [[ "${direct_rows}" -eq 100000 && "${gwama_rows}" -eq 100000 && "${easy_rows}" -eq 100000 ]] \
    || die "One or more comparator paths did not retain all 100,000 fixture rows."

  {
    printf 'tool\tversion\tinput_fixture\tcommand\trows_input\trows_output\telapsed_seconds\tpeak_rss_kb\tstatus\tnotes\n'
    printf 'MultiGWAS-Explorer\tlocal rerun\tfixture_all_100k/multigwas.tsv\tperl DiffGWASDeps/benchmark_diff_methods.pl --input fixture_all_100k/multigwas.tsv --prefix ALL --rho 0 --limit 100000\t100000\t%s\t%s\t%s\tPASS\tDirect formula audit measured with /usr/bin/time -v\n' \
      "${direct_rows}" "${direct_elapsed}" "${direct_rss}"
    printf 'GWAMA\t2.2.2\tfixture_all_100k/female.gwama.tsv + male.gwama.tsv\tGWAMA --filelist gwama.in --output gwama --quantitative --sex\t100000\t%s\t%s\t%s\tPASS\tOfficial Ubuntu package measured with /usr/bin/time -v\n' \
      "${gwama_rows}" "${gwama_elapsed}" "${gwama_rss}"
    printf 'EasyStrata\t8.6\tfixture_all_100k/easystrata.tsv\tRscript run_easystrata.R (CALCPDIFF option 1)\t100000\t%s\t%s\t%s\tPASS\tOfficial package measured with /usr/bin/time -v\n' \
      "${easy_rows}" "${easy_elapsed}" "${easy_rss}"
  } > "${RUN_ROOT}/official_tool_runs.tsv"

  awk -F '\t' '
    NR==1 {next}
    $4!=100000 || $5!=0 || $10!=0 || $11!=0 || $12!=0 {bad=1}
    END {exit bad}
  ' "${FIXTURE_DIR}/official_comparison.tsv" \
    || die "Comparator acceptance criteria failed; inspect official_comparison.tsv."
  log "Comparator benchmark PASS: ${FIXTURE_DIR}/official_comparison.tsv"
}

run_full_qc() {
  require_command perl
  require_file "${MERGED_LONG_GZ}"
  local qc_input="${QC_STANDARDIZED_LONG_GZ}"
  local generated_qc_input="${RUN_ROOT}/qc_source/unfiltered_pairwise.stdized.tsv.gz"
  if [[ "${QC_INPUT_EXPLICIT}" -eq 0 && -s "${generated_qc_input}" ]]; then
    qc_input="${generated_qc_input}"
  fi
  require_file "${qc_input}"
  require_file "${GNU_TIME}"
  mkdir -p "${RUN_ROOT}/harmonization_audit" "${RUN_ROOT}/unfiltered_qc" "${RUN_ROOT}/unfiltered_qc_harmonized"

  log "Auditing all merged-long rows and pair/allele categories."
  "${GNU_TIME}" -v -o "${RUN_ROOT}/harmonization_audit.time.txt" \
    perl "${DEPS_DIR}/audit_allele_harmonization.pl" \
      --input "${MERGED_LONG_GZ}" \
      --outdir "${RUN_ROOT}/harmonization_audit"

  log "Running the pre-exclusion unfiltered QC sensitivity scan."
  "${GNU_TIME}" -v -o "${RUN_ROOT}/unfiltered_qc.time.txt" \
    perl "${DEPS_DIR}/qc_long_report.pl" \
      --input "${qc_input}" \
      --outdir "${RUN_ROOT}/unfiltered_qc" \
      --candidate-p 1e-5 \
      --no-exclude-strand-ambiguous \
      --max-eaf-abs-diff 0.20

  log "Running the harmonization-filtered unfiltered QC scan used in the revision."
  "${GNU_TIME}" -v -o "${RUN_ROOT}/unfiltered_qc_harmonized.time.txt" \
    perl "${DEPS_DIR}/qc_long_report.pl" \
      --input "${qc_input}" \
      --outdir "${RUN_ROOT}/unfiltered_qc_harmonized" \
      --candidate-p 1e-5 \
      --exclude-strand-ambiguous \
      --max-eaf-abs-diff 0.20

  local pre_seconds audit_seconds filtered_seconds
  pre_seconds="$(elapsed_from_gnu_time "${RUN_ROOT}/unfiltered_qc.time.txt")"
  audit_seconds="$(elapsed_from_gnu_time "${RUN_ROOT}/harmonization_audit.time.txt")"
  filtered_seconds="$(elapsed_from_gnu_time "${RUN_ROOT}/unfiltered_qc_harmonized.time.txt")"
  {
    printf 'TASK\tINPUT_ROWS\tELAPSED_SECONDS\tSTATUS\tNOTES\n'
    printf 'pre_exclusion_unfiltered_qc\t22292345\t%s\tPASS\tFull-stream raw-Z lambda and candidate sensitivity analysis\n' "${pre_seconds}"
    printf 'full_source_harmonization_audit\t44930352\t%s\tPASS\tPair-presence, allele category, ambiguity, duplicate, and frequency audit\n' "${audit_seconds}"
    printf 'harmonization_filtered_unfiltered_qc\t22292345\t%s\tPASS\tRevised ambiguity and EAF checks before lambda and candidate enumeration\n' "${filtered_seconds}"
  } > "${RUN_ROOT}/full_scan_runtime.tsv"

  awk -F '\t' '
    NR>1 {both+=$3; ambiguous+=$11; eligible+=$16}
    END {exit !(both==22292345 && ambiguous==3386012 && eligible==18906333)}
  ' "${RUN_ROOT}/harmonization_audit/harmonization_audit.tsv" \
    || die "Full-source harmonization counts differ from the report."
  log "Full-source and QC scans completed under ${RUN_ROOT}."
}

run_power() {
  require_file "${SCRIPT_DIR}/top_variant_strata_input.tsv"
  mkdir -p "${RUN_ROOT}"
  if command -v Rscript >/dev/null 2>&1; then
    Rscript "${SCRIPT_DIR}/calculate_detectable_effects.R" \
      "${SCRIPT_DIR}/top_variant_strata_input.tsv" \
      "${RUN_ROOT}/top_variant_detectable_effects.tsv"
  else
    require_command docker
    if ! docker image inspect "${COMPARATOR_IMAGE}" >/dev/null 2>&1; then
      docker build \
        --tag "${COMPARATOR_IMAGE}" \
        --file "$(docker_path "${SCRIPT_DIR}/Dockerfile.comparators")" \
        "$(docker_path "${SCRIPT_DIR}")"
    fi
    docker run --rm \
      --volume "$(docker_path "${PROJECT_ROOT}"):/project:ro" \
      --volume "$(docker_path "${RUN_ROOT}"):/results" \
      "${COMPARATOR_IMAGE}" \
      Rscript /project/benchmark/calculate_detectable_effects.R \
        /project/benchmark/top_variant_strata_input.tsv \
        /results/top_variant_detectable_effects.tsv
  fi
  log "Power results: ${RUN_ROOT}/top_variant_detectable_effects.tsv"
}

run_agent() {
  mkdir -p "${RUN_ROOT}/agent_interface"
  run_powershell_file "${SCRIPT_DIR}/run_agent_interface_benchmark.ps1" \
    -Runs 10 \
    -RecoveryRuns 3 \
    -Port "${AGENT_BENCHMARK_PORT:-18081}" \
    -OutputDir "$(powershell_path "${RUN_ROOT}/agent_interface")"
  if command -v Rscript >/dev/null 2>&1; then
    (
      cd "${RUN_ROOT}/agent_interface"
      Rscript "${SCRIPT_DIR}/agent_interface/summarize_agent_benchmark.R"
    )
  else
    require_command docker
    if ! docker image inspect "${COMPARATOR_IMAGE}" >/dev/null 2>&1; then
      docker build \
        --tag "${COMPARATOR_IMAGE}" \
        --file "$(docker_path "${SCRIPT_DIR}/Dockerfile.comparators")" \
        "$(docker_path "${SCRIPT_DIR}")"
    fi
    docker run --rm \
      --volume "$(docker_path "${PROJECT_ROOT}"):/project:ro" \
      --volume "$(docker_path "${RUN_ROOT}"):/results" \
      --workdir /results/agent_interface \
      "${COMPARATOR_IMAGE}" \
      Rscript /project/benchmark/agent_interface/summarize_agent_benchmark.R
  fi
  log "Agent summary tables: ${RUN_ROOT}/agent_interface"
}

generate_common_candidates() {
  require_command perl
  require_file "${WIDE_GZ}"
  require_file "${RUNNER_CONFIG}"
  mkdir -p "${RUN_ROOT}"
  perl "${DEPS_DIR}/generate_requested_top_hits_csv.pl" \
    --input "${WIDE_GZ}" \
    --output "${RUN_ROOT}/pgc_common_ld_candidates.csv" \
    --runner-config "${RUNNER_CONFIG}" \
    --top-hit-mode common_association \
    --top-hit-signal-thrshds 5e-8 \
    --top-hit-dist-bp 0 \
    --maf-threshold 0.01
}

generate_differential_candidates() {
  require_command perl
  require_file "${WIDE_GZ}"
  require_file "${RUNNER_CONFIG}"
  perl "${DEPS_DIR}/generate_requested_top_hits_csv.pl" \
    --input "${WIDE_GZ}" \
    --output "${RUN_ROOT}/PGC_SCZ_SAS_local_top_hits_manhattan_common_top_hits.csv" \
    --runner-config "${RUNNER_CONFIG}" \
    --top-hit-mode differential \
    --top-hit-focus-pvar ALL_DIFF_P \
    --top-hit-signal-thrshds 1e-6 \
    --top-hit-dist-bp 0 \
    --maf-threshold 0.01
}

write_archive_coverage() {
  local eur_json="${RUN_ROOT}/pgc_common_haploreg_candidate_ld.tsv.EUR.summary.json"
  local asn_json="${RUN_ROOT}/pgc_common_haploreg_candidate_ld.tsv.ASN.summary.json"
  local combined="${RUN_ROOT}/pgc_common_haploreg_candidate_ld.tsv"
  require_file "${eur_json}"
  require_file "${asn_json}"
  require_file "${combined}"
  EUR_JSON="${eur_json}" ASN_JSON="${asn_json}" COMBINED_CACHE="${combined}" \
  perl -MJSON::PP -e '
    sub loadj { open my $f,"<",$_[0] or die $!; local $/; decode_json(<$f>) }
    my $e=loadj($ENV{EUR_JSON}); my $a=loadj($ENV{ASN_JSON});
    open my $c,"<",$ENV{COMBINED_CACHE} or die $!; my $lines=-1; $lines++ while <$c>;
    print join("\t",qw(POPULATION ARCHIVE_ROWS_SCANNED UNIQUE_CANDIDATE_RSIDS MATCHED_QUERY_RSIDS MISSING_QUERY_RSIDS RETAINED_CANDIDATE_EDGES MINIMUM_R2 SOURCE_ARCHIVE)),"\n";
    for my $x ($e,$a) {
      print join("\t",$x->{population},$x->{archive_rows_scanned},$x->{candidate_count},$x->{matched_query_count},$x->{missing_query_count},$x->{retained_candidate_edges},$x->{minimum_r2},$x->{archive_source}),"\n";
    }
    print join("\t","COMBINED",$e->{archive_rows_scanned}+$a->{archive_rows_scanned},$e->{candidate_count},"NA","NA",$lines,$e->{minimum_r2},"Candidate-restricted_EUR_plus_ASN_cache"),"\n";
  ' > "${RUN_ROOT}/pgc_common_haploreg_archive_coverage.tsv"
}

run_archive_ld() {
  require_command perl
  require_command sha256sum
  require_file "${GNU_TIME}"
  generate_common_candidates
  local candidates="${RUN_ROOT}/pgc_common_ld_candidates.csv"
  local cache="${RUN_ROOT}/pgc_common_haploreg_candidate_ld.tsv"
  local sqlite="${RUN_ROOT}/pgc_common_haploreg_candidate_ld.sqlite"
  local signals="ALL_GROUP1_P ALL_GROUP2_P EUR_GROUP1_P EUR_GROUP2_P ASN_GROUP1_P ASN_GROUP2_P"

  log "Streaming EUR and ASN HaploReg archives in parallel (or using HAPLOREG_LD_ARCHIVE_DIR)."
  "${GNU_TIME}" -v -o "${RUN_ROOT}/haploreg_archive_extract.time.txt" \
    bash "${DEPS_DIR}/build_haploreg_ld_candidate_cache.sh" \
      --candidates "${candidates}" \
      --output "${cache}" \
      --populations "EUR ASN" \
      --min-r2 0.2
  write_archive_coverage

  log "Building SQLite and running the first complete archive-only clump."
  "${GNU_TIME}" -v -o "${RUN_ROOT}/archive_clump_first.time.txt" \
    perl "${DEPS_DIR}/clump_top_hits_with_haploreg_cache.pl" \
      --candidates "${candidates}" \
      --cache "${cache}" \
      --sqlite "${sqlite}" \
      --signal-columns "${signals}" \
      --signal-threshold 5e-8 \
      --populations "EUR ASN" \
      --r2-threshold 0.2 \
      --population-rule ANY \
      --fallback-distance-bp 1000000 \
      --rebuild-sqlite \
      --output-leads "${RUN_ROOT}/pgc_common_ld_cached_r2ge0p2_leads_local.csv" \
      --output-audit "${RUN_ROOT}/pgc_common_ld_cached_r2ge0p2_audit_local.tsv" \
      > "${RUN_ROOT}/archive_clump_first.stdout.txt"

  local repeat
  for repeat in 1 2; do
    "${GNU_TIME}" -v -o "${RUN_ROOT}/archive_clump_reuse${repeat}.time.txt" \
      perl "${DEPS_DIR}/clump_top_hits_with_haploreg_cache.pl" \
        --candidates "${candidates}" \
        --cache "${cache}" \
        --sqlite "${sqlite}" \
        --signal-columns "${signals}" \
        --signal-threshold 5e-8 \
        --populations "EUR ASN" \
        --r2-threshold 0.2 \
        --population-rule ANY \
        --fallback-distance-bp 1000000 \
        --output-leads "${RUN_ROOT}/pgc_common_ld_cached_r2ge0p2_leads_local_reuse${repeat}.csv" \
        --output-audit "${RUN_ROOT}/pgc_common_ld_cached_r2ge0p2_audit_local_reuse${repeat}.tsv" \
        > "${RUN_ROOT}/archive_clump_reuse${repeat}.stdout.txt"
  done

  sha256sum \
    "${RUN_ROOT}/pgc_common_ld_cached_r2ge0p2_leads_local"*.csv \
    "${RUN_ROOT}/pgc_common_ld_cached_r2ge0p2_audit_local"*.tsv \
    > "${RUN_ROOT}/archive_clump_sha256.txt"
  cmp -s \
    "${RUN_ROOT}/pgc_common_ld_cached_r2ge0p2_leads_local.csv" \
    "${RUN_ROOT}/pgc_common_ld_cached_r2ge0p2_leads_local_reuse1.csv" \
    || die "First indexed lead rerun was not byte-identical."
  cmp -s \
    "${RUN_ROOT}/pgc_common_ld_cached_r2ge0p2_audit_local.tsv" \
    "${RUN_ROOT}/pgc_common_ld_cached_r2ge0p2_audit_local_reuse1.tsv" \
    || die "First indexed audit rerun was not byte-identical."
  log "Archive/cache benchmark complete."
}

oda() {
  require_command perl
  require_file "${SAS_RUNNER}"
  perl "${SAS_RUNNER}" "$@"
}

upload_sas_dependencies() {
  require_file "${QUERY_LD_MACRO}"
  require_file "${GET_TOP_SIGNAL_LD_MACRO}"
  require_file "${GET_TOP_SIGNAL_DISTANCE_MACRO}"
  mkdir -p "${RUN_ROOT}/sas_dependency_upload"
  oda \
    --upload-file "${QUERY_LD_MACRO}" \
    --upload-file "${GET_TOP_SIGNAL_LD_MACRO}" \
    --upload-file "${GET_TOP_SIGNAL_DISTANCE_MACRO}" \
    --persistent --session-id "${SAS_SESSION_ID}" \
    --output-prefix "${RUN_ROOT}/sas_dependency_upload/output"
}

run_sas_ld() {
  require_file "${QUERY_MULTI_LD_MACRO}"
  mkdir -p "${RUN_ROOT}/sas_haploreg_high_ld" "${RUN_ROOT}/sas_haploreg_pairwise" \
           "${RUN_ROOT}/sas_fixture_ld" "${RUN_ROOT}/sas_real_diff_ld" \
           "${RUN_ROOT}/sas_dependencies" "${RUN_ROOT}/sas_query_multi_upload"
  upload_sas_dependencies
  cp "${QUERY_MULTI_LD_MACRO}" "${RUN_ROOT}/sas_dependencies/QueryMulti_LD_SNPs_at_Haploreg4.sas"
  oda --upload-file "${RUN_ROOT}/sas_dependencies/QueryMulti_LD_SNPs_at_Haploreg4.sas" \
    --persistent --session-id "${SAS_SESSION_ID}" \
    --output-prefix "${RUN_ROOT}/sas_query_multi_upload/output"

  oda --file "${SCRIPT_DIR}/run_haploreg_ld_benchmark.sas" \
    --persistent --session-id "${SAS_SESSION_ID}" --run-timeout-seconds 3600 \
    --output-prefix "${RUN_ROOT}/sas_haploreg_high_ld/output"
  oda \
    --download-file '~/haploreg_ld_eur_r2ge0p5.tsv' \
    --download-file '~/haploreg_ld_asn_r2ge0p5.tsv' \
    --download-local-path "${RUN_ROOT}/haploreg_ld_eur_r2ge0p5.tsv" \
    --download-local-path "${RUN_ROOT}/haploreg_ld_asn_r2ge0p5.tsv" \
    --persistent --session-id "${SAS_SESSION_ID}" \
    --output-prefix "${RUN_ROOT}/sas_haploreg_high_ld/download"

  oda --file "${SCRIPT_DIR}/run_haploreg_pairwise_ld_benchmark.sas" \
    --persistent --session-id "${SAS_SESSION_ID}" --run-timeout-seconds 3600 \
    --output-prefix "${RUN_ROOT}/sas_haploreg_pairwise/output"
  oda \
    --download-file '~/haploreg_pairwise_ld.tsv' \
    --download-local-path "${RUN_ROOT}/haploreg_pairwise_ld.tsv" \
    --persistent --session-id "${SAS_SESSION_ID}" \
    --output-prefix "${RUN_ROOT}/sas_haploreg_pairwise/download"

  oda --file "${SCRIPT_DIR}/test_ld_independent_top_hits.sas" \
    --persistent --session-id "${SAS_SESSION_ID}" --run-timeout-seconds 3600 \
    --output-prefix "${RUN_ROOT}/sas_fixture_ld/output"
  oda \
    --download-file '~/ld_fixture_independent_leads.tsv' \
    --download-file '~/ld_fixture_audit.tsv' \
    --download-file '~/ld_fixture_distance_leads.tsv' \
    --download-local-path "${RUN_ROOT}/ld_fixture_independent_leads.tsv" \
    --download-local-path "${RUN_ROOT}/ld_fixture_audit.tsv" \
    --download-local-path "${RUN_ROOT}/ld_fixture_distance_leads.tsv" \
    --persistent --session-id "${SAS_SESSION_ID}" \
    --output-prefix "${RUN_ROOT}/sas_fixture_ld/download"

  generate_differential_candidates
  oda --upload-file "${RUN_ROOT}/PGC_SCZ_SAS_local_top_hits_manhattan_common_top_hits.csv" \
    --persistent --session-id "${SAS_SESSION_ID}" \
    --output-prefix "${RUN_ROOT}/sas_real_diff_ld/upload"
  oda --file "${SCRIPT_DIR}/run_real_pgc_diff_ld_independence.sas" \
    --persistent --session-id "${SAS_SESSION_ID}" --run-timeout-seconds 3600 \
    --output-prefix "${RUN_ROOT}/sas_real_diff_ld/output"
  oda \
    --download-file '~/pgc_diff_ld_independent_leads.tsv' \
    --download-file '~/pgc_diff_ld_audit.tsv' \
    --download-local-path "${RUN_ROOT}/pgc_diff_ld_independent_leads.tsv" \
    --download-local-path "${RUN_ROOT}/pgc_diff_ld_audit.tsv" \
    --persistent --session-id "${SAS_SESSION_ID}" \
    --output-prefix "${RUN_ROOT}/sas_real_diff_ld/download"
  log "Nominated-SNP and differential SAS LD stages completed."
}

run_sas_common() {
  mkdir -p "${RUN_ROOT}/sas_common_hybrid"
  [[ -s "${RUN_ROOT}/pgc_common_ld_candidates.csv" ]] || generate_common_candidates
  require_file "${RUN_ROOT}/pgc_common_haploreg_candidate_ld.tsv"
  upload_sas_dependencies
  oda \
    --upload-file "${RUN_ROOT}/pgc_common_ld_candidates.csv" \
    --upload-file "${RUN_ROOT}/pgc_common_haploreg_candidate_ld.tsv" \
    --persistent --session-id "${SAS_SESSION_ID}" \
    --run-timeout-seconds 7200 \
    --output-prefix "${RUN_ROOT}/sas_common_hybrid/upload"
  oda --file "${SCRIPT_DIR}/run_real_pgc_common_ld_cached_r2ge0p2.sas" \
    --persistent --session-id "${SAS_SESSION_ID}" --run-timeout-seconds 7200 \
    --output-prefix "${RUN_ROOT}/sas_common_hybrid/output"
  oda \
    --download-file '~/pgc_common_ld_cached_r2ge0p2_leads.tsv' \
    --download-file '~/pgc_common_ld_cached_r2ge0p2_audit.tsv' \
    --download-local-path "${RUN_ROOT}/pgc_common_ld_cached_r2ge0p2_leads.tsv" \
    --download-local-path "${RUN_ROOT}/pgc_common_ld_cached_r2ge0p2_audit.tsv" \
    --persistent --session-id "${SAS_SESSION_ID}" \
    --output-prefix "${RUN_ROOT}/sas_common_hybrid/download"
  log "SAS archive-cache plus live-web fallback stage completed."
}

audit_action_count() {
  local file="$1" action="$2"
  awk -F '\t' -v wanted="${action}" '
    NR==1 {for(i=1;i<=NF;i++){gsub(/\r/,"",$i); h[$i]=i}; next}
    {v=$h["selection_action"]; gsub(/\r/,"",$v); if($v==wanted)n++}
    END {print n+0}
  ' "${file}"
}

audit_status_count() {
  local file="$1" pattern="$2"
  awk -F '\t' -v pattern="${pattern}" '
    NR==1 {for(i=1;i<=NF;i++){gsub(/\r/,"",$i); h[$i]=i}; next}
    {
      action=$h["selection_action"]; status=$h["query_status"];
      gsub(/\r/,"",action); gsub(/\r/,"",status);
      if(action=="SELECTED_LEAD" && status ~ pattern)n++
    }
    END {print n+0}
  ' "${file}"
}

write_ld_summary() {
  local fixture_audit="${RUN_ROOT}/ld_fixture_audit.tsv"
  local fixture_distance="${RUN_ROOT}/ld_fixture_distance_leads.tsv"
  local diff_audit="${RUN_ROOT}/pgc_diff_ld_audit.tsv"
  local archive_audit="${RUN_ROOT}/pgc_common_ld_cached_r2ge0p2_audit_local.tsv"
  local hybrid_audit="${RUN_ROOT}/pgc_common_ld_cached_r2ge0p2_audit.tsv"
  local required
  for required in "${fixture_audit}" "${fixture_distance}" "${diff_audit}" "${archive_audit}" "${hybrid_audit}"; do
    require_file "${required}"
  done

  local fixture_selected fixture_pruned fixture_complete distance_selected
  local diff_candidates diff_selected diff_pruned diff_complete
  local archive_candidates archive_selected archive_pruned archive_fallback archive_complete archive_partial archive_failed
  local hybrid_candidates hybrid_selected hybrid_pruned hybrid_fallback hybrid_complete hybrid_partial hybrid_failed
  fixture_selected="$(audit_action_count "${fixture_audit}" SELECTED_LEAD)"
  fixture_pruned="$(audit_action_count "${fixture_audit}" PRUNED_LD)"
  fixture_complete="$(audit_status_count "${fixture_audit}" '^OK')"
  distance_selected="$(data_rows "${fixture_distance}")"
  diff_candidates="$(data_rows "${diff_audit}")"
  diff_selected="$(audit_action_count "${diff_audit}" SELECTED_LEAD)"
  diff_pruned="$(audit_action_count "${diff_audit}" PRUNED_LD)"
  diff_complete="$(audit_status_count "${diff_audit}" '^OK')"
  archive_candidates="$(data_rows "${archive_audit}")"
  archive_selected="$(audit_action_count "${archive_audit}" SELECTED_LEAD)"
  archive_pruned="$(audit_action_count "${archive_audit}" PRUNED_LD)"
  archive_fallback="$(audit_action_count "${archive_audit}" PRUNED_DISTANCE_FALLBACK)"
  archive_complete="$(audit_status_count "${archive_audit}" '^OK_LOCAL_CACHE$')"
  archive_partial="$(audit_status_count "${archive_audit}" '^PARTIAL_LOCAL_CACHE$')"
  archive_failed="$(audit_status_count "${archive_audit}" '^NO_LD_RESPONSE$')"
  hybrid_candidates="$(data_rows "${hybrid_audit}")"
  hybrid_selected="$(audit_action_count "${hybrid_audit}" SELECTED_LEAD)"
  hybrid_pruned="$(audit_action_count "${hybrid_audit}" PRUNED_LD)"
  hybrid_fallback="$(audit_action_count "${hybrid_audit}" PRUNED_DISTANCE_FALLBACK)"
  hybrid_complete="$(audit_status_count "${hybrid_audit}" '^OK')"
  hybrid_partial="$(audit_status_count "${hybrid_audit}" '^PARTIAL')"
  hybrid_failed="$(audit_status_count "${hybrid_audit}" '^(NO_LD_RESPONSE|FAILED)')"

  {
    printf 'BENCHMARK\tMODE\tSIGNAL_THRESHOLD\tLD_POPULATIONS\tLD_R2_THRESHOLD\tLD_POPULATION_RULE\tCANDIDATES\tSELECTED_LEADS\tPRUNED_LD\tPRUNED_DISTANCE_FALLBACK\tQUERY_COMPLETE_LEADS\tQUERY_PARTIAL_LEADS\tQUERY_FAILED_LEADS\tNOTES\n'
    printf 'four_snp_fixture_distance\tdifferential\t1e-6\tNA\tNA\tNA\t4\t%s\t0\t%s\tNA\tNA\tNA\tLegacy 1-Mb physical-distance selector\n' "${distance_selected}" "$((4-distance_selected))"
    printf 'four_snp_fixture_ld\tdifferential\t1e-6\tEUR\t0.5\tANY\t4\t%s\t%s\t0\t%s\t0\t0\tFour-SNP HaploReg regression fixture\n' "${fixture_selected}" "${fixture_pruned}" "${fixture_complete}"
    printf 'pgc_real_differential_ld\tdifferential\t1e-6\tEUR ASN\t0.1\tANY\t%s\t%s\t%s\t0\t%s\t0\t0\tMAF-passing real PGC differential candidates\n' "${diff_candidates}" "${diff_selected}" "${diff_pruned}" "${diff_complete}"
    printf 'pgc_real_common_archive_cache_local\tcommon_association\t5e-8\tEUR ASN\t0.2\tANY\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tArchive-only SQLite clump; timing and hashes are retained beside this table\n' \
      "${archive_candidates}" "${archive_selected}" "${archive_pruned}" "${archive_fallback}" "${archive_complete}" "${archive_partial}" "${archive_failed}"
    printf 'pgc_real_common_archive_cache_sas_hybrid\tcommon_association\t5e-8\tEUR ASN\t0.2\tANY\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tSAS ODA cache plus live-web fallback; timing is retained in wrapper logs\n' \
      "${hybrid_candidates}" "${hybrid_selected}" "${hybrid_pruned}" "${hybrid_fallback}" "${hybrid_complete}" "${hybrid_partial}" "${hybrid_failed}"
  } > "${RUN_ROOT}/ld_top_hit_benchmark_summary.tsv"
  log "Consolidated LD summary: ${RUN_ROOT}/ld_top_hit_benchmark_summary.tsv"
}

run_web_only_control() {
  mkdir -p "${RUN_ROOT}/sas_common_web_only"
  [[ -s "${RUN_ROOT}/pgc_common_ld_candidates.csv" ]] || generate_common_candidates
  upload_sas_dependencies
  oda --upload-file "${RUN_ROOT}/pgc_common_ld_candidates.csv" \
    --persistent --session-id "${SAS_SESSION_ID}" \
    --output-prefix "${RUN_ROOT}/sas_common_web_only/upload"
  log "Starting the documented negative control; non-zero/timeout is expected."
  set +e
  oda --file "${SCRIPT_DIR}/run_real_pgc_common_ld_independence.sas" \
    --persistent --session-id "${SAS_SESSION_ID}" \
    --run-timeout-seconds "${WEB_ONLY_TIMEOUT_SECONDS:-3458}" \
    --output-prefix "${RUN_ROOT}/sas_common_web_only/output"
  local status=$?
  set -e
  printf 'exit_status\t%s\nexpected_result\tnoncompletion_or_timeout\ntimeout_seconds\t%s\n' \
    "${status}" "${WEB_ONLY_TIMEOUT_SECONDS:-3458}" \
    > "${RUN_ROOT}/sas_common_web_only/negative_control_status.tsv"
  if [[ "${status}" -eq 0 ]]; then
    log "WARNING: web-only control completed in this rerun; inspect and update the report rather than assuming the historical timeout."
  else
    log "Web-only negative control ended non-zero as expected; status=${status}."
  fi
}

validate_results() {
  "${BASH:-bash}" "${SCRIPT_DIR}/validate_revision_benchmark.sh"
}

target="${1:-help}"
mkdir -p "${RUN_ROOT}"
case "${target}" in
  help|-h|--help) usage ;;
  show-inputs) show_inputs ;;
  prepare-pipeline) prepare_pipeline ;;
  prepare-qc-source) prepare_qc_source ;;
  comparators) run_comparators ;;
  full-qc) run_full_qc ;;
  power) run_power ;;
  agent) run_agent ;;
  ld-archive) run_archive_ld ;;
  ld-sas) run_sas_ld ;;
  ld-sas-common) run_sas_common ;;
  summarize-ld) write_ld_summary ;;
  web-only-control) run_web_only_control ;;
  validate-results|validate-published) validate_results ;;
  local-report)
    run_comparators
    run_full_qc
    run_power
    run_agent
    run_archive_ld
    ;;
  all)
    run_comparators
    run_full_qc
    run_power
    run_agent
    run_archive_ld
    run_sas_ld
    run_sas_common
    write_ld_summary
    ;;
  *) die "Unknown target: ${target}. Run with --help." ;;
esac
