#!/usr/bin/env bash
set -euo pipefail
check_dir="$(cd "$(/usr/bin/dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
cd "${check_dir}/.."
runtime_python="${PIPELINE_PYTHON_BIN:-python3}"
export PIPELINE_PYTHON_BIN="${runtime_python}"
perl_probe="$(perl -e 'print "EXPLORER_PERL_RUNTIME_READY\n"')"
if [ "${perl_probe}" != EXPLORER_PERL_RUNTIME_READY ]; then
  printf 'ERROR: Perl did not execute the runtime probe. Check PATH; use Perl from the same Cygwin installation as Bash.\n' >&2
  exit 1
fi

# Offline regressions: no ODA credentials, Java session or scientific downloads.
for test in \
  test_saspy_python_discovery test_poc_input_boundaries test_runner_preflight_parity \
  test_persistent_bulk_dispatch test_one_shot_status_path \
  test_python_path_boundary test_sas_supplied_autocall \
  test_literal_execution_contexts test_deferred_include_paths \
  test_deferred_macro_literals test_double_quoted_deferred_includes \
  test_include_abbreviation test_include_file_lists test_include_source_order \
  test_resolve_macro_dependencies test_unquote_macro_dependencies \
  test_macro_quoting_preflight test_inline_macro_preflight \
  test_open_code_macro_diagnostics test_explicit_sas_target test_local_sas_home test_sas_com_threading \
  test_macro_span_scaling test_nested_include_graph test_copilot_runner_review \
  test_quoted_filename_upload test_sas_log_diagnostics \
  test_expected_failure_evaluation test_java_setup_help test_mcp_bounded_completion \
  test_runner_alias_options test_runtime_port_manifest; do
  perl "tests/${test}.pl"
done
"${runtime_python}" tests/test_persistent_archive.py
"${runtime_python}" tests/test_result_browser.py
for test in test_cygwin_python_handoff test_macroparas_autoload test_macroparas_gtex test_macroparas_html; do
  "${runtime_python}" "tests/${test}.py"
done
perl tests/test_macroparas_autoload.pl
perl tests/test_macroparas_gtex.pl
"${runtime_python}" tools/sync_macroparas_bundle.py --check
"${runtime_python}" tests/test_macroquery_port_manifest.py
"${runtime_python}" tests/test_efficient_sas_submit.py
perl DiffGWASDeps/test_sas_oda_macro_dir_fallback.pl
perl DiffGWASDeps/test_sas_oda_include_lookup.pl
"${runtime_python}" - <<'PY'
from pathlib import Path
# Apptainer runs this check against a read-only SIF; do not write __pycache__.
for name in ('DiffGWASDeps/sas_oda_session_server.py', 'DiffGWASDeps/sas_archive_transfer.py', 'DiffGWASDeps/result_browser.py'):
    compile(Path(name).read_bytes(), name, 'exec')
print('Python runtime syntax: PASS (no cache writes)')
PY
perl -c tools/profile_sas_source.pl
for runner in run_sas_codes_or_script_in_ODA.pl run_sas_codes_or_files_in_ODA.pl \
  run_sas_codes_or_file_in_ODA.pl run_sas_codes_or_files_in_SAS_ODA.pl; do
  perl -c "$runner"
done
printf 'SAS runtime offline regressions: PASS\n'
