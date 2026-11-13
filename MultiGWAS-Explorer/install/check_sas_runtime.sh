#!/usr/bin/env bash
set -euo pipefail
check_dir="$(cd "$(/usr/bin/dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
cd "${check_dir}/.."
runtime_python="${PIPELINE_PYTHON_BIN:-python3}"
export PIPELINE_PYTHON_BIN="${runtime_python}"

# Offline regressions: no ODA credentials, Java session or scientific downloads.
for test in \
  test_poc_input_boundaries test_runner_preflight_parity \
  test_persistent_bulk_dispatch test_one_shot_status_path \
  test_python_path_boundary test_sas_supplied_autocall \
  test_literal_execution_contexts test_deferred_include_paths \
  test_deferred_macro_literals test_double_quoted_deferred_includes \
  test_include_abbreviation test_include_file_lists test_include_source_order \
  test_resolve_macro_dependencies test_unquote_macro_dependencies \
  test_macro_quoting_preflight test_inline_macro_preflight \
  test_open_code_macro_diagnostics test_local_sas_home test_sas_com_threading; do
  perl "tests/${test}.pl"
done
"${runtime_python}" tests/test_persistent_archive.py
"${runtime_python}" -m py_compile DiffGWASDeps/sas_oda_session_server.py DiffGWASDeps/sas_archive_transfer.py
perl -c tools/profile_sas_source.pl
printf 'SAS runtime offline regressions: PASS\n'
