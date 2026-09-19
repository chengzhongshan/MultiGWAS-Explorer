#!/usr/bin/env bash
# Repair repo-local modules after a Cygwin Perl upgrade, then run local tests.
set -euo pipefail

SCRIPT_DIR="$(cd "$(/usr/bin/dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=install/common.sh
. "${SCRIPT_DIR}/common.sh"

command_exists cygpath \
  || die "install/repair_and_test_cygwin.sh must run inside Cygwin"

log_file="${PIPELINE_ROOT}/cygwin-repair-test.log"
exec > >(tee "${log_file}") 2>&1

log "Cygwin repair log: ${log_file}"
log "Perl version: $(perl -e 'print $^V')"
log "Perl archname: $(perl -MConfig -e 'print $Config{archname}')"

activate_perl_env
activate_python_env
ensure_cpanm

# Repair modules that can prevent cpanm or the main cpanfile pass from loading.
install_cygwin_legacy_perl_deps
install_pdl_perl_deps

# Fill in all remaining pipeline and server dependencies, then validate them.
install_perl_deps
ensure_local_hts_tools
bash "${SCRIPT_DIR}/check_pipeline_install.sh"
bash "${SCRIPT_DIR}/run_plotting_example.sh" \
  "${PIPELINE_ROOT}/example-output-cygwin-repair"

log "Cygwin repair and local pipeline tests completed successfully"
