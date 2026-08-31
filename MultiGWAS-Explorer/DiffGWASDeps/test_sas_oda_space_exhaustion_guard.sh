#!/usr/bin/env bash
set -euo pipefail

DEPS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
# shellcheck source=sas_oda_failure_guard.sh
source "${DEPS_DIR}/sas_oda_failure_guard.sh"

TEST_TMP="$(mktemp -d)"
cleanup() {
  rm -f \
    "${TEST_TMP}/space.log" \
    "${TEST_TMP}/space.log.non_retryable_space_failure.txt" \
    "${TEST_TMP}/quota.log" \
    "${TEST_TMP}/compile.log" \
    "${TEST_TMP}/source_echo.log" \
    "${TEST_TMP}/timeout.log" \
    "${TEST_TMP}/classify_space.json" \
    "${TEST_TMP}/classify_timeout.json" \
    "${TEST_TMP}/native.run.status.json.non_retryable_infrastructure_abort.txt" \
    "${TEST_TMP}/native_warning.txt" \
    "${TEST_TMP}/warning.txt"
  rmdir "${TEST_TMP}" 2>/dev/null || true
}
trap cleanup EXIT

printf '%s\n' 'ERROR: Insufficient space in file WORK.FINAL.DATA.' > "${TEST_TMP}/space.log"
printf '%s\n' 'ERROR: UNIX errno = 122 (Disk quota exceeded).' > "${TEST_TMP}/quota.log"
printf '%s\n' 'ERROR 180-322: Statement is not valid or it is used out of proper order.' > "${TEST_TMP}/compile.log"
printf '%s\n' '123  ERROR: Insufficient space in file WORK.EXAMPLE.DATA.' > "${TEST_TMP}/source_echo.log"
printf '%s\n' 'SAS submit timed out' > "${TEST_TMP}/timeout.log"

sas_oda_log_has_space_exhaustion "${TEST_TMP}/space.log"
sas_oda_log_has_space_exhaustion "${TEST_TMP}/quota.log"
! sas_oda_log_has_space_exhaustion "${TEST_TMP}/source_echo.log"
! sas_oda_log_has_space_exhaustion "${TEST_TMP}/timeout.log"
sas_oda_log_has_terminal_sas_error "${TEST_TMP}/space.log"
sas_oda_log_has_terminal_sas_error "${TEST_TMP}/compile.log"
! sas_oda_log_has_terminal_sas_error "${TEST_TMP}/source_echo.log"

sas_oda_report_space_exhaustion \
  "${TEST_TMP}/space.log" \
  'synthetic local-GTF test' 2> "${TEST_TMP}/warning.txt"
test -s "${TEST_TMP}/space.log.non_retryable_space_failure.txt"
grep -q 'retryable=false' "${TEST_TMP}/space.log.non_retryable_space_failure.txt"
grep -q 'NON-RETRYABLE' "${TEST_TMP}/warning.txt"

! sas_oda_report_native_abort 1 "${TEST_TMP}/native.run.status.json" 'synthetic native-abort test'
sas_oda_report_native_abort \
  134 \
  "${TEST_TMP}/native.run.status.json" \
  'synthetic native-abort test' 2> "${TEST_TMP}/native_warning.txt"
test -s "${TEST_TMP}/native.run.status.json.non_retryable_infrastructure_abort.txt"
grep -q 'failure_class=sas_oda_infrastructure_abort' "${TEST_TMP}/native.run.status.json.non_retryable_infrastructure_abort.txt"
grep -q 'retryable=false' "${TEST_TMP}/native.run.status.json.non_retryable_infrastructure_abort.txt"
grep -q 'SIGABRT' "${TEST_TMP}/native_warning.txt"

set +e
perl "${DEPS_DIR}/../run_sas_codes_or_script_in_ODA.pl" \
  --classify-sas-log "${TEST_TMP}/space.log" > "${TEST_TMP}/classify_space.json"
CLASSIFY_SPACE_RC=$?
set -e
test "${CLASSIFY_SPACE_RC}" -eq 73
grep -q '"failure_class":"sas_oda_space_exhaustion"' "${TEST_TMP}/classify_space.json"
grep -q '"retryable":false' "${TEST_TMP}/classify_space.json"

perl "${DEPS_DIR}/../run_sas_codes_or_script_in_ODA.pl" \
  --classify-sas-log "${TEST_TMP}/timeout.log" > "${TEST_TMP}/classify_timeout.json"
grep -q '"retryable":true' "${TEST_TMP}/classify_timeout.json"

echo 'SAS_ODA_SPACE_EXHAUSTION_GUARD_TEST_PASSED'
