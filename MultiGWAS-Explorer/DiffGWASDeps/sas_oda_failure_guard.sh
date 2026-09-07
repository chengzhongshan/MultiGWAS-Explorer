#!/usr/bin/env bash

# Shared, side-effect-free SAS ODA failure classification helpers.  Callers
# decide which exit code to use after a non-retryable failure is reported.

sas_oda_log_has_terminal_sas_error() {
  local logfile="${1:-}"
  [[ -n "${logfile}" && -s "${logfile}" ]] || return 1

  # Match SAS ERROR: and numbered ERROR nnn-nnn: records at the start of a
  # real log line. Numbered submitted-source echoes such as "123 ERROR: ..."
  # deliberately do not match.
  grep -Eq '^[[:space:]]*ERROR([[:space:]]+[0-9]+-[0-9]+)?:' "${logfile}"
}

sas_oda_log_has_space_exhaustion() {
  local logfile="${1:-}"
  [[ -n "${logfile}" && -s "${logfile}" ]] || return 1

  # Anchor SAS ERROR records so an ERROR example echoed as submitted source
  # code does not trigger the guard.  UNIX errno 28 is ENOSPC and 122 is the
  # Linux disk-quota error used by some ODA hosts.
  grep -Eiq \
    '^[[:space:]]*ERROR:[[:space:]]+(Insufficient space in|.*(No space left on device|Disk quota exceeded|file system is full|disk is full|quota[^[:cntrl:]]*exceed))|^[[:space:]]*(ERROR:[[:space:]]*)?UNIX errno[[:space:]]*=[[:space:]]*(28|122)([^0-9]|$)' \
    "${logfile}"
}

sas_oda_write_space_exhaustion_marker() {
  local logfile="${1:?SAS log path is required}"
  local context="${2:-SAS ODA job}"
  local marker="${logfile}.non_retryable_space_failure.txt"
  local detected_at
  detected_at="$(date -u '+%Y-%m-%dT%H:%M:%SZ' 2>/dev/null || date)"

  {
    printf 'failure_class=sas_oda_space_exhaustion\n'
    printf 'retryable=false\n'
    printf 'detected_at=%s\n' "${detected_at}"
    printf 'context=%s\n' "${context}"
    printf 'sas_log=%s\n' "${logfile}"
    printf 'recommended_action=Free ODA storage or reduce WORK data before starting a new run.\n'
  } > "${marker}"
  printf '%s\n' "${marker}"
}

sas_oda_report_space_exhaustion() {
  local logfile="${1:?SAS log path is required}"
  local context="${2:-SAS ODA job}"
  local marker

  sas_oda_log_has_space_exhaustion "${logfile}" || return 1
  marker="$(sas_oda_write_space_exhaustion_marker "${logfile}" "${context}")"

  printf '%s\n' \
    'ERROR: SAS ODA exhausted its WORK/storage space (for example, "ERROR: Insufficient space in file WORK....").' \
    "ERROR: ${context} is classified as NON-RETRYABLE with the same script and inputs." \
    'ERROR: Stopping immediately instead of submitting the failed SAS script again.' \
    "ERROR: SAS diagnostic log preserved at: ${logfile}" \
    "ERROR: Failure marker written to: ${marker}" \
    'ERROR: Free ODA storage, reduce the uploaded/local subset, or use the low-WORK implementation before rerunning.' >&2
  return 0
}

sas_oda_report_native_abort() {
  local exit_code="${1:-0}"
  local status_file="${2:-}"
  local context="${3:-SAS ODA job}"
  [[ "${exit_code}" -eq 134 ]] || return 1

  local marker
  if [[ -n "${status_file}" ]]; then
    marker="${status_file}.non_retryable_infrastructure_abort.txt"
  else
    marker="sas_oda.non_retryable_infrastructure_abort.txt"
  fi
  {
    printf 'failure_class=sas_oda_infrastructure_abort\n'
    printf 'retryable=false\n'
    printf 'exit_code=134\n'
    printf 'signal=SIGABRT\n'
    printf 'context=%s\n' "${context}"
    [[ -n "${status_file}" ]] && printf 'last_status=%s\n' "${status_file}"
    printf 'recommended_action=Do not resubmit the incomplete ODA job automatically; use the local gnuplot fallback or start a fresh ODA session manually.\n'
  } > "${marker}"

  printf '%s\n' \
    'ERROR: The local SAS ODA/SASPy process terminated with exit 134 (SIGABRT).' \
    "ERROR: ${context} is classified as a NON-RETRYABLE infrastructure abort for this invocation." \
    'ERROR: Stopping immediately instead of blindly resubmitting a possibly still-running remote SAS job.' \
    "ERROR: Failure marker written to: ${marker}" >&2
  return 0
}

# A remote SAS process can disappear without returning a SAS ERROR line.  This
# is commonly seen when ODA kills a session for WORK/quota or memory pressure,
# but the client cannot prove the exact server-side cause.  Keep this distinct
# from the explicit disk-full code so callers can diagnose it without blindly
# resubmitting the same large job.
sas_oda_log_has_remote_termination() {
  local logfile="${1:-}"
  [[ -n "${logfile}" && -s "${logfile}" ]] || return 1
  grep -Eiq 'No SAS process attached|SAS process has terminated unexpectedly|SAS submit returned empty output and the SAS session was no longer usable|session server.*(terminated|closed)|remote SAS.*(terminated|disconnected)' "${logfile}"
}

sas_oda_write_remote_termination_marker() {
  local logfile="${1:-}"
  local status_file="${2:-}"
  local context="${3:-SAS ODA job}"
  local marker
  if [[ -n "${status_file}" ]]; then
    marker="${status_file}.non_retryable_remote_termination.txt"
  elif [[ -n "${logfile}" ]]; then
    marker="${logfile}.non_retryable_remote_termination.txt"
  else
    marker="sas_oda.non_retryable_remote_termination.txt"
  fi
  {
    printf 'failure_class=sas_oda_remote_session_termination\n'
    printf 'error_code=SAS_ODA_REMOTE_SESSION_TERMINATED\n'
    printf 'exit_code=%s\n' "${SAS_ODA_REMOTE_TERMINATION_EXIT_CODE:-74}"
    printf 'retryable=false\n'
    printf 'context=%s\n' "${context}"
    [[ -n "${logfile}" ]] && printf 'sas_log=%s\n' "${logfile}"
    [[ -n "${status_file}" ]] && printf 'sas_status=%s\n' "${status_file}"
    printf 'diagnosis=ODA terminated the remote SAS session without a definitive SAS log; WORK/quota or memory exhaustion is possible but unconfirmed.\n'
    printf 'recommended_action=Preserve diagnostics, reduce the input or WORK footprint, and start a fresh ODA session before rerunning.\n'
  } > "${marker}"
  printf '%s\n' "${marker}"
}

sas_oda_report_remote_termination() {
  local logfile="${1:-}"
  local status_file="${2:-}"
  local context="${3:-SAS ODA job}"
  local exit_code="${4:-0}"
  local detected=0
  [[ "${exit_code}" =~ ^[0-9]+$ && "${exit_code}" -eq "${SAS_ODA_REMOTE_TERMINATION_EXIT_CODE:-74}" ]] && detected=1
  sas_oda_log_has_remote_termination "${logfile}" && detected=1
  if [[ "${detected}" -eq 0 && -n "${status_file}" && -s "${status_file}" ]]; then
    grep -Eiq 'sas_oda_remote_session_termination|No SAS process attached|SAS process has terminated unexpectedly|session was no longer usable' "${status_file}" && detected=1
  fi
  [[ "${detected}" -eq 1 ]] || return 1
  local marker
  marker="$(sas_oda_write_remote_termination_marker "${logfile}" "${status_file}" "${context}")"
  printf '%s\n' \
    'ERROR: The remote SAS ODA session terminated without a definitive SAS log.' \
    "ERROR: ${context} is classified as NON-RETRYABLE (possible WORK/quota or memory exhaustion; unconfirmed)." \
    'ERROR: The same large submission will not be retried automatically.' \
    "ERROR: Failure marker written to: ${marker}" >&2
  return 0
}
