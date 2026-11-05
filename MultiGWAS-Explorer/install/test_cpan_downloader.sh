#!/usr/bin/env bash
set -euo pipefail
export PATH=/usr/bin:$PATH
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/common.sh"
temporary_root="$(cd "${TMPDIR:-/tmp}" && pwd -P)"
fixture="$(mktemp -d "$temporary_root/pipeline-cpan-test.XXXXXX")"
cleanup() {
  [[ -d "$fixture" && "${fixture##*/}" == pipeline-cpan-test.* ]] || return 0
  [[ "$(cd "$fixture/.." && pwd -P)" == "$temporary_root" ]] || return 1
  rm -rf -- "$fixture"
}
trap cleanup EXIT
args_log="$fixture/args"
perl() {
  printf '%s\n' "$@" > "$args_log"
  return "${test_exit:-0}"
}
PIPELINE_CPANM_BIN="$fixture/local script/cpanm"
PIPELINE_PERL_LOCAL_DIR="$fixture/local library/perl5"
PIPELINE_CPAN_MIRROR='https://cpan.metacpan.org'
count=0
assert() {
  count=$((count+1))
  if "$@"; then printf 'ok %s\n' "$count";
  else printf 'not ok %s\n' "$count"; exit 1; fi
}
read_args() {
  args=()
  while IFS= read -r arg; do args+=("$arg"); done < "$args_log"
}
PIPELINE_PLATFORM_TAG=cygwin
pipeline_cpanm --local-lib "$PIPELINE_PERL_LOCAL_DIR" --notest Mojolicious
read_args
assert test "${args[0]}" = "$PIPELINE_CPANM_BIN"
assert test "${args[1]}" = --from
assert test "${args[2]}" = "$PIPELINE_CPAN_MIRROR"
assert test "${args[3]}" = --no-lwp
assert test "${args[4]}" = --no-wget
assert test "${args[5]}" = --curl
assert test "${args[6]}" = --verbose
assert test "${args[7]}" = --local-lib
assert test "${args[8]}" = "$PIPELINE_PERL_LOCAL_DIR"
assert test "${args[9]}" = --notest
assert test "${args[10]}" = Mojolicious
for platform in linux darwin; do
  PIPELINE_PLATFORM_TAG="$platform"
  pipeline_cpanm --local-lib "$PIPELINE_PERL_LOCAL_DIR" --notest Mojolicious
  read_args
  assert test "${#args[@]}" = 7
  assert test "${args[1]}" = --from
  assert test "${args[2]}" = "$PIPELINE_CPAN_MIRROR"
  assert test "${args[4]}" = "$PIPELINE_PERL_LOCAL_DIR"
  assert test "${args[6]}" = Mojolicious
done
PIPELINE_CPAN_MIRROR='https://cpan.example.invalid/custom/'
pipeline_cpanm --reinstall --configure-timeout 900 --build-timeout 7200 PDL
read_args
assert test "${args[2]}" = "$PIPELINE_CPAN_MIRROR"
assert test "${args[3]}" = --reinstall
assert test "${args[4]}" = --configure-timeout
assert test "${args[5]}" = 900
assert test "${args[6]}" = --build-timeout
assert test "${args[7]}" = 7200
assert test "${args[8]}" = PDL
test_exit=23
set +e
pipeline_cpanm Mojolicious
rc=$?
set -e
assert test "$rc" = 23
printf '1..%s\n' "$count"
