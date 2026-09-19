#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${SCRIPT_DIR}/common.sh"
test_root="$(mktemp -d)"
trap 'rm -rf "$test_root"' EXIT

# Reproduce first-install state: no library directory exists yet. A module
# written after activation must be loadable by a child Perl build process.
PIPELINE_LOCAL_DIR="$test_root/local"
PIPELINE_PERL_LOCAL_DIR="$PIPELINE_LOCAL_DIR/perl5"
PIPELINE_PERL_ABI_STAMP="$PIPELINE_PERL_LOCAL_DIR/.perl-abi"
PIPELINE_VENDOR_PERL_DIR="$test_root/vendor"
unset PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MM_OPT PERL_MB_OPT
ensure_perl_abi_compatible
[ -s "$PIPELINE_PERL_ABI_STAMP" ] || die "Perl ABI stamp was not created"
activate_perl_env
printf 'package InstallProbe; sub value { 42 } 1;\n' > "$PIPELINE_PERL_LOCAL_DIR/lib/perl5/InstallProbe.pm"
perl -MInstallProbe -e 'exit(InstallProbe::value() == 42 ? 0 : 1)'
first_lib_root="$PERL_LOCAL_LIB_ROOT"
activate_perl_env
[ "$PERL_LOCAL_LIB_ROOT" = "$first_lib_root" ] || die "Repeated activation duplicates PERL_LOCAL_LIB_ROOT"

# A changed ABI stamp must rotate the entire local tree, not repair individual
# XS modules. This fixture is isolated under mktemp and contains no real CPAN
# installation.
printf 'intentionally incompatible\n' > "$PIPELINE_PERL_ABI_STAMP"
if (activate_perl_env) >/dev/null 2>&1; then
  die "Perl environment activation accepted an incompatible ABI stamp"
fi
ensure_perl_abi_compatible
perl_local_abi_matches_current || die "Perl ABI stamp was not refreshed"
[ ! -f "$PIPELINE_PERL_LOCAL_DIR/lib/perl5/InstallProbe.pm" ] \
  || die "Incompatible Perl tree was not rotated"
find "$PIPELINE_LOCAL_DIR" -maxdepth 1 -type d -name 'perl5.incompatible-*' \
  -print -quit | grep -q . || die "Incompatible Perl tree backup was not retained"
activate_perl_env

# A foreign cpanm on PATH must never be selected. Avoid a network request by
# placing a standalone script at the repository-local bootstrap destination.
mkdir -p "$test_root/foreign" "$PIPELINE_LOCAL_DIR/bin"
printf '#!/bin/sh\nexit 91\n' > "$test_root/foreign/cpanm"
chmod +x "$test_root/foreign/cpanm"
printf '# standalone test fixture\n' > "$PIPELINE_LOCAL_DIR/bin/cpanm"
export PATH="$test_root/foreign:$PATH"
ensure_cpanm
[ "$PIPELINE_CPANM_BIN" = "$PIPELINE_LOCAL_DIR/bin/cpanm" ] || die "Foreign cpanm selected"

# A JAVA_HOME path with spaces should be returned as a single executable path.
mkdir -p "$test_root/JDK with spaces/bin"
printf '#!/usr/bin/env sh\nexit 0\n' > "$test_root/JDK with spaces/bin/java"
chmod +x "$test_root/JDK with spaces/bin/java"
unset SASPY_JAVA
export JAVA_HOME="$test_root/JDK with spaces"
resolved_java="$(resolve_unix_java_for_saspy)"
[ "$resolved_java" = "$JAVA_HOME/bin/java" ] || die "JAVA_HOME was not respected"
"$resolved_java" -version
log "First-install Perl environment and Java path checks passed"
