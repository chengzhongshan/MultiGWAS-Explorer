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
PIPELINE_VENDOR_PERL_DIR="$test_root/vendor"
unset PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MM_OPT PERL_MB_OPT
activate_perl_env
printf 'package InstallProbe; sub value { 42 } 1;\n' > "$PIPELINE_PERL_LOCAL_DIR/lib/perl5/InstallProbe.pm"
perl -MInstallProbe -e 'exit(InstallProbe::value() == 42 ? 0 : 1)'
first_lib_root="$PERL_LOCAL_LIB_ROOT"
activate_perl_env
[ "$PERL_LOCAL_LIB_ROOT" = "$first_lib_root" ] || die "Repeated activation duplicates PERL_LOCAL_LIB_ROOT"

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
