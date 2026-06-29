#!/usr/bin/env bash
set -euo pipefail

_install_cygwin_source="${BASH_SOURCE[0]:-$0}"
SCRIPT_DIR="$(cd "$(/usr/bin/dirname "${_install_cygwin_source}")" && pwd)"
# shellcheck source=install/common.sh
. "${SCRIPT_DIR}/common.sh"

resolve_cygwin_root_windows() {
  if command_exists mount; then
    local mounted_root
    mounted_root="$(mount | awk '$3 == "/" { print $1; exit }' 2>/dev/null || true)"
    if [[ -n "${mounted_root}" ]]; then
      printf '%s\n' "${mounted_root}"
      return 0
    fi
  fi
  if command_exists cygpath; then
    cygpath -w /
    return 0
  fi
  return 1
}

command_exists cygpath || die "install/install_cygwin.sh must be run inside a Cygwin shell; on Windows use install/install_windows_portable_cygwin.ps1 first"

CYGWIN_SETUP_URL="${CYGWIN_SETUP_URL:-https://cygwin.com/setup-x86_64.exe}"
CYGWIN_SETUP_MIRROR="${CYGWIN_SETUP_MIRROR:-https://mirrors.kernel.org/sourceware/cygwin/}"
CYGWIN_SETUP_EXE="${CYGWIN_SETUP_EXE:-${PIPELINE_INSTALL_DIR}/cache/setup-x86_64.exe}"
CYGWIN_ROOT_WINDOWS="${CYGWIN_ROOT_WINDOWS:-$(resolve_cygwin_root_windows)}"
CYGWIN_PACKAGES="${CYGWIN_PACKAGES:-bash,curl,cygwin,gcc-core,gcc-g++,gnuplot-base,ImageMagick,libgd-devel,make,perl,perl-File-Which,perl-GD,perl-JSON,perl-JSON-MaybeXS,perl-Mojolicious,pkg-config,python3,python312,python312-devel,python312-imaging,python312-pip,python312-setuptools,python312-wheel,unzip,wget,which,zip}"
CYGWIN_SKIP_PACKAGE_UPDATE="${CYGWIN_SKIP_PACKAGE_UPDATE:-0}"
CYGWIN_ONLY_PACKAGE_UPDATE="${CYGWIN_ONLY_PACKAGE_UPDATE:-0}"

if [[ "${CYGWIN_SKIP_PACKAGE_UPDATE}" != "1" ]]; then
if [[ -s "${CYGWIN_SETUP_EXE}" ]]; then
  log "Reusing existing Cygwin setup helper at ${CYGWIN_SETUP_EXE}"
else
  log "Downloading Cygwin setup helper into ${CYGWIN_SETUP_EXE}"
  download_url "${CYGWIN_SETUP_URL}" "${CYGWIN_SETUP_EXE}"
fi
chmod +x "${CYGWIN_SETUP_EXE}"

  log "Installing or updating required Cygwin packages"
  "${CYGWIN_SETUP_EXE}" \
    -q \
    -B \
    -g \
    -n \
    -N \
    -d \
    --no-write-registry \
    -R "${CYGWIN_ROOT_WINDOWS}" \
    -s "${CYGWIN_SETUP_MIRROR}" \
    -P "${CYGWIN_PACKAGES}"
else
  log "Skipping Cygwin package refresh because CYGWIN_SKIP_PACKAGE_UPDATE=1"
fi

if [[ "${CYGWIN_ONLY_PACKAGE_UPDATE}" == "1" ]]; then
  log "Stopping after package refresh because CYGWIN_ONLY_PACKAGE_UPDATE=1"
  exit 0
fi

create_python_venv
install_perl_deps
ensure_local_hts_tools
run_pipeline_check

log "Cygwin installation completed"
