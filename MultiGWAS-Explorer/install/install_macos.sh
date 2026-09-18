#!/usr/bin/env bash
set -euo pipefail

_install_macos_source="${BASH_SOURCE[0]:-$0}"
SCRIPT_DIR="$(cd "$(/usr/bin/dirname "${_install_macos_source}")" && pwd)"
# shellcheck source=install/common.sh
. "${SCRIPT_DIR}/common.sh"

ensure_xcode_clt() {
  if xcode-select -p >/dev/null 2>&1; then
    return 0
  fi
  warn "Xcode Command Line Tools are required for Perl module compilation"
  xcode-select --install || true
  die "Install the Xcode Command Line Tools and rerun install/install_macos.sh"
}

ensure_homebrew() {
  if ! command_exists brew; then
    log "Installing Homebrew"
    NONINTERACTIVE=1 /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
  fi
  if [ -x /opt/homebrew/bin/brew ]; then
    eval "$(/opt/homebrew/bin/brew shellenv)"
  elif [ -x /usr/local/bin/brew ]; then
    eval "$(/usr/local/bin/brew shellenv)"
  fi
}

ensure_macports() {
  local version="2.12.5"
  local package="MacPorts-${version}-15-Sequoia.pkg"
  local expected_sha256="10a048e235ba252eb31ca5030dbf560409ea73cacb3267a00e3983205e9b0e36"
  local package_path="${PIPELINE_ROOT}/tools/${package}"
  local actual_sha256=""

  if [ -x /opt/local/bin/port ]; then
    prepend_path /opt/local/sbin
    prepend_path /opt/local/bin
    return 0
  fi
  [ "$(sw_vers -productVersion | cut -d. -f1)" = "15" ] \
    || die "The pinned MacPorts fallback supports Intel macOS 15; set PIPELINE_MACOS_PACKAGE_MANAGER=homebrew to override"
  log "Installing verified MacPorts ${version} for Intel macOS 15"
  download_url \
    "https://distfiles.macports.org/MacPorts/${package}" \
    "${package_path}"
  if command_exists sha256sum; then
    actual_sha256="$(sha256sum "${package_path}" | awk '{print $1}')"
  else
    actual_sha256="$(shasum -a 256 "${package_path}" | awk '{print $1}')"
  fi
  [ "${actual_sha256}" = "${expected_sha256}" ] \
    || die "MacPorts installer checksum mismatch: ${package_path}"
  sudo /usr/sbin/installer -pkg "${package_path}" -target /
  prepend_path /opt/local/sbin
  prepend_path /opt/local/bin
}

persist_github_actions_path() {
  local dir=""
  [ -n "${GITHUB_PATH:-}" ] || return 0
  for dir in "$@"; do
    [ -d "${dir}" ] || continue
    printf '%s\n' "${dir}" >> "${GITHUB_PATH}"
  done
}

macports_cmd() {
  sudo /opt/local/bin/port -N "$@"
}

brew_cmd() {
  if [ -x /opt/homebrew/bin/brew ]; then
    /usr/bin/arch -arm64 /opt/homebrew/bin/brew "$@"
  elif [ "$(uname -m)" = "arm64" ] && [ -x /usr/local/bin/brew ]; then
    warn "Using Intel Homebrew under /usr/local on Apple Silicon; prefer installing ARM Homebrew under /opt/homebrew"
    /usr/bin/arch -x86_64 /usr/local/bin/brew "$@"
  else
    brew "$@"
  fi
}

binary_supports_current_arch() {
  local bin="$1"
  local current_arch=""
  local file_info=""
  [ -n "${bin}" ] || return 1
  [ -x "${bin}" ] || return 1
  current_arch="$(uname -m)"
  file_info="$(file -L "${bin}" 2>/dev/null || true)"
  case "${current_arch}:${file_info}" in
    arm64:*arm64*|arm64:*arm64e*|x86_64:*x86_64*) return 0 ;;
  esac
  return 1
}

select_macos_python() {
  local cand=""
  for cand in \
    /opt/local/bin/python3.12 \
    /opt/local/bin/python3 \
    /opt/homebrew/bin/python3 \
    /opt/homebrew/opt/python@3.14/bin/python3 \
    /opt/homebrew/opt/python@3.13/bin/python3 \
    /opt/homebrew/opt/python@3.12/bin/python3 \
    /usr/local/bin/python3 \
    /usr/bin/python3 \
    python3; do
    if command_exists "${cand}" || [ -x "${cand}" ]; then
      cand="$(command -v "${cand}" 2>/dev/null || printf '%s\n' "${cand}")"
      binary_supports_current_arch "${cand}" || continue
      if "${cand}" - <<'PY' >/dev/null 2>&1
import sysconfig
import pathlib

inc = sysconfig.get_config_var("INCLUDEPY")
raise SystemExit(0 if inc and pathlib.Path(inc, "Python.h").exists() else 1)
PY
      then
        printf '%s\n' "${cand}"
        return 0
      fi
    fi
  done
  return 1
}

ensure_xcode_clt
macos_package_manager="${PIPELINE_MACOS_PACKAGE_MANAGER:-}"
if [ -z "${macos_package_manager}" ]; then
  if [ "$(uname -m)" = "x86_64" ] && [ "$(sw_vers -productVersion | cut -d. -f1)" -ge 15 ]; then
    macos_package_manager="macports"
  else
    macos_package_manager="homebrew"
  fi
fi

case "${macos_package_manager}" in
  macports)
    ensure_macports
    # GitHub Actions starts every workflow step in a fresh non-login shell.
    # Publish the MacPorts paths so the following rendering step can find the
    # dependencies that this step just installed.
    persist_github_actions_path /opt/local/bin /opt/local/sbin
    log "Installing Intel macOS packages with MacPorts"
    macports_cmd selfupdate
    macports_cmd install \
      bash curl gd2 htslib ImageMagick openjdk21 openssl pkgconfig \
      python312 py312-pip wget
    macports_cmd install gnuplot \
      +pangocairo -aquaterm -luaterm -qt -qt5 -wxwidgets -x11
    macports_cmd select --set python3 python312
    macports_cmd select --set pip3 pip312
    export OPENSSL_PREFIX="/opt/local"
    export PKG_CONFIG_PATH="/opt/local/lib/pkgconfig:/opt/local/share/pkgconfig${PKG_CONFIG_PATH:+:${PKG_CONFIG_PATH}}"
    export CPPFLAGS="-I/opt/local/include${CPPFLAGS:+ ${CPPFLAGS}}"
    export LDFLAGS="-L/opt/local/lib${LDFLAGS:+ ${LDFLAGS}}"
    if [ -d /Library/Java/JavaVirtualMachines/jdk-21-macports.jdk/Contents/Home ]; then
      export JAVA_HOME=/Library/Java/JavaVirtualMachines/jdk-21-macports.jdk/Contents/Home
      prepend_path "${JAVA_HOME}/bin"
    fi
    ;;
  homebrew)
    ensure_homebrew
    log "Installing macOS packages with Homebrew"
    brew_cmd update
    brew_cmd install bash curl gd htslib imagemagick openjdk openssl@3 pkg-config python wget
    export OPENSSL_PREFIX="$(brew_cmd --prefix openssl@3)"
    ;;
  *)
    die "Unsupported PIPELINE_MACOS_PACKAGE_MANAGER '${macos_package_manager}'; use homebrew or macports"
    ;;
esac

prepend_path "${PIPELINE_LOCAL_DIR}/bin"
if ! command_exists gnuplot || ! gnuplot -e 'set terminal pngcairo' >/dev/null 2>&1; then
  if [ "${macos_package_manager}" = "macports" ]; then
    die "MacPorts gnuplot is installed but does not provide the pngcairo terminal"
  elif [ "${PIPELINE_MACOS_GNUPLOT:-headless}" = brew ]; then
    brew_cmd install gnuplot
  else
    brew_cmd install cairo pango
    bash "${SCRIPT_DIR}/build_local_gnuplot.sh"
  fi
fi

make_project_scripts_executable

PIPELINE_MACOS_PYTHON_BIN="$(select_macos_python || true)"
[ -n "${PIPELINE_MACOS_PYTHON_BIN}" ] || die "Could not find a Python with headers for the current CPU architecture"
log "Using ${PIPELINE_MACOS_PYTHON_BIN} for Python packages"

create_python_venv "${PIPELINE_MACOS_PYTHON_BIN}"
install_perl_deps
ensure_local_hts_tools
run_pipeline_check

log "macOS installation completed"
