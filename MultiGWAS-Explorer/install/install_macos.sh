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
  file_info="$(file "${bin}" 2>/dev/null || true)"
  case "${current_arch}:${file_info}" in
    arm64:*arm64*|arm64:*arm64e*|x86_64:*x86_64*) return 0 ;;
  esac
  return 1
}

select_macos_python() {
  local cand=""
  for cand in \
    /opt/homebrew/bin/python3 \
    /opt/homebrew/opt/python@3.14/bin/python3 \
    /opt/homebrew/opt/python@3.13/bin/python3 \
    /opt/homebrew/opt/python@3.12/bin/python3 \
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
ensure_homebrew

log "Installing macOS packages with Homebrew"
brew_cmd update
brew_cmd install bash curl cpanminus gd gnuplot htslib imagemagick pkg-config python wget

make_project_scripts_executable

PIPELINE_INLINE_PYTHON_BIN="$(select_macos_python || true)"
[ -n "${PIPELINE_INLINE_PYTHON_BIN}" ] || die "Could not find a Python with headers for the current CPU architecture"
export PIPELINE_INLINE_PYTHON_BIN
log "Using ${PIPELINE_INLINE_PYTHON_BIN} for Python packages and Inline::Python"

create_python_venv "${PIPELINE_INLINE_PYTHON_BIN}"
install_perl_deps
ensure_local_hts_tools
run_pipeline_check

log "macOS installation completed"
