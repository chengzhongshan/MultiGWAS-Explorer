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

ensure_xcode_clt
ensure_homebrew

log "Installing macOS packages with Homebrew"
brew update
brew install bash curl gd gnuplot htslib imagemagick pkg-config python wget

create_python_venv "$(command -v python3)"
install_perl_deps
ensure_local_hts_tools
run_pipeline_check

log "macOS installation completed"
