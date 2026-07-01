#!/usr/bin/env bash
set -euo pipefail

_install_ubuntu_source="${BASH_SOURCE[0]:-$0}"
SCRIPT_DIR="$(cd "$(/usr/bin/dirname "${_install_ubuntu_source}")" && pwd)"
# shellcheck source=install/common.sh
. "${SCRIPT_DIR}/common.sh"

APT_GET="${APT_GET:-apt-get}"
SKIP_APT="${PIPELINE_SKIP_APT:-0}"
if [[ "${SKIP_APT}" =~ ^(1|true|yes|y|on)$ ]]; then
  log "Skipping Ubuntu package installation because PIPELINE_SKIP_APT=${SKIP_APT}"
else
  if [ "$(id -u)" -eq 0 ]; then
    SUDO=()
  else
    SUDO=(sudo)
    if ! sudo -n true >/dev/null 2>&1; then
      if [ -t 0 ]; then
        log "Ubuntu package installation requires sudo; requesting credentials"
        sudo -v || die "sudo authentication failed. If packages are already installed, rerun with PIPELINE_SKIP_APT=1."
      else
        die "Ubuntu package installation requires sudo, but this shell cannot accept a password prompt. Run 'sudo bash install/install_ubuntu.sh' from a normal terminal, or run 'sudo -v' there first and rerun this installer. If apt packages are already installed, rerun with PIPELINE_SKIP_APT=1."
      fi
    fi
  fi

  log "Installing Ubuntu packages"
  DEBIAN_FRONTEND=noninteractive "${SUDO[@]}" "${APT_GET}" update
  DEBIAN_FRONTEND=noninteractive "${SUDO[@]}" "${APT_GET}" install -y \
    bash \
    build-essential \
    curl \
    default-jre-headless \
    git \
    gnuplot-nox \
    imagemagick \
    libbz2-dev \
    libcurl4-openssl-dev \
    libgd-dev \
    liblzma-dev \
    make \
    perl \
    pkg-config \
    python3 \
    python3-dev \
    python3-pip \
    python3-venv \
    tabix \
    unzip \
    wget \
    zip \
    zlib1g-dev
fi

UBUNTU_PYTHON_BIN="${PIPELINE_UBUNTU_PYTHON_BIN:-/usr/bin/python3}"
[ -x "${UBUNTU_PYTHON_BIN}" ] || die "Ubuntu Python not found at ${UBUNTU_PYTHON_BIN}; install python3 or set PIPELINE_UBUNTU_PYTHON_BIN"
if ! "${UBUNTU_PYTHON_BIN}" - <<'PY'
import sys
raise SystemExit(0 if sys.version_info >= (3, 8) else 1)
PY
then
  die "Python at ${UBUNTU_PYTHON_BIN} is too old; Pillow>=10 and current saspy require Python >=3.8. Set PIPELINE_UBUNTU_PYTHON_BIN to a newer Python."
fi
log "Using ${UBUNTU_PYTHON_BIN} for repo-local Python environment"
create_python_venv "${UBUNTU_PYTHON_BIN}"
install_perl_deps
ensure_local_hts_tools
run_pipeline_check

log "Ubuntu installation completed"
