#!/usr/bin/env bash
set -euo pipefail

_install_ubuntu_source="${BASH_SOURCE[0]:-$0}"
SCRIPT_DIR="$(cd "$(/usr/bin/dirname "${_install_ubuntu_source}")" && pwd)"
# shellcheck source=install/common.sh
. "${SCRIPT_DIR}/common.sh"

APT_GET="${APT_GET:-apt-get}"
if [ "$(id -u)" -eq 0 ]; then
  SUDO=()
else
  SUDO=(sudo)
fi

log "Installing Ubuntu packages"
DEBIAN_FRONTEND=noninteractive "${SUDO[@]}" "${APT_GET}" update
DEBIAN_FRONTEND=noninteractive "${SUDO[@]}" "${APT_GET}" install -y \
  bash \
  build-essential \
  curl \
  default-jre-headless \
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

create_python_venv "$(command -v python3)"
install_perl_deps
ensure_local_hts_tools
run_pipeline_check

log "Ubuntu installation completed"
