#!/usr/bin/env bash
set -euo pipefail

_build_local_htslib_source="${BASH_SOURCE[0]:-$0}"
SCRIPT_DIR="$(cd "$(/usr/bin/dirname "${_build_local_htslib_source}")" && pwd)"
# shellcheck source=install/common.sh
. "${SCRIPT_DIR}/common.sh"

SRC_DIR="${PIPELINE_ROOT}/tools/htslib-${PIPELINE_HTSLIB_VERSION}"
TARBALL="${PIPELINE_ROOT}/tools/htslib-${PIPELINE_HTSLIB_VERSION}.tar.bz2"

if command_exists bgzip && command_exists tabix; then
  log "System bgzip/tabix already available; skipping local htslib build"
  exit 0
fi

if [ ! -d "${SRC_DIR}" ]; then
  if [ ! -f "${TARBALL}" ]; then
    log "Downloading htslib ${PIPELINE_HTSLIB_VERSION} into ${TARBALL}"
    download_url "${PIPELINE_HTSLIB_URL}" "${TARBALL}"
  fi
  log "Extracting ${TARBALL}"
  mkdir -p "${PIPELINE_ROOT}/tools"
  tar -xjf "${TARBALL}" -C "${PIPELINE_ROOT}/tools"
fi

log "Building htslib from ${SRC_DIR}"
cd "${SRC_DIR}"
make clean >/dev/null 2>&1 || true
./configure --prefix="${PIPELINE_LOCAL_DIR}"
make -j "$(num_cpus)"
make install

if [ ! -x "${PIPELINE_LOCAL_DIR}/bin/bgzip" ] || [ ! -x "${PIPELINE_LOCAL_DIR}/bin/tabix" ]; then
  die "htslib build completed but bgzip/tabix were not installed under ${PIPELINE_LOCAL_DIR}/bin"
fi

log "Installed repo-local bgzip/tabix into ${PIPELINE_LOCAL_DIR}/bin"
