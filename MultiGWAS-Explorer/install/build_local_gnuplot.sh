#!/usr/bin/env bash
# Build the PNG renderer without pulling in Qt or other GUI toolkits.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${SCRIPT_DIR}/common.sh"
version=6.0.4
archive="${PIPELINE_ROOT}/tools/gnuplot-${version}.tar.gz"
source_dir="${PIPELINE_ROOT}/tools/gnuplot-${version}"
expected_sha256=458d94769625e73d5f6232500f49cbadcb2b183380d43d2266a0f9701aeb9c5b
download_url "https://downloads.sourceforge.net/project/gnuplot/gnuplot/${version}/gnuplot-${version}.tar.gz" "$archive"
if command_exists sha256sum; then
  actual_sha256="$(sha256sum "$archive" | awk '{print $1}')"
else
  actual_sha256="$(shasum -a 256 "$archive" | awk '{print $1}')"
fi
[ "$actual_sha256" = "$expected_sha256" ] || die "gnuplot archive checksum mismatch: ${archive}"
tar -xzf "$archive" -C "${PIPELINE_ROOT}/tools"
(
  cd "$source_dir"
  ./configure --prefix="${PIPELINE_LOCAL_DIR}" \
    --without-qt --disable-wxwidgets --without-x --without-lua \
    --without-readline --without-latex
  make -j"$(num_cpus)"
  make install
)
"${PIPELINE_LOCAL_DIR}/bin/gnuplot" -e 'set terminal pngcairo' \
  || die "Built gnuplot lacks pngcairo; install cairo and pango development packages"
log "Installed headless gnuplot ${version} under ${PIPELINE_LOCAL_DIR}"
