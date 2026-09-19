#!/usr/bin/env bash
set -euo pipefail

_install_check_source="${BASH_SOURCE[0]:-$0}"
SCRIPT_DIR="$(cd "$(/usr/bin/dirname "${_install_check_source}")" && pwd)"
# shellcheck source=install/common.sh
. "${SCRIPT_DIR}/common.sh"
cd "${PIPELINE_ROOT}"

need_cmd() {
  local cmd="$1"
  command_exists "$cmd" || die "Required command not found: ${cmd}"
}

activate_perl_env
activate_python_env
prepend_path "${PIPELINE_LOCAL_DIR}/bin"
prepend_path "${PIPELINE_ROOT}"
prepend_path "${PIPELINE_ROOT}/DiffGWASDeps"

need_cmd bash
need_cmd perl
need_cmd gnuplot
log "perl archname: $(perl -MConfig -e 'print $Config{archname}')"
log "gnuplot on PATH: $(command -v gnuplot)"
log "gnuplot version: $(gnuplot --version | head -n 1)"

if ! command_exists bgzip; then
  die "bgzip not found in DiffGWASDeps, the pipeline root, local/bin, or PATH"
fi
if ! command_exists tabix; then
  die "tabix not found in DiffGWASDeps, the pipeline root, local/bin, or PATH"
fi
if ! command_exists magick && ! command_exists convert; then
  die "ImageMagick executable not found as magick or convert"
fi

[ -n "${PIPELINE_PYTHON_BIN}" ] || die "Repo-local Python environment not found; run an install script first"

"${PIPELINE_PYTHON_BIN}" - <<'PY'
import PIL
import saspy
import os
import shutil
import subprocess
import sys
print("python imports ok")
cfg_names = None
try:
    import saspy.sascfg_personal as personal_cfg
    print(f"saspy personal config: {getattr(personal_cfg, '__file__', 'unknown')}")
    cfg_names = getattr(personal_cfg, 'SAS_config_names', None)
    print(f"saspy config names: {cfg_names}")
except Exception:
    personal_cfg = None
if not cfg_names or 'oda' not in cfg_names:
    raise SystemExit("saspy ODA config was not provisioned in the repo-local install")
oda_cfg = getattr(personal_cfg, 'oda', None)
if not isinstance(oda_cfg, dict):
    raise SystemExit("saspy ODA profile is missing its oda configuration dictionary")
if not oda_cfg.get('java'):
    raise SystemExit("saspy ODA profile has no Java executable; rerun the platform installer")
if 'saspyiom.jar' not in (oda_cfg.get('classpath') or ''):
    raise SystemExit("saspy ODA profile has no usable IOM classpath; rerun the platform installer")
def local_path(path):
    if sys.platform == 'cygwin' and len(path) > 2 and path[1] == ':':
        return subprocess.check_output(['cygpath', '-u', path], text=True).strip()
    return path

java_path = local_path(oda_cfg['java'])
if not os.path.isfile(java_path) and not shutil.which(java_path):
    raise SystemExit("saspy ODA Java path does not exist; rerun the platform installer")
separator = ';' if sys.platform in ('cygwin', 'win32') else ':'
for jar in oda_cfg['classpath'].split(separator):
    if jar and not os.path.isfile(local_path(jar)):
        raise SystemExit("saspy ODA classpath contains a missing JAR; rerun the platform installer: " + jar)
PY

if command_exists cygpath; then
  java_bin="$(resolve_windows_java_for_saspy || true)"
  [ -n "$java_bin" ] && java_bin="$(cygpath -u "$java_bin")"
else
  java_bin="$(resolve_unix_java_for_saspy || true)"
fi
[ -n "$java_bin" ] || die "Java runtime not found. Install a JDK and set JAVA_HOME, SASPY_JAVA (Unix), or SASPY_JAVA_WIN (Windows)."
"$java_bin" -version || die "Java could not start: ${java_bin}"

gd_version="$(perl -MGD -e 'print $GD::VERSION')" || die "GD cannot load. Use the Perl/Cygwin installation that built local modules, or reinstall dependencies with the current Perl."
pdl_version="$(perl -MPDL -e 'my $x = sequence(3); die "PDL arithmetic failed\n" unless $x->at(2) == 2; print $PDL::VERSION')" \
  || die "PDL cannot load and run arithmetic. Reinstall dependencies with the current Perl/Cygwin installation."
[ -n "${pdl_version}" ] || die "PDL returned an empty version; reinstall dependencies with the current Perl/Cygwin installation."
log "GD version: $gd_version"
log "PDL version: $pdl_version"
perl -e "require JSON::PP; require JSON::MaybeXS; require File::Which; require DBI; require DBD::SQLite; require GD; require Mojolicious::Lite; require MCP::Server; require PDL; require Text::CSV; 1;" >/dev/null
if command_exists uname && uname -s | grep -qi '^CYGWIN'; then
  perl -e "require JSON; require Inline::Python; require Compress::Raw::Zlib; require IO::Uncompress::Gunzip; require Compress::Raw::Bzip2; require IO::Uncompress::Bunzip2; 1;" >/dev/null
fi
perl -I DiffGWASDeps -MSAS_ODA_Runner -e "print qq{SAS_ODA_Runner ok\n};"
perl -MIO::Socket::SSL -MNet::SSLeay -MHTTP::Tiny -e 'my ($ok, $why) = HTTP::Tiny->can_ssl; die "Perl HTTPS unavailable: $why\n" unless $ok; print "Perl HTTPS support ok\n";'
perl DiffGWASDeps/test_sas_oda_debug_macro_guard.pl >/dev/null
perl DiffGWASDeps/test_sas_oda_connection_lifecycle.pl >/dev/null
perl -I DiffGWASDeps -c auto_prepare_and_run_diff_gwas.pl >/dev/null
perl -c auto_prepare_and_run_diff_gwas_with_gunplot.pl >/dev/null
perl -I DiffGWASDeps -c server.pl >/dev/null
perl -I DiffGWASDeps -c run_sas_codes_or_script_in_ODA.pl >/dev/null
perl -c DiffGWASDeps/gnuplot/pdl_gunplot_manhattan.pl >/dev/null
perl -c DiffGWASDeps/gnuplot/pdl_gunplot_forest.pl >/dev/null
perl -c DiffGWASDeps/gnuplot/pdl_gunplot_local_locus.pl >/dev/null
perl DiffGWASDeps/test_gnuplot_directory_layout.pl >/dev/null
perl DiffGWASDeps/test_ld_heatmap_rendering.pl >/dev/null
perl DiffGWASDeps/test_ld_heatmap_contract.pl >/dev/null
perl DiffGWASDeps/test_ld_cache_queries.pl >/dev/null
perl install/test_sort_long_gwas.pl >/dev/null
perl install/test_precomputed_paths.pl >/dev/null
perl install/test_forest_text.pl
perl install/test_requested_hit_genes.pl
"${BASH:-bash}" -n DiffGWASDeps/run_sas_oda_manhattan4diffgwas_download_png.sh
"${BASH:-bash}" -n DiffGWASDeps/run_sas_oda_local_top_hits_manhattan_download_png.sh
"${BASH:-bash}" -n DiffGWASDeps/run_sas_oda_local_top_hits_with_gtf_download_html.sh
"${BASH:-bash}" -n DiffGWASDeps/run_sas_oda_top_hits_forest_plot_download_html.sh

log "Pipeline dependency smoke test completed successfully"
log "Local dependencies and synthetic plotting passed; SAS ODA login was not tested"
