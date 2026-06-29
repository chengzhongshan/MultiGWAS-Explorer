#!/usr/bin/env bash
set -euo pipefail

_install_check_source="${BASH_SOURCE[0]:-$0}"
SCRIPT_DIR="$(cd "$(/usr/bin/dirname "${_install_check_source}")" && pwd)"
# shellcheck source=install/common.sh
. "${SCRIPT_DIR}/common.sh"

need_cmd() {
  local cmd="$1"
  command_exists "$cmd" || die "Required command not found: ${cmd}"
}

activate_perl_env
activate_python_env
prepend_path "${PIPELINE_LOCAL_DIR}/bin"

need_cmd bash
need_cmd perl
need_cmd gnuplot
log "perl archname: $(perl -MConfig -e 'print $Config{archname}')"
log "gnuplot on PATH: $(command -v gnuplot)"
log "gnuplot version: $(gnuplot --version | head -n 1)"

if ! command_exists bgzip && [ ! -x "${PIPELINE_LOCAL_DIR}/bin/bgzip.exe" ]; then
  die "bgzip not found on PATH or under local/bin"
fi
if ! command_exists tabix && [ ! -x "${PIPELINE_LOCAL_DIR}/bin/tabix.exe" ]; then
  die "tabix not found on PATH or under local/bin"
fi
if ! command_exists magick && ! command_exists convert; then
  die "ImageMagick executable not found as magick or convert"
fi

[ -n "${PIPELINE_PYTHON_BIN}" ] || die "Repo-local Python environment not found; run an install script first"

"${PIPELINE_PYTHON_BIN}" - <<'PY'
import PIL
import saspy
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
PY

log "GD version: $(perl -MGD -e 'print $GD::VERSION')"
log "PDL version: $(perl -MPDL -e 'print $PDL::VERSION')"
perl -e "require JSON; require JSON::MaybeXS; require File::Which; require GD; require Mojolicious::Lite; require MCP::Server; require PDL; 1;" >/dev/null
perl -I DiffGWASDeps -MSAS_ODA_Runner -e "print qq{SAS_ODA_Runner ok\n};"
perl -I DiffGWASDeps -c auto_prepare_and_run_diff_gwas.pl >/dev/null
perl -c auto_prepare_and_run_diff_gwas_with_gunplot.pl >/dev/null
perl -I DiffGWASDeps -c server.pl >/dev/null
perl -I DiffGWASDeps -c run_sas_codes_or_script_in_ODA.pl >/dev/null
perl -c DiffGWASDeps/gunplot/pdl_gunplot_manhattan.pl >/dev/null
perl -c DiffGWASDeps/gunplot/pdl_gunplot_forest.pl >/dev/null
perl -c DiffGWASDeps/gunplot/pdl_gunplot_local_locus.pl >/dev/null
bash -n DiffGWASDeps/run_sas_oda_manhattan4diffgwas_download_png.sh
bash -n DiffGWASDeps/run_sas_oda_local_top_hits_manhattan_download_png.sh
bash -n DiffGWASDeps/run_sas_oda_local_top_hits_with_gtf_download_html.sh

log "Pipeline dependency smoke test completed successfully"
