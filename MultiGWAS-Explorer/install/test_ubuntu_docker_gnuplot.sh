#!/usr/bin/env bash
set -euo pipefail

_smoke_source="${BASH_SOURCE[0]:-$0}"
SCRIPT_DIR="$(cd "$(/usr/bin/dirname "${_smoke_source}")" && pwd)"
# shellcheck source=install/common.sh
. "${SCRIPT_DIR}/common.sh"

IMAGE_NAME="${IMAGE_NAME:-multigwas-explorer-pipeline:latest}"
SPEC_PATH="configs/spec_pgc_scz_sex_common_automation.json"
TARGET_SNP="rs185665940"
PLOTS="local_manhattan,local_gtf"
DISPLAY_GWAS=""
OUTPUT_STEM="PGC_SCZ_GUNPLOT_DOCKER_TEST"
HOST_DATA_ROOT=""
INCLUDE_MANHATTAN=0
SKIP_BUILD=0
FORCE=0

usage() {
  cat <<'EOF'
Usage:
  bash install/test_ubuntu_docker_gnuplot.sh [options]

Options:
  --spec PATH              Spec JSON to run inside the Ubuntu Docker image.
                           Default: configs/spec_pgc_scz_sex_common_automation.json
  --target-snp RSID        Inquiry SNP to reproduce. Default: rs185665940
  --plots LIST             Comma list from manhattan,local_manhattan,local_gtf.
                           Default: local_manhattan,local_gtf
  --include-manhattan      Prepend the slower genome-wide Manhattan stage.
  --display-gwas LIST      Optional displayed GWAS override.
  --output-stem STEM       Output basename stem for the Docker smoke test.
                           Default: PGC_SCZ_GUNPLOT_DOCKER_TEST
  --data-root-host PATH    Host directory to mount onto the spec's /mnt/<drive> root.
                           If omitted, the script auto-detects /cygdrive/<drive> or /mnt/<drive>.
  --image NAME             Docker image name. Default: multigwas-explorer-pipeline:latest
  --skip-build             Reuse the existing image instead of rebuilding it.
  --force                  Forward --force to the gunplot wrapper.
  --help                   Show this help text.

Notes:
  - This smoke test keeps PIPELINE_WORKDIR=/opt/MultiGWAS-Explorer so the container
    uses the Linux-installed repo inside the image.
  - The default quick path validates one top differential SNP with
    local_manhattan + local_gtf. Add --include-manhattan for the slower
    genome-wide renderer.
EOF
}

while [ $# -gt 0 ]; do
  case "$1" in
    --spec)
      [ $# -ge 2 ] || die "--spec requires a value"
      SPEC_PATH="$2"
      shift 2
      ;;
    --target-snp)
      [ $# -ge 2 ] || die "--target-snp requires a value"
      TARGET_SNP="$2"
      shift 2
      ;;
    --plots)
      [ $# -ge 2 ] || die "--plots requires a value"
      PLOTS="$2"
      shift 2
      ;;
    --include-manhattan)
      INCLUDE_MANHATTAN=1
      shift
      ;;
    --display-gwas)
      [ $# -ge 2 ] || die "--display-gwas requires a value"
      DISPLAY_GWAS="$2"
      shift 2
      ;;
    --output-stem)
      [ $# -ge 2 ] || die "--output-stem requires a value"
      OUTPUT_STEM="$2"
      shift 2
      ;;
    --data-root-host)
      [ $# -ge 2 ] || die "--data-root-host requires a value"
      HOST_DATA_ROOT="$2"
      shift 2
      ;;
    --image)
      [ $# -ge 2 ] || die "--image requires a value"
      IMAGE_NAME="$2"
      shift 2
      ;;
    --skip-build)
      SKIP_BUILD=1
      shift
      ;;
    --force)
      FORCE=1
      shift
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      die "Unknown argument: $1"
      ;;
  esac
done

command_exists docker || die "docker is required for this smoke test"
docker version >/dev/null 2>&1 || die "docker is not ready; start Docker Desktop or the Docker daemon first"

abs_path_portable() {
  perl -MCwd=abs_path -e 'my $p = shift; defined $p or exit 1; my $abs = abs_path($p); exit 1 unless defined $abs; print $abs;' "$1"
}

docker_mount_path() {
  local path="$1"
  if command_exists cygpath; then
    cygpath -m "$path"
  else
    abs_path_portable "$path"
  fi
}

json_string() {
  local file="$1"
  local key="$2"
  perl -MJSON::PP -e '
    local $/;
    my ($file, $key) = @ARGV;
    open my $fh, "<", $file or die "Cannot read $file: $!";
    my $data = decode_json(<$fh>);
    my $value = $data->{$key};
    print defined $value ? $value : q{};
  ' "$file" "$key"
}

container_mount_root_from_path() {
  local path="$1"
  case "$path" in
    /mnt/[A-Za-z]/*|/mnt/[A-Za-z])
      printf '/mnt/%s\n' "$(printf '%s' "$path" | sed -E 's#^/mnt/([A-Za-z]).*#\L\1#')"
      ;;
    [A-Za-z]:/*|[A-Za-z]:\\*)
      printf '/mnt/%s\n' "$(printf '%s' "$path" | sed -E 's#^([A-Za-z]):.*#\L\1#')"
      ;;
    *)
      return 1
      ;;
  esac
}

default_host_root_for_container_mount() {
  local container_root="$1"
  local drive=""
  case "$container_root" in
    /mnt/[A-Za-z])
      drive="${container_root#/mnt/}"
      if [ -d "/cygdrive/${drive}" ]; then
        printf '/cygdrive/%s\n' "$drive"
        return 0
      fi
      if [ -d "/mnt/${drive}" ]; then
        printf '/mnt/%s\n' "$drive"
        return 0
      fi
      ;;
  esac
  return 1
}

map_container_path_to_host() {
  local container_path="$1"
  local container_root="$2"
  local host_root="$3"
  case "$container_path" in
    "${container_root}"/*)
      printf '%s/%s\n' "${host_root%/}" "${container_path#${container_root}/}"
      ;;
    "${container_root}")
      printf '%s\n' "${host_root%/}"
      ;;
    *)
      return 1
      ;;
  esac
}

csv_contains() {
  local csv="$1"
  local needle="$2"
  local item=""
  OLD_IFS="${IFS}"
  IFS=','
  for item in $csv; do
    item="$(printf '%s' "$item" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')"
    if [ "$item" = "$needle" ]; then
      IFS="${OLD_IFS}"
      return 0
    fi
  done
  IFS="${OLD_IFS}"
  return 1
}

docker_image_exists() {
  local image="$1"
  local repo="$image"
  local tag="latest"
  local line=""
  if [[ "$image" == *:* ]]; then
    repo="${image%:*}"
    tag="${image##*:}"
  fi
  while IFS= read -r line; do
    [ -n "$line" ] || continue
    if [ "$line" = "${repo}:${tag}" ]; then
      return 0
    fi
  done <<EOF
$(docker images --format '{{.Repository}}:{{.Tag}}')
EOF
  return 1
}

if [ "${INCLUDE_MANHATTAN}" -eq 1 ] && ! csv_contains "$PLOTS" "manhattan"; then
  PLOTS="manhattan,${PLOTS}"
fi

SPEC_ABS="$(abs_path_portable "${SPEC_PATH}")" || die "Spec file not found: ${SPEC_PATH}"
SPEC_DIR="$(cd "$(/usr/bin/dirname "${SPEC_ABS}")" && pwd)"
SPEC_BASENAME="$(basename "${SPEC_ABS}")"
SPEC_INPUT_DIR="$(json_string "${SPEC_ABS}" "input_dir")"
SPEC_OUTPUT_DIR="$(json_string "${SPEC_ABS}" "output_dir")"
[ -n "${SPEC_INPUT_DIR}" ] || die "Spec is missing input_dir: ${SPEC_ABS}"
[ -n "${SPEC_OUTPUT_DIR}" ] || die "Spec is missing output_dir: ${SPEC_ABS}"

CONTAINER_INPUT_ROOT="$(container_mount_root_from_path "${SPEC_INPUT_DIR}")" || die "Could not derive a /mnt/<drive> input root from ${SPEC_INPUT_DIR}"
CONTAINER_OUTPUT_ROOT="$(container_mount_root_from_path "${SPEC_OUTPUT_DIR}")" || die "Could not derive a /mnt/<drive> output root from ${SPEC_OUTPUT_DIR}"
[ "${CONTAINER_INPUT_ROOT}" = "${CONTAINER_OUTPUT_ROOT}" ] || die "This smoke test expects input_dir and output_dir to share one mounted root, but saw ${CONTAINER_INPUT_ROOT} vs ${CONTAINER_OUTPUT_ROOT}"
CONTAINER_DATA_ROOT="${CONTAINER_INPUT_ROOT}"

if [ -n "${HOST_DATA_ROOT}" ]; then
  HOST_DATA_ROOT_ABS="$(abs_path_portable "${HOST_DATA_ROOT}")" || die "Host data root not found: ${HOST_DATA_ROOT}"
else
  HOST_DATA_ROOT_ABS="$(default_host_root_for_container_mount "${CONTAINER_DATA_ROOT}" || true)"
  [ -n "${HOST_DATA_ROOT_ABS}" ] || die "Could not auto-detect a host mount for ${CONTAINER_DATA_ROOT}; rerun with --data-root-host"
fi

[ -d "${HOST_DATA_ROOT_ABS}" ] || die "Host data root directory not found: ${HOST_DATA_ROOT_ABS}"
HOST_OUTPUT_DIR="$(map_container_path_to_host "${SPEC_OUTPUT_DIR}" "${CONTAINER_DATA_ROOT}" "${HOST_DATA_ROOT_ABS}")" || die "Could not map output_dir ${SPEC_OUTPUT_DIR} onto host root ${HOST_DATA_ROOT_ABS}"
[ -d "${HOST_OUTPUT_DIR}" ] || die "Host output directory not found: ${HOST_OUTPUT_DIR}"

REPO_CACHE_DIR="${PIPELINE_ROOT}/cache/docker_ubuntu_validation"
/usr/bin/mkdir -p "${REPO_CACHE_DIR}"
timestamp="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${REPO_CACHE_DIR}/docker_gnuplot_smoke_${timestamp}.log"
CONTAINER_HELPER="${REPO_CACHE_DIR}/test_ubuntu_docker_gnuplot.container.sh"
DOCKER_PIPELINE_ROOT="$(docker_mount_path "${PIPELINE_ROOT}")"

exec > >(tee -a "${LOG_FILE}") 2>&1

log "Docker smoke-test log: ${LOG_FILE}"
log "Spec: ${SPEC_ABS}"
log "Docker image: ${IMAGE_NAME}"
log "Target SNP: ${TARGET_SNP}"
log "Plots: ${PLOTS}"
log "Container data root: ${CONTAINER_DATA_ROOT}"
log "Host data root: ${HOST_DATA_ROOT_ABS}"
log "Host output dir: ${HOST_OUTPUT_DIR}"

if [ "${SKIP_BUILD}" -eq 0 ]; then
  build_started="$(date +%s)"
  log "Building Docker image ${IMAGE_NAME}"
  docker build -t "${IMAGE_NAME}" "${DOCKER_PIPELINE_ROOT}"
  log "Docker image build finished in $(($(date +%s) - build_started))s"
else
  docker_image_exists "${IMAGE_NAME}" || die "Docker image not found: ${IMAGE_NAME}"
  log "Reusing existing Docker image ${IMAGE_NAME}"
fi

MANHATTAN_PREFIX="${OUTPUT_STEM}_manhattan"
LOCAL_PREFIX="${OUTPUT_STEM}_local_top_hits_manhattan"
LOCAL_TOP_HITS_CSV="${OUTPUT_STEM}_top_hits.csv"
LOCAL_GTF_HTML="${OUTPUT_STEM}_local_top_hits_with_gtf.html"
LOCAL_GTF_BASE="${LOCAL_GTF_HTML%.html}"

{
  printf '#!/usr/bin/env bash\n'
  printf 'set -euo pipefail\n'
  printf 'export SMOKE_SPEC_SRC=%q\n' "/hostspec/${SPEC_BASENAME}"
  printf 'export SMOKE_SPEC_DST=%q\n' "/tmp/${SPEC_BASENAME%.json}.docker_gnuplot.json"
  printf 'export SMOKE_OUTPUT_PREFIX=%q\n' "${MANHATTAN_PREFIX}"
  printf 'export SMOKE_LOCAL_OUTPUT_PREFIX=%q\n' "${LOCAL_PREFIX}"
  printf 'export SMOKE_LOCAL_TOP_HITS_CSV=%q\n' "${LOCAL_TOP_HITS_CSV}"
  printf 'export SMOKE_OUTPUT_HTML_BASENAME=%q\n' "${LOCAL_GTF_HTML}"
  printf 'export SMOKE_HTML_TITLE_SUFFIX=%q\n' " (Ubuntu Docker smoke test)"
  printf 'export SMOKE_LOCAL_HTML_TITLE=%q\n' "PGC schizophrenia local plots (Ubuntu Docker smoke test)"
  printf 'export SMOKE_TARGET_SNP=%q\n' "${TARGET_SNP}"
  printf 'export SMOKE_PLOTS=%q\n' "${PLOTS}"
  printf 'export SMOKE_DISPLAY_GWAS=%q\n' "${DISPLAY_GWAS}"
  printf 'export SMOKE_FORCE=%q\n' "${FORCE}"
  cat <<'EOF'
cd /opt/MultiGWAS-Explorer

python3 - <<'PY'
import json
import os
from pathlib import Path

src = Path(os.environ["SMOKE_SPEC_SRC"])
dst = Path(os.environ["SMOKE_SPEC_DST"])
data = json.loads(src.read_text())
data["workdir"] = "/opt/MultiGWAS-Explorer"
data["output_prefix"] = os.environ["SMOKE_OUTPUT_PREFIX"]
data["local_output_prefix"] = os.environ["SMOKE_LOCAL_OUTPUT_PREFIX"]
data["local_top_hits_csv_basename"] = os.environ["SMOKE_LOCAL_TOP_HITS_CSV"]
data["output_html_basename"] = os.environ["SMOKE_OUTPUT_HTML_BASENAME"]
data["html_title"] = (data.get("html_title") or "Gunplot Manhattan Plot") + os.environ["SMOKE_HTML_TITLE_SUFFIX"]
data["local_html_title"] = os.environ["SMOKE_LOCAL_HTML_TITLE"]
data["open_result"] = 0
dst.write_text(json.dumps(data, indent=2) + "\n")
print(dst)
PY

cmd=(
  perl ./auto_prepare_and_run_diff_gwas_with_gunplot.pl
  --spec "${SMOKE_SPEC_DST}"
  --plots "${SMOKE_PLOTS}"
  --target-snps "${SMOKE_TARGET_SNP}"
)

if [ -n "${SMOKE_DISPLAY_GWAS}" ]; then
  cmd+=(--display-gwas "${SMOKE_DISPLAY_GWAS}")
fi
if [ "${SMOKE_FORCE}" = "1" ]; then
  cmd+=(--force)
fi

printf '[install] Running Docker gunplot smoke test:'
printf ' %q' "${cmd[@]}"
printf '\n'
"${cmd[@]}"
EOF
} > "${CONTAINER_HELPER}"

chmod +x "${CONTAINER_HELPER}"

DOCKER_SPEC_DIR="$(docker_mount_path "${SPEC_DIR}")"
DOCKER_DATA_ROOT="$(docker_mount_path "${HOST_DATA_ROOT_ABS}")"
DOCKER_HELPER="$(docker_mount_path "${CONTAINER_HELPER}")"

run_started="$(date +%s)"
docker run --rm \
  -e PIPELINE_WORKDIR=/opt/MultiGWAS-Explorer \
  -v "${DOCKER_SPEC_DIR}:/hostspec:ro" \
  -v "${DOCKER_DATA_ROOT}:${CONTAINER_DATA_ROOT}" \
  -v "${DOCKER_HELPER}:/tmp/test_ubuntu_docker_gnuplot.sh:ro" \
  "${IMAGE_NAME}" \
  bash /tmp/test_ubuntu_docker_gnuplot.sh
log "Docker gunplot smoke test finished in $(($(date +%s) - run_started))s"

require_file() {
  local path="$1"
  [ -s "$path" ] || die "Expected smoke-test output was not created or is empty: $path"
  log "Verified output: $path"
}

if csv_contains "${PLOTS}" "manhattan"; then
  require_file "${HOST_OUTPUT_DIR}/${MANHATTAN_PREFIX}.png"
  require_file "${HOST_OUTPUT_DIR}/${MANHATTAN_PREFIX}.html"
fi

if csv_contains "${PLOTS}" "local_manhattan"; then
  require_file "${HOST_OUTPUT_DIR}/${LOCAL_PREFIX}.png"
  require_file "${HOST_OUTPUT_DIR}/${LOCAL_PREFIX}.html"
  require_file "${HOST_OUTPUT_DIR}/${LOCAL_PREFIX}_${TARGET_SNP}.png"
fi

if csv_contains "${PLOTS}" "local_gtf"; then
  require_file "${HOST_OUTPUT_DIR}/${LOCAL_GTF_BASE}.html"
  require_file "${HOST_OUTPUT_DIR}/${LOCAL_GTF_BASE}_${TARGET_SNP}.png"
fi

if csv_contains "${PLOTS}" "local_manhattan" || csv_contains "${PLOTS}" "local_gtf"; then
  require_file "${HOST_OUTPUT_DIR}/${LOCAL_TOP_HITS_CSV}"
fi

log "Ubuntu Docker gunplot smoke test completed successfully"
