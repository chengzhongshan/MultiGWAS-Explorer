#!/usr/bin/env bash
set -euo pipefail

_pipeline_install_source="${BASH_SOURCE[0]:-$0}"
PIPELINE_INSTALL_DIR="$(cd "$(/usr/bin/dirname "${_pipeline_install_source}")" && pwd)"
PIPELINE_ROOT="$(cd "${PIPELINE_INSTALL_DIR}/.." && pwd)"
PIPELINE_LOCAL_DIR="${PIPELINE_ROOT}/local"
PIPELINE_VENDOR_PERL_DIR="${PIPELINE_ROOT}/vendor/perl5"
PIPELINE_PLATFORM_TAG="${PIPELINE_PLATFORM_TAG:-}"
if [ -z "${PIPELINE_PLATFORM_TAG}" ]; then
  case "$(uname -s 2>/dev/null || printf '%s' unknown)" in
    CYGWIN*) PIPELINE_PLATFORM_TAG="cygwin" ;;
    Linux)   PIPELINE_PLATFORM_TAG="linux" ;;
    Darwin)  PIPELINE_PLATFORM_TAG="darwin" ;;
    MSYS*|MINGW*|Windows_NT) PIPELINE_PLATFORM_TAG="mswin32" ;;
    *)       PIPELINE_PLATFORM_TAG="" ;;
  esac
fi
PIPELINE_PERL_LOCAL_DIR="${PIPELINE_PERL_LOCAL_DIR:-${PIPELINE_LOCAL_DIR}/perl5${PIPELINE_PLATFORM_TAG:+-${PIPELINE_PLATFORM_TAG}}}"
PIPELINE_VENV_DIR="${PIPELINE_ROOT}/.venv-pipeline"
PIPELINE_PYTHON_RECORD_FILE="${PIPELINE_VENV_DIR}/.python-bin"
PIPELINE_REQUIREMENTS_FILE="${PIPELINE_INSTALL_DIR}/requirements-pipeline.txt"
PIPELINE_CPANFILE="${PIPELINE_ROOT}/cpanfile"
PIPELINE_HTSLIB_VERSION="${PIPELINE_HTSLIB_VERSION:-1.20}"
PIPELINE_HTSLIB_URL="${PIPELINE_HTSLIB_URL:-https://github.com/samtools/htslib/releases/download/${PIPELINE_HTSLIB_VERSION}/htslib-${PIPELINE_HTSLIB_VERSION}.tar.bz2}"
PIPELINE_CPANM_BIN=""
PIPELINE_PYTHON_BIN="${PIPELINE_PYTHON_BIN:-}"

if [[ "${PIPELINE_INSTALL_DEBUG:-0}" =~ ^(1|true|yes|y|on)$ ]]; then
  export PS4='+ ${BASH_SOURCE[0]:-bash}:${LINENO}:${FUNCNAME[0]:-main}: '
  if [ -n "${PIPELINE_INSTALL_DEBUG_LOG:-}" ]; then
    mkdir -p "$(/usr/bin/dirname "${PIPELINE_INSTALL_DEBUG_LOG}")"
    exec 9>>"${PIPELINE_INSTALL_DEBUG_LOG}"
    export BASH_XTRACEFD=9
  fi
  set -x
fi

log() {
  printf '[install] %s\n' "$*"
}

warn() {
  printf '[install] WARNING: %s\n' "$*" >&2
}

die() {
  printf '[install] ERROR: %s\n' "$*" >&2
  exit 1
}

command_exists() {
  command -v "$1" >/dev/null 2>&1
}

make_project_scripts_executable() {
  log "Ensuring project shell and Perl scripts are executable"
  find "${PIPELINE_ROOT}" -type f \( -name '*.sh' -o -name '*.pl' \) -exec chmod a+x {} +
}

num_cpus() {
  if command_exists nproc; then
    nproc
  elif command_exists sysctl; then
    sysctl -n hw.ncpu 2>/dev/null || echo 2
  else
    echo 2
  fi
}

prepend_path() {
  local dir="$1"
  [ -d "$dir" ] || return 0
  case ":${PATH:-}:" in
    *":${dir}:"*) ;;
    *) PATH="${dir}${PATH:+:${PATH}}" ;;
  esac
  export PATH
}

prepend_env_list() {
  local var="$1"
  local dir="$2"
  local current="${!var:-}"
  [ -d "$dir" ] || return 0
  case ":${current}:" in
    *":${dir}:"*) ;;
    *)
      if [ -n "$current" ]; then
        export "${var}=${dir}:${current}"
      else
        export "${var}=${dir}"
      fi
      ;;
  esac
}

is_perl_arch_dir() {
  local dir="$1"
  local name=""
  [ -d "$dir" ] || return 1
  name="$(basename "$dir")"
  case "$name" in
    *-thread-multi|*linux*|*gnu*|*darwin*|*MSWin32*|*cygwin*|x86_64*|aarch64*|arm64*|i[3-6]86*)
      return 0
      ;;
  esac
  [ -d "${dir}/auto" ] && return 0
  return 1
}

current_perl_archname() {
  if [ -n "${PIPELINE_PERL_ARCHNAME:-}" ]; then
    printf '%s\n' "${PIPELINE_PERL_ARCHNAME}"
    return 0
  fi
  PIPELINE_PERL_ARCHNAME="$(perl -MConfig -e 'print $Config{archname}' 2>/dev/null || true)"
  printf '%s\n' "${PIPELINE_PERL_ARCHNAME}"
}

perl_arch_matches_current() {
  local dir="$1"
  local name="" current_arch="" current_os=""
  [ -d "$dir" ] || return 1
  name="$(basename "$dir")"
  current_arch="$(current_perl_archname)"
  current_os="$(uname -s 2>/dev/null | tr '[:upper:]' '[:lower:]')"

  if [ -n "$current_arch" ]; then
    case "$(printf '%s' "$name" | tr '[:upper:]' '[:lower:]')" in
      "$(printf '%s' "$current_arch" | tr '[:upper:]' '[:lower:]')" )
        return 0
        ;;
    esac
    case "$(printf '%s' "$current_arch" | tr '[:upper:]' '[:lower:]')" in
      *"$(printf '%s' "$name" | tr '[:upper:]' '[:lower:]')"* )
        return 0
        ;;
    esac
    case "$(printf '%s' "$name" | tr '[:upper:]' '[:lower:]')" in
      *"$(printf '%s' "$current_arch" | tr '[:upper:]' '[:lower:]')"* )
        return 0
        ;;
    esac
  fi

  case "$current_os" in
    cygwin*)
      [[ "$name" == *cygwin* ]]
      return
      ;;
    linux*)
      [[ "$name" == *linux* || "$name" == *gnu* ]]
      return
      ;;
    darwin*)
      [[ "$name" == *darwin* ]]
      return
      ;;
    msys*|mingw*|windows_nt*)
      [[ "$name" == *MSWin32* ]]
      return
      ;;
  esac
  return 1
}

download_url() {
  local url="$1"
  local dest="$2"
  mkdir -p "$(dirname "$dest")"
  if command_exists curl; then
    if [ "${PIPELINE_CURL_INSECURE:-0}" = "1" ]; then
      curl -kLfsS "$url" -o "$dest"
    else
      curl -LfsS "$url" -o "$dest"
    fi
  elif command_exists wget; then
    if [ "${PIPELINE_CURL_INSECURE:-0}" = "1" ]; then
      wget --no-check-certificate -qO "$dest" "$url"
    else
      wget -qO "$dest" "$url"
    fi
  else
    die "Need curl or wget to download ${url}"
  fi
}

find_system_python() {
  local cand
  for cand in \
    /usr/bin/python3 \
    /usr/bin/python \
    /bin/python3 \
    /bin/python \
    /usr/bin/python3.* \
    /bin/python3.* \
    /usr/local/bin/python3 \
    /usr/local/bin/python \
    python3 \
    python; do
    if [ -x "$cand" ] || command_exists "$cand"; then
      if "$cand" -c "import sys; print(sys.executable)" >/dev/null 2>&1; then
        command -v "$cand"
        return 0
      fi
    fi
  done
  return 1
}

resolve_venv_python() {
  local cand
  if [ -f "${PIPELINE_PYTHON_RECORD_FILE}" ]; then
    if IFS= read -r cand < "${PIPELINE_PYTHON_RECORD_FILE}"; then
      if [ -n "$cand" ] && { [ -x "$cand" ] || command_exists "$cand"; }; then
        printf '%s\n' "$cand"
        return 0
      fi
    fi
  fi
  for cand in \
    "${PIPELINE_VENV_DIR}/bin/python" \
    "${PIPELINE_VENV_DIR}/bin/python3" \
    "${PIPELINE_VENV_DIR}/Scripts/python.exe" \
    "${PIPELINE_VENV_DIR}/Scripts/python"; do
    if [ -x "$cand" ]; then
      printf '%s\n' "$cand"
      return 0
    fi
  done
  return 1
}

resolve_venv_site_packages() {
  local site
  for site in \
    "${PIPELINE_VENV_DIR}/Lib/site-packages" \
    "${PIPELINE_VENV_DIR}/lib/site-packages" \
    "${PIPELINE_VENV_DIR}"/lib/python*/site-packages; do
    if [ -d "$site" ]; then
      printf '%s\n' "$site"
      return 0
    fi
  done
  return 1
}

resolve_saspy_site_packages() {
  local site_packages=""
  site_packages="$(resolve_venv_site_packages || true)"
  if [ -n "$site_packages" ] && [ -d "${site_packages}/saspy" ]; then
    printf '%s\n' "$site_packages"
    return 0
  fi
  return 1
}

resolve_windows_java_for_saspy() {
  local cand="" unix_cand=""
  if command_exists cygpath; then
    for unix_cand in \
      /usr/local/jdk/bin/java.exe \
      /usr/bin/java.exe \
      /bin/java.exe; do
      if [ -f "$unix_cand" ]; then
        cygpath -w "$unix_cand"
        return 0
      fi
    done
  fi
  for cand in \
    "${SASPY_JAVA_WIN:-}" \
    'C:\Program Files (x86)\Common Files\Oracle\Java\java8path\java.exe' \
    'C:\Program Files\Common Files\Oracle\Java\java8path\java.exe' \
    'C:\Program Files (x86)\Common Files\Oracle\Java\javapath\java.exe' \
    'C:\Program Files\Common Files\Oracle\Java\javapath\java.exe'; do
    [ -n "$cand" ] || continue
    if command_exists cygpath; then
      unix_cand="$(cygpath -u "$cand" 2>/dev/null || true)"
      if [ -n "$unix_cand" ] && [ -f "$unix_cand" ]; then
        printf '%s\n' "$cand"
        return 0
      fi
    elif [ -f "$cand" ]; then
      printf '%s\n' "$cand"
      return 0
    fi
  done
  return 1
}

resolve_saspy_reference_dir() {
  local candidate="" unix_candidate=""
  local -a candidates=()
  if [ -n "${PIPELINE_SASPY_REFERENCE_DIR:-}" ]; then
    candidates+=("${PIPELINE_SASPY_REFERENCE_DIR}")
  fi
  if [ -n "${PIPELINE_SASPY_REFERENCE_SITEPKG:-}" ]; then
    candidates+=("${PIPELINE_SASPY_REFERENCE_SITEPKG}/saspy" "${PIPELINE_SASPY_REFERENCE_SITEPKG}")
  fi
  if [ -n "${SASPY_REFERENCE_SITEPKG:-}" ]; then
    candidates+=("${SASPY_REFERENCE_SITEPKG}/saspy" "${SASPY_REFERENCE_SITEPKG}")
  fi
  candidates+=("${PIPELINE_INSTALL_DIR}/saspy-java-supplement")
  if [ -n "${USERPROFILE:-}" ]; then
    candidates+=(
      "${USERPROFILE}\\Downloads\\cygwin-portable-20210411\\cygwin-portable\\App\\cygwin\\usr\\local\\lib\\python3.9\\site-packages\\saspy"
      "${USERPROFILE}\\Downloads\\cygwin-portable-20210411\\cygwin-portable\\App\\cygwin\\usr\\local\\lib\\python3.9\\site-packages"
    )
  fi
  if [ -n "${HOME:-}" ]; then
    candidates+=(
      "/usr/local/anaconda3/lib/python3.12/site-packages/saspy"
      "/usr/local/anaconda3/lib/python3.12/site-packages"
      "${HOME}/Desktop/shared/SynchronizationVersions/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/.venv-pipeline/Lib/site-packages/saspy"
      "${HOME}/Desktop/shared/SynchronizationVersions/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/.venv-pipeline/Lib/site-packages"
    )
  fi
  for candidate in ${candidates[@]+"${candidates[@]}"}; do
    [ -n "$candidate" ] || continue
    unix_candidate="$candidate"
    if command_exists cygpath && [[ "$candidate" == [A-Za-z]:\\* ]]; then
      unix_candidate="$(cygpath -u "$candidate" 2>/dev/null || true)"
    fi
    [ -n "$unix_candidate" ] || continue
    if [ -d "${unix_candidate}/java" ]; then
      printf '%s\n' "${unix_candidate}"
      return 0
    fi
    if [ -d "${unix_candidate}/saspy/java" ]; then
      printf '%s\n' "${unix_candidate}/saspy"
      return 0
    fi
  done
  return 1
}

sync_saspy_reference_java_assets() {
  local saspy_dir="$1"
  local ref_dir=""
  ref_dir="$(resolve_saspy_reference_dir || true)"
  [ -n "$ref_dir" ] || return 0
  [ -d "${ref_dir}/java" ] || return 0
  mkdir -p "${saspy_dir}/java"
  cp -Rf "${ref_dir}/java/." "${saspy_dir}/java/"
  log "Synced reference SASPy Java assets from ${ref_dir} into ${saspy_dir}/java"
}

build_saspy_windows_classpath() {
  local saspy_dir="$1"
  local jar="" classpath_win=""
  local -a ordered_jars=() cpw_parts=()
  local rel=""
  local preferred_rel=(
    "java/saspyiom.jar"
    "java/iomclient/log4j-1.2-api-2.12.4.jar"
    "java/iomclient/log4j-api-2.12.4.jar"
    "java/iomclient/log4j-core-2.12.4.jar"
    "java/iomclient/sas.security.sspi.jar"
    "java/iomclient/sas.core.jar"
    "java/iomclient/sas.svc.connection.jar"
    "java/iomclient/sas.rutil.jar"
    "java/iomclient/sas.rutil.nls.jar"
    "java/iomclient/sastpj.rutil.jar"
    "java/thirdparty/glassfish-corba-internal-api.jar"
    "java/thirdparty/glassfish-corba-omgapi.jar"
    "java/thirdparty/glassfish-corba-orb.jar"
    "java/thirdparty/pfl-basic.jar"
    "java/thirdparty/pfl-tf.jar"
  )

  jar_already_listed() {
    local candidate="$1" item=""
    for item in ${ordered_jars[@]+"${ordered_jars[@]}"}; do
      [ "$item" = "$candidate" ] && return 0
    done
    return 1
  }

  for rel in "${preferred_rel[@]}"; do
    jar="${saspy_dir}/${rel}"
    [ -f "$jar" ] || continue
    ordered_jars+=("$jar")
  done

  for jar in \
    "${saspy_dir}"/java/iomclient/*.jar \
    "${saspy_dir}"/java/thirdparty/*.jar; do
    [ -f "$jar" ] || continue
    jar_already_listed "$jar" && continue
    ordered_jars+=("$jar")
  done

  if command_exists cygpath; then
    for jar in ${ordered_jars[@]+"${ordered_jars[@]}"}; do
      [ -f "$jar" ] || continue
      cpw_parts+=("$(cygpath -w "$jar")")
    done
  fi

  if [ "${#cpw_parts[@]}" -gt 0 ]; then
    local IFS=';'
    classpath_win="${cpw_parts[*]}"
  fi
  printf '%s\n' "${classpath_win}"
}

ensure_authinfo_from_reference_source() {
  [ -f "${HOME}/.authinfo" ] || [ -f "${HOME}/_authinfo" ] && return 0
  local unix_auth="" source_label=""
  local -a auth_candidates=()
  if [ -n "${PIPELINE_ODA_AUTHINFO_SOURCE:-}" ]; then
    auth_candidates+=("${PIPELINE_ODA_AUTHINFO_SOURCE}")
  fi
  if [ -n "${SASPY_AUTHINFO_SOURCE:-}" ]; then
    auth_candidates+=("${SASPY_AUTHINFO_SOURCE}")
  fi
  if [ -n "${HOME:-}" ]; then
    auth_candidates+=("${HOME}/.authinfo" "${HOME}/_authinfo")
  fi
  if [ -n "${USERPROFILE:-}" ]; then
    auth_candidates+=("${USERPROFILE}\\.authinfo" "${USERPROFILE}\\_authinfo")
  fi
  for source_label in ${auth_candidates[@]+"${auth_candidates[@]}"}; do
    [ -n "$source_label" ] || continue
    if [[ "$source_label" == [A-Za-z]:\\* ]]; then
      command_exists cygpath || continue
      unix_auth="$(cygpath -u "$source_label" 2>/dev/null || true)"
    else
      unix_auth="$source_label"
    fi
    [ -n "$unix_auth" ] || continue
    [ -f "$unix_auth" ] || continue
    cp "$unix_auth" "${HOME}/.authinfo"
    chmod 600 "${HOME}/.authinfo" || true
    log "Copied ${source_label} into ${HOME}/.authinfo for SAS ODA authkey reuse"
    return 0
  done
}

build_saspy_unix_classpath() {
  local saspy_dir="$1"
  local jar="" classpath_unix=""
  local -a ordered_jars=()
  local rel=""
  local preferred_rel=(
    "java/saspyiom.jar"
    "java/iomclient/log4j-1.2-api-2.12.4.jar"
    "java/iomclient/log4j-api-2.12.4.jar"
    "java/iomclient/log4j-core-2.12.4.jar"
    "java/iomclient/sas.security.sspi.jar"
    "java/iomclient/sas.core.jar"
    "java/iomclient/sas.svc.connection.jar"
    "java/iomclient/sas.rutil.jar"
    "java/iomclient/sas.rutil.nls.jar"
    "java/iomclient/sastpj.rutil.jar"
    "java/thirdparty/glassfish-corba-internal-api.jar"
    "java/thirdparty/glassfish-corba-omgapi.jar"
    "java/thirdparty/glassfish-corba-orb.jar"
    "java/thirdparty/pfl-basic.jar"
    "java/thirdparty/pfl-tf.jar"
  )
  for rel in "${preferred_rel[@]}"; do
    jar="${saspy_dir}/${rel}"
    [ -f "$jar" ] || continue
    ordered_jars+=("$jar")
  done
  for jar in \
    "${saspy_dir}"/java/iomclient/*.jar \
    "${saspy_dir}"/java/thirdparty/*.jar; do
    [ -f "$jar" ] || continue
    case ":${classpath_unix}:" in
      *":${jar}:"*) continue ;;
    esac
    ordered_jars+=("$jar")
  done
  if [ "${#ordered_jars[@]}" -gt 0 ]; then
    local IFS=':'
    classpath_unix="${ordered_jars[*]}"
  fi
  printf '%s\n' "${classpath_unix}"
}

resolve_unix_java_for_saspy() {
  local cand=""
  for cand in \
    "${SASPY_JAVA:-}" \
    /usr/bin/java \
    /usr/local/bin/java \
    java; do
    [ -n "$cand" ] || continue
    if [ -x "$cand" ]; then
      printf '%s\n' "$cand"
      return 0
    fi
    if command_exists "$cand"; then
      command -v "$cand"
      return 0
    fi
  done
  return 1
}

configure_saspy_oda_profile() {
  activate_python_env
  local site_packages="" saspy_dir="" cfg_file="" java_bin="" classpath_value="" path_sep=":" is_cygwin=0
  site_packages="$(resolve_saspy_site_packages || true)"
  [ -n "$site_packages" ] || return 0
  saspy_dir="${site_packages}/saspy"
  [ -d "$saspy_dir" ] || return 0
  cfg_file="${saspy_dir}/sascfg_personal.py"
  if command_exists uname && uname -s | grep -qi '^CYGWIN'; then
    is_cygwin=1
  fi
  sync_saspy_reference_java_assets "${saspy_dir}"
  if [ "${is_cygwin}" -eq 1 ]; then
    java_bin="$(resolve_windows_java_for_saspy || true)"
    [ -n "$java_bin" ] || warn "Could not resolve a Windows java.exe for SASPy ODA config; SAS ODA sessions may still fail until SASPY_JAVA_WIN is set"
    classpath_value="$(build_saspy_windows_classpath "${saspy_dir}")"
    path_sep=";"
  else
    java_bin="$(resolve_unix_java_for_saspy || true)"
    [ -n "$java_bin" ] || warn "Could not resolve a Java runtime for SASPy ODA config; SAS ODA sessions may still fail until SASPY_JAVA is set"
    classpath_value="$(build_saspy_unix_classpath "${saspy_dir}")"
    path_sep=":"
  fi
  cat > "${cfg_file}" <<PY
import os

_DEFAULT_CLASSPATH = r'${classpath_value}'
cpW = _DEFAULT_CLASSPATH
cpL = _DEFAULT_CLASSPATH

SAS_config_names = ['oda', 'default']
SAS_config_options = {'lock_down': False, 'verbose': True, 'prompt': True}

_ODA_BASE = {
    'java': os.environ.get('SASPY_JAVA_WIN', os.environ.get('SASPY_JAVA', r'${java_bin}')),
    'iomhost': os.environ.get('SASPY_ODA_HOST', 'odaws01-usw2.oda.sas.com'),
    'iomport': int(os.environ.get('SASPY_ODA_PORT', '8591') or '8591'),
    'authkey': os.environ.get('SASPY_ODA_AUTHKEY', 'oda'),
    'encoding': 'utf-8',
    'classpath': os.environ.get('SASPY_ODA_CLASSPATH', _DEFAULT_CLASSPATH),
}

oda = dict(_ODA_BASE)
default = dict(_ODA_BASE)
PY
  ensure_authinfo_from_reference_source
  log "Configured repo-local SASPy ODA profile at ${cfg_file}"
}

python_module_available() {
  local module_name="$1"
  [ -n "${PIPELINE_PYTHON_BIN}" ] || return 1
  "${PIPELINE_PYTHON_BIN}" - <<PY >/dev/null 2>&1
import importlib.util
import sys
sys.exit(0 if importlib.util.find_spec("${module_name}") else 1)
PY
}

prepare_python_requirements() {
  local req_file="${PIPELINE_REQUIREMENTS_FILE}"
  local filtered_req=""

  if command_exists uname && uname -s | grep -qi '^CYGWIN' && python_module_available PIL; then
    mkdir -p "${PIPELINE_VENV_DIR}"
    filtered_req="${PIPELINE_VENV_DIR}/requirements-cygwin.txt"
    grep -viE '^[[:space:]]*Pillow([[:space:]]*([<>=!~].*)?)?$' "${req_file}" > "${filtered_req}"
    printf '%s\n' "${filtered_req}"
    return 0
  fi

  printf '%s\n' "${req_file}"
}

prepare_perl_cpanfile() {
  local cpanfile="${PIPELINE_CPANFILE}"
  local filtered_cpanfile=""

  if command_exists uname && uname -s | grep -qi '^CYGWIN'; then
    mkdir -p "${PIPELINE_LOCAL_DIR}"
    filtered_cpanfile="${PIPELINE_LOCAL_DIR}/cpanfile-cygwin"
    grep -vE "requires '(File::Which|GD|JSON|JSON::MaybeXS|Mojolicious)';" "${cpanfile}" > "${filtered_cpanfile}"
    printf '%s\n' "${filtered_cpanfile}"
    return 0
  fi

  printf '%s\n' "${cpanfile}"
}

activate_python_env() {
  local site_packages=""
  local venv_python=""
  PIPELINE_PYTHON_BIN=""
  venv_python="$(resolve_venv_python || true)"
  if [ -n "$venv_python" ]; then
    PIPELINE_PYTHON_BIN="$venv_python"
    prepend_path "$(/usr/bin/dirname "${PIPELINE_PYTHON_BIN}")"
    site_packages="$("${PIPELINE_PYTHON_BIN}" - <<'PY'
import sysconfig
print(sysconfig.get_path("purelib"))
PY
)"
    if [ -z "$site_packages" ] || [ ! -d "$site_packages" ]; then
      site_packages="$(resolve_venv_site_packages || true)"
    fi
    if [ -n "$site_packages" ] && [ -d "$site_packages" ]; then
      prepend_env_list PYTHONPATH "$site_packages"
    fi
  fi
  export PIPELINE_PYTHON_BIN
}

create_python_venv() {
  local sys_python="${1:-}"
  local target_site="${PIPELINE_VENV_DIR}/Lib/site-packages"
  local target_mode=0
  local req_file=""
  if [ -z "$sys_python" ]; then
    sys_python="$(find_system_python || true)"
  fi
  [ -n "$sys_python" ] || die "Could not find python3 or python"

  if command_exists uname && uname -s | grep -qi '^CYGWIN'; then
    log "Using repo-local Python site-packages under ${target_site} for portable Cygwin"
    rm -rf "${PIPELINE_VENV_DIR}"
    mkdir -p "${target_site}"
    printf '%s\n' "$sys_python" > "${PIPELINE_PYTHON_RECORD_FILE}"
    target_mode=1
  elif [ -z "$(resolve_venv_python || true)" ]; then
    log "Creating repo-local Python environment under ${PIPELINE_VENV_DIR}"
    if ! "$sys_python" -m venv "${PIPELINE_VENV_DIR}"; then
      warn "Python venv creation failed; falling back to repo-local site-packages under ${target_site}"
      rm -rf "${PIPELINE_VENV_DIR}"
      mkdir -p "${target_site}"
      printf '%s\n' "$sys_python" > "${PIPELINE_PYTHON_RECORD_FILE}"
      target_mode=1
    fi
  elif [ -f "${PIPELINE_PYTHON_RECORD_FILE}" ]; then
    target_mode=1
  fi

  activate_python_env
  req_file="$(prepare_python_requirements)"
  if [ "${target_mode}" -eq 0 ] && ! "${PIPELINE_PYTHON_BIN}" -m pip --version >/dev/null 2>&1; then
    warn "The created Python environment does not contain pip; falling back to repo-local site-packages under ${target_site}"
    rm -rf "${PIPELINE_VENV_DIR}"
    mkdir -p "${target_site}"
    printf '%s\n' "$sys_python" > "${PIPELINE_PYTHON_RECORD_FILE}"
    target_mode=1
    activate_python_env
  fi
  if [ "${target_mode}" -eq 1 ]; then
    "${PIPELINE_PYTHON_BIN}" -m pip install --upgrade -r "${req_file}" --target "${target_site}"
  else
    "${PIPELINE_PYTHON_BIN}" -m pip install --upgrade pip setuptools wheel
    "${PIPELINE_PYTHON_BIN}" -m pip install -r "${req_file}"
  fi
  configure_saspy_oda_profile
  if [ -n "${req_file}" ] && [ "${req_file}" != "${PIPELINE_REQUIREMENTS_FILE}" ]; then
    rm -f "${req_file}"
  fi
}

activate_perl_env() {
  local base="${PIPELINE_PERL_LOCAL_DIR}/lib/perl5"
  local arch
  mkdir -p "${PIPELINE_PERL_LOCAL_DIR}"
  prepend_path "${PIPELINE_PERL_LOCAL_DIR}/bin"
  prepend_path "${PIPELINE_LOCAL_DIR}/bin"
  if [ -d "${PIPELINE_VENDOR_PERL_DIR}" ]; then
    prepend_env_list PERL5LIB "${PIPELINE_VENDOR_PERL_DIR}"
  fi
  if [ -d "$base" ]; then
    prepend_env_list PERL5LIB "$base"
    for arch in "$base"/*; do
      is_perl_arch_dir "$arch" || continue
      perl_arch_matches_current "$arch" || continue
      prepend_env_list PERL5LIB "$arch"
    done
  fi
  export PERL_LOCAL_LIB_ROOT="${PIPELINE_PERL_LOCAL_DIR}${PERL_LOCAL_LIB_ROOT:+:${PERL_LOCAL_LIB_ROOT}}"
  export PERL_MB_OPT="--install_base ${PIPELINE_PERL_LOCAL_DIR}"
  export PERL_MM_OPT="INSTALL_BASE=${PIPELINE_PERL_LOCAL_DIR}"
}

ensure_cpanm() {
  activate_perl_env
  if command_exists cpanm; then
    PIPELINE_CPANM_BIN="$(command -v cpanm)"
    return 0
  fi
  PIPELINE_CPANM_BIN="${PIPELINE_LOCAL_DIR}/bin/cpanm"
  if [ ! -f "${PIPELINE_CPANM_BIN}" ]; then
    log "Bootstrapping cpanminus into ${PIPELINE_CPANM_BIN}"
    download_url "https://cpanmin.us/" "${PIPELINE_CPANM_BIN}"
    chmod +x "${PIPELINE_CPANM_BIN}"
  fi
}

install_perl_deps() {
  local cpanfile_to_use=""
  local module_name=""
  local modules=()
  local regular_modules=()
  local needs_pdl=0
  activate_perl_env
  activate_python_env
  ensure_cpanm
  cpanfile_to_use="$(prepare_perl_cpanfile)"
  log "Installing repo-local Perl dependencies from ${cpanfile_to_use}"
  while IFS= read -r module_name; do
    [ -n "${module_name}" ] || continue
    modules+=("${module_name}")
  done < <(awk -F"'" '/^[[:space:]]*requires[[:space:]]+/ { print $2 }' "${cpanfile_to_use}")
  [ "${#modules[@]}" -gt 0 ] || die "No Perl modules were parsed from ${cpanfile_to_use}"
  for module_name in "${modules[@]}"; do
    case "${module_name}" in
      PDL) needs_pdl=1 ;;
      *) regular_modules+=("${module_name}") ;;
    esac
  done
  if [ "${#regular_modules[@]}" -gt 0 ]; then
    perl "${PIPELINE_CPANM_BIN}" \
      --local-lib-contained "${PIPELINE_PERL_LOCAL_DIR}" \
      --notest \
      "${regular_modules[@]}"
  fi
  if [ "${needs_pdl}" -eq 1 ]; then
    install_pdl_perl_deps
  fi
  if [ -n "${cpanfile_to_use}" ] && [ "${cpanfile_to_use}" != "${PIPELINE_CPANFILE}" ]; then
    rm -f "${cpanfile_to_use}"
  fi
  activate_perl_env
}

install_pdl_perl_deps() {
  if perl -MPDL -e1 >/dev/null 2>&1; then
    log "PDL is already installed and loadable"
    return 0
  fi
  log "Installing PDL with extended Cygwin-friendly build timeouts"
  MAKEFLAGS="${MAKEFLAGS:--j$(num_cpus)}" perl "${PIPELINE_CPANM_BIN}" \
    --local-lib-contained "${PIPELINE_PERL_LOCAL_DIR}" \
    --notest \
    --configure-timeout 900 \
    --build-timeout 7200 \
    PDL
}

ensure_local_hts_tools() {
  prepend_path "${PIPELINE_LOCAL_DIR}/bin"
  if command_exists bgzip && command_exists tabix; then
    log "Using bgzip/tabix from PATH"
    return 0
  fi
  if [ -x "${PIPELINE_LOCAL_DIR}/bin/bgzip.exe" ] && [ -x "${PIPELINE_LOCAL_DIR}/bin/tabix.exe" ]; then
    log "Using bundled Windows bgzip/tabix under ${PIPELINE_LOCAL_DIR}/bin"
    return 0
  fi
  log "bgzip/tabix not found; building a repo-local htslib copy"
  "${PIPELINE_INSTALL_DIR}/build_local_htslib.sh"
  prepend_path "${PIPELINE_LOCAL_DIR}/bin"
}

run_pipeline_check() {
  (cd "${PIPELINE_ROOT}" && "${PIPELINE_INSTALL_DIR}/check_pipeline_install.sh")
}
