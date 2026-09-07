#!/usr/bin/env bash
set -euo pipefail

# Download and prepare the official PLINK2 1000 Genomes phase-3 reference.
# The files are intentionally kept outside git.  Use --whole to prepare the
# complete reference and --make-bed to create a biallelic-SNP BED fileset.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
OUT_DIR="${PLINK2_1KG_OUT_DIR:-${SCRIPT_DIR}/../cache/plink2_1kg_phase3}"
CHR="${PLINK2_1KG_CHR:-2}"
WHOLE=0
MAKE_BED=0
PLINK2_BIN="${PLINK2:-plink2}"
FORCE=0
while (($#)); do
  case "$1" in
    --chr) CHR="$2"; shift 2 ;;
    --whole) WHOLE=1; shift ;;
    --make-bed) MAKE_BED=1; shift ;;
    --output-dir) OUT_DIR="$2"; shift 2 ;;
    --plink2) PLINK2_BIN="$2"; shift 2 ;;
    --force) FORCE=1; shift ;;
    -h|--help) echo 'Usage: prepare_plink2_1kg_phase3_reference.sh [--whole|--chr 2] [--make-bed] [--output-dir DIR] [--plink2 EXE] [--force]'; exit 0 ;;
    *) echo "ERROR: unknown option: $1" >&2; exit 2 ;;
  esac
done
mkdir -p "$OUT_DIR"
if [[ "$WHOLE" == 1 ]]; then
  prefix="$OUT_DIR/all_phase3"
  pgen_url='https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst?dl=1'
  pvar_url='https://www.dropbox.com/s/odlexvo8fummcvt/all_phase3.pvar.zst?dl=1'
else
  [[ "$CHR" == 2 ]] || { echo "ERROR: chromosome bootstrap currently pins chr2 URLs" >&2; exit 2; }
  prefix="$OUT_DIR/chr2_phase3"
  pgen_url='https://www.dropbox.com/s/6qlhq2mdawa27f2/chr2_phase3.pgen.zst?dl=1'
  pvar_url='https://www.dropbox.com/s/bof8v3odxtd8ihm/chr2_phase3.pvar.zst?dl=1'
fi
psam_url='https://www.dropbox.com/scl/fi/haqvrumpuzfutklstazwk/phase3_corrected.psam?dl=1&rlkey=0yyifzj2fb863ddbmsv4jkeq6'

download() {
  local url="$1" dest="$2"
  if [[ "$FORCE" != 1 && -s "$dest" ]]; then echo "[reuse] $dest"; return; fi
  echo "[download] $dest"
  if command -v curl.exe >/dev/null 2>&1; then
    curl.exe -L --fail --retry 3 --retry-delay 3 -o "$dest" "$url"
  elif command -v curl >/dev/null 2>&1; then
    curl -L --fail --retry 3 --retry-delay 3 -o "$dest" "$url"
  else
    echo 'ERROR: curl is required to download the PLINK2 reference.' >&2; exit 1
  fi
}
download "$pgen_url" "${prefix}.pgen.zst"
download "$pvar_url" "${prefix}.pvar.zst"
download "$psam_url" "${prefix}.psam"
if [[ ! -s "${prefix}.pgen" ]]; then
  echo "[prepare] Decompressing ${prefix}.pgen.zst"
  "$PLINK2_BIN" --zst-decompress "${prefix}.pgen.zst" "${prefix}.pgen"
fi
[[ -s "${prefix}.pgen" && -s "${prefix}.pvar.zst" && -s "${prefix}.psam" ]] || {
  echo "ERROR: Prepared PLINK2 fileset is incomplete: ${prefix}" >&2; exit 1;
}
if [[ "$MAKE_BED" == 1 ]]; then
  bed_prefix="${prefix}_biallelic"
  if [[ "$FORCE" == 1 || ! -s "${bed_prefix}.bed" || ! -s "${bed_prefix}.bim" || ! -s "${bed_prefix}.fam" ]]; then
    echo "[prepare] Creating biallelic ACGT PLINK BED fileset: ${bed_prefix}"
    "$PLINK2_BIN" --pfile "$prefix" vzs --snps-only just-acgt --max-alleles 2 --make-bed --out "$bed_prefix"
  fi
  echo "PLINK2_1KG_BFILE=${bed_prefix}"
fi
echo "PLINK2_1KG_PFILE=${prefix}"
echo 'PLINK2_1KG_BUILD=GRCh37_hg19'
echo "PLINK2_1KG_CHR=$([[ "$WHOLE" == 1 ]] && echo whole || echo "$CHR")"
