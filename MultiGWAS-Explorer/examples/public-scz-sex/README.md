# Public schizophrenia GWAS by sex: completed example

This directory contains the presentation outputs from the September 2026
end-to-end validation with the public PGC European female and male
schizophrenia summary statistics (GRCh37/hg19). It is a software validation
record rather than a claim about biological significance.

Open `results.html` after cloning the repository. It contains 11 gnuplot and 11
SAS OnDemand for Academics images: genome-wide Manhattan, local Manhattan,
local gene-track, and forest plots. `results_gallery_manifest.json` records the
SHA-256 checksum of each copied image.

## Validation summary

- 15,758,068 public source records were processed.
- 6,650,636 female-versus-male comparisons passed independent numerical
  validation, including 194,288 on chromosome X.
- Differential effects were checked independently as female minus male;
  standard errors, Z scores, and two-sided P values were recalculated.
- All four plot families decoded successfully for both backends.
- Inquiry targets were `rs2232429`, `rs185665940`, and `rs62604261`.
- Gnuplot local Manhattan and GTF heatmaps use phased EUR 1000 Genomes Phase 3
  threshold-0 LD caches. Their complete sidecars contain 18,917, 15,367, and
  7,519 non-reference R² values; the 0.1 cutoff controls marker symbols only.

The machine-readable row counts and targets are in `numeric_validation.json`
and `targets.txt`. See the
[full validation report](../../install/PUBLIC_GWAS_RESULTS.md) for scope and
limitations.

## Reproduce the run

From the inner `MultiGWAS-Explorer` directory:

```bash
. install/common.sh
activate_perl_env
activate_python_env
export PATH="$PIPELINE_ROOT/local/bin:$PATH"
bash install/check_pipeline_install.sh
perl install/run_public_sex_gwas.pl \
  --output-dir "$PWD/public-gwas-test" \
  --backend both
```

The complete command downloads about 820 MB and needs several GB for generated
tables. The SAS backend requires an existing SAS ODA/SASPy configuration. The
[real-data test guide](../../install/PUBLIC_GWAS_TEST.md) explains phased runs,
reuse of existing downloads, output validation, and troubleshooting.
