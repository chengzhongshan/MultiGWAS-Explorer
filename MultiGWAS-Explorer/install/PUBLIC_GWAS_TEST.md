# End-to-end test with public female and male GWAS

This test uses the European female and male schizophrenia GWAS from
[Trubetskoy et al., Nature 2022](https://doi.org/10.1038/s41586-022-04434-5),
distributed in the [PGC scz2022 release](https://figshare.com/articles/dataset/scz2022/19426775).
It includes autosomes and chromosome X, in GRCh37/hg19 coordinates. Four
compressed source files total approximately 820 MB. Keep several GB free for
intermediate tables. A full run can take hours on a slow computer.

The downloader, orchestration, and numerical validation are Perl scripts.
SAS plotting still uses the pipeline's SASPy dependency. Read the source
release documentation before interpreting results; this is a software
integration test, including chromosome X handling, rather than validation of
every biological assumption in a sex comparison.

## Ubuntu or Windows portable Cygwin

Open the same Cygwin installation that installed the dependencies on Windows.
Perl modules built with a different Perl version can fail to load even when
both installations are Cygwin. On Ubuntu use a terminal.

Start in the inner `MultiGWAS-Explorer` directory containing `install/`:

```bash
. install/common.sh
activate_perl_env
activate_python_env
export PATH="$PIPELINE_ROOT/local/bin:$PIPELINE_ROOT/DiffGWASDeps:$PATH"
bash install/check_pipeline_install.sh

perl install/run_public_sex_gwas.pl \
  --output-dir "$PWD/public-gwas-test" \
  --backend gnuplot
```

To test both plotting backends with an already configured SAS ODA account:

```bash
perl install/run_public_sex_gwas.pl \
  --output-dir "$PWD/public-gwas-test" \
  --backend both
```

Do not put SAS passwords in command lines, scripts, or Git. Use the existing
SASPy configuration/authinfo setup supported by the installer.

For existing downloads, add `--input-dir /path/to/PGC_SCZ_Sex_Stratified_GWASs`.
The four exact filenames and published MD5 checksums are in
`prepare_public_sex_gwas.pl`; existing files are verified before use.
An existing GENCODE cache can be supplied with `--gtf-cache-dir /path/to/cache/gtf`.

## Run and inspect each phase

```bash
perl install/run_public_sex_gwas.pl --output-dir "$PWD/public-gwas-test" --phase prepare
perl install/run_public_sex_gwas.pl --output-dir "$PWD/public-gwas-test" --phase preprocess
perl install/run_public_sex_gwas.pl --output-dir "$PWD/public-gwas-test" --phase validate
perl install/run_public_sex_gwas.pl --output-dir "$PWD/public-gwas-test" --phase plots --backend both
```

| Phase | What it checks |
| --- | --- |
| Prepare | Public downloads, MD5 integrity, explicit sex groups/build, separate run configuration directory |
| Preprocess | Raw OR-to-log-effect conversion, merge, coordinate sort, bgzip/tabix, conservative allele/frequency filtering, differential calculation, standardized and wide tables |
| Validate | Every differential row: finite values, positive SEs, female-minus-male effects, independent SE/Z/P calculations, ambiguous-allele exclusion, frequency discordance, presence of all autosomes and X |
| Plots | Genome-wide Manhattan, local Manhattan, GENCODE gene tracks, forest plots; inquiry loci selected from real output for common association, differential association, and X |

The numerical validator uses `POSIX::erfc` independently of the production
P-value approximation. It assumes `rho=0` for disjoint female/male strata.
The chosen inquiry SNPs are test targets; selecting a minimum P value does
not establish statistical significance. Inquiry plots do not test automatic
LD clumping. The generated spec retains MAF/LD settings for separate automatic
lead-selection runs; report LD service/reference failures separately.

Inspect `source_manifest.json`, the pipeline's row-count/harmonization
manifests, `numeric_validation.json`, `targets.txt`, and `test_run_*.json`.
The run records report process exit status; inspect generated PNG/HTML files
and SAS logs as well. A successful login alone does not establish that SAS
plots succeeded. Keep source data and generated results outside Git.

The small `test_sort_long_gwas.pl` regression is also part of the installation
check. It tests a genuine tabix index and autosomal/X interval queries using a
directory with spaces. It is a synthetic edge-case test, separate from the
public-data workflow above.
