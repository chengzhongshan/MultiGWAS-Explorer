# MultiGWAS-Explorer

MultiGWAS-Explorer compares GWAS summary statistics and produces genome-wide
Manhattan plots, local Manhattan plots, gene-track plots, and forest plots.
It supports local **gnuplot** rendering and **SAS OnDemand for Academics**,
with Perl scripts for preprocessing and workflow automation.

<img width="2400" height="2510" alt="MultiGWAS-Explorer pipeline overview" src="https://github.com/user-attachments/assets/a56b8675-2a48-41c9-ba2c-1d1e793603ee" />

## Features

- Compare sex-, population-, or dataset-stratified GWAS summary statistics.
- Prepare coordinate-sorted, indexed tables and pairwise differential effects.
- Plot common-association signals, differential signals, or selected inquiry SNPs.
- Run from the command line or through the Perl MCP interface with different AI agents.

## Install

Supported installation paths include Windows portable Cygwin, Ubuntu/Linux,
macOS Intel and Apple Silicon, Docker, and Apptainer.

```bash
git clone https://github.com/chengzhongshan/MultiGWAS-Explorer.git
cd MultiGWAS-Explorer/MultiGWAS-Explorer
```

Choose your platform in the [installation and usage guide](docs/README.md#installation).
It includes prerequisites, commands, container setup, SAS configuration, and
troubleshooting links. Installation commands run from the inner pipeline
directory shown above.

## Try a local example

After installing:

```bash
bash install/check_pipeline_install.sh
bash install/run_plotting_example.sh
```

Open `example-output/example.png`. This synthetic example needs no GWAS
download or SAS account.

## Test with public GWAS data

The [Perl real-data test guide](MultiGWAS-Explorer/install/PUBLIC_GWAS_TEST.md)
provides a complete female/male schizophrenia example with autosomes and X.
It covers download checksums, numerical validation, and both plotting backends.
Run it from the inner pipeline directory after installation:

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

The command downloads the public PGC schizophrenia female/male files, verifies
their published checksums, validates all differential-effect calculations, and
renders gnuplot and SAS ODA figures. SAS ODA must already be configured. The
four downloads are about 820 MB and the complete test needs several GB of free
space.

### Completed real-data gallery

The following 18 figures come from the full 6,650,636-comparison Windows
portable-Cygwin validation of the sex-stratified schizophrenia GWAS. Inquiry
targets are `rs2232429` (common association), `rs185665940` (autosomal
differential signal), and `rs62604261` (chromosome-X differential signal).

The local Manhattan panels show -log10(P) by Z-score color in gnuplot and by chromosome color in SAS. Gene-track panels use signed LD r² × sign(Z) from the EUR 1000 Genomes Phase 3 reference. Each figure opens at full resolution when clicked.

#### SAS OnDemand for Academics outputs

**Genome-wide Manhattan: chromosome-colored association P**

<a href="MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_manhattan.png"><img src="MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_manhattan.png" alt="Genome-wide Manhattan: chromosome-colored association P" width="1100"></a>

**Local Manhattan: chromosome colors**

<a href="MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_local_top_hits_manhattan.png"><img src="MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_local_top_hits_manhattan.png" alt="Local Manhattan: chromosome colors" width="1000"></a>

**Local Manhattan with gene tracks: LD r² × sign(Z) — rs2232429**

<a href="MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_local_top_hits_with_gtf_rs2232429.png"><img src="MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_local_top_hits_with_gtf_rs2232429.png" alt="Local Manhattan with gene tracks: LD r² × sign(Z) — rs2232429" width="760"></a>

<details><summary>View the remaining 4 SAS ODA figures</summary>

**Local Manhattan with gene tracks: LD r² × sign(Z) — rs185665940**

<a href="MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_local_top_hits_with_gtf_rs185665940.png"><img src="MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_local_top_hits_with_gtf_rs185665940.png" alt="Local Manhattan with gene tracks: LD r² × sign(Z) — rs185665940" width="760"></a>

**Local Manhattan with gene tracks: LD r² × sign(Z) — rs62604261**

<a href="MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_local_top_hits_with_gtf_rs62604261.png"><img src="MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_local_top_hits_with_gtf_rs62604261.png" alt="Local Manhattan with gene tracks: LD r² × sign(Z) — rs62604261" width="760"></a>

**Forest plot: female**

<a href="MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_top_hits_forest_EUR_FEMALE.png"><img src="MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_top_hits_forest_EUR_FEMALE.png" alt="Forest plot: female" width="800"></a>

**Forest plot: male**

<a href="MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_top_hits_forest_EUR_MALE.png"><img src="MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_top_hits_forest_EUR_MALE.png" alt="Forest plot: male" width="800"></a>

</details>

#### gnuplot outputs

**Genome-wide Manhattan: chromosome-colored association P**

<a href="MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_manhattan.png"><img src="MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_manhattan.png" alt="Genome-wide Manhattan: chromosome-colored association P" width="1100"></a>

**Local Manhattan: Z-score colors**

<a href="MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_local_top_hits_manhattan.png"><img src="MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_local_top_hits_manhattan.png" alt="Local Manhattan: Z-score colors" width="1000"></a>

**Local Manhattan with gene tracks: LD r² × sign(Z) — rs2232429**

<a href="MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_local_top_hits_with_gtf_rs2232429.png"><img src="MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_local_top_hits_with_gtf_rs2232429.png" alt="Local Manhattan with gene tracks: LD r² × sign(Z) — rs2232429" width="760"></a>

<details><summary>View the remaining 8 gnuplot figures</summary>

**Local Manhattan: Z-score colors — rs185665940**

<a href="MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_local_top_hits_manhattan_rs185665940.png"><img src="MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_local_top_hits_manhattan_rs185665940.png" alt="Local Manhattan: Z-score colors — rs185665940" width="1000"></a>

**Local Manhattan: Z-score colors — rs2232429**

<a href="MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_local_top_hits_manhattan_rs2232429.png"><img src="MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_local_top_hits_manhattan_rs2232429.png" alt="Local Manhattan: Z-score colors — rs2232429" width="1000"></a>

**Local Manhattan: Z-score colors — rs62604261**

<a href="MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_local_top_hits_manhattan_rs62604261.png"><img src="MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_local_top_hits_manhattan_rs62604261.png" alt="Local Manhattan: Z-score colors — rs62604261" width="1000"></a>

**Local Manhattan with gene tracks: LD r² × sign(Z) — rs185665940**

<a href="MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_local_top_hits_with_gtf_rs185665940.png"><img src="MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_local_top_hits_with_gtf_rs185665940.png" alt="Local Manhattan with gene tracks: LD r² × sign(Z) — rs185665940" width="760"></a>

**Local Manhattan with gene tracks: LD r² × sign(Z) — rs62604261**

<a href="MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_local_top_hits_with_gtf_rs62604261.png"><img src="MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_local_top_hits_with_gtf_rs62604261.png" alt="Local Manhattan with gene tracks: LD r² × sign(Z) — rs62604261" width="760"></a>

**Forest plot: female**

<a href="MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_top_hits_forest_EUR_FEMALE.png"><img src="MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_top_hits_forest_EUR_FEMALE.png" alt="Forest plot: female" width="800"></a>

**Forest plot: male**

<a href="MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_top_hits_forest_EUR_MALE.png"><img src="MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_top_hits_forest_EUR_MALE.png" alt="Forest plot: male" width="800"></a>

**Forest plot: combined**

<a href="MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_top_hits_forest_combined.png"><img src="MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_top_hits_forest_combined.png" alt="Forest plot: combined" width="800"></a>

</details>

Clone the repository and open
[`examples/public-scz-sex/results.html`](MultiGWAS-Explorer/examples/public-scz-sex/results.html)
for the standalone gallery and its
[validation record](MultiGWAS-Explorer/examples/public-scz-sex/README.md).

See the [test results and limitations](MultiGWAS-Explorer/install/PUBLIC_GWAS_RESULTS.md)
and [platform installation checks](https://github.com/chengzhongshan/MultiGWAS-Explorer/actions/workflows/installation.yml).

## Documentation

- [Installation, usage, validation, and revision notes](docs/README.md)
- [Detailed pipeline reference](MultiGWAS-Explorer/README.md): input formats,
  configuration, filtering, plotting options, SAS utilities, and MCP integration
- [Benchmark documentation](MultiGWAS-Explorer/benchmark/README.md)

The main entry points are `auto_prepare_and_run_diff_gwas_with_gnuplot.pl`
for local gnuplot and `auto_prepare_and_run_diff_gwas.pl` for SAS ODA. Their
historical filenames are retained for compatibility.
