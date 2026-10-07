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

The following 22 figures come from the full 6,650,636-comparison Windows
portable-Cygwin validation of the sex-stratified schizophrenia GWAS. Inquiry
targets are `rs2232429` (common association), `rs185665940` (autosomal
differential signal), and `rs62604261` (chromosome-X differential signal).
Click any figure to open its full-resolution PNG.

#### SAS OnDemand for Academics outputs

| Genome-wide Manhattan | Combined local Manhattan |
| --- | --- |
| [![SAS ODA genome-wide Manhattan](MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_manhattan.png)](MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_manhattan.png) | [![SAS ODA combined local Manhattan](MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_local_top_hits_manhattan.png)](MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_local_top_hits_manhattan.png) |
| Combined signed-LD gene-track output | Gene-track export: panel 1 |
| [![SAS ODA combined signed LD gene track](MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_local_top_hits_with_gtf.png)](MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_local_top_hits_with_gtf.png) | [![SAS ODA gene track panel 1](MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_local_top_hits_with_gtf_part1.png)](MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_local_top_hits_with_gtf_part1.png) |
| Gene-track export: panel 2 | Gene-track export: panel 3 |
| [![SAS ODA gene track panel 2](MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_local_top_hits_with_gtf_part2.png)](MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_local_top_hits_with_gtf_part2.png) | [![SAS ODA gene track panel 3](MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_local_top_hits_with_gtf_part3.png)](MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_local_top_hits_with_gtf_part3.png) |
| Signed-LD gene track: `rs2232429` | Signed-LD gene track: `rs185665940` |
| [![SAS ODA signed LD gene track rs2232429](MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_local_top_hits_with_gtf_rs2232429.png)](MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_local_top_hits_with_gtf_rs2232429.png) | [![SAS ODA signed LD gene track rs185665940](MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_local_top_hits_with_gtf_rs185665940.png)](MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_local_top_hits_with_gtf_rs185665940.png) |
| Signed-LD gene track: `rs62604261` | Female forest plot |
| [![SAS ODA signed LD gene track rs62604261](MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_local_top_hits_with_gtf_rs62604261.png)](MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_local_top_hits_with_gtf_rs62604261.png) | [![SAS ODA female forest plot](MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_top_hits_forest_EUR_FEMALE.png)](MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_top_hits_forest_EUR_FEMALE.png) |
| Male forest plot | |
| [![SAS ODA male forest plot](MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_top_hits_forest_EUR_MALE.png)](MultiGWAS-Explorer/examples/public-scz-sex/figures/sas/PUBLIC_SCZ_EUR_SEX_SAS_top_hits_forest_EUR_MALE.png) | |

#### gnuplot outputs

| Genome-wide Manhattan | Combined signed-LD local Manhattan |
| --- | --- |
| [![gnuplot genome-wide Manhattan](MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_manhattan.png)](MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_manhattan.png) | [![gnuplot combined signed LD local Manhattan](MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_local_top_hits_manhattan.png)](MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_local_top_hits_manhattan.png) |
| Signed-LD local Manhattan: `rs2232429` | Signed-LD local Manhattan: `rs185665940` |
| [![gnuplot signed LD local Manhattan rs2232429](MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_local_top_hits_manhattan_rs2232429.png)](MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_local_top_hits_manhattan_rs2232429.png) | [![gnuplot signed LD local Manhattan rs185665940](MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_local_top_hits_manhattan_rs185665940.png)](MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_local_top_hits_manhattan_rs185665940.png) |
| Signed-LD local Manhattan: `rs62604261` | Signed-LD gene track: `rs2232429` |
| [![gnuplot signed LD local Manhattan rs62604261](MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_local_top_hits_manhattan_rs62604261.png)](MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_local_top_hits_manhattan_rs62604261.png) | [![gnuplot signed LD gene track rs2232429](MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_local_top_hits_with_gtf_rs2232429.png)](MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_local_top_hits_with_gtf_rs2232429.png) |
| Signed-LD gene track: `rs185665940` | Signed-LD gene track: `rs62604261` |
| [![gnuplot signed LD gene track rs185665940](MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_local_top_hits_with_gtf_rs185665940.png)](MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_local_top_hits_with_gtf_rs185665940.png) | [![gnuplot signed LD gene track rs62604261](MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_local_top_hits_with_gtf_rs62604261.png)](MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_local_top_hits_with_gtf_rs62604261.png) |
| Female forest plot | Male forest plot |
| [![gnuplot female forest plot](MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_top_hits_forest_EUR_FEMALE.png)](MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_top_hits_forest_EUR_FEMALE.png) | [![gnuplot male forest plot](MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_top_hits_forest_EUR_MALE.png)](MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_top_hits_forest_EUR_MALE.png) |
| Combined female/male forest plot | |
| [![gnuplot combined forest plot](MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_top_hits_forest_combined.png)](MultiGWAS-Explorer/examples/public-scz-sex/figures/gnuplot/PUBLIC_SCZ_EUR_SEX_GNUPLOT_top_hits_forest_combined.png) | |

Clone the repository and open
[`examples/public-scz-sex/results.html`](MultiGWAS-Explorer/examples/public-scz-sex/results.html)
for the standalone gallery. The
[example record](MultiGWAS-Explorer/examples/public-scz-sex/README.md).

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
