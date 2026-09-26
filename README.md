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
- Run from the command line or through the Perl MCP interface.

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

## Update your GitHub repository

From the outer clone directory (the one containing `.git`), preview and then
upload source changes:

```bash
./update_github_upload.sh --dry-run
./update_github_upload.sh "Describe the pipeline change"
```

The script stages tracked edits and new source or documentation files. It
excludes SAS ODA and gnuplot run results, even if they were staged earlier;
the local files remain available. The preview does not change the Git index
or contact GitHub. If GitHub has newer commits, the upload stops before
committing so you can integrate those changes first.

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
Its `results.html` page groups gnuplot and SAS figures in one place.

See the [test results and limitations](MultiGWAS-Explorer/install/PUBLIC_GWAS_RESULTS.md)
and [platform installation checks](https://github.com/chengzhongshan/MultiGWAS-Explorer/actions/workflows/installation.yml).

## Documentation

- [Installation, usage, validation, and revision notes](docs/README.md)
- [Detailed pipeline reference](MultiGWAS-Explorer/README.md): input formats,
  configuration, filtering, plotting options, SAS utilities, and MCP integration
- [Benchmark documentation](MultiGWAS-Explorer/benchmark/README.md)

The main entry points are `auto_prepare_and_run_diff_gwas_with_gunplot.pl`
for local gnuplot and `auto_prepare_and_run_diff_gwas.pl` for SAS ODA. Their
historical filenames are retained for compatibility.
