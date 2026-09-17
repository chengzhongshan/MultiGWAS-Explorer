# Installation, usage, and revision notes

MultiGWAS-Explorer is a cross-platform GWAS plotting pipeline for differential
and common-association analyses. It combines local preprocessing scripts,
SAS OnDemand for Academics plotting wrappers, and an alternative gunplot / PDL
backend for genome-wide Manhattan plots, local Manhattan plots, local GTF
gene-track plots, and forest plots.
<img width="2400" height="2510" alt="Figure1_pipeline_overview_Zhongshan" src="https://github.com/user-attachments/assets/a56b8675-2a48-41c9-ba2c-1d1e793603ee" />
For the full project guide, advanced troubleshooting, and validation notes,
see [MultiGWAS-Explorer/README.md](../MultiGWAS-Explorer/README.md).

## Reviewer-driven scientific revision

The 2026 reviewer revision adds conservative allele harmonization, raw-P-value
inference, MAF-aware lead filtering, and greedy HaploReg LD clumping as the
default top-hit method. Physical distance is retained only as an explicit
legacy method or a labeled fallback when an LD query is unresolved. The
repository also contains executable GWAMA 2.2.2 and EasyStrata 8.6 comparison
commands, full-source QC scripts, LD audit tables, and the deterministic MCP
interface evaluation.

```bash
cd MultiGWAS-Explorer
bash benchmark/run_differential_gwas_revision_benchmarks.sh show-inputs
bash benchmark/run_differential_gwas_revision_benchmarks.sh validate-results
```

See [benchmark/README.md](../MultiGWAS-Explorer/benchmark/README.md). The interface
timing used exact JSON-RPC tool calls—not free-form prompts—and its complete
request records are published under `benchmark/agent_interface/`.

## Main Scripts

- `auto_prepare_and_run_diff_gwas.pl`
  Main automation entry point for the SAS ODA workflow.
- `auto_prepare_and_run_diff_gwas_with_gunplot.pl`
  Main automation entry point for the non-SAS gunplot workflow.
- `run_sas_codes_or_script_in_ODA.pl`
  Low-level helper for SAS ODA submit, upload, download, delete, and session
  reuse. Repeated upload/download arguments are transferred through one SASPy
  session.

## Installation

[Installation CI](https://github.com/chengzhongshan/MultiGWAS-Explorer/actions/workflows/installation.yml)
tests Ubuntu, Windows portable Cygwin, macOS ARM/Intel, Docker, and Apptainer. Check the
individual job results for the revision you are installing; SAS ODA login and
remote plotting are separate integration checks.
See [validation scope and commands](../MultiGWAS-Explorer/install/VALIDATION.md).

For an end-to-end test using public female and male schizophrenia GWAS,
including chromosome X and optional SAS ODA plotting, see the
[Perl real-data test guide](../MultiGWAS-Explorer/install/PUBLIC_GWAS_TEST.md).

First clone the repository and enter the pipeline directory. All installation
and container commands below start in this directory:

```bash
git clone https://github.com/chengzhongshan/MultiGWAS-Explorer.git
cd MultiGWAS-Explorer/MultiGWAS-Explorer
```

If you downloaded a ZIP, enter its `MultiGWAS-Explorer` subdirectory instead.
Do not repeat the directory change before each command below.

### Windows

Recommended for portable Cygwin:

Install 7-Zip (`7z` on PATH) and a Java JDK first. Set `JAVA_HOME` to the JDK
directory if `java.exe` is not on PATH. The bootstrap installs Cygwin packages;
it does not install a Windows JDK.

```powershell
powershell -NoProfile -ExecutionPolicy Bypass `
  -File .\install\install_windows_portable_cygwin.ps1
```

The portable bootstrap defaults to
`%USERPROFILE%\CygwinPortablePipeline`. If Cygwin `curl` fails with a
self-signed certificate chain while bootstrapping repo-local dependencies, rerun
with the explicit opt-in:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass `
  -File .\install\install_windows_portable_cygwin.ps1 `
  -AllowInsecureDownloads
```

The pipeline also ships Windows `bgzip.exe` / `tabix.exe` under
`DiffGWASDeps/` and places that directory on the portable-Cygwin runtime path.
The installer can still download htslib 1.20 into `tools/` and build newer
repo-local copies when needed.

After installation, open the portable shell with:

```text
C:\Users\<username>\CygwinPortablePipeline\CygwinPortable.exe
```

If needed, start it directly from PowerShell:

```powershell
& "$env:USERPROFILE\CygwinPortablePipeline\App\Runtime\Cygwin\bin\mintty.exe" -
```

Inside Cygwin, change into the project by using `/mnt/c/...` paths, then run
the smoke test or pipeline:

```bash
cd /mnt/c/Users/<username>/Desktop/MultiGWAS-Explorer-main/MultiGWAS-Explorer
bash install/check_pipeline_install.sh
perl ./auto_prepare_and_run_diff_gwas_with_gunplot.pl \
  --spec configs/spec_pgc_scz_sex_common_automation.json \
  --plots local_manhattan,local_gtf \
  --target-snps rs185665940
```

If you are already inside a supported Cygwin shell:

```bash
bash install/install_cygwin.sh
```

### Linux

Recommended hosts are supported Ubuntu LTS/current releases, Debian-like Linux
systems with equivalent packages, or the Docker/Singularity paths below for old
or locked-down machines. On Ubuntu, run:

```bash
sudo bash install/install_ubuntu.sh
bash install/check_pipeline_install.sh
```

For non-Ubuntu Linux, install equivalent system packages first, then run the
repo-local phase with `PIPELINE_SKIP_APT=1 bash install/install_ubuntu.sh`.
See [MultiGWAS-Explorer/README.md](../MultiGWAS-Explorer/README.md) for package
details and legacy Ubuntu troubleshooting.

Perl packages are installed into the repository using the HTTPS mirror
`https://cpan.metacpan.org`. To use an institutional mirror, set
`PIPELINE_CPAN_MIRROR` when running the installer. Existing system Perl modules
can satisfy dependencies; newly installed modules stay under `local/perl5-<platform>`.

### macOS

The installer supports Intel and Apple Silicon and provisions OpenJDK. Apple
Silicon uses Homebrew; Intel macOS 15 uses a pinned, checksum-verified MacPorts
installer because current Homebrew releases no longer support Intel macOS.
Set `PIPELINE_MACOS_PACKAGE_MANAGER=homebrew` only to override this selection
for an Intel machine with a separately maintained Homebrew installation.

```bash
bash install/install_macos.sh
```

### Post-install check

Run this on any host install:

```bash
bash install/check_pipeline_install.sh
```

## Containers

Run these commands from the pipeline directory containing `Dockerfile` and
`install/`, not the outer Git checkout directory.

### Docker

Build:

```bash
docker build -t multigwas-explorer-pipeline:latest .
```

Smoke test:

```bash
docker run --rm -it multigwas-explorer-pipeline:latest \
  bash -lc "cd /opt/MultiGWAS-Explorer && bash install/check_pipeline_install.sh"
```

Interactive container with mounted GWAS data and SAS ODA authinfo:

```bash
docker run --rm -it \
  -e PIPELINE_WORKDIR=/opt/MultiGWAS-Explorer \
  -v /path/to/_authinfo:/root/_authinfo:ro \
  -v /path/to/gwas_drive_e:/mnt/e \
  -v /path/to/gwas_drive_g:/mnt/g \
  multigwas-explorer-pipeline:latest bash
```

### Singularity / Apptainer

Build and test:

```bash
bash install/singularity/build_apptainer_image.sh
apptainer exec MultiGWAS-Explorer_pipeline.sif \
  bash -lc "cd /opt/MultiGWAS-Explorer && bash install/check_pipeline_install.sh"
```

## Quick Start

### 1. Check SAS ODA login

Useful before running any SAS-based plotting step:

```bash
perl ./run_sas_codes_or_script_in_ODA.pl --check-sas-oda-login-only
```

For a tiny direct submit test, prefer the repo-local runtime and a short
timeout:

```bash
. install/common.sh
activate_perl_env
activate_python_env
SAS_ODA_RUN_TIMEOUT_SECONDS=90 \
./run_sas_codes_or_script_in_ODA.pl --code "proc print data=sashelp.class;run;"
```

If that test prints repeated `Waiting for SAS ODA session server response while
reading response header...` messages, the SAS code usually has not started yet;
the local SASPy Java/IOM bridge is still creating or answering through an ODA
session. Stop stale bridge processes with
`./run_sas_codes_or_script_in_ODA.pl --kill-saspy-sessions`, then rerun from
the repo-local environment. More details are in
[MultiGWAS-Explorer/README.md](../MultiGWAS-Explorer/README.md).

Persistent SAS ODA submits are bounded by
`SAS_ODA_SESSION_SUBMIT_TIMEOUT_SECONDS` so a wedged macro bootstrap returns a
clear timeout instead of hanging forever. When global macro autoload is enabled,
remote macro calls are resolved by the normal `~/Macros/importallmacros_ue`
bootstrap rather than by injecting individual targeted loaders for submacros.

On Linux, result HTML opens through a real browser binary such as Chrome or
Firefox before falling back to `xdg-open`. Override the browser with
`OPEN_RESULT_BROWSER=google-chrome-stable`, override a remote desktop display
with `OPEN_RESULT_DISPLAY=:20`, or set `OPEN_RESULT=0` to only print the saved
HTML path.

### 2. Self-contained local plotting example

After installation, create a synthetic local Manhattan plot without downloading
GWAS data or using SAS ODA credentials:

```bash
bash install/run_plotting_example.sh
```

Open `example-output/example.png`. Supply an output directory as the first
argument to save the example elsewhere. The dependency smoke test also checks
synthetic LD heatmap rendering. It checks that Java starts, but does not log in
to SAS ODA.

### PGC data example

The following example requires your PGC summary statistics and reference files.
Edit `input_dir`, `output_dir`, `workdir`, and PLINK/reference paths in the spec
for your machine before running it; the bundled spec contains author-specific
paths and is not a fresh-install test:

```bash
perl ./auto_prepare_and_run_diff_gwas_with_gunplot.pl \
  --spec configs/spec_pgc_scz_sex_common_automation.json \
  --plots local_manhattan,local_gtf \
  --target-snps rs185665940
```

To include the slower genome-wide Manhattan panel:

```bash
perl ./auto_prepare_and_run_diff_gwas_with_gunplot.pl \
  --spec configs/spec_pgc_scz_sex_common_automation.json \
  --plots manhattan,local_manhattan,local_gtf \
  --target-snps rs185665940
```

### 3. Quick SAS ODA run

Example SAS ODA workflow using the bundled schizophrenia spec:

```bash
perl ./auto_prepare_and_run_diff_gwas.pl \
  --spec configs/spec_pgc_scz_sex_common_automation.json \
  --step plot_manhattan \
  --step plot_local_manhattan \
  --step plot_local_gtf \
  --target-snps rs185665940
```

Explicit target-SNP local-GTF runs create their requested-target CSV
before indexed GTF extraction, use bundled tabix/bgzip by default on Windows,
and prepare all known ODA uploads as one manifest. Existing local-GTF output is
reused only when its request key matches the current target/configuration, so a
new target no longer requires `--force` merely to avoid an unrelated result.

### 4. Forest plot examples

Gunplot:

```bash
perl ./auto_prepare_and_run_diff_gwas_with_gunplot.pl \
  --spec configs/spec_pgc_scz_sex_common_automation.json \
  --plots forest \
  --target-snps rs185665940
```

SAS ODA:

```bash
perl ./auto_prepare_and_run_diff_gwas.pl \
  --spec configs/spec_pgc_scz_sex_common_automation.json \
  --step plot_forest \
  --target-snps rs185665940
```

## Example Configs

- `configs/spec_pgc_scz_sex_common_automation.json`
  Example sex-stratified schizophrenia workflow.
- `configs/spec_pgc_scz_ancestry_diff_automation.json`
  Example ancestry-differential schizophrenia workflow.

## Notes

- Use the gunplot workflow first if you want a fast end-to-end validation.
- Use the SAS ODA workflow when you want the SAS-rendered plot outputs.
- On container runs, keep `PIPELINE_WORKDIR=/opt/MultiGWAS-Explorer` so the
  wrappers use the Linux-installed environment inside the image.
- For detailed options, top-hit filtering behavior, troubleshooting, and file
  management commands, see [MultiGWAS-Explorer/README.md](../MultiGWAS-Explorer/README.md).
