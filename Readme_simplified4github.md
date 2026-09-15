# MultiGWAS-Explorer

MultiGWAS-Explorer is a cross-platform GWAS plotting pipeline for differential
and common-association analyses. It combines local preprocessing scripts,
SAS OnDemand for Academics plotting wrappers, and an alternative gunplot / PDL
backend for genome-wide Manhattan plots, local Manhattan plots, local GTF
gene-track plots, and forest plots.
<img width="2400" height="2510" alt="Figure1_pipeline_overview_Zhongshan" src="https://github.com/user-attachments/assets/a56b8675-2a48-41c9-ba2c-1d1e793603ee" />
For the full project guide, advanced troubleshooting, and validation notes,
see [MultiGWAS-Explorer/README.md](MultiGWAS-Explorer/README.md).
The [main README](README.md) is the canonical installation guide.

## Main Scripts

- `auto_prepare_and_run_diff_gwas.pl`
  Main automation entry point for the SAS ODA workflow.
- `auto_prepare_and_run_diff_gwas_with_gunplot.pl`
  Main automation entry point for the non-SAS gunplot workflow.
- `run_sas_codes_or_script_in_ODA.pl`
  Low-level helper for SAS ODA submit, upload, download, delete, and session
  reuse.

## Installation

Clone the repository and enter its nested pipeline directory before running
any of the installation or container commands below:

```bash
git clone https://github.com/chengzhongshan/MultiGWAS-Explorer.git
cd MultiGWAS-Explorer/MultiGWAS-Explorer
```

### Windows

Recommended for portable Cygwin:

Install 7-Zip (`7z` on PATH) and a Java JDK first. Set `JAVA_HOME` to the JDK
directory if `java.exe` is not on PATH.

```powershell
powershell -NoProfile -ExecutionPolicy Bypass `
  -File .\install\install_windows_portable_cygwin.ps1
```

If you are already inside a supported Cygwin shell:

```bash
bash install/install_cygwin.sh
```

### Linux

Ubuntu:

```bash
sudo bash install/install_ubuntu.sh
```

### macOS

```bash
bash install/install_macos.sh
```

### Post-install check

Run this on any host install:

```bash
bash install/check_pipeline_install.sh
```

## Containers

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

### 2. Quick gunplot validation

For a self-contained test without SAS ODA credentials or external GWAS data:

```bash
bash install/run_plotting_example.sh
```

Open `example-output/example.png`. The PGC examples below require the external
data and reference files; edit the bundled spec's input, output, work-directory,
and PLINK/reference paths to match your machine first:

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
  management commands, see [MultiGWAS-Explorer/README.md](MultiGWAS-Explorer/README.md).
