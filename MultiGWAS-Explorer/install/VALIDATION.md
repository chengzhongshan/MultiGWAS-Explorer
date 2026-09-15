# Installation validation

The [installation workflow](https://github.com/chengzhongshan/MultiGWAS-Explorer/actions/workflows/installation.yml)
installs dependencies and renders a synthetic plot on Ubuntu 24.04, Windows
Server 2022/Cygwin, macOS 15 ARM, macOS 15 Intel, Docker, and Apptainer.
Inspect the results for the commit you intend to install. A workflow definition
alone is not evidence that a platform passed.

## Checks performed

- `bash install/test_install_environment.sh`: first-install Perl library
  visibility, repeatable environment activation, foreign cpanm exclusion, and
  Java paths containing spaces.
- Platform installer: system packages, repository-local Python and Perl
  dependencies, SASPy configuration, and dependency smoke check.
- `bash install/check_pipeline_install.sh`: required modules, Java startup,
  entry-point syntax, synthetic LD rendering, and plot-contract tests. This can
  be invoked by absolute path from another directory.
- `bash install/run_plotting_example.sh`: synthetic GWAS data and a persistent
  PNG under `example-output/`, with no user GWAS data or SAS credentials.

CI uploads the installer log and the example output where available. Container
jobs upload build logs. Repository-local dependency caches are not restored in
CI; the jobs exercise installation against each runner's supplied base image.

## Initial findings and fixes (September 2026)

Default CPAN mirror downloads failed repeatedly on Ubuntu. The installer now
selects an HTTPS mirror explicitly (`PIPELINE_CPAN_MIRROR` overrides it).
PDL configuration could not see newly installed prerequisites. The installer
now creates and activates the local Perl library directory before installation
and uses `--local-lib` for CPAN builds.

The first CI run of the fix branch passed Ubuntu, macOS ARM, and Docker:
[run 34919847861](https://github.com/chengzhongshan/MultiGWAS-Explorer/actions/runs/34919847861).
That run also exposed a foreign Strawberry Perl cpanm on Windows and a
90-minute Qt dependency build on Intel macOS. The follow-up uses a standalone
repository-local cpanm and builds PNG-capable gnuplot without Qt when needed.

Local testing uses Windows 10 x64 with portable Cygwin and Ubuntu 24.04 under
WSL1. Those local machines already have system software; they are not fresh
operating-system images. The Ubuntu test uses fresh repository-local Perl and
Python dependencies. Windows bootstrap and smoke checks passed locally.

## Limits

These checks do not authenticate to SAS OnDemand, upload data, validate full
GWAS analyses, or exercise a site's HPC scheduler. Apptainer CI builds with
root on a disposable runner; an HPC site's unprivileged/fakeroot policy may
require different launch instructions. Other Linux distributions and macOS
versions require their own validation.
