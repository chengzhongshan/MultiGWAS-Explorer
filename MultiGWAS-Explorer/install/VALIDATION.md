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

### Verified result

All six jobs passed at source commit
`f073931eb08ae019106a9595ab016a3a66b4740e`:
[run 34967910284](https://github.com/chengzhongshan/MultiGWAS-Explorer/actions/runs/34967910284).

| Environment | Installer/build and example |
| --- | --- |
| Ubuntu 24.04 x64 | Passed |
| Windows Server 2022 / portable Cygwin | Passed |
| macOS 15 Apple Silicon | Passed |
| macOS 15 Intel | Passed |
| Docker / Ubuntu 24.04 | Passed |
| Apptainer / Ubuntu 24.04 | Passed |

Intel macOS used the local gnuplot build without Qt. The Windows run selected
the repository-local cpanm rather than Strawberry Perl's copy. Later
documentation-only edits do not change the tested installer source.

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

## Local SAS ODA validation (2026-09-19)

The portable Cygwin installation was tested with Perl 5.44.0, PDL 2.106,
GD 2.91, Python 3.12, and Windows Java 8. Inline::Python resolves
`cygperl5_44.dll`. The preceding local public schizophrenia female/male GWAS
run validated 6,650,636 differential rows across autosomes and chromosome X.

The real ODA login probe (`proc setinit;run;`) and the public sex-GWAS SAS
plotting phase passed using saved local credentials. The plotting run disabled
both gnuplot fallback options and tested `rs2232429`, `rs185665940`, and
`rs62604261` with 500 kb local half-windows. Seven downloaded PNGs decoded
successfully: genome-wide Manhattan, local Manhattan, three local GTF panels,
and female/male forest panels. Generated outputs and credentials are not
committed. The local results directory contains `test_run_plots.json`,
`image_validation_sas.json`, and `results.html`.

Two configuration issues were corrected during this validation:

- The copied SASPy profile had an empty Java setting and absolute JAR paths
  pointing at the previous checkout. Regenerating the profile restored login;
  the Cygwin repair harness now regenerates it. The installation check rejects
  an empty Java setting or missing JAR before opening an ODA connection. Both
  rejection paths and the valid profile were tested.
- The SAS local-GTF wrapper omitted `--reference-build` when calling the
  annotation extractor, causing an hg38 status label for the explicit hg19
  run. It now passes the resolved build. Fresh extraction reported hg19 and
  returned 3,456 annotation rows from the GRCh37/lift37 source; its SHA-256
  matched the earlier subset, confirming the annotation data had not changed.

The initial six SAS job logs had no `ERROR:` or `FATAL` lines. Nonfatal warnings
remain for SQL table rewrites, gene-layout/label handling, and Arial font
substitution. This validation does not establish that every SAS option or
external LD service works; it covers the explicit-target plotting workflow
above. The dependency smoke test also passed after the profile-check changes.

## Forest gene-label regression (2026-09-19)

Visual review exposed a failure that PNG decoding alone did not detect:
explicit SNP requests bypassed GTF annotation, and the SAS macro treated the
gene strings as independent categorical coordinates. Three missing labels
therefore collapsed into one centered `NA`. The CSV helper now annotates
explicit requests, and a right-side YAXISTABLE uses the SNP row coordinates.
Forest CSV reuse is opt-in, preventing reruns from silently retaining old
missing labels. The per-panel separator variable is also initialized when
there is no hit-class boundary.

The new `install/test_requested_hit_genes.pl` passed all 12 checks covering
overlap, nearest gene, repeated gene names, chr23/chrX matching, unavailable
chromosomes, GTF provenance, and manual overrides. It runs in the installation
smoke test, which passed on Cygwin Perl 5.44. A real SAS render using the old
all-NA CSV confirmed visually that all three labels occupy their own SNP rows.

The full cached GENCODE v49lift37 lookup returned:

| SNP | hg19 position | Gene | Positional relationship |
| --- | --- | --- | --- |
| rs2232429 | chr6:28359632 | ZSCAN12 | Overlaps gene interval |
| rs185665940 | chr2:72269028 | CYP26B1 | Nearest selected gene, 87,339 bp away |
| rs62604261 | chrX:123643668 | TENM1 | Overlaps gene interval |

These mappings were cross-checked against the extracted hg19 GTF gene
intervals. Positional assignment does not establish causal gene involvement.

A forced end-to-end SAS forest rerun regenerated the previously all-NA CSV
without enabling reuse. Both female and male PNGs were inspected: all three
gene symbols are present, unclipped, and aligned with the correct SNP rows.
The final forest SAS log had no ERROR, FATAL, or WARNING diagnostics. The
earlier all-NA rendering was retained locally as a repeated-label regression
artifact; test images and public GWAS data are not committed.

## Gnuplot forest follow-up (2026-09-19)

The gnuplot panels still used a pre-fix CSV with empty gene fields. Unlike the
SAS macro, the gnuplot renderer already assigned right-axis labels to numeric
SNP row coordinates; its remaining problem was unconditional cache reuse.
The gnuplot wrapper now regenerates the CSV by default, with intentional reuse
available through `REUSE_FOREST_TOP_HITS_CSV=1`. Regenerated CSVs trigger panel
and combined-image rendering explicitly, even if file timestamps coincide.

A public sex-GWAS forest rerun without `--force` refreshed the old CSV and both
panels. The combined image was inspected: ZSCAN12, CYP26B1, and TENM1 are
visible and aligned with the correct SNPs in both sexes. The extended
`install/test_forest_text.pl` passed eight checks, including repeated gene
names at distinct row coordinates and nonblank text in all three right-margin
row regions.

## Limits

The CI smoke checks do not authenticate to SAS OnDemand, upload data, validate full
GWAS analyses, or exercise a site's HPC scheduler. Apptainer CI builds with
root on a disposable runner; an HPC site's unprivileged/fakeroot policy may
require different launch instructions. Other Linux distributions and macOS
versions require their own validation.
