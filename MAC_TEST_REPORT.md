# Mac installation and schizophrenia example test — 2026-09-20

Repository: https://github.com/chengzhongshan/MultiGWAS-Explorer  
Tested branch: `main`, including the macOS, PLINK2, and SAS ODA fixes recorded below
Host: Apple Silicon (arm64), macOS 26.6.2, 18 GiB RAM

Status: **PASS** for a fresh native macOS installation and the complete public
PGC schizophrenia female-versus-male GWAS example with both gnuplot and SAS
ODA backends.

## Installation

- Installed and activated repository-local Perl and Python environments.
- Used Homebrew Perl 5.44.0 to avoid a reproducible PDL crash under the macOS
  system Perl 5.34.
- Used `/usr/bin/curl` while bootstrapping Homebrew to avoid an incompatible
  Anaconda TLS certificate configuration.
- `install/check_pipeline_install.sh`: PASS, including native bgzip/tabix,
  PNG rendering, PLINK alias checks, and local-GTF regressions.
- Synthetic plotting, installation-environment, and coordinate-sort
  regressions: PASS. The coordinate-sort regression crosses multiple BGZF
  members and validates target lookup after the first member.

## Public schizophrenia workflow

- All four PGC source downloads passed their pinned checksums.
- Preprocessing completed for autosomes 1–22 and chromosome X.
- Independent numerical validation: PASS for all 6,650,636 standardized rows,
  including beta, standard error, Z-score, and P-value calculations.
- The final wide plotting subset contains 1,251,499 rows.
- Target loci: `rs2232429` (ZSCAN12), `rs185665940` (CYP26B1), and
  `rs62604261` (TENM1).
- Genome-wide Manhattan, local Manhattan, local GTF, and forest plot families
  all rendered and decoded successfully with both plotting backends.

## SAS ODA validation

- SAS ODA authentication and remote execution passed from this Mac.
- The genome-wide Manhattan plot, local Manhattan plot, three target-specific
  local GTF plots, and two forest panels downloaded and decoded successfully.
- Each target-specific GTF panel contains adjacent genes and its own signed
  `r2 * sign(Z)` colorbar from `-1` to `1`, referenced to the SNP named in the
  panel: `rs2232429`, `rs185665940`, or `rs62604261`.
- The gallery manifest contains seven validated SAS figures and the generated
  gallery displays them alongside the gnuplot results.

## Phase 3 LD validation

- PLINK v2.0.0-a.6.39 M1 (19 Sep 2026) ran natively on Apple Silicon.
- Prepared the official phased 1000 Genomes Phase 3 GRCh37 reference (84.8
  million variants; 9.8 GiB local cache).
- PLINK2 restricted LD calculations to 503 EUR samples.
- All three local-GTF manifests passed direct-source and provenance checks for
  phased EUR GRCh37 LD and signed `r2 * sign(Z)` heatmap coloring.
- Rendered LD point counts were 689 for `rs2232429`, 9 for `rs185665940`, and
  22 for `rs62604261`.

## Fixes made during the Mac test

- Hardened macOS curl/Perl installation and Perl ABI selection (`fe3c5ca`).
- Fixed plotting extractors to read every member of BGZF files and added a
  multi-member regression (`4fe1085`).
- Added cross-platform auto-detection of `cache/plink2_bin/plink2` on
  macOS/Linux and `plink2.exe` on Windows; updated setup documentation.
- Replaced Bash 4-only case conversion and `mapfile` use in SAS wrappers so
  they run under the Bash 3.2 supplied with macOS; added a portable lock when
  the Linux `flock` command is unavailable.
- Changed multi-target SAS heatmaps to build and use one PLINK2 LD cache per
  target. Each target now runs through the single-locus GTF path, preserving
  both the target-specific signed-LD scale and the adjacent gene track.

## Results

- Gallery: `MultiGWAS-Explorer/public-gwas-test/results.html`
- Numerical report: `MultiGWAS-Explorer/public-gwas-test/numeric_validation.json`
- Image report: `MultiGWAS-Explorer/public-gwas-test/image_validation_gunplot.json`
- SAS image report: `MultiGWAS-Explorer/public-gwas-test/image_validation_sas.json`
- Figures and signed-LD manifests:
  `MultiGWAS-Explorer/public-gwas-test/PUBLIC_SCZ_EUR_SEX_GUNPLOT_*`

The SAS ODA and native gnuplot workflows both completed without fallback.
