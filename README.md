# MultiGWAS-Explorer
MultiGWAS-Explorer: an AI-assisted pipeline for comparative analysis and visualization of shared and differential loci across multiple GWASs
<img width="2400" height="2510" alt="Figure1_pipeline_overview_Zhongshan" src="https://github.com/user-attachments/assets/a56b8675-2a48-41c9-ba2c-1d1e793603ee" />

# MultiGWAS-Explorer Workflow

This repository combines a Perl MCP server, GWAS preprocessing utilities, and
SAS OnDemand for Academics plotting wrappers into one workflow for differential
GWAS analysis. It is designed to take one or more GWAS summary-statistics
tables, normalize them into a common layout, compute pairwise differential
effects, and generate shareable genome-wide and local visualizations.

## What This Repository Does

- Merges raw GWAS summary-statistics files into a long comparison table.
- Sorts and prepares coordinate-aware outputs for downstream extraction.
- Computes pairwise differential effect sizes and standardized Z/P signals.
- Builds wide-format GWAS subsets for plotting.
- Produces genome-wide Manhattan plots, local Manhattan plots, and local GTF
  gene-track plots through SAS ODA.
- Produces the same genome-wide Manhattan, local Manhattan, and local GTF plot
  families through the alternative gunplot / PDL backend.
- Produces top-hit forest plots through both the SAS ODA and gunplot paths for
  either one inquiry SNP or multiple common / differential top hits.
- Lets users render either the default multi-GWAS comparison view or a custom
  displayed GWAS subset, including a single selected GWAS track.
- Lets users render those same plot families for explicit inquiry SNPs through
  `--target-snps`.
- Supports both direct command-line use and MCP-driven orchestration.

## Main Entry Points

- `auto_prepare_and_run_diff_gwas.pl`
  High-level automation entry point for config generation, preprocessing, and
  plot submission. It now also accepts `--gwas-dir` for auto-detecting either
  raw multi-file GWAS inputs or a single merged-wide GWAS table and can emit
  an inferred spec through `--preview-spec` / `--generate-spec-only`.
- `auto_prepare_and_run_diff_gwas_with_gunplot.pl`
  Parallel non-SAS visualization pipeline for genome-wide Manhattan, combined
  local Manhattan, local GTF-style plots, and forest plots rendered through
  gnuplot / PDL, with the same displayed-GWAS and inquiry-SNP flexibility as
  the SAS ODA path.
- `server.pl`
  Local Perl MCP server that exposes reusable tools, including the automated
  differential-GWAS runner.
- `DiffGWASDeps/`
  Helper Perl, Bash, and SAS scripts used by the automation layer.

## Merged-Wide GWAS Support

The automation entrypoint now recognizes a merged-wide GWAS table as its own
source mode:

- `source_mode=merged_gwas_table`
- one shared locus table with columns such as `CHR`, `BP`, `SNP`, `A1`, `A2`
- one association block per cohort or GWAS, such as:
  - `BETA_DS_ALL`, `SE_DS_ALL`, `P_DS_ALL`
  - `BETA_MP2PRT`, `SE_MP2PRT`, `P_MP2PRT`
- optional extra association tracks such as meta-analysis `P` / `Z` columns

The merged-table path normalizes that input into the same plotting-wide schema
used elsewhere in the repository, then drives:

- genome-wide Manhattan plots
- local Manhattan panels
- local GTF plots
- top-hit forest plots

The auto-detector also skips generated artifacts such as prior
`*.merged_plotwide.tsv.gz`, `gunplot`, `png`, and `html` outputs so a rerun
against a directory like `AOA_GWAS_Data/` keeps pointing back to the original
merged input table instead of recycling derived files.

## Top-Hit MAF Safeguard

Top common or differential SNP selection now applies a minor-allele-frequency
safeguard before the local Manhattan and local GTF plotting stages:

- default rule: keep only top-hit candidates with `MAF > 0.01`
- first choice: use the matching GWAS frequency columns already carried through
  the standardized and wide differential tables, such as
  `*_GROUP1_FRQ_A`, `*_GROUP1_FRQ_U`, `*_GROUP2_FRQ_A`, and
  `*_GROUP2_FRQ_U`
- fallback: if those GWAS frequencies are absent, use an optional local
  gnomAD lookup table provided through:
  - spec key `gnomad_freq_file`
  - spec key `gnomad_population_map`
- conservative edge case: if both GWAS frequencies and a configured gnomAD
  lookup are unavailable for a candidate, the selector marks it
  `maf_source=UNKNOWN` and filters it out at the default safeguard threshold
- override: change the cutoff through spec key `top_hit_maf_threshold`

The SAS ODA local-top-hit wrappers now generate the requested top-hit CSV
locally with `DiffGWASDeps/generate_requested_top_hits_csv.pl` before the SAS
submit, so the MAF-aware locus list is reused consistently by both the SAS ODA
and gunplot paths instead of being recomputed differently inside SAS.

Those requested top-hit CSVs now preserve the MAF provenance used for each
retained locus:

- `selected_maf`
- `maf_source`
- `gwas_group1_maf`
- `gwas_group2_maf`
- `gwas_pair_maf_min`
- `gnomad_maf`
- `gnomad_pops`
- `maf_filter_decision`
- `maf_filter_reason`

Important refresh note:

- older wide plotting subsets created before this change may not contain the
  `*_FRQ_A` / `*_FRQ_U` columns yet
- older auto-generated runner configs may also predate the new pair metadata
  used by the common-association MAF filter
- if so, rerun at least:

```bash
perl ./auto_prepare_and_run_diff_gwas.pl \
  --spec configs/spec_pgc_scz_sex_common_automation.json \
  --step extract_wide_subset \
  --force
```

Current schizophrenia validation note:

- on 2026-06-12, the real PGC sex-differential schizophrenia validation
  retained lead SNP `rs185665940` with GWAS-derived `selected_maf = 0.013`
  (`gwas_group1_maf = 0.0135`, `gwas_group2_maf = 0.013`)
- the default differential top-hit selector now uses a two-step ladder of
  `1e-6` then `1e-5`; if no MAF-passing loci survive `1e-6`, the selector
  automatically retries at `1e-5` unless the caller already provided an
  explicit multi-threshold ladder
- on the same validation run, the real common-association selection retained
  `15` loci and every retained locus had `selected_maf > 0.01`
- because the current PGC runner config leaves `TOP_HIT_GNOMAD_FREQ_FILE`
  empty, candidates without usable GWAS frequency fields are conservatively
  filtered instead of being kept without an MAF check
- the rare-variant rejection branch and the gnomAD-fallback branch are covered
  by the regression harness:

```bash
perl DiffGWASDeps/test_top_hit_maf_filter.pl --no-real --keep-workdir
```

For a full real-data rerun, use the same helper without `--no-real`, or run the
two real checks separately when you want clearer timing and logs:

```bash
perl DiffGWASDeps/generate_requested_top_hits_csv.pl \
  --input /mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects.stdized.wide_beta_se_p_p_lt_0p05.final.tsv.gz \
  --output tmp_maf_validation_20260612/real_pgc_differential_top_hits.csv \
  --runner-config configs/auto_PGC_SCZ_female_vs_male_diff_effects_runner.json \
  --maf-threshold 0.01

perl DiffGWASDeps/generate_requested_top_hits_csv.pl \
  --input /mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects.stdized.wide_beta_se_p_p_lt_0p05.final.tsv.gz \
  --output tmp_maf_validation_20260612/real_pgc_common_top_hits.csv \
  --runner-config configs/auto_PGC_SCZ_female_vs_male_diff_effects_runner.json \
  --top-hit-mode common_association \
  --top-hit-signal-thrshds "5e-8 1e-6 1e-5" \
  --maf-threshold 0.01
```

## Forest Plot Support

Top-hit forest plots are now supported in both rendering backends:

- SAS ODA forest plots are driven by the existing
  `DiffGWASDeps/beta2OR_forest_plot.sas` macro.
- gunplot forest plots are rendered through
  `DiffGWASDeps/gunplot/pdl_gunplot_forest.pl` and are intentionally styled to
  stay close to the SAS ODA forest output.

Current behavior:

- `--plots forest` in the gunplot wrapper, or `--step plot_forest` in the SAS
  ODA wrapper, renders a forest plot artifact family instead of Manhattan /
  local-GTF families.
- For a single SNP, the forest plot uses cohort names on the y-axis and
  `OR and 95% CI` on the x-axis.
- For multiple SNPs, the forest output uses one panel per displayed GWAS /
  cohort track, puts SNP IDs on the left y-axis, puts adjacent gene labels on
  the right y-axis, separates differential and common hits with a horizontal
  dashed divider, and adds a significance star for genome-wide significant
  points (`P < 5e-8`).
- The same top-hit CSV generation and MAF-aware requested-hit selection are
  reused by both the SAS ODA and gunplot forest paths.

Example gunplot single-SNP forest rerun:

```bash
perl ./auto_prepare_and_run_diff_gwas_with_gunplot.pl \
  --spec configs/spec_pgc_scz_sex_common_automation.json \
  --plots forest \
  --target-snps rs185665940
```

Example gunplot multi-SNP forest rerun:

```bash
perl ./auto_prepare_and_run_diff_gwas_with_gunplot.pl \
  --spec configs/spec_pgc_scz_sex_common_automation.json \
  --plots forest \
  --target-snps rs185665940,rs4950119
```

Example SAS ODA forest rerun:

```bash
perl ./auto_prepare_and_run_diff_gwas.pl \
  --spec ./configs/spec_pgc_scz_sex_common_automation.json \
  --step plot_forest \
  --target-snps rs185665940
```

## Supplementary Table Regeneration

The manuscript-table helper now treats the supplementary common-hit and
differential-hit tables as full-strata exports by default.

Canonical outputs:

- `manuscript_assets/tables/Table_S1_all_common_association_loci.csv`
- `manuscript_assets/tables/Table_S2_differential_loci.csv`

Those standard filenames now include:

- association `P`, `BETA`, and `SE` for all pooled and ancestry-specific GWAS
  strata currently carried by the standardized wide table
- pairwise differential `P`, `BETA`, and `SE` for `ALL`, `EUR`, and `ASN`
- standardized differential `P` columns such as
  `ALL_STD_DIFF_P`, `EUR_STD_DIFF_P`, and `ASN_STD_DIFF_P`
- retained MAF QC fields such as `selected_maf` and `maf_source`

For compatibility with downstream manuscript editing, the same regeneration
step also writes mirrored convenience copies:

- `manuscript_assets/tables/Table_S1_all_common_association_loci_full_strata.csv`
- `manuscript_assets/tables/Table_S2_differential_loci_full_strata.csv`

The main-text tables remain intentionally narrower:

- `manuscript_assets/tables/Table_2_representative_common_loci.csv`
- `manuscript_assets/tables/Table_1_top_differential_locus.csv`

Recommended rerun:

```bash
perl DiffGWASDeps/regenerate_manuscript_hit_tables.pl \
  --config configs/spec_pgc_scz_sex_common_automation.json \
  --wide /mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects.stdized.wide_beta_se_p_p_lt_0p05.final.tsv.gz \
  --common-loci tmp_common_verify_postfix_1e6.tsv \
  --gtf cache/gtf/gencode.v49lift37.annotation.gtf.gz \
  --output-dir manuscript_assets/tables
```

Legacy-output recovery helpers:

- `DiffGWASDeps/augment_common_hits_table_s1.pl`
  Backfills an older compact `Table_S1_all_common_association_loci.csv` into a
  richer full-strata CSV when the original regeneration step was not rerun.
- `DiffGWASDeps/export_augmented_table_s1_excel.ps1`
  Writes a Windows `.xlsx` workbook from a regenerated or backfilled full
  strata `Table_S1` CSV when a manually editable Excel version is needed.

## One-Command Installation

The repository now ships platform-specific bootstrap scripts plus shared Perl
and Python dependency manifests:

- `install/install_windows_portable_cygwin.ps1`
- `install/install_windows_portable_cygwin.cmd`
- `install/install_cygwin.sh`
- `install/install_ubuntu.sh`
- `install/test_ubuntu_docker_gnuplot.sh`
- `install/install_macos.sh`
- `Dockerfile`
- `install/singularity/MultiGWAS-Explorer_pipeline.def`
- `install/singularity/build_apptainer_image.sh`
- `cpanfile`
- `install/requirements-pipeline.txt`

These installers set up a repo-local runtime instead of depending on a user's
personal `PERL5LIB` or a globally preconfigured Python:

- Perl modules are installed under platform-specific repo-local trees such as
  `local/perl5-cygwin/`, `local/perl5-linux/`, or `local/perl5-darwin/`
- Python packages such as `saspy` and `Pillow` are installed under
  `.venv-pipeline/`
- `bgzip` / `tabix` are taken from system packages when available, or built
  into `local/bin/` through `install/build_local_htslib.sh`

Recommended entry points:

Windows portable Cygwin bootstrap from PowerShell:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass `
  -File .\install\install_windows_portable_cygwin.ps1
```

This Windows bootstrap is built around
`MachinaCore/CygwinPortable`. By default it creates an isolated portable root
under `H:\TMP4SAS\CygwinPortablePipeline`, refreshes the required Cygwin
packages there, and then runs the same repo-local phase-2 installer that the
normal Cygwin path uses. During development, this portable path was exercised
in isolated `H:\TMP4SAS\...` directories so the pipeline could be validated
without depending on the user's preexisting global Cygwin installation.

If you are already inside that portable shell, or inside another supported
Cygwin shell, run the phase-2 installer directly:

```bash
bash install/install_cygwin.sh
```

The Windows/Cygwin install path now intentionally keeps its repo-local Perl
modules separate from Linux and macOS builds. During cross-platform validation,
the main failure was not simply "missing GD headers"; it was accidental reuse
of Linux-built repo-local Perl modules from a Cygwin shell, which caused
`GD.pm` and `Compress::Raw::Zlib` version mismatches. The installer/runtime
stack now avoids that by preferring `local/perl5-cygwin/` on portable Cygwin
instead of sharing one generic `local/perl5/` tree across operating systems.

Ubuntu:

```bash
bash install/install_ubuntu.sh
```

If an Ubuntu smoke test cannot be launched through the bundled Vagrant harness,
the same installer can be validated through Docker Desktop with an isolated
`ubuntu:24.04` container. During the current validation cycle, that Docker
fallback confirmed:

- `bash install/install_ubuntu.sh`
- `bash install/check_pipeline_install.sh`
- the top-level gunplot wrapper for `manhattan`, `local_manhattan`, and
  `local_gtf`

During the same 2026-06-08 Docker validation, the SAS ODA login probe inside
the image did not complete successfully, so the containerized SAS ODA path
should still be treated as unresolved until the SASPy/IOM startup issue is
fixed:

- `perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl --check-sas-oda-login-only`
- current observed failure: `SASIOConnectionTerminated: No SAS process attached.
  SAS process has terminated unexpectedly.`

Practical Ubuntu/Docker note:

- keep a repo-local Gencode cache at:
  - `cache/gtf/gencode.v49.annotation.gtf.gz`
  so local-GTF reruns are not blocked by transient EBI download errors
- under Docker, SAS ODA housekeeping calls such as remote upload verification
  and remote file deletion can be slower than on portable Cygwin, so a valid
  local-GTF figure may already exist even if a host-side Docker command later
  times out during cleanup
- for validation reruns, `KEEP_REMOTE_PLOT_DATA=1` is recommended, and any
  final `PGC_SCZ_SAS_local_top_hits_with_gtf_common.html` /
  `PGC_SCZ_SAS_local_top_hits_with_gtf_common.png` output should be checked
  before treating the whole run as failed

Why Ubuntu/Docker validation can feel slow:

- the first uncached image build installs Ubuntu packages, Java, gnuplot,
  ImageMagick, Python development headers, and repo-local Perl/Python
  dependencies inside the container; on the 2026-06-08 validation of this
  repository, `docker build -t multigwas-explorer-pipeline:latest .` took about
  4.2 minutes end to end
- once the image already exists, the post-build validation check is quick; the
  same 2026-06-08 run completed
  `bash install/check_pipeline_install.sh` inside the built image in about
  10 seconds
- the slowest scientific stage in Ubuntu Docker is usually the genome-wide
  gunplot Manhattan render rather than the local plots; on the same date, a
  one-SNP containerized gunplot validation of
  `manhattan,local_manhattan,local_gtf` took about 9 minutes 43 seconds total,
  and about 9 minutes 10 seconds of that time came from `plot_manhattan`
- the reason is data volume and bind-mounted I/O: the bundled schizophrenia
  wide subset contains 2,417,954 rows, the Docker Manhattan run scanned all of
  them, kept 2,396,261 rows across 9 plotted P columns, and wrote a
  ~110 MB intermediate `.plot.tsv` before gnuplot rendered the final PNG
- by contrast, the same Docker validation completed the one-locus
  `local_manhattan` and `local_gtf` stages in about 12 seconds and 11 seconds,
  respectively

macOS:

```bash
bash install/install_macos.sh
```

Post-install smoke test:

```bash
bash install/check_pipeline_install.sh
```

## Container Images

The repository also now ships a saved containerized installation path for users
who prefer Docker or Singularity / Apptainer instead of a host-native install.

Docker:

```bash
docker build -t multigwas-explorer-pipeline:latest .
docker run --rm -it multigwas-explorer-pipeline:latest \
  bash -lc "cd /opt/MultiGWAS-Explorer && bash install/check_pipeline_install.sh"
```

For real plotting runs, mount your GWAS data directories plus the SAS ODA
authinfo file into the container, for example:

```bash
docker run --rm -it \
  -e PIPELINE_WORKDIR=/opt/MultiGWAS-Explorer \
  -v /path/to/_authinfo:/root/_authinfo:ro \
  -v /path/to/gwas_drive_e:/mnt/e \
  -v /path/to/gwas_drive_g:/mnt/g \
  multigwas-explorer-pipeline:latest bash
```

Keep `PIPELINE_WORKDIR=/opt/MultiGWAS-Explorer` when running from the image. This
forces the wrappers to use the Linux-installed repository inside the container
even if you also mount a Windows `G:` drive for data. Avoid mounting a host
repository checkout over `/opt/MultiGWAS-Explorer`, because that would hide the
image's Linux `.venv-pipeline/` and repo-local Perl modules.

Singularity / Apptainer:

```bash
bash install/singularity/build_apptainer_image.sh
apptainer exec MultiGWAS-Explorer_pipeline.sif \
  bash -lc "cd /opt/MultiGWAS-Explorer && bash install/check_pipeline_install.sh"
```

The Apptainer definition is stored in:

- `install/singularity/MultiGWAS-Explorer_pipeline.def`

and should be built from the repository root because it copies the current repo
tree into `/opt/MultiGWAS-Explorer` inside the image.

Validation status:

- the Docker image build and post-build runtime smoke test were both completed
  successfully on the current workstation
- the Singularity / Apptainer definition is saved for users who need that HPC
  path, but it was not executed on the current workstation in this validation
  cycle

Portable Cygwin validation note:

- the top-level `auto_prepare_and_run_diff_gwas_with_gunplot.pl` wrapper was
  validated from an isolated portable Cygwin install under `H:\TMP4SAS\...`
- `bash install/check_pipeline_install.sh` now also reports the active Perl
  archname and GD version, which helped confirm that the portable shell was
  using `x86_64-cygwin-threads-multi` plus the Cygwin GD build rather than a
  leaked Linux repo-local module
- the wrapper now prefers `gnuplot` from the active shell `PATH` before trying
  older Windows-specific fallback locations
- in a healthy portable-Cygwin run, you should see:
  - `Using gnuplot executable: gnuplot`
- a practical wrapper-level validation command after installation is:

```bash
perl ./auto_prepare_and_run_diff_gwas_with_gunplot.pl \
  --spec configs/spec_pgc_scz_sex_common_automation.json \
  --plots manhattan,local_manhattan,local_gtf
```

Ubuntu Docker validation note:

- the same gunplot wrapper command above was also validated from an isolated
  Ubuntu Docker image after `bash install/install_ubuntu.sh`
- the repository now also ships a convenience smoke-test entry point for that
  path:

```bash
bash install/test_ubuntu_docker_gnuplot.sh
```

- by default, that script reproduces the top differential schizophrenia SNP
  `rs185665940` with `local_manhattan,local_gtf`; add
  `--include-manhattan` when you also want the slower genome-wide gnuplot
  Manhattan panel in the same Docker validation run
- on 2026-06-08, the rebuilt Ubuntu Docker image also validated the updated
  gunplot layout defaults:
  - genome-wide gunplot Manhattan now prints the GWAS track labels at the top
    of each subplot, matching the SAS ODA-style multi-track layout more
    closely
  - combined gunplot local Manhattan now defaults to a bottom gene-track view
    instead of the older vertical SNP/gene text labels
  - if you still prefer the old compact text-label style for local Manhattan,
    override it with:
    - `--local-manhattan-annotation labels`
- in the same 2026-06-08 Docker validation cycle, the repo-local SAS ODA
  helper bootstrap did not complete successfully inside the image:

```bash
perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl \
  --check-sas-oda-login-only
```

- current observed failure: `SASIOConnectionTerminated: No SAS process attached.
  SAS process has terminated unexpectedly.`

- for a quick Docker gnuplot verification, prefer a one-SNP local run first:

```bash
cp configs/spec_pgc_scz_sex_common_automation.json \
  /tmp/spec_pgc_scz_sex_common_automation.docker_gnuplot.json
perl ./auto_prepare_and_run_diff_gwas_with_gunplot.pl \
  --spec /tmp/spec_pgc_scz_sex_common_automation.docker_gnuplot.json \
  --plots local_manhattan,local_gtf \
  --target-snps rs185665940
```

- when you also include `manhattan`, expect that genome-wide stage to dominate
  runtime; in the 2026-06-08 Ubuntu Docker validation, a full
  `manhattan,local_manhattan,local_gtf` run for `rs185665940` finished
  successfully, but `plot_manhattan` alone took about 9 minutes 10 seconds
- after the 2026-06-08 layout update and image rebuild, the focused smoke
  tests completed in:
  - about 13 seconds for `--plots local_manhattan`
  - about 9 minutes 6 seconds for `--plots manhattan`
- after the 2026-06-11 Manhattan style refresh, a rebuilt Docker image again
  validated the genome-wide gnuplot panel against the schizophrenia dataset:
  - the repeating chromosome palette now follows the same SAS ODA-style color
    family across chromosomes instead of coloring by track index
  - the per-track GWAS sublabels remain at the top of each subplot
  - the PNG no longer injects an extra plot title above the panel, which makes
    the default output closer to the manuscript SAS figure
  - measured timing for that rerun was about 231 seconds for the image rebuild
    and about 8 minutes 46 seconds for the `--plots manhattan` smoke stage
- if you do not want to overwrite existing figure basenames in the mounted data
  directory, make a temporary copy of the spec inside the container and change
  `output_prefix`, `local_output_prefix`, `local_top_hits_csv_basename`, and
  `output_html_basename` before running the wrapper
- when that Docker smoke-test script is launched from a Windows portable-Cygwin
  shell, it now converts the repository path into a Docker-friendly host build
  context automatically; this avoids the older `path "/mnt/g/.../perlMCP4Gemini_Paper" not found`
  build failure under Docker Desktop

Container validation note:

- the Docker image path now uses the same `install/install_ubuntu.sh` logic as
  the host Ubuntu installer rather than a separate custom package recipe
- the current repository also saves a matching Singularity / Apptainer
  definition so the same Linux stack can be reproduced on HPC-style systems

Genome-wide gunplot cache-safety note:

- during validation, a stale interrupted rerun left a corrupted cached
  genome-wide wide subset on disk, which made the gunplot Manhattan figure look
  drastically different from the SAS ODA figure
- the gunplot wrapper now validates the cached wide subset against its manifest
  before reuse instead of trusting file existence alone
- `DiffGWASDeps/extract_significant_diff_gwas.pl` now writes the wide subset
  and its manifest atomically, which reduces the chance of leaving a
  half-overwritten genome-wide input behind after an interrupted rerun
- if a future genome-wide gunplot Manhattan plot suddenly looks truncated or
  structurally wrong, rerun:

```bash
perl ./auto_prepare_and_run_diff_gwas.pl \
  --spec configs/spec_pgc_scz_sex_common_automation.json \
  --step extract_wide_subset \
  --force
```

Reference-build-aware local GTF note:

- both the SAS ODA and gunplot entry scripts now accept:
  - `--reference-build hg19|hg38|t2t`
- the same setting can also be stored in the spec JSON as:
  - `reference_build`
- when no explicit override is given, the pipeline tries to infer the build in
  this order:
  - explicit spec fields such as `reference_build`
  - header tokens such as `POS_HG38` or `BP_HG19`
  - filename/path tokens such as `hg19`, `grch38`, `hs1`, or `chm13`
  - final fallback to `hg38`
- the built-in local-GTF profiles currently map to:
  - hg19 / GRCh37: GENCODE v49 lift37
  - hg38 / GRCh38: GENCODE v49
  - T2T / hs1 / CHM13v2.0: UCSC hs1 RefSeq GTF
- during validation on 2026-06-11, temporary spec reruns verified that
  explicit `reference_build=hg19` selected `FM.GTF_HG19` / `gtf_hg19`, and
  `reference_build=t2t` selected `FM.GTF_T2T` / `gtf_t2t`
- the bundled PGC schizophrenia sex example is now pinned explicitly to
  `hg19` because the older DANER-style source files do not carry an explicit
  build token in their headers
- the bundled PGC schizophrenia ancestry example is also treated as `hg19`;
  its raw PGCsumstatsVCF headers directly report `##genomeReference="GRCh37"`
- if your filenames and headers do not clearly encode the build, set
  `reference_build` explicitly instead of relying on the fallback

Genome-wide SAS-style gunplot note:

- the genome-wide gunplot Manhattan renderer now follows the same repeated
  chromosome palette family and top-of-panel GWAS label placement used by the
  SAS ODA multi-track figure
- small visual differences can still remain because SAS ODA and gnuplot do not
  rasterize points identically
- in practice, the remaining differences after the 2026-06-11 refresh were
  mostly limited to backend-specific point packing and antialiasing rather than
  different chromosome color logic or misplaced subplot labels

First SAS ODA login behavior:

- on the first SAS-backed run, if no saved SASPy authinfo entry for authkey
  `oda` is present, `run_sas_codes_or_script_in_ODA.pl` now prompts for the
  SAS ODA account/email and password
- the helper validates those credentials immediately with:
  - `proc setinit;run;`
- if validation fails, the helper warns that the supplied account/password may
  be wrong and does not keep the failed credentials
- if validation succeeds, the helper saves the working credentials into the
  SASPy authinfo file and later SAS ODA runs reuse them automatically
- for noninteractive use, the same bootstrap can be supplied explicitly with:
  - `--sas-oda-account EMAIL`
  - `--sas-oda-password PASS`
- to force a credential refresh even when a saved authinfo entry already
  exists, use:
  - `--prompt-sas-oda-auth`
- to validate the saved login directly without running a plot, use:

```bash
perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl \
  --check-sas-oda-login-only
```

## VS Code + Codex Workflow Across Operating Systems

This repository can be used directly from VS Code with Codex while keeping the
same repo-local scripts and MCP tools underneath.

Shared pattern:

1. open this repository as the VS Code workspace;
2. finish the matching platform installer first;
3. open an integrated terminal rooted in the repository;
4. start `server.pl` in one terminal and keep it running;
5. register the local MCP endpoint with Codex:

```bash
codex mcp add perl-bio --url http://127.0.0.1:8080/mcp
codex mcp list
```

6. start a fresh Codex session in that same workspace; and
7. ask Codex to run either the SAS ODA workflow or the alternative local
   gunplot workflow through the repository scripts.

Windows through portable Cygwin:

- preferred bootstrap:
  - `powershell -NoProfile -ExecutionPolicy Bypass -File .\install\install_windows_portable_cygwin.ps1`
- inside the Cygwin terminal, use repository paths such as:
  - `/cygdrive/g/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper`
- you can either launch VS Code from that shell or point the VS Code
  integrated terminal profile to the portable Cygwin `bash.exe`

macOS:

- install with:
  - `bash install/install_macos.sh`
- use the VS Code integrated `zsh` or `bash` terminal for both `server.pl`
  and any direct command-line validation

Ubuntu Linux:

- install with:
  - `bash install/install_ubuntu.sh`
- use the VS Code integrated `bash` terminal in the same way

This editor-centered pattern is also summarized in the manuscript supplement so
that the same VS Code plus Codex workflow can be reproduced across Windows,
macOS, and Ubuntu Linux.

The same repo-local runtime is then reused by:

- `auto_prepare_and_run_diff_gwas.pl`
- `auto_prepare_and_run_diff_gwas_with_gunplot.pl`
- `server.pl`
- `run_sas_codes_or_script_in_ODA.pl`

When these repo-local environments exist, the main entry scripts, MCP server,
and SAS ODA helper stack now prefer them automatically.

For the gunplot wrapper specifically, the runtime now prefers a real `gnuplot`
command from the active shell `PATH` before trying older Windows-specific
fallback binaries. This keeps portable Cygwin, Ubuntu, and macOS installs
self-contained instead of silently borrowing a developer-specific
`gnuplot.exe`. It now also validates cached genome-wide wide subsets against
their manifest row counts before reuse, so a stale interrupted rerun is less
likely to poison later Manhattan plots.

## Pipeline Overview

The default workflow is:

1. Start from raw GWAS files or a precomputed differential GWAS table.
2. Validate headers and infer or load the expected column mappings.
3. Merge compatible raw inputs by group when needed.
4. Sort and standardize the differential GWAS signals.
5. Extract a smaller wide-format subset for plotting.
6. Upload only the compact plotting subset to SAS ODA.
7. Generate one or more plots:
   genome-wide Manhattan
   local top-hit Manhattan
   local top-hit GTF
8. Download HTML/PNG artifacts and optionally clean remote SAS ODA inputs.

## Parallel Gunplot Pipeline

The repository also keeps a separate gunplot workflow in parallel with the SAS
ODA pipeline:

- `auto_prepare_and_run_diff_gwas_with_gunplot.pl`

This pipeline should stay independent from `auto_prepare_and_run_diff_gwas.pl`.
Use it when you want local rendering without SAS ODA, or when you want to
prototype layout logic without modifying the SAS ODA workflow.

Both pipelines now share one track-selection model:

- use pair prefixes such as `ALL`, `EUR`, or `ASN` to display differential
  standardized-difference tracks
- use GWAS labels such as `ALL_FEMALE`, `ALL_MALE`, `EUR_FEMALE`, or
  `EUR_MALE` to display single-GWAS association tracks
- use `--display-gwas` to request any subset or ordering of those tracks
- combine `--display-gwas` with `--target-snps` when you want inquiry-SNP local
  Manhattan or local GTF panels for specific rsIDs

Examples:

SAS ODA single-GWAS plot set:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --spec ./configs/spec_pgc_scz_sex_common_automation.json \
  --display-gwas ALL_FEMALE \
  --step plot_manhattan \
  --step plot_local_manhattan \
  --step plot_local_gtf
```

Noninteractive first-run SAS ODA credential bootstrap:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --spec ./configs/spec_pgc_scz_sex_common_automation.json \
  --step plot_manhattan \
  --sas-oda-account your_email@example.com \
  --sas-oda-password 'your_password'
```

gunplot single-GWAS plot set:

```bash
perl auto_prepare_and_run_diff_gwas_with_gunplot.pl \
  --spec ./configs/spec_pgc_scz_sex_common_automation.json \
  --display-gwas ALL_FEMALE \
  --plots manhattan,local_manhattan,local_gtf
```

gunplot forest plot for one inquiry SNP:

```bash
perl auto_prepare_and_run_diff_gwas_with_gunplot.pl \
  --spec ./configs/spec_pgc_scz_sex_common_automation.json \
  --plots forest \
  --target-snps rs185665940
```

Fast Docker-first gunplot validation:

```bash
perl auto_prepare_and_run_diff_gwas_with_gunplot.pl \
  --spec ./configs/spec_pgc_scz_sex_common_automation.json \
  --plots local_manhattan,local_gtf \
  --target-snps rs185665940
```

If that succeeds and you want the full genome-wide figure too, add
`manhattan`. On Ubuntu Docker, the genome-wide step is usually much slower than
the single-locus steps because it has to scan and thin millions of rows before
rendering.

Inquiry-SNP local panels with a custom displayed GWAS subset:

```bash
perl auto_prepare_and_run_diff_gwas_with_gunplot.pl \
  --spec ./configs/spec_pgc_scz_sex_common_automation.json \
  --display-gwas EUR,EUR_FEMALE,EUR_MALE \
  --target-snps rs185665940 \
  --plots local_manhattan,local_gtf
```

Example:

```bash
perl auto_prepare_and_run_diff_gwas_with_gunplot.pl \
  --spec ./configs/spec_pgc_scz_sex_common_automation.json \
  --target-snps rs17425819,rs185665940 \
  --plots local_manhattan \
  --local-manhattan-columns 2 \
  --local-manhattan-annotation gtf
```

Important gunplot local-Manhattan implementation detail:

- the default combined local-Manhattan annotation mode is now `gtf`, so the
  final combined figure prefers a bottom gene track instead of the older
  vertical SNP/gene label pair
- if you prefer the previous compact label mode, override it with:
  - `--local-manhattan-annotation labels`
- when `--local-manhattan-annotation gtf` is used, the combined local
  Manhattan panel is now driven by one unified scaled dataset rather than by
  stitching together unrelated x coordinates afterward
- the renderer first combines association rows from all loci, then combines
  the matching GTF rows, tags both with locus membership, and rescales genomic
  position into one shared panel x system
- the resulting artifact is:
  - `*.combined_scaled.tsv`
- that table contains both association and gene/exon rows and is the preferred
  debugging artifact for the combined gunplot local-Manhattan GTF mode
- older artifacts such as `*.combined.tsv` or `*.combined.gp` may still exist
  from the label-only renderer, but for the GTF annotation mode the
  `*.combined_scaled.tsv` table is now the source of truth

Both pipelines now also support optional user-designated adjacent-gene labels
for target SNPs:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --spec ./configs/spec_pgc_scz_sex_common_automation.json \
  --target-snps rs17425819,rs185665940,rs2564978 \
  --target-snp-genes rs17425819:JAK2,rs185665940:FANCL,rs2564978:CR1 \
  --step plot_local_manhattan
```

The same `target_snp_genes` mapping can also be stored in the spec JSON. In
both the SAS ODA and gunplot paths, these user-supplied labels now override the
automatic HaploReg / nearest-GTF fallback and are exported as `gene_source=USER`
in the top-hit CSV.

## PGC SCZ chrX Notes

For the PGC schizophrenia sex-stratified example in this repository, missing
chromosome X data can come from two different causes, and they should be
distinguished carefully.

Expected source-data behavior:

- `ALL_FEMALE` and `ALL_MALE` are autosome-only in the current input bundle.
- `EUR_FEMALE`, `EUR_MALE`, `ASN_FEMALE`, and `ASN_MALE` each have a separate
  chrX supplement file in addition to their autosomal file.

Concrete evidence in this repository:

- [configs/spec_pgc_scz_sex_common_automation.json](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/configs/spec_pgc_scz_sex_common_automation.json)
  wires `ALL_*` to one autosome file each, while `EUR_*` and `ASN_*` each use
  an extra `chrX` file.
- [merge_scz_sex_stratified_long.pl](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/merge_scz_sex_stratified_long.pl)
  labels the pooled strata as `*_AUTOSOME`, while EUR/ASN groups include the
  chrX supplements explicitly.
- [PGC_SCZ_female_vs_male_diff_effects_merged_long.manifest.tsv](</e:/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects_merged_long.manifest.tsv>)
  records the same source-file composition.

Observed row counts from the raw source files:

- `daner_PGC_SCZ_w3_81_0618a_all_female.gz`: `0` chrX rows
- `daner_PGC_SCZ_w3_81_0618a_all_male.gz`: `0` chrX rows
- `daner_scz_w3_HRC_chrX_eur_fem_deduped_0518e.gz`: `231,908` chrX rows
- `daner_scz_w3_HRC_chrX_eur_mal_deduped_0518e.gz`: `232,899` chrX rows
- `daner_scz_w3_HRC_chrX_asn_fem_run2.gz`: `193,305` chrX rows
- `daner_scz_w3_HRC_chrX_asn_mal_run2.gz`: `195,213` chrX rows

Important debugging implication:

- if chrX is missing from `ALL_*`, that is expected from the current source
  data bundle
- if chrX is missing from `EUR_*` or `ASN_*` downstream, suspect a pipeline
  artifact rather than the raw inputs

For this project, one concrete failure mode has already been seen:

- the sorted merged long file still contained EUR/ASN chrX rows
- but a later `PGC_SCZ_female_vs_male_diff_effects.tsv.gz` on disk was
  inconsistent with its own manifest and had been rewritten later with fewer
  rows than expected
- in that state, downstream standardized and wide outputs can appear to have
  lost chrX even though the raw and merged-long inputs still contain it

Recommended checks when chrX seems to disappear:

1. inspect the input composition in the spec and merged-long manifest first
2. verify whether chrX exists in the merged long sorted table before blaming
   plotting
3. compare the differential output row count against its manifest; if the file
   was partially regenerated or overwritten later, rebuild the differential,
   standardized, and wide descendants from the sorted merged long input

One plotting note:

- the gunplot pipeline now removes chrX from final figures by default unless
  `--no-remove-X-chr` is supplied
- that affects plotting only, not the upstream merged or differential data

## Supported Input Modes

- `raw_pgc_vcf_sumstats`
  Use this when starting from raw GWAS summary-statistics files that need to be
  merged, contrasted, and standardized.
- `precomputed_diff`
  Use this when pairwise differential effects already exist but still need
  standardization and plotting.
- `precomputed_diff_stdized`
  Use this when the standardized long differential GWAS table already exists
  and only extraction plus plotting remain.

## Header Detection and New GWAS Formats

The raw-input workflow can recognize several common header families and now
accepts GEMMA-style outputs such as `chr`, `rs`, `ps`, `allele1`, `allele0`,
`beta`, `se`, and `p_lrt`.

When automatic header guessing is not enough, you can provide extra aliases:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --gwas-dir /mnt/e/path/to/gwas_dir \
  --raw-column-alias-config ./raw_aliases.json \
  --generate-spec-only
```

Example alias file:

```json
{
  "ID": ["MARKER_ID"],
  "POS": ["GENOMIC_POS"],
  "PVAL": ["PVALUE_LRT"]
}
```

## Quick Start

Generate a draft spec from a GWAS directory:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --gwas-dir /mnt/e/path/to/gwas_dir \
  --generate-spec-only
```

Preview the inferred merged-wide spec for the AOA dataset without writing it:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --gwas-dir AOA_GWAS_Data/ \
  --preview-spec \
  --generate-spec-only
```

Write the inferred AOA merged-wide spec to disk:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --gwas-dir AOA_GWAS_Data/ \
  --spec-out configs/auto_aoa_merged.spec.json \
  --generate-spec-only
```

Run the merged-wide AOA plotting workflow from the same directory scan:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --gwas-dir AOA_GWAS_Data/ \
  --spec-out configs/auto_aoa_merged.spec.json \
  --plots manhattan,local_manhattan,local_gtf,forest
```

When the default top-hit threshold leaves no eligible loci in merged-wide mode,
prefer rerunning a targeted stage with explicit SNPs instead of assuming the
merged table was parsed incorrectly:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --spec configs/auto_aoa_merged.spec.json \
  --target-snps 1_790112_GATTT_G \
  --step plot_forest
```

Preview the inferred spec without writing it:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --gwas-dir /mnt/e/path/to/gwas_dir \
  --preview-spec \
  --generate-spec-only
```

Run the full pipeline from a spec:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --spec ./configs/your_spec.json
```

Generate only the derived config files:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --spec ./configs/your_spec.json \
  --mode configs
```

Rerun only one plotting stage:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --spec ./configs/your_spec.json \
  --step plot_local_gtf
```

## SAS ODA Wrapper Utilities

The low-level SAS ODA helper:

- `run_sas_codes_or_script_in_ODA.pl`

The helper now owns first-run SAS ODA credential bootstrap too:

- it checks whether the SASPy authinfo file already has an `authkey=oda` entry
- if not, it prompts interactively for the SAS ODA account/password
- it validates that login with `proc setinit;run;`
- only a successful login is written back to the authinfo file
- later SAS ODA wrapper calls then reuse that saved entry automatically

now supports bulk file operations in a single invocation. This is useful both
for manual debugging and for reducing repeated wrapper startup overhead in
automation.

Practical remote-path lessons from the recent SAS ODA debugging:

- quote remote home-directory paths like `'~/a.txt'` in shell examples so the
  literal `~` reaches the helper unchanged
- the helper now normalizes both `~/...` and absolute SAS home paths such as
  `/home/...` consistently for download, delete, and file-info operations
- delete requests now verify that the target path no longer resolves after the
  helper reports success

Examples:

```bash
perl -S run_sas_codes_or_script_in_ODA.pl \
  --upload-file a.tsv.gz \
  --upload-file b.tsv.gz
```

```bash
perl -S run_sas_codes_or_script_in_ODA.pl \
  --download-file '~/a.txt' \
  --download-file '~/b.txt' \
  --download-local-path ./a.txt \
  --download-local-path ./b.txt
```

```bash
perl -S run_sas_codes_or_script_in_ODA.pl \
  --delete-file old_a.tsv.gz \
  --delete-file old_b.tsv.gz
```

```bash
perl -S run_sas_codes_or_script_in_ODA.pl \
  --delete-file-rgx '.*\.png$'
```

Regex-based deletes now match both:

- the bare basename, such as `plot.png`
- the resolved remote path, such as `~/plot.png`

So both of these are valid:

```bash
perl -S run_sas_codes_or_script_in_ODA.pl --delete-file-rgx '.*\.png$'
perl -S run_sas_codes_or_script_in_ODA.pl --delete-file-rgx '~\/.*\.png'
```

The low-level SAS ODA helper now also supports hard submit timeouts through:

- `SAS_ODA_RUN_TIMEOUT_SECONDS`
- `SAS_ODA_RUN_TIMEOUT_GRACE_SECONDS`

Example:

```bash
SAS_ODA_RUN_TIMEOUT_SECONDS=300 \
SAS_ODA_RUN_TIMEOUT_GRACE_SECONDS=20 \
perl -S run_sas_codes_or_script_in_ODA.pl \
  --file ./run_sas_oda_local_top_hits_with_gtf.sas
```

Important implementation detail from debugging this project:

- the timeout path in `run_sas_codes_or_script_in_ODA.pl` now runs the real
  SAS submit in a fresh worker process instead of forking an already-loaded
  `Inline::Python` / `saspy` session
- this matters because the older in-process fork path could return an empty
  SAS result payload even for trivial code like `%put HELLO;`
- if a direct `SAS_ODA_Runner->run_file(...)` replay succeeds but the helper
  reports empty log/output, suspect the helper transport layer first, not the
  SAS macro stack
- after the worker-process fix, timeout-enabled submits and direct replays
  should agree much more closely

For `%include`-heavy SAS debugging, the same helper now also auto-detects
submissions that contain `%include` and switches into a safer debug path:

- it tries to refresh remote `~/macro.sas` files from matching local project
  copies before running the include
- it can run the resolved included file as a standalone debug submit before the
  parent script, so failures can be isolated to the included SAS file itself
- it runs a local preflight scan for likely compile blockers such as unmatched
  `/* */` comment structure or unterminated quotes
- it records line-numbered source context for suspicious lines in
  `output.html.info.txt`
- it attempts remote `PROC PRINTTO` log capture for the included file
- it fails fast instead of hanging silently when the include looks broken
- it keeps the `=== Submitted SAS Codes or file ===` section clean: helper
  `%include` diagnostics now live under `=== Dependency Logs ===` instead of
  being appended as bracketed banners to the submitted file/code display

Example:

```bash
SAS_ODA_RUN_TIMEOUT_SECONDS=60 \
SAS_ODA_RUN_TIMEOUT_GRACE_SECONDS=10 \
perl -S run_sas_codes_or_script_in_ODA.pl \
  --output-prefix lattice_compile_check \
  --code "%include '~/Lattice_gscatter_over_bed_track.sas'"
```

The resulting debug bundle is written to:

- `./lattice_compile_check/output.html.info.txt`

and now includes:

- local include-preflight findings
- refreshed-remote-file notes when a local replacement was uploaded first
- remote source excerpts around flagged lines when available
- the remote `PROC PRINTTO` log status if SAS ODA managed to create one

The include-debug wrapper now also uses a valid short SAS fileref for
`PROC PRINTTO`. This avoids a false compile failure from overlong fileref names
inside the helper-generated debug prologue.

Important SAS macro debugging lesson from this project:

- do not try to disable a macro abort with `*%abort 255;`
- in SAS macro code, that form can still execute `%abort`
- use either `/*%abort 255;*/` or `%* %abort 255;` instead
- this matters especially for `%include`-driven ODA runs, because a hidden
  macro abort can make the SAS ODA job fail while returning no useful inline
  log text through saspy
- when that happens, `output.html.info.txt` may show:
  - normal include-preflight notes
  - an empty `=== SAS Log ===` section
  - an empty `=== Output ===` section
- treat that pattern as "SAS ODA returned no usable log artifact", not as proof
  that the submitted SAS script ran successfully
- the right next steps are:
  - inspect the included macro for macro-statement comments like `*%abort ...;`
  - rerun through `run_sas_codes_or_script_in_ODA.pl` so include preflight can
    refresh the remote macro copy and save a debug bundle
  - preserve rendered runner inputs with `KEEP_RENDERED_DEBUG_FILES=1` when the
    failure happened inside the full plotting pipeline

Another practical lesson from the local-GTF reruns:

- a top-level automation command can time out in the outer terminal even after
  the actual SAS ODA plot has already completed and written local artifacts
- before concluding the plot failed, inspect the newest:
  - `run_local_hits_with_gtf_*/output.html.info.txt`
  - downloaded HTML or PNG outputs
  - selected top-hit CSV
- if the saved SAS log contains repeated markers like:
  - `The final figure is put here`
  - `Lattice gscatter plot is completed!`
  then treat the local GTF stage itself as successful, even if the parent
  automation process is still busy in a later post-run step

Scientific lesson from the merged-AOA genome-wide Manhattan rerun:

- when building a compact custom SAS plotting subset, normalize chromosome
  labels before the final sort and SAS import
- strip an optional `chr` prefix, map `X -> 23`, map `Y -> 24`, and drop rows
  whose chromosome or base-pair position still cannot be parsed
- sort on those normalized coordinate fields, not on the raw chromosome text
- otherwise SAS can import a literal `CHR='X'` as numeric missing `.`, which
  sorts before `chr1` and creates a false extra chromosome block at the far
  left of the genome-wide Manhattan plot
- for this repository, the safer pattern is to keep the genome-wide upload
  compact, retain only required plotting columns, and either normalize the
  subset upstream or do the `CHR_RAW -> CHR` conversion explicitly in the SAS
  import data step before plotting

Recommended `run_sas_codes_or_script_in_ODA.pl` file-handling debug process:

1. Validate the saved SAS ODA login first:

```bash
perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl \
  --check-sas-oda-login-only
```

2. If the issue is listing or remote file existence, stay in file-management
   mode and probe that path before any SAS code submit:

```bash
perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl \
  --output-prefix oda_file_probe \
  --dir4listing '~' \
  --file-info '~/some_remote_file.txt'
```

3. If upload/download behavior looks suspicious, test transfer-only mode with a
   small file and an explicit output prefix:

```bash
perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl \
  --output-prefix oda_transfer_probe \
  --upload-file ./small_test.txt \
  --download-file '~/small_test.txt' \
  --download-local-path ./small_test.roundtrip.txt
```

4. If stale remote files may be confusing the rerun, explicitly check both the
   basename and the resolved remote path, then remove only the exact target you
   intend to refresh:

```bash
perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl \
  --output-prefix oda_stale_remote_probe \
  --file-info 'Lattice_gscatter_over_bed_track.sas' \
  --file-info '~/Lattice_gscatter_over_bed_track.sas'
```

5. Only after listing, file-info, upload, and download succeed should you move
   on to SAS code execution probes:

```bash
perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl \
  --output-prefix oda_code_probe \
  --code "%put HELLO_FROM_ODA;"
```

6. If the real failure involves `%include`, rerun the include target or the
   parent code with a short timeout and a dedicated output prefix so the debug
   bundle is easy to find:

```bash
SAS_ODA_RUN_TIMEOUT_SECONDS=60 \
SAS_ODA_RUN_TIMEOUT_GRACE_SECONDS=10 \
perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl \
  --output-prefix include_debug_probe \
  --code "%include '~/Lattice_gscatter_over_bed_track.sas';"
```

7. Inspect the generated debug bundle, especially:
   - `<output-prefix>/output.html.info.txt`
   - `=== Dependency Logs ===`
   - `=== SAS Log ===`
   - `=== Output ===`

What to look for:

- listing or `--file-info` failures before any submit
  - likely SAS ODA home-directory / control-plane / session-side file-handling
    problem, not a plotting-macro bug
- upload succeeds but download or `--file-info` for the same path fails
  - likely remote-path mismatch, stale-home collision, or transfer/control-plane
    issue
- empty `=== SAS Log ===` plus empty `=== Output ===`
  - likely helper transport failure, early `%abort`, or early SAS ODA abort
- populated `=== Dependency Logs ===`
  - include-preflight findings, refreshed remote include notes, or remote
    compile-log clues
- remote file-info failures after submit
  - likely transfer/control-plane issue rather than a plotting macro bug

For file-handling-only incidents, do not start by debugging `%include`,
graphics macros, or SAS compile logs. First confirm:

- login works
- `--dir4listing` works
- `--file-info` works for the exact remote path you care about
- upload and download both work on a small round-trip file
- quoted `~/...` paths and absolute `/home/...` paths both resolve the same
  way for the helper action you are testing

Pay attention to SAS macro loading as a separate layer:

- file-only helper actions such as listing, `--file-info`, upload, download,
  and delete should not be diagnosed first as SAS macro problems
- when a new SAS ODA session is created, the runner auto-loads macros from
  `~/Macros` once via `importallmacros_ue(...)`
- when a persistent session is reused, that global macro bootstrap is not run
  again
- when the submitted SAS program already contains self-contained `%include`
  usage, the helper disables the global `importallmacros_ue` bootstrap for that
  submit and relies on the included files instead
- if logs mention `Macro load may have failed`, first decide whether the
  failing operation was actually a file-handling step or a true SAS-code step
- if the failure is a true SAS-code step, also inspect for open-code macro
  control statements such as `%if`, `%do`, `%end`, `%else`, `%goto`, or
  `%return` outside a `%macro/%mend` block, because the helper now flags those
  patterns during preflight
- do not add redundant `%importallmacros_ue` calls inside pipeline-generated
  SAS scripts unless you are deliberately debugging macro bootstrap behavior

Helpful extras:

- `--dry-run`
  - print resolved session settings without connecting
- `--prompt-sas-oda-auth`
  - force credential refresh
- `KEEP_RENDERED_DEBUG_FILES=1`
  - preserve rendered runner inputs when the failure occurred inside a higher
    level plot wrapper
- `output.run.status.json`
  - live status sidecar for the current helper run
- `<output-prefix>/output.macro_bootstrap.log.txt`
  - bootstrap-only trace/log file for the global `~/Macros` autoload path
  - it is created as soon as macro bootstrap starts, so even a timeout should
    still leave a partial trace file behind

Quick AI Debug Ladder for `run_sas_codes_or_script_in_ODA.pl`
-------------------------------------------------------------

For future AI-driven debugging, start with the smallest possible probe and only
move to `%include` or plotting code after the lower layer succeeds.

1. login-only probe:

```bash
perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl \
  --check-sas-oda-login-only
```

2. file-only probe:

```bash
perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl \
  --output-prefix oda_file_probe \
  --dir4listing '~' \
  --file-info '~/importallmacros_ue.sas'
```

3. small upload/download round trip:

```bash
perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl \
  --output-prefix oda_transfer_probe \
  --upload-file ./small_test.txt \
  --download-file '~/small_test.txt' \
  --download-local-path ./small_test.roundtrip.txt
```

4. minimal no-macro SAS submit:

```bash
perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl \
  --output-prefix oda_code_probe \
  --code "%put HELLO_FROM_ODA;"
```

5. minimal HTML-producing SAS submit:

```bash
cat > codex_proc_print_smoke.sas <<'EOF'
proc print data=sashelp.class;
run;
EOF

perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl \
  --output-prefix oda_proc_print_probe \
  --file codex_proc_print_smoke.sas
```

6. forced default `~/Macros` bootstrap timing probe:

```bash
cat > codex_put_smoke.sas <<'EOF'
%put HELLO_FROM_ODA;
EOF

SAS_ODA_AUTOLOAD_MACROS=1 \
perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl \
  --output-prefix oda_macro_bootstrap_probe \
  --file codex_put_smoke.sas \
  --run-timeout-seconds 30
```

7. manual in-SAS `~/Macros` load timing probe:

```bash
cat > codex_manual_importallmacros_smoke.sas <<'EOF'
%let _pipeline_start=%sysfunc(datetime());
%let _home=%sysfunc(pathname(HOME));
%include "&_home/importallmacros_ue.sas";
%importallmacros_ue(MacroDir=&_home/Macros,fileRgx=.,verbose=0);
%let _pipeline_end=%sysfunc(datetime());
%put MACRO_LOAD_ELAPSED_SECONDS=%sysfunc(round(%sysevalf(&_pipeline_end-&_pipeline_start),0.01));
%put HELLO_AFTER_MACRO_LOAD;
EOF

SAS_ODA_AUTOLOAD_MACROS=0 \
perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl \
  --output-prefix oda_manual_macro_probe \
  --file codex_manual_importallmacros_smoke.sas
```

How to read the new bootstrap diagnostics:

- `Upload step: macro bootstrap helper: importallmacros_ue.sas ...`
  - the tiny helper file is being uploaded
  - this is not the actual macro bootstrap submit
- `SAS ODA macro bootstrap started at ...`
  - the wrapper has entered the real global `~/Macros` autoload submit
- `output.run.status.json`
  - inspect `bootstrap_started_at`, `bootstrap_finished_at`,
    `bootstrap_elapsed_seconds`, `bootstrap_ok`, and `bootstrap_log_path`
- `output.macro_bootstrap.log.txt`
  - inspect this first when the run times out during macro bootstrap
  - if only `Bootstrap Start:` is filled and the file still says
    `Status: running`, the wrapper reached bootstrap entry but never got a
    completed `sess.submit(...)` return
- `output.html.info.txt`
  - inspect `=== Status Snapshot ===` and `=== Dependency Logs ===` together

Practical interpretation:

- if the manual in-SAS macro timing probe is fast but the forced default
  bootstrap probe stalls, the problem is in the helper bootstrap path rather
  than in the `~/Macros` library itself
- if `%put HELLO_FROM_ODA;` fails while macro autoload is skipped, fix the base
  SAS ODA submit path before touching macro bootstrap logic
- when debugging from PowerShell on Windows, prefer `--file` for multi-token
  SAS code because quoting through PowerShell into Cygwin can silently truncate
  inline `--code` payloads
- `SAS_ODA_RUN_TIMEOUT_SECONDS` and `SAS_ODA_RUN_TIMEOUT_GRACE_SECONDS`
  - keep failed probes short and reproducible

## Plot Types

- Genome-wide Manhattan
  Uses standardized differential P values plus per-GWAS association P tracks.
- Local top-hit Manhattan
  Selects lead loci and renders focused multi-track local views.
- Local GTF
  Adds a bottom gene track around the selected local region using a region-
  limited Gencode subset generated inside SAS ODA.

## Common-Association Top-Hit Mode

The pipeline supports a second local-hit selection mode focused on loci with
replicable single-GWAS association:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --spec ./configs/your_spec.json \
  --step plot_local_manhattan \
  --get-common-associations
```

This mode keeps loci where:

- one GWAS provides the strongest association signal
- another GWAS shows at least nominal association
- both GWASs point in the same effect direction

Important selector details:

- `top_hit_dist_bp` is the total exclusion span, not the half-window. The
  selector internally prunes within `BP +/- 0.5 * top_hit_dist_bp`.
- The common-association selector now defaults to `top_hit_dist_bp=1e6`
  unless your spec overrides it.
- The local-top-hit wrappers no longer hard-cap the helper at `15` loci.
  Use `top_hit_max_loci` in the spec or `TOP_HIT_MAX_LOCI` in the runner
  environment when you want an explicit cap; `0` means "no hard limit".
- In the bundled PGC schizophrenia sex-stratified validation on
  2026-06-14, the shared-association verifier found `11,703` candidates at
  `COMMON_ASSOC_P < 5e-8` and retained `103` loci after MAF QC and
  `1e6`-span pruning. The earlier undercount came from an overly broad
  `1e8` pruning span plus a hidden wrapper cap of `15` loci.

## Local GTF Gene-Track Options

Local GTF plots are protein-coding-only by default. The plotting scripts build
a region-limited Gencode subset so local gene tracks can still be made more
informative even when the original SAS ODA gene table is sparse or uses a
different schema. If you want non-coding genes too, enable them explicitly in
the config/spec.

The local GTF wrappers now also prefer the in-repo plotting macro:

- `DiffGWASDeps/SNP_Local_Manhattan_With_GTF.sas`

instead of relying only on whatever copy may already exist inside SAS ODA.
This makes local GTF reruns more reproducible and gives the pipeline a single
versioned place to patch local-axis and label behavior.

For large windows, the local top-hit GTF path now pre-extracts the requested
Gencode region locally and uploads only that subset to SAS ODA. This avoids the
older failure mode where SAS `WORK` had to materialize an oversized GTF table
for the selected locus window.
That uploaded subset is now gzip-compressed before transfer, and the SAS import
path reads it explicitly through:

```sas
filename gtfdata zip "~/local_gtf_subset_*.tsv.gz" gzip;
```

The local SNP/GTF plotting path now also separates:

- gene-search expansion distance
- final displayed x-axis bounds

For local GTF plots, the displayed x-axis is now forced back to the observed
association-signal span from the uploaded GWAS subset, so a large gene-search
window no longer causes the final chr axis to start at `0` just because the
expanded gene-track window crosses the chromosome start.

The local GTF wrapper is also more resilient for long SAS ODA runs:

- it defaults to one-shot ODA submits instead of relying on the slower local
  persistent-session relay
- if the first submit returns an incomplete control-plane result, it retries
  once automatically
- if SAS already finished remotely but the expected final HTML download is
  flaky, the wrapper can now reuse the helper-saved `sas_res_*.html` artifact
  and separately download the final PNG path reported in the SAS log
- when that PNG is available, the wrapper now opens a figure-first HTML page
  that embeds the completed plot and keeps the raw SAS HTML beside it as a
  `.sasraw.html` sidecar, instead of opening the sparse helper HTML directly

For target-SNP extraction, the standardized long differential GWAS output is
now bgzip-indexed with `tabix` when `bgzip/tabix` are available locally. The
single-SNP wide extractor then uses that index for the region/window pass after
it resolves the target SNP location, instead of doing a second full-file scan.

The SAS ODA plotting wrappers now use the same timeout/retry pattern across:

- `DiffGWASDeps/run_sas_oda_local_top_hits_with_gtf_download_html.sh`
- `DiffGWASDeps/run_sas_oda_local_top_hits_manhattan_download_png.sh`
- `DiffGWASDeps/run_sas_oda_single_snp_with_gtf_download_html.sh`

Shared wrapper controls:

- `ODA_HELPER_TIMEOUT_SECONDS`
- `ODA_HELPER_TIMEOUT_GRACE_SECONDS`

Stage-specific submit controls:

- local top-hit GTF:
  - `GTF_SUBMIT_TIMEOUT_SECONDS`
  - `GTF_SUBMIT_TIMEOUT_GRACE_SECONDS`
  - `GTF_SUBMIT_MAX_ATTEMPTS`
- local top-hit Manhattan:
  - `LOCAL_MH_SUBMIT_TIMEOUT_SECONDS`
  - `LOCAL_MH_SUBMIT_TIMEOUT_GRACE_SECONDS`
  - `LOCAL_MH_SUBMIT_MAX_ATTEMPTS`
- single-SNP local GTF:
  - `SINGLE_SNP_GTF_SUBMIT_TIMEOUT_SECONDS`
  - `SINGLE_SNP_GTF_SUBMIT_TIMEOUT_GRACE_SECONDS`
  - `SINGLE_SNP_GTF_SUBMIT_MAX_ATTEMPTS`

These wrappers now fail fast or retry when SAS ODA stalls during:

- upload/download/list/file-info helper calls
- the actual SAS submit

Local GTF tracks are protein-coding-only by default:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --spec ./configs/your_spec.json \
  --step plot_local_gtf
```

If you explicitly want non-coding genes too, set it in the spec:

```json
{
  "include_non_protein_coding_genes_in_local_gtf": 1
}
```

When a locus has many overlapping genes, the local GTF stage now keeps the
figure width and height unchanged and instead slightly increases the lower
gene-track share by auto-tuning the SAS `pct4neg_y` ratio.

The current SAS local-GTF wrappers now start from the same larger lower-track
base before that auto-tuning runs:

- `DiffGWASDeps/run_sas_oda_local_top_hits_with_gtf_download_html.sh`
- `DiffGWASDeps/run_sas_oda_single_snp_with_gtf_download_html.sh`

Both wrappers now default `GTF_PCT4NEG_Y` to `1.4`, which keeps the bottom
gene/exon track visibly larger in manuscript-scale figures and one-locus debug
reruns alike.

Use a dedicated local-GTF window without changing the local Manhattan window:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --spec ./configs/your_spec.json \
  --step plot_local_gtf \
  --local-gtf-window-bp 2e7
```

`local_gtf_window_bp` now controls two things together:

- the genomic half-window used to extract the local GTF subset
- the actual displayed x-axis half-window passed into the final
  `SNP_Local_Manhattan_With_GTF` call

So if you request `--local-gtf-window-bp 1e9`, the final local GTF plot is
expected to render approximately `+/-1e9 bp` around the selected top hit,
rather than silently staying near the older default display distance.

For practical pipeline testing, a relaxed large window such as `1e8` is usually
a better first rerun target than `1e9`. It is large enough to stress the local
GTF path while still keeping the uploaded subset and SAS ODA rendering more
manageable.

Example:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --spec ./configs/your_spec.json \
  --step plot_local_gtf \
  --get-common-associations \
  --local-gtf-window-bp 1e8
```

## Local Manhattan Layout Controls

Local top-hit Manhattan panels now support manual label-placement overrides for
cases where SNP-gene labels are hard to optimize automatically.

Available CLI overrides:

- `--local-max-hits-per-fig N`
  Sets the requested maximum number of local top-hit columns per figure. The
  pipeline now allows up to `30` columns per panel.
- `--local-manhattan-angle4xaxis-label N`
  Rotates the SNP-gene x-axis labels.
- `--local-manhattan-xgrp-y-pos N`
  Moves the SNP-gene label row vertically.
- `--local-manhattan-yoffset-top N`
  Adjusts the label offset for the upper local track.
- `--local-manhattan-yoffset-bottom N`
  Adjusts the label offset for the lower local track.
- `--local-manhattan-fontsize N`
  Changes local SNP-gene label font size.
- `--local-manhattan-y-axis-label-size N`
  Changes local Manhattan y-axis title size.
- `--local-manhattan-y-axis-value-size N`
  Changes local Manhattan y-axis tick-label size.

Example rerun with manual label tuning:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --spec ./configs/your_spec.json \
  --step plot_local_manhattan \
  --local-max-hits-per-fig 30 \
  --local-manhattan-xgrp-y-pos -2.5 \
  --local-manhattan-yoffset-top 14 \
  --local-manhattan-yoffset-bottom 0.5 \
  --local-manhattan-angle4xaxis-label 60 \
  --local-manhattan-fontsize 3.0
```

The same values can be stored in the spec JSON:

```json
{
  "local_max_hits_per_fig": 30,
  "local_manhattan_xgrp_y_pos": "-2.5",
  "local_manhattan_yoffset_top": "14",
  "local_manhattan_yoffset_bottom": "0.5",
  "local_manhattan_angle4xaxis_label": "60",
  "local_manhattan_fontsize": "3.0",
  "local_manhattan_y_axis_label_size": "2.8",
  "local_manhattan_y_axis_value_size": "2.2"
}
```

## SAS ODA Strategy

This project intentionally keeps most heavy GWAS processing local and reserves
SAS ODA for compact plotting inputs. That reduces upload time, avoids oversized
remote jobs, and makes reruns more practical.

The automation layer can also reuse one shared wide-format plotting subset
across multiple plot stages in the same run. When `keep_remote_plot_data` is
enabled in the spec, later reruns can skip re-uploading the same large plot
input file.

For SAS local-Manhattan reruns with `target_snps`, the wrapper is now more
aggressive about staying compact:

- it can build a union local-window wide subset directly from `SOURCE_LONG_GZ`
  using `extract_single_snp_wide_diff_gwas.pl`
- it caches that compact subset under `cache/local_manhattan_reuse/`
- it uploads the compact target-window subset instead of the full wide plot
  input whenever that targeted path is available

This is especially helpful when the default wide plot input is hundreds of
megabytes and the run only needs a few local target SNP windows.

For local GTF plots, the safer default is now intentionally stricter than the
local Manhattan default:

- `local_max_hits_per_fig` continues to control local Manhattan batching
- `local_gtf_max_hits_per_fig` defaults to `1`
- if you set `--local-max-hits-per-fig N`, that explicit CLI override still
  applies to both local Manhattan and local GTF reruns
- if you want a local-GTF-only override in JSON without changing local
  Manhattan batching, use:

```json
{
  "local_gtf_max_hits_per_fig": 1
}
```

The local GTF shell runner now also auto-scales the default
`GTF_DESIGN_HEIGHT` from the number of stacked scatter tracks in
`GTF_ASSOC_PVARS`. If you leave `GTF_DESIGN_HEIGHT` unset, taller multi-track
plots get a modest height increase automatically; explicit user height settings
still win.

For the current manuscript-style local-GTF reruns, both SAS wrappers also now
share a practical baseline height of `1000` pixels before any explicit user
override. In other words, the default single-SNP and batched local-GTF paths
now start from the same larger canvas plus the same larger bottom gene-track
share, which produces more readable nearby-gene labels in the final PNGs.

For top SNP labels in the local GTF headroom, keep these two SAS macro
parameters in mind:

- `Yoffset4textlabels=2.5`
  Controls the SNP-label position within the top headroom. It moves the target
  SNP labels up or down in y-axis-value units, and the default `2.5` works for
  many cases.
- `yoffset4max_drawmarkersontop=0.15`
  Controls the height of the top headroom itself when SNP labels are drawn on
  top. In that mode, it overrides the normal `yaxis_offset4max` path and is
  internally adjusted again from the number of scatter tracks before the final
  `offsetmax` is assigned to the y-axis.

Important caveat:

- `Yoffset4textlabels` is also auto-tuned internally, so in single-label cases
  manual changes may appear weaker than expected.
- if one top SNP still does not sit in the middle of the headroom, the
  lower-level fallback is to edit the lattice macro directly near the single-
  label branch around the `1000/&track_height` ratio and rerun the local GTF
  script for validation.

The SAS local-Manhattan shell runner is now also more informative when a step
fails before SAS submission:

- the outer automation wrapper captures shell-wrapper stderr
- it now reports the latest likely local log path, such as
  `run_local_hits_manhattan_png_*/output.html.info.txt`
- mid-upload failures are therefore easier to distinguish from SAS macro or ODA
  submit failures

## Local GTF Overflow Debugging

When a large-window local GTF rerun looks like it "finished" but the HTML
contains blank regions, check the saved `run_local_hits_with_gtf_*/output.html.info.txt`
before trusting the figure.

Important failure signature:

- `ERROR: Insufficient space in file WORK.FINAL.DATA.`
- related `_DOCTMP...` / insufficient-space / damaged-dataset messages
- very large scatter-group counts in the lattice macro path

Important interpretation:

- later markers such as `The final figure is put here` or
  `Lattice gscatter plot is completed!` do not prove the figure is valid after
  `WORK` overflow has already damaged the intermediate SAS datasets
- in that case, reduce at least one of:
  - local GTF window size
  - local GTF loci per run or per figure
  - stacked GTF track count

The new safer defaults were added specifically to reduce this failure mode:

- local GTF defaults to one locus per figure or batch
- default local GTF height increases with stacked-track count

## MCP Server Integration

`server.pl` exposes the automation entry point as the MCP tool
`auto_prepare_and_run_diff_gwas`.

That tool supports:

- background execution with PID polling
- full pipeline runs or config-only runs
- targeted reruns by exact step name
- common-association local-hit mode
- local GTF window overrides that control both the extracted gene subset and
  the displayed local GTF plot range
- local-GTF-specific batching through `local_gtf_max_hits_per_fig`, which now
  defaults to `1` to reduce SAS ODA `WORK` overflow on large windows
- automatic default local GTF height scaling from the number of stacked
  `GTF_ASSOC_PVARS` tracks when `GTF_DESIGN_HEIGHT` is left unset
- gzipped local GTF subset upload for large-window local gene-track reruns
- local Manhattan label-placement overrides
- local top-hit column cap up to 30 panels per figure
- `target_snps` and optional `target_snp_genes` overrides for user-designated
  adjacent-gene labels in both SAS local Manhattan and local GTF outputs
- the local GTF protein-coding-only option

This keeps the MCP server focused on stable orchestration while scientific
logic stays in versioned project scripts.

## Repository Layout

- `auto_prepare_and_run_diff_gwas.pl`
  Main automation wrapper.
- `DiffGWASDeps/`
  Differential GWAS helpers, SAS templates, and shell wrappers.
- `configs/`
  Presets, generated configs, and example specs.
- `gwas-sas-oda-workflow/`
  Scientific skill documentation for Codex-assisted use.
- `server.pl`
  Perl MCP server.

## Example Outputs

Typical outputs include:

- merged long differential GWAS tables
- standardized long GWAS tables
- wide-format plot subsets
- genome-wide Manhattan PNG and HTML wrappers
- local top-hit Manhattan HTML wrappers
- local GTF HTML outputs
- local top-hit CSV exports
- manifest files describing generated artifacts

## Additional Documentation

- [GENERALIZE_DIFF_GWAS_PIPELINE.md](./GENERALIZE_DIFF_GWAS_PIPELINE.md)
  Deeper notes on the generalized pipeline contract, plotting behavior, and MCP
  wrapper behavior.
- [gwas-sas-oda-workflow/SKILL.md](./gwas-sas-oda-workflow/SKILL.md)
  Workflow guidance for Codex or other agents using this repository.
