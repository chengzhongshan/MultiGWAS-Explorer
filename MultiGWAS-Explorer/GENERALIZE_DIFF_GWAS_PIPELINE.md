# Generalizing the Differential GWAS Pipeline

This pipeline now has a rollback snapshot in:

- `backup_generalize_diff_gwas_20260424_224826/`

## What changed

The two key Perl extractors are no longer tied only to the PGC schizophrenia female-vs-male layout:

- `DiffGWASDeps/extract_significant_diff_gwas.pl`
- `DiffGWASDeps/extract_single_snp_wide_diff_gwas.pl`

They now accept a reusable JSON config via `--config`, plus command-line overrides.

## Config contract

Use the template:

- [configs/diff_gwas_config_schema.json](/G:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/configs/diff_gwas_config_schema.json)

Important fields:

- `project_tag`: used in auto-generated output names
- `input`: long standardized differential GWAS input
- `pair_col`: grouping column, usually `PAIR_TAG`
- `base_cols`: key columns used to define unique SNP rows
- `value_fields`: measures to pivot into wide format
  - when top-hit MAF filtering should prefer GWAS frequencies, also keep
    `GROUP1_FRQ_A`, `GROUP1_FRQ_U`, `GROUP2_FRQ_A`, `GROUP2_FRQ_U`,
    `GROUP1_INFO`, and `GROUP2_INFO`
- `filter_fields`: P-like fields that drive subset retention
- `pair_map`: map from raw pair tag to output prefix
- `prefix_order`: preferred output prefix order
- `window_bp`: used by the single-locus extractor
- `char_lengths`: character-column lengths for generated SAS input blocks
- `alias_map`: derived SAS aliases such as `ALL_STD_P = ALL_STD_DIFF_P`

Top-hit filtering controls now live one level up in the comparison spec /
runner-config layer:

- `top_hit_maf_threshold`
  - default: `0.01`
- `gnomad_freq_file`
  - optional local TSV or TSV.GZ lookup used only when GWAS allele frequencies
    are unavailable for a top-hit candidate
- `gnomad_population_map`
  - optional token map such as `EUR:NFE,ASN:EAS,AFR:AFR`

## Presets

Current working preset for the original project:

- [configs/preset_pgc_scz_sex_diff.json](/G:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/configs/preset_pgc_scz_sex_diff.json)

Template for ancestry or any other pairwise/groupwise differential GWAS:

- [configs/preset_generic_ancestry_diff.json](/G:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/configs/preset_generic_ancestry_diff.json)
- [configs/preset_pgc_scz_ancestry_diff.json](/G:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/configs/preset_pgc_scz_ancestry_diff.json)

## Automation layer

The pipeline now also has a higher-level automation entry point:

- [auto_prepare_and_run_diff_gwas.pl](/G:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/auto_prepare_and_run_diff_gwas.pl)

This script takes one comparison spec JSON and will:

- validate the input comparison definition
- auto-generate the matching merge, diff, preset, and runner config JSON files
- reuse existing outputs when present, unless `--force` is requested
- run the needed local prep stages
- trigger the requested SAS ODA runners for genome-wide Manhattan, local
  top-hit Manhattan, local top-hit GTF plots, and top-hit forest plots
- optionally narrow those plot families to a user-selected displayed GWAS subset
  through `--display-gwas`, including true single-GWAS rendering mode when only
  one GWAS association track is selected
- optionally render those same plot families for explicit inquiry SNPs through
  `--target-snps`

The repository also keeps a parallel non-SAS plotting entry point:

- [auto_prepare_and_run_diff_gwas_with_gunplot.pl](/G:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/auto_prepare_and_run_diff_gwas_with_gunplot.pl)

Keep that gunplot pipeline separate from the SAS ODA workflow. It is intended
for independent rendering and layout experimentation, not as an automatic
fallback inside the SAS ODA pipeline.

Forest plotting is now part of both entry points:

- SAS ODA path:
  - `--step plot_forest`
- gunplot path:
  - `--plots forest`

Shared forest-plot behavior:

- single-SNP mode renders one forest plot with cohort names on the y-axis and
  `OR and 95% CI` on the x-axis
- multi-SNP mode renders one panel per cohort / displayed GWAS track, with SNP
  IDs on the left y-axis and nearby gene labels on the right y-axis
- common and differential hits are separated by a horizontal dashed divider
- genome-wide significant points (`P < 5e-8`) receive a star marker
- both backends reuse the same requested top-hit CSV generation logic, so the
  same MAF-aware top-hit filtering feeds the forest plot and the local locus
  panels

## Manuscript Table Outputs

The manuscript-table refresh helper now writes the supplementary hit tables in
full-strata form by default instead of emitting a compact pooled-only layout.

Canonical supplementary outputs:

- `manuscript_assets/tables/Table_S1_all_common_association_loci.csv`
- `manuscript_assets/tables/Table_S2_differential_loci.csv`

These standard filenames now retain:

- pooled and ancestry-specific association `P`, `BETA`, and `SE`
- pairwise differential `P`, `BETA`, and `SE` for `ALL`, `EUR`, and `ASN`
- standardized differential `P` values such as:
  - `ALL_STD_DIFF_P`
  - `EUR_STD_DIFF_P`
  - `ASN_STD_DIFF_P`
- top-hit QC provenance such as:
  - `selected_maf`
  - `maf_source`
  - `gwas_group1_maf`
  - `gwas_group2_maf`

For manuscript-editing convenience, the same helper also mirrors those richer
exports to:

- `manuscript_assets/tables/Table_S1_all_common_association_loci_full_strata.csv`
- `manuscript_assets/tables/Table_S2_differential_loci_full_strata.csv`

The narrower main-text tables are still emitted separately:

- `manuscript_assets/tables/Table_2_representative_common_loci.csv`
- `manuscript_assets/tables/Table_1_top_differential_locus.csv`

Recommended regeneration command:

```bash
perl DiffGWASDeps/regenerate_manuscript_hit_tables.pl \
  --config configs/spec_pgc_scz_sex_common_automation.json \
  --wide /mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects.stdized.wide_beta_se_p_p_lt_0p05.final.tsv.gz \
  --common-loci tmp_common_verify_postfix_1e6.tsv \
  --gtf cache/gtf/gencode.v49lift37.annotation.gtf.gz \
  --output-dir manuscript_assets/tables
```

Legacy recovery helpers:

- `DiffGWASDeps/augment_common_hits_table_s1.pl`
  Enriches an older pooled-only `Table_S1` CSV from the wide GWAS table when a
  fresh regeneration was not performed.
- `DiffGWASDeps/export_augmented_table_s1_excel.ps1`
  Converts a full-strata `Table_S1` CSV into a Windows `.xlsx` workbook with
  grouped headers for manuscript editing.

## Top-Hit MAF Safeguard

Local top-hit selection now applies the same MAF safeguard before either the
gunplot or SAS ODA local-panel renderers consume the selected loci:

- prefer GWAS-derived MAF when the wide plotting subset carries
  `*_GROUP1_FRQ_A`, `*_GROUP1_FRQ_U`, `*_GROUP2_FRQ_A`, and `*_GROUP2_FRQ_U`
- convert effect-allele frequency to MAF as `min(AF, 1-AF)`
- for pairwise differential top hits, use the smaller of the two compared
  GWAS-specific MAF values
- only fall back to gnomAD if those GWAS frequencies are absent for the
  candidate SNP
- if neither GWAS nor configured gnomAD MAF is available, mark the candidate
  as `maf_source = UNKNOWN` and conservatively filter it out at the default
  safeguard threshold
- filter candidates at `MAF <= top_hit_maf_threshold`

Implementation details:

- `DiffGWASDeps/TopHitMAF.pm`
  - shared GWAS / gnomAD MAF parsing helpers
- `DiffGWASDeps/gunplot/select_top_hits_from_wide.pl`
  - differential top-hit selector with GWAS-first and gnomAD-fallback MAF
    filtering
- `DiffGWASDeps/verify_common_association_loci.pl`
  - common-association selector with the same MAF policy
- `DiffGWASDeps/generate_requested_top_hits_csv.pl`
  - local helper now called by the SAS ODA local-top-hit wrappers before SAS
    submit so SAS reuses the same filtered requested-hit CSV instead of
    recomputing a different locus list internally
  - exported requested-top-hit CSVs now preserve
    `gwas_group1_maf`, `gwas_group2_maf`, `gwas_pair_maf_min`,
    `maf_filter_decision`, and `maf_filter_reason` alongside the selected SNP

Validation summary from the current schizophrenia project:

- the bundled regression helper
  `DiffGWASDeps/test_top_hit_maf_filter.pl` now checks:
  - GWAS-first differential filtering
  - gnomAD-fallback differential filtering
  - common-association filtering
- on 2026-06-12, the real PGC sex-differential lead SNP `rs185665940` passed
  with GWAS-derived `selected_maf = 0.013`
- on the same run, the real common-association selection retained `15` loci
  and every retained locus had `selected_maf > 0.01`
- the same real runner config did not define a gnomAD lookup file, so rows
  missing GWAS frequency fields were conservatively filtered as
  `maf_source = UNKNOWN` rather than being kept without an MAF check

Convenient regression commands:

```bash
perl DiffGWASDeps/test_top_hit_maf_filter.pl --no-real --keep-workdir
perl DiffGWASDeps/test_top_hit_maf_filter.pl --keep-workdir
```

## Cross-platform bootstrap

The repository now includes a repo-native installer layer so users do not need
to rely on personal `PERL5LIB` state or a preconfigured Python environment.

Shipped installer entry points:

- `install/install_windows_portable_cygwin.ps1`
- `install/install_windows_portable_cygwin.cmd`
- `install/install_cygwin.sh`
- `install/install_ubuntu.sh`
- `install/install_macos.sh`
- `Dockerfile`
- `install/singularity/MultiGWAS-Explorer_pipeline.def`
- `install/singularity/build_apptainer_image.sh`

Shared dependency manifests:

- `cpanfile`
- `install/requirements-pipeline.txt`

The installers provision:

- repo-local Perl modules under platform-specific trees such as
  `local/perl5-cygwin/`, `local/perl5-linux/`, and `local/perl5-darwin/`
- repo-local Python packages such as `saspy` and `Pillow` under
  `.venv-pipeline/`
- `bgzip` / `tabix` from the system when present, or a repo-local htslib build
  through `install/build_local_htslib.sh`
- native plotting prerequisites such as `gnuplot` and `ImageMagick`

Build-aware local GTF selection is now also part of the top-level contract:

- the spec JSON can carry:
  - `reference_build`
- both top-level entry scripts can take:
  - `--reference-build hg19|hg38|t2t`
- build resolution currently follows:
  - explicit override
  - header token detection such as `POS_HG38` or `BP_HG19`
  - filename/path token detection such as `hg19`, `grch38`, `hs1`, or `chm13`
  - fallback default `hg38`
- the built-in GTF profiles currently map to:
  - hg19 / GRCh37 via GENCODE v49 lift37
  - hg38 / GRCh38 via GENCODE v49
  - T2T / hs1 / CHM13v2.0 via the UCSC hs1 RefSeq GTF
- if the input headers and filenames do not carry trustworthy build tokens,
  set `reference_build` explicitly in the spec instead of relying on fallback

Post-install validation:

```bash
bash install/check_pipeline_install.sh
```

### Containerized deployment

The same repo-local Ubuntu installation path is also saved as a container build
definition for users who prefer Docker or Singularity / Apptainer.

Docker build and smoke test:

```bash
docker build -t multigwas-explorer-pipeline:latest .
docker run --rm -it multigwas-explorer-pipeline:latest \
  bash -lc "cd /opt/MultiGWAS-Explorer && bash install/check_pipeline_install.sh"
```

Convenience Ubuntu Docker gnuplot smoke test:

```bash
bash install/test_ubuntu_docker_gnuplot.sh
```

By default, that helper reproduces the top differential schizophrenia SNP
`rs185665940` with `local_manhattan,local_gtf`. Add `--include-manhattan` when
you also want the slower genome-wide gnuplot Manhattan panel in the same run.
After the 2026-06-08 layout refresh, that Docker smoke-test path also verified
that:

- genome-wide gunplot Manhattan now places the GWAS sublabels at the top of
  each subplot, which is closer to the SAS ODA genome-wide layout
- combined gunplot local Manhattan now defaults to a bottom gene-track view
  instead of the older vertical SNP/gene text labels
- users who still prefer the older compact label mode can restore it with:
  - `--local-manhattan-annotation labels`

Measured Ubuntu Docker timing on 2026-06-08 for this repository:

- first uncached `docker build` completed in about 4.2 minutes
- post-build `bash install/check_pipeline_install.sh` inside the image
  completed in about 10 seconds
- a one-SNP gunplot validation of
  `manhattan,local_manhattan,local_gtf` completed successfully in about
  9 minutes 43 seconds total, with `plot_manhattan` alone taking about
  9 minutes 10 seconds
- after rebuilding the image with the 2026-06-08 layout changes, focused smoke
  tests completed in about:
  - 13 seconds for `--plots local_manhattan`
  - 9 minutes 6 seconds for `--plots manhattan`
- in that run, the genome-wide gunplot Manhattan stage scanned 2,417,954 wide
  rows and wrote an intermediate `.plot.tsv` of about 110 MB, which explains
  why Docker validation can feel slow even when the container is healthy
- after the 2026-06-11 Manhattan style refresh, another rebuilt Ubuntu Docker
  run revalidated the schizophrenia genome-wide panel with the SAS-like
  repeated chromosome palette plus top-of-panel GWAS labels still in place
  while removing the extra default plot title from the PNG
- measured timing for that 2026-06-11 rerun was about:
  - 231 seconds for the image rebuild
  - 8 minutes 46 seconds for the focused `--plots manhattan` smoke stage

Example interactive Docker runtime with mounted data and SAS ODA authinfo:

```bash
docker run --rm -it \
  -e PIPELINE_WORKDIR=/opt/MultiGWAS-Explorer \
  -v /path/to/_authinfo:/root/_authinfo:ro \
  -v /path/to/gwas_drive_e:/mnt/e \
  -v /path/to/gwas_drive_g:/mnt/g \
  multigwas-explorer-pipeline:latest bash
```

Keep `PIPELINE_WORKDIR=/opt/MultiGWAS-Explorer` so the wrapper uses the Linux
runtime inside the image even if the host also mounts a Windows `G:` drive.
Avoid overlaying a host checkout onto `/opt/MultiGWAS-Explorer`, because that can
hide the image's Linux-installed `.venv-pipeline/` and repo-local Perl modules.
When the convenience smoke-test helper is launched from a Windows portable
Cygwin shell, it now converts the repository checkout path into a
Docker-friendly host build context automatically; this avoids the earlier
Docker Desktop error that reported a missing build context for
`/mnt/g/.../perlMCP4Gemini_Paper`.

Singularity / Apptainer build and smoke test:

```bash
bash install/singularity/build_apptainer_image.sh
apptainer exec MultiGWAS-Explorer_pipeline.sif \
  bash -lc "cd /opt/MultiGWAS-Explorer && bash install/check_pipeline_install.sh"
```

Important container note:

- the Dockerfile and the Apptainer definition both call the same
  `install/install_ubuntu.sh` installer, so the Linux container path stays
  aligned with the validated host-Ubuntu path instead of becoming a second,
  separate installation recipe
- for local-GTF workflows, keeping `cache/gtf/gencode.v49.annotation.gtf.gz`
  inside the repo or mounted into the container avoids treating a transient
  EBI download error as a pipeline failure
- in the current validation cycle, the Docker image build and the post-build
  runtime smoke test both succeeded
- the Singularity / Apptainer definition is saved for downstream users, but it
  was not executed on the current workstation in this validation pass

Ubuntu Docker fallback validation:

- when the bundled Vagrant harness could not be completed from a
  non-administrative Windows shell, the Ubuntu path was validated instead
  through Docker Desktop with an isolated `ubuntu:24.04` image
- that Docker fallback confirmed:
  - `bash install/install_ubuntu.sh`
  - `bash install/check_pipeline_install.sh`
  - the top-level gunplot wrapper for `manhattan`, `local_manhattan`, and
    `local_gtf`
- in the same 2026-06-08 Docker validation cycle, the containerized SAS ODA
  login probe did not complete successfully:
  - `perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl --check-sas-oda-login-only`
  - current observed failure: `SASIOConnectionTerminated: No SAS process attached. SAS process has terminated unexpectedly.`
- practical lesson from that Docker validation:
  - keep `cache/gtf/gencode.v49.annotation.gtf.gz` available locally so
    local-GTF runs do not depend on a fresh EBI download
  - on Ubuntu/Docker, SAS ODA housekeeping calls can outlast conservative
    host-side automation timeouts, so final local-GTF artifacts should be
    checked before concluding that the scientific plot itself failed
  - for the gunplot path, prefer a quick one-SNP `local_manhattan,local_gtf`
    validation before adding `manhattan`, because the genome-wide stage is the
    dominant runtime cost in Docker

Portable Cygwin wrapper validation:

- the top-level gunplot wrapper was validated from an isolated portable Cygwin
  repo copy under `H:\TMP4SAS\...`
- during that validation, the real Windows/Cygwin GD problem turned out to be
  mixed-platform repo-local Perl libraries rather than only missing native
  headers; the runtime now prefers `local/perl5-cygwin/` so portable Cygwin
  does not accidentally load Linux `GD.pm` or zlib XS modules
- the wrapper now prefers `gnuplot` from the active shell `PATH` before trying
  older Windows-specific fallback locations
- in a healthy portable-Cygwin run, you should see:
  - `Using gnuplot executable: gnuplot`
- practical validation command:

```bash
perl ./auto_prepare_and_run_diff_gwas_with_gunplot.pl \
  --spec configs/spec_pgc_scz_sex_common_automation.json \
  --plots manhattan,local_manhattan,local_gtf
```

Quick Ubuntu Docker gunplot validation:

```bash
perl ./auto_prepare_and_run_diff_gwas_with_gunplot.pl \
  --spec configs/spec_pgc_scz_sex_common_automation.json \
  --plots local_manhattan,local_gtf \
  --target-snps rs185665940
```

Quick gunplot forest validation:

```bash
perl ./auto_prepare_and_run_diff_gwas_with_gunplot.pl \
  --spec configs/spec_pgc_scz_sex_common_automation.json \
  --plots forest \
  --target-snps rs185665940
```

Genome-wide gunplot cache-safety note:

- during validation, an interrupted rerun left a corrupted cached
  genome-wide wide subset that made the gunplot Manhattan figure diverge
  sharply from the SAS ODA figure
- the gunplot wrapper now validates the cached wide subset against its
  manifest before reuse instead of trusting file existence alone
- `DiffGWASDeps/extract_significant_diff_gwas.pl` now writes the wide subset
  and manifest atomically, which reduces the chance of leaving a partially
  overwritten genome-wide input behind
- if a future gunplot genome-wide Manhattan plot looks structurally truncated
  or obviously wrong, rebuild the plotting subset with:

```bash
perl ./auto_prepare_and_run_diff_gwas.pl \
  --spec configs/spec_pgc_scz_sex_common_automation.json \
  --step extract_wide_subset \
  --force
```

Reference-build validation note:

- on 2026-06-11, temporary config-generation reruns were used to verify that:
  - `reference_build=hg19` emitted `REFERENCE_BUILD=hg19`,
    `GTF_DSD=FM.GTF_HG19`, and `GTF_LOCAL_DSD=gtf_hg19`
  - `reference_build=t2t` emitted `REFERENCE_BUILD=t2t`,
    `GTF_DSD=FM.GTF_T2T`, and `GTF_LOCAL_DSD=gtf_t2t`
- the bundled PGC schizophrenia sex spec is now pinned explicitly to `hg19`
  because the older DANER-style files do not expose a build token directly
- the bundled PGC schizophrenia ancestry inputs also resolve to `hg19`; their
  raw PGCsumstatsVCF headers report `##genomeReference="GRCh37"`

Genome-wide visual-parity note:

- the gunplot genome-wide Manhattan renderer now colors points by chromosome
  palette index rather than track index, which restores the same repeating
  chromosome-color logic used by the SAS ODA figure
- the gunplot PNG also now omits the extra default top title so the figure
  frame is closer to the manuscript SAS layout
- the main residual differences after the 2026-06-11 refresh were backend
  rasterization details such as point packing and antialiasing, not different
  chromosome ordering or mislabeled subplot headers

First-run SAS ODA authentication is now built into the vendored helper layer:

- if no saved SASPy authinfo entry exists for authkey `oda`, the first
  SAS-backed run prompts for the SAS ODA account/email and password
- the helper validates the supplied login by submitting:
  - `proc setinit;run;`
- failed credentials are rejected with a warning and are not kept
- successful credentials are saved into the SASPy authinfo file so later SAS
  ODA runs do not need to prompt again
- noninteractive bootstrap is also available with:
  - `--sas-oda-account`
  - `--sas-oda-password`
- a saved entry can be refreshed intentionally with:
  - `--prompt-sas-oda-auth`

The SAS ODA helper stack and the top-level entry scripts now prefer these
repo-local Perl and Python runtimes automatically when they exist.

### Recommended Windows bootstrap

For Windows, the preferred entry point is now the portable Cygwin bootstrap:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass `
  -File .\install\install_windows_portable_cygwin.ps1
```

This wrapper uses `MachinaCore/CygwinPortable` as the base runtime, downloads
or reuses the portable archive, refreshes the required Cygwin packages inside
an isolated portable root, writes an isolated `fstab`, and then runs the same
`install/install_cygwin.sh` phase-2 installer used inside a normal Cygwin
shell.

Default portable root:

- `H:\TMP4SAS\CygwinPortablePipeline`

If the user is already inside that portable shell, or already inside another
compatible Cygwin shell, the phase-2 installer can still be called directly:

```bash
bash install/install_cygwin.sh
```

During development, the Windows bootstrap was repeatedly exercised in isolated
`H:\TMP4SAS\...` directories so package refresh, repo-local Python setup, and
portable path handling could be checked without relying on the user's global
Cygwin tree. The gunplot wrapper now also prefers `gnuplot` from the active
portable shell `PATH`, which prevents the installed workflow from silently
borrowing a host-specific Windows `gnuplot.exe`. It also validates cached
genome-wide wide subsets against their manifest row counts before reuse, so a
stale interrupted rerun is less likely to poison later Manhattan plots.

The gunplot path now mirrors the same displayed-GWAS flexibility as the SAS ODA
path:

- `--display-gwas ALL` for one differential track
- `--display-gwas ALL_FEMALE` for one single-GWAS association track
- `--display-gwas EUR,EUR_FEMALE,EUR_MALE` for a mixed custom track set
- `--target-snps rs123` to render local Manhattan and local GTF context around
  an inquiry SNP with the selected displayed tracks

The same inquiry-SNP flexibility now also applies to forest plots:

- one inquiry SNP:
  - one forest plot with cohorts on the y-axis
- multiple inquiry SNPs or multi-hit exports:
  - one forest panel per displayed GWAS / cohort track
  - left y-axis = SNP ID
  - right y-axis = nearby gene label
  - styling intentionally aligned to the SAS ODA
    `beta2OR_forest_plot.sas` output

Example combined local-Manhattan gunplot rerun:

```bash
perl auto_prepare_and_run_diff_gwas_with_gunplot.pl \
  --spec ./configs/spec_pgc_scz_sex_common_automation.json \
  --target-snps rs17425819,rs185665940 \
  --plots local_manhattan \
  --local-manhattan-columns 2 \
  --local-manhattan-annotation gtf
```

Example single-GWAS plot set in the gunplot path:

```bash
perl auto_prepare_and_run_diff_gwas_with_gunplot.pl \
  --spec ./configs/spec_pgc_scz_sex_common_automation.json \
  --display-gwas ALL_FEMALE \
  --plots manhattan,local_manhattan,local_gtf
```

Example gunplot forest rerun for one target SNP:

```bash
perl auto_prepare_and_run_diff_gwas_with_gunplot.pl \
  --spec ./configs/spec_pgc_scz_sex_common_automation.json \
  --plots forest \
  --target-snps rs185665940
```

Example gunplot forest rerun for multiple target SNPs:

```bash
perl auto_prepare_and_run_diff_gwas_with_gunplot.pl \
  --spec ./configs/spec_pgc_scz_sex_common_automation.json \
  --plots forest \
  --target-snps rs185665940,rs4950119
```

Example inquiry-SNP local plot set with a custom displayed GWAS subset:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --spec ./configs/spec_pgc_scz_sex_common_automation.json \
  --display-gwas EUR,EUR_FEMALE,EUR_MALE \
  --target-snps rs185665940 \
  --step plot_local_manhattan \
  --step plot_local_gtf
```

### VS Code plus Codex usage across operating systems

The same repository can be driven from VS Code with Codex while preserving the
same local scripts and MCP server underneath.

Recommended pattern:

1. open the repository as the VS Code workspace;
2. finish the OS-specific installer first;
3. use an integrated terminal rooted in the repository;
4. start `server.pl` in one terminal;
5. register the MCP endpoint with:

```bash
codex mcp add perl-bio --url http://127.0.0.1:8080/mcp
codex mcp list
```

6. start a fresh Codex session in that same workspace; and
7. ask Codex to run either the SAS ODA or gunplot entry script.

Platform-specific terminal notes:

- Windows:
  - prefer `install/install_windows_portable_cygwin.ps1`
  - use the portable Cygwin `bash` terminal in VS Code
  - repository paths inside that shell should look like
    `/cygdrive/g/.../perlMCP4Gemini_Paper`
- macOS:
  - use the integrated `zsh` or `bash` terminal after
    `bash install/install_macos.sh`
- Ubuntu Linux:
  - use the integrated `bash` terminal after
    `bash install/install_ubuntu.sh`

This setup keeps Codex, the local MCP server, and the repo checkout inside one
editor workspace instead of splitting execution across unrelated shells.

Example noninteractive first-run SAS ODA login bootstrap:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --spec ./configs/spec_pgc_scz_sex_common_automation.json \
  --step plot_manhattan \
  --sas-oda-account your_email@example.com \
  --sas-oda-password 'your_password'
```

Important gunplot note:

- the default combined local-Manhattan annotation mode is now `gtf`, so the
  final combined figure prefers a bottom gene track instead of the older
  vertical SNP/gene label pair
- if you prefer the previous compact label mode, override it with:
  - `--local-manhattan-annotation labels`
- for combined local Manhattan with `--local-manhattan-annotation gtf`, the
  renderer now combines association rows across loci, combines the matching
  GTF rows, tags both by top-SNP membership, and rescales genomic position into
  one shared x coordinate system before plotting
- the key debug artifact for that mode is:
  - `*.combined_scaled.tsv`
- that unified table is now the preferred source of truth for the combined
  gunplot local-Manhattan GTF path

Both the SAS ODA and gunplot workflows now also support optional
user-designated adjacent-gene labels for explicit target SNPs.

Example:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --spec ./configs/spec_pgc_scz_sex_common_automation.json \
  --target-snps rs17425819,rs185665940,rs2564978 \
  --target-snp-genes rs17425819:JAK2,rs185665940:FANCL,rs2564978:CR1 \
  --step plot_local_manhattan
```

You can also store the same mapping in the spec JSON as `target_snp_genes`.
These user-supplied labels now override the automatic HaploReg / nearest-GTF
fallback and propagate into the exported SAS-style top-hit CSV with
`gene_source=USER`.

## PGC Sex-Stratified chrX Troubleshooting

For the PGC schizophrenia sex-stratified example, absence of chrX can reflect
either the original source-file design or a later downstream artifact.

Expected source-file design in this repository:

- `ALL_FEMALE` and `ALL_MALE` are autosome-only
- `EUR_FEMALE`, `EUR_MALE`, `ASN_FEMALE`, and `ASN_MALE` each add a separate
  chrX supplement file

Repository evidence:

- [configs/spec_pgc_scz_sex_common_automation.json](/G:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/configs/spec_pgc_scz_sex_common_automation.json)
  shows `ALL_*` using autosomal files only, while `EUR_*` and `ASN_*` each use
  an additional `chrX` file
- [merge_scz_sex_stratified_long.pl](/G:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/merge_scz_sex_stratified_long.pl)
  encodes the same logic with `*_AUTOSOME` tags for the pooled strata
- [PGC_SCZ_female_vs_male_diff_effects_merged_long.manifest.tsv](</E:/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_SCZ_Sex_Stratified_GWASs/PGC_SCZ_female_vs_male_diff_effects_merged_long.manifest.tsv>)
  confirms which raw files were actually merged

Observed raw-file chrX counts:

- `all_female.gz`: `0`
- `all_male.gz`: `0`
- `chrX_eur_fem`: `231,908`
- `chrX_eur_mal`: `232,899`
- `chrX_asn_fem`: `193,305`
- `chrX_asn_mal`: `195,213`

Practical interpretation:

- missing chrX in `ALL_*` is expected from the input bundle
- missing chrX in `EUR_*` or `ASN_*` should trigger a downstream integrity
  check

Important downstream failure mode already observed:

- the sorted merged long file still contained EUR/ASN chrX rows
- but a later `PGC_SCZ_female_vs_male_diff_effects.tsv.gz` on disk had fewer
  rows than its manifest claimed and no chrX rows
- its manifest timestamp and data-file timestamp no longer matched, which
  strongly suggested the differential output had been partially regenerated or
  overwritten after the manifest was written

Recommended recovery sequence when that pattern appears:

1. trust the raw-file composition and merged-long manifest first
2. verify chrX presence in the sorted merged long table
3. compare the differential output row count against its manifest
4. if inconsistent, regenerate:
   - `*_diff_effects.tsv.gz`
   - `*.stdized.tsv.gz`
   - `*.wide_beta_se_p_p_lt_0p05*.tsv.gz`
5. then rerun SAS ODA or gunplot figures

Plotting reminder:

- the gunplot workflow now removes chrX from final figures by default unless
  `--no-remove-X-chr` is set
- that behavior does not mean chrX was absent upstream

The main tested ancestry example spec is:

- [configs/spec_pgc_scz_ancestry_diff_automation.json](/G:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/configs/spec_pgc_scz_ancestry_diff_automation.json)

Example:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --spec ./configs/spec_pgc_scz_ancestry_diff_automation.json
```

To only generate configs:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --spec ./configs/spec_pgc_scz_ancestry_diff_automation.json \
  --mode configs
```

To rerun only selected plot stages:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --spec ./configs/spec_pgc_scz_ancestry_diff_automation.json \
  --plots local_manhattan,local_gtf
```

For SAS local-Manhattan target-SNP reruns, the wrapper is now optimized to
avoid uploading the full wide subset when only a few loci are needed:

- it builds per-target local windows directly from `SOURCE_LONG_GZ`
  with `DiffGWASDeps/extract_single_snp_wide_diff_gwas.pl`
- it unions those windows into one compact wide subset
- it caches the result under `cache/local_manhattan_reuse/`
- when available, it uploads that compact subset instead of the much larger
  global wide plot input

This makes local-Manhattan reruns much less sensitive to interrupted or flaky
large-file uploads to SAS ODA.

To keep the default protein-coding-focused local GTF bottom track, no extra
flag is needed. Local GTF bottom tracks are protein-coding-only by default:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --spec ./configs/spec_pgc_scz_ancestry_diff_automation.json \
  --step plot_local_gtf
```

To explicitly include non-coding genes again, set
`"include_non_protein_coding_genes_in_local_gtf": 1` in the spec JSON.

When a local GTF locus has many overlapping genes, the pipeline now keeps the
overall figure width and height fixed and slightly increases the lower
gene-track share by auto-tuning the SAS `pct4neg_y` ratio.

The current SAS local-GTF wrappers now start from the same larger lower-track
base before that auto-tuning path is applied:

- `DiffGWASDeps/run_sas_oda_local_top_hits_with_gtf_download_html.sh`
- `DiffGWASDeps/run_sas_oda_single_snp_with_gtf_download_html.sh`

Both wrappers now default `GTF_PCT4NEG_Y` to `1.4`, which makes the bottom
gene/exon track more readable in manuscript-scale figures and single-locus
debug reruns.

## MCP wrapper

The same automation script is now exposed through the local Perl MCP server as:

- `auto_prepare_and_run_diff_gwas`

That MCP tool accepts:

- `spec_file`
- `mode`
- `plots`
- `skip_plots`
- `force`
- `output_file`
- `pid`

Recommended MCP pattern:

1. call it once with `spec_file`
2. save the returned PID
3. poll the same tool with `spec_file` plus `pid`

Tested through MCP with:

- `./configs/spec_pgc_scz_ancestry_diff_automation.json`

Verified MCP outcomes:

- `mode=configs` completed correctly
- full mode with `skip_plots=true` completed correctly and reused the
  ancestry artifacts already on disk

Implementation notes learned during MCP testing:

- the wrapper must launch `auto_prepare_and_run_diff_gwas.pl` by explicit path,
  not `perl -S`, because the daemon PATH does not reliably include the current
  directory
- when a custom `output_file` is supplied, the wrapper must save that path with
  the PID so later status polling reads the correct log
- an interrupted `force=true` ancestry rerun can leave
  `PGC_SCZ_ancestry_diff_effects_merged_long.manifest.tsv` at zero bytes,
  which causes later runs to think the merge stage must rerun; restore a valid
  manifest or rerun the merge stage before trusting reuse behavior

## Shared ODA data upload

When one automation call requests multiple plot stages that use the same wide
gz subset, the pipeline now uploads that data file to SAS ODA once at the
automation layer and reuses it across the requested plot runners.

The implemented pattern is:

1. upload shared wide gz once
2. run downstream plot wrappers with:
   - `SKIP_DATA_UPLOAD=1`
   - `CLEAN_ODA_INPUT=0`
3. after the last plot stage, delete the shared remote gz once if cleanup is
   enabled

This avoids paying the large gz upload cost multiple times in one command such
as:

```bash
perl auto_prepare_and_run_diff_gwas.pl --spec x.json --plots manhattan,local_gtf
```

For repeated reruns of the same comparison, add this to the spec:

```json
"keep_remote_plot_data": 1
```

That tells the automation to:

- keep the uploaded wide gz in SAS ODA after plotting
- check whether the remote gz already exists on later runs
- skip re-upload when it is already present

The individual shell runners still support standalone behavior; when run by
themselves they upload the gz normally unless `SKIP_DATA_UPLOAD=1` is supplied.

## SAS ODA helper improvements

The shared low-level wrapper:

- `run_sas_codes_or_script_in_ODA.pl`

now supports bulk remote file management in one call:

- repeated `--upload-file`
- repeated `--download-file`
- repeated `--download-local-path`
- repeated `--delete-file`
- repeated `--file-info`
- repeated `--delete-file-rgx`

Practical remote-path lessons from the recent SAS ODA debugging:

- quote remote home-directory paths like `'~/plot.png'` in shell examples so
  the literal `~` reaches the helper unchanged
- the helper now normalizes both `~/...` and absolute SAS home paths such as
  `/home/...` consistently for download, delete, and file-info operations
- delete requests now verify that the target path no longer resolves after the
  helper reports success

Regex-based deletion first lists the selected remote directory and then matches
each entry against both:

- the basename, such as `plot.png`
- the resolved remote path, such as `~/plot.png`

That means either of these styles works:

```bash
perl -S run_sas_codes_or_script_in_ODA.pl --delete-file-rgx '.*\.png$'
perl -S run_sas_codes_or_script_in_ODA.pl --delete-file-rgx '~\/.*\.png'
```

This is especially useful when cleaning up batches of uploaded PNG, HTML, SAS,
or subset files after repeated SAS ODA plot reruns.

The same low-level helper is now more defensive for `%include` debugging too.
When a submitted SAS code block or SAS file contains `%include`, the helper can:

- auto-refresh a remote include target such as `~/Lattice_gscatter_over_bed_track.sas`
  from a matching local project copy before the real submit
- run a local preflight scan for likely compile blockers such as malformed
  nested block comments or unterminated quotes
- record line-numbered source context for flagged lines in the saved
  `output.html.info.txt`
- attempt a remote `PROC PRINTTO` compile log for the include target
- fail fast with an explicit error when the include looks broken, instead of
  waiting for a long SAS ODA timeout and returning an empty log

The helper timeout path itself was also hardened after a real local-GTF replay
debugging session:

- older timeout handling forked a live Perl process that had already loaded
  `Inline::Python` / `saspy`
- in that state, the child could return an empty SAS result payload even for a
  simple `%put HELLO;`
- the timeout path now launches a fresh worker process to run the real
  `SAS_ODA_Runner` submit and serialize the result back to the parent
- if a preserved rendered SAS script succeeds through direct
  `SAS_ODA_Runner->run_file(...)` replay but the helper reports empty log/output,
  treat that as a helper transport bug first, not automatically as a macro bug

This is particularly helpful when the local macro has already been corrected
but SAS ODA is still using a stale uploaded copy in the remote home directory.
The helper now tries to refresh that remote copy first when it can resolve a
matching local file path.

Scientific lesson from the merged-AOA genome-wide Manhattan rerun:

- compact custom SAS plotting subsets must already be coordinate-safe before
  plotting
- strip an optional `chr` prefix, map `X -> 23`, map `Y -> 24`, and drop rows
  whose chromosome or base-pair position still cannot be parsed
- sort on the normalized chromosome/position fields that SAS will actually
  read, not on the raw chromosome text
- otherwise a literal `CHR='X'` can be imported by SAS as numeric missing `.`
  and sort before `chr1`, which creates a false extra chromosome block at the
  far left of the genome-wide Manhattan figure
- for custom merged study tables, prefer a compact upload with only the
  required plotting columns and apply the `CHR_RAW -> CHR` normalization
  either upstream or explicitly in the SAS import data step

Merged-wide GWAS quickstart for AOA:

- the generalized automation entrypoint now recognizes a single merged study
  table as `source_mode=merged_gwas_table`
- this is intended for tables with one shared locus block such as
  `CHR/BP/SNP/A1/A2` and repeated cohort association blocks such as
  `BETA_DS_ALL/SE_DS_ALL/P_DS_ALL` and
  `BETA_MP2PRT/SE_MP2PRT/P_MP2PRT`
- optional extra tracks such as meta-analysis `P` / `Z` columns can stay in
  the same merged input and will be carried through when configured

Preview the inferred spec without writing it:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --gwas-dir AOA_GWAS_Data/ \
  --preview-spec \
  --generate-spec-only
```

Write the inferred spec:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --gwas-dir AOA_GWAS_Data/ \
  --spec-out configs/auto_aoa_merged.spec.json \
  --generate-spec-only
```

Run the merged-wide AOA plotting workflow:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --gwas-dir AOA_GWAS_Data/ \
  --spec-out configs/auto_aoa_merged.spec.json \
  --plots manhattan,local_manhattan,local_gtf,forest
```

Merged-wide implementation notes:

- the directory scan now ignores generated artifacts such as
  `*.merged_plotwide.tsv.gz`, `gunplot`, `png`, `html`, and manifest files so
  reruns continue to target the original merged input
- the normalization path uses:
  - `DiffGWASDeps/convert_merged_gwas_to_plotwide.pl`
  - `DiffGWASDeps/generate_sas_wide_import_include.pl`
- when the default auto top-hit threshold yields no retained loci, prefer an
  explicit rerun with `--target-snps` or a relaxed threshold instead of
  assuming the merged table was mis-parsed

For `%include`-driven debugging, the helper behavior is now split more cleanly:

- it can run a resolved included file as a standalone debug submit before the
  parent script
- it keeps helper `%include` diagnostics under `=== Dependency Logs ===`
  instead of appending bracketed helper banners to the saved
  `=== Submitted SAS Codes or file ===` section

One concrete SAS macro debugging lesson from this work:

- `*%abort 255;` is not a safe way to disable a macro abort inside SAS macro
  code
- that `* ... ;` form is a statement comment, and the macro processor can still
  execute `%abort`
- use `/*%abort 255;*/` or `%* %abort 255;` instead

Why this matters for SAS ODA:

- a hidden `%abort` inside an included macro can terminate the ODA run before
  saspy returns a usable inline log
- in that failure mode, the saved `output.html.info.txt` can contain only:
  - include-preflight notes
  - an empty `=== SAS Log ===`
  - an empty `=== Output ===`
- treat that pattern as a control-plane failure or early compile/runtime abort,
  not as evidence that the plotting script completed

Recommended debugging response:

1. inspect included macros for macro statements that were "commented out" with
   `*...;`, especially `%abort`
2. rerun the failing `%include` through `run_sas_codes_or_script_in_ODA.pl`
   so the helper can refresh stale remote copies and save preflight findings
3. rerun the pipeline with `KEEP_RENDERED_DEBUG_FILES=1` so the exact rendered
   SAS script and local input subsets are preserved for replay

Recommended `run_sas_codes_or_script_in_ODA.pl` file-handling debug process:

1. validate the saved SAS ODA login first:

```bash
perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl \
  --check-sas-oda-login-only
```

2. if the issue is listing or remote file existence, stay in file-management
   mode and probe that path before any SAS code submit:

```bash
perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl \
  --output-prefix oda_file_probe \
  --dir4listing '~' \
  --file-info '~/some_remote_file.txt'
```

3. if upload/download behavior looks suspicious, test transfer-only mode with a
   small file and an explicit output prefix:

```bash
perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl \
  --output-prefix oda_transfer_probe \
  --upload-file ./small_test.txt \
  --download-file '~/small_test.txt' \
  --download-local-path ./small_test.roundtrip.txt
```

4. if stale remote files may be confusing the rerun, explicitly check both the
   basename and the resolved remote path, then remove only the exact target you
   intend to refresh:

```bash
perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl \
  --output-prefix oda_stale_remote_probe \
  --file-info 'Lattice_gscatter_over_bed_track.sas' \
  --file-info '~/Lattice_gscatter_over_bed_track.sas'
```

5. only after listing, file-info, upload, and download succeed should you move
   on to SAS code execution probes:

```bash
perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl \
  --output-prefix oda_code_probe \
  --code "%put HELLO_FROM_ODA;"
```

6. if the real failure involves `%include`, rerun the include target or the
   parent code with a short timeout and a dedicated output prefix so the debug
   bundle is easy to find:

```bash
SAS_ODA_RUN_TIMEOUT_SECONDS=60 \
SAS_ODA_RUN_TIMEOUT_GRACE_SECONDS=10 \
perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl \
  --output-prefix include_debug_probe \
  --code "%include '~/Lattice_gscatter_over_bed_track.sas';"
```

7. inspect the generated debug bundle, especially:
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
- `SAS_ODA_RUN_TIMEOUT_SECONDS` and `SAS_ODA_RUN_TIMEOUT_GRACE_SECONDS`
  - keep failed probes short and reproducible

Quick AI Debug Ladder for `run_sas_codes_or_script_in_ODA.pl`
-------------------------------------------------------------

When debugging this helper in the future, keep the probes in this order so AI
does not jump directly into plotting macros:

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

3. transfer-only probe:

```bash
perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl \
  --output-prefix oda_transfer_probe \
  --upload-file ./small_test.txt \
  --download-file '~/small_test.txt' \
  --download-local-path ./small_test.roundtrip.txt
```

4. smallest SAS submit with macro autoload skipped automatically:

```bash
perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl \
  --output-prefix oda_code_probe \
  --code "%put HELLO_FROM_ODA;"
```

5. smallest HTML-producing SAS submit:

```bash
cat > codex_proc_print_smoke.sas <<'EOF'
proc print data=sashelp.class;
run;
EOF

perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl \
  --output-prefix oda_proc_print_probe \
  --file codex_proc_print_smoke.sas
```

6. forced default `~/Macros` bootstrap probe:

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

7. manual in-SAS `~/Macros` timing probe:

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

New helper artifacts to inspect:

- `<output-prefix>/output.run.status.json`
  - live sidecar with fields such as `bootstrap_started_at`,
    `bootstrap_finished_at`, `bootstrap_elapsed_seconds`, `bootstrap_ok`, and
    `bootstrap_log_path`
- `<output-prefix>/output.macro_bootstrap.log.txt`
  - bootstrap-only log/trace file for the default `~/Macros` autoload path
  - created as soon as bootstrap starts, so timeout cases still leave behind a
    partial trace
- `<output-prefix>/output.html.info.txt`
  - now includes `=== Status Snapshot ===` in addition to dependency logs,
    SAS log, and output summary

Interpret the progress lines carefully:

- `Upload step: macro bootstrap helper: importallmacros_ue.sas ...`
  - only the small helper file is being uploaded
- `SAS ODA macro bootstrap started at ...`
  - the actual bootstrap submit has begun
- if `output.macro_bootstrap.log.txt` still shows `Status: running` with only
  `Bootstrap Start:` populated after timeout, the wrapper entered bootstrap but
  never got a completed `sess.submit(...)` return

Windows note:

- when debugging from PowerShell into Cygwin, prefer `--file` for anything
  beyond the trivial `%put HELLO_FROM_ODA;` probe because inline `--code`
  quoting can be silently mangled before SAS ever sees it

One more operational rule from the local-GTF reruns:

- the outer automation command can time out after the local plot stage has
  already completed and written fresh HTML / CSV artifacts
- before treating that as a plotting failure, inspect:
  - the newest `run_local_hits_with_gtf_*/output.html.info.txt`
  - the downloaded final HTML or PNG
  - the selected-hit CSV
- if the saved log contains repeated completion markers such as
  `The final figure is put here` and `Lattice gscatter plot is completed!`,
  then the local GTF rendering itself succeeded and the remaining delay is
  likely in a later automation step or wrapper cleanup path

Critical `WORK`-overflow lesson from the larger local-GTF reruns:

- if the saved SAS log reports `ERROR: Insufficient space in file WORK.FINAL.DATA.`
  or related `_DOCTMP...` / insufficient-space / damaged-dataset failures,
  treat that as the primary root cause even if the log later prints
  figure-completion markers
- after that failure, SAS can still emit a remote PNG path and
  `Lattice gscatter plot is completed!`, but the figure may contain blank
  regions because the intermediate plotting datasets were already corrupted
- the first mitigations should be:
  - reduce the local GTF window size
  - reduce local GTF loci per figure or per run
  - reduce stacked `GTF_ASSOC_PVARS` track count for exploratory reruns

## Plotting specification

Current preferred plotting rule for differential-comparison workflows:

- genome-wide Manhattan plots:
  - include standardized differential P for all pairwise comparisons
  - include raw association P for each individual GWAS
  - do not prioritize raw differential P as a default displayed track

- local top-hit Manhattan plots:

## Local GTF rendering details

The local GTF path now prefers the project copy of:

- `DiffGWASDeps/SNP_Local_Manhattan_With_GTF.sas`

instead of assuming the correct macro is already present inside SAS ODA.
That macro is uploaded together with the patched lattice helper so local GTF
reruns use a versioned plotting stack from this repository.

For the local GTF display window, the pipeline now distinguishes between:

- the larger GTF search/extraction region
- the final displayed association-signal x-axis span

This matters for large windows such as `1e8`, where gene-track expansion can
cross the left chromosome boundary. The patched local macro now forces the
final displayed x-axis back to the min/max signal positions from the uploaded
GWAS subset, rather than leaving the final chr axis anchored at `0`.

The shell wrapper for local top-hit GTF plots is also more recovery-oriented:

- it defaults to one-shot ODA submits because that path has been more reliable
  for large local GTF jobs
- it retries once when the first submit returns an incomplete control-plane log
- if SAS already produced results remotely, it can reuse the helper-saved
  `sas_res_*.html` artifact and separately download the final PNG path reported
  in the SAS log
- when that PNG is present, the wrapper now writes a small figure-first final
  HTML that embeds the completed PNG and keeps the recovered/raw SAS HTML as a
  sidecar `*.sasraw.html`, so the user-opened result is the informative plot
  - follow the same track mix as the genome-wide Manhattan plot when driven by
    the generalized runner config
  - now also emit a CSV file containing the exact top hits selected for those
    local panels
  - the CSV is ordered by plotted top-hit rank and starts with:
    `hit_order`, `panel_index`, `CHR`, `BP`, `SNP`, `EFFECT_ALLELE`,
    `OTHER_ALLELE`, `REFERENCE_ALLELE`, `ALTERNATIVE_ALLELE`, `gene`,
    `snp_gene`, the focus signal used for ranking, `selected_maf`,
    `maf_source`, `gwas_group1_maf`, `gwas_group2_maf`,
    `gwas_pair_maf_min`, `gnomad_maf`, `gnomad_pops`,
    `maf_filter_decision`, and `maf_filter_reason`
  - after those leading columns, it appends the remaining wide-table signal
    columns for the selected hits, so users can inspect per-GWAS association
    P-values as well as the pairwise differential and standardized-differential
    signals in one shareable CSV
  - local top-hit selection now supports two modes:
    - differential mode: choose loci from significant pairwise differential
      signals
    - common-association mode: choose loci with strong single-GWAS association
      and at least one nominal association in another GWAS with the same
      effect direction
  - local panels now accept up to `30` columns per figure through either:
    - CLI: `--local-max-hits-per-fig 30`
    - spec: `"local_max_hits_per_fig": 30`
  - difficult SNP-gene label layouts can now be adjusted manually through:
    - CLI:
      - `--local-manhattan-angle4xaxis-label`
      - `--local-manhattan-xgrp-y-pos`
      - `--local-manhattan-yoffset-top`
      - `--local-manhattan-yoffset-bottom`
      - `--local-manhattan-fontsize`
      - `--local-manhattan-y-axis-label-size`
      - `--local-manhattan-y-axis-value-size`
    - spec:
      - `local_manhattan_angle4xaxis_label`
      - `local_manhattan_xgrp_y_pos`
      - `local_manhattan_yoffset_top`
      - `local_manhattan_yoffset_bottom`
      - `local_manhattan_fontsize`
      - `local_manhattan_y_axis_label_size`
      - `local_manhattan_y_axis_value_size`
  - when HaploReg does not provide a usable nearby gene for a plotted top hit,
    the local Manhattan path now falls back to a region-limited Gencode subset
    and assigns the nearest adjacent gene before plotting

The same SAS ODA timeout hardening now applies across the sibling plotting
wrappers too:

- `DiffGWASDeps/run_sas_oda_local_top_hits_manhattan_download_png.sh`
- `DiffGWASDeps/run_sas_oda_single_snp_with_gtf_download_html.sh`

That shared hardening includes:

- a hard timeout around the low-level ODA helper calls used for upload,
  download, delete, directory listing, and file-info checks
- a hard timeout around the actual SAS submit stage
- a small retry loop when the helper log is clearly incomplete or shows a
  known SAS ODA transport failure

Relevant environment controls are now:

- shared helper layer:
  - `ODA_HELPER_TIMEOUT_SECONDS`
  - `ODA_HELPER_TIMEOUT_GRACE_SECONDS`
- low-level submit helper:
  - `SAS_ODA_RUN_TIMEOUT_SECONDS`
  - `SAS_ODA_RUN_TIMEOUT_GRACE_SECONDS`
- local top-hit Manhattan:
  - `LOCAL_MH_SUBMIT_TIMEOUT_SECONDS`
  - `LOCAL_MH_SUBMIT_TIMEOUT_GRACE_SECONDS`
  - `LOCAL_MH_SUBMIT_MAX_ATTEMPTS`
- single-SNP local GTF:
  - `SINGLE_SNP_GTF_SUBMIT_TIMEOUT_SECONDS`
  - `SINGLE_SNP_GTF_SUBMIT_TIMEOUT_GRACE_SECONDS`
  - `SINGLE_SNP_GTF_SUBMIT_MAX_ATTEMPTS`
- local top-hit GTF:
  - `GTF_SUBMIT_TIMEOUT_SECONDS`
  - `GTF_SUBMIT_TIMEOUT_GRACE_SECONDS`
  - `GTF_SUBMIT_MAX_ATTEMPTS`

- local GTF plots:
  - build a region-limited Gencode subset around the selected locus or loci
  - keep the bottom track protein-coding-focused by default
  - allow non-protein-coding inclusion through:
    - spec: `"include_non_protein_coding_genes_in_local_gtf": 1`
  - if you want to enforce the default explicitly from the CLI, use:
    - `--exclude-non-protein-coding-genes-in-local-gtf`
  - use a safer default per-figure or per-run batching rule than local
    Manhattan:
    - `local_max_hits_per_fig` continues to control local Manhattan batching
    - spec: `"local_gtf_max_hits_per_fig": 1`
    - CLI: `--local-max-hits-per-fig N` still overrides both local Manhattan
      and local GTF when you explicitly request it
  - if `GTF_DESIGN_HEIGHT` is left unset, the local GTF wrapper now raises the
    default figure height modestly from the number of stacked
    `GTF_ASSOC_PVARS` tracks, which helps multi-track plots avoid cramped
    panels before you need manual tuning
  - for the current shared SAS defaults, both the single-SNP and batched
    local-GTF wrappers now also begin from a practical baseline height of
    `1000` pixels before any explicit user override

For local GTF top-SNP labels, document both the label-position control and the
headroom-height control from the lattice macro:

- `Yoffset4textlabels=2.5`
  - adjusts the SNP-label position within the top headroom
  - acts in y-axis-value units rather than figure-fraction units
  - is often sufficient as a starting value, but may still be changed
    internally by the macro auto-tuning path
- `yoffset4max_drawmarkersontop=0.15`
  - controls the top headroom height when SNP labels are drawn above the
    scatter tracks
  - replaces the ordinary `yaxis_offset4max` path in that top-label mode
  - is then adjusted internally again from scatter-track count before the
    final `offsetmax` value is assigned to the y-axis

Practical debugging rule:

- if there is only one SNP label on top and changing `Yoffset4textlabels`
  barely moves it, inspect the single-label branch inside
  `Lattice_gscatter_over_bed_track.sas`
- the fallback manual control is the ratio around `1000/&track_height`, which
  can be reduced or increased to move that single top label lower or higher in
  the reserved headroom

For target-SNP extraction, the pipeline can now take advantage of a tabix index
on the standardized long differential GWAS:

- `DiffGWASDeps/standardize_diff_gwas_zscore.pl` now writes bgzip-compatible
  output and creates `output.tbi` automatically when local `bgzip/tabix` are
  available
- `DiffGWASDeps/extract_single_snp_wide_diff_gwas.pl` still has a streaming
  fallback, but when the input has a `.tbi` it now uses `tabix` for the
  target-window extraction pass
- if the caller already knows the target coordinate, the extractor also accepts:
  - `--target-chr`
  - `--target-bp`
  which skips the initial full-file SNP-location scan too

For common-association mode, the automation entrypoint now supports:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --spec your_spec.json \
  --step plot_local_manhattan \
  --get-common-associations
```

You can also provide a starting threshold directly:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --spec your_spec.json \
  --step plot_local_manhattan \
  --get-common-associations=5e-8
```

The common-association threshold ladder defaults to:

- `5e-8`
- `1e-6`
- `1e-5`

The differential top-hit selector now defaults to:

- `1e-6`
- fallback `1e-5` only when no MAF-passing differential loci survive `1e-6`

If you provide an explicit multi-threshold ladder yourself, that explicit
ladder is respected as-is.

Important details for common-association mode:

- the `THR` in `--get-common-associations=THR` is a single-GWAS association
  P-value threshold, not a differential-association threshold
- the pipeline uses the minimum per-GWAS association P across the stacked raw
  GWAS tracks as `COMMON_ASSOC_P`
- `top_hit_dist_bp` is the total exclusion span. The selector keeps one lead
  signal within `BP +/- 0.5 * top_hit_dist_bp`, so `1e6` means one lead per
  1 Mb span, not a 1 Mb half-window on both sides
- the default common-association pruning span is now `1e6` unless overridden
  in the spec or runner config
- the local common-hit wrappers no longer force `--max-hits 15`; use
  `top_hit_max_loci` / `TOP_HIT_MAX_LOCI` when you want a deliberate cap,
  and leave it at `0` for an uncapped export
- a locus is then retained only if:
  - one GWAS provides that strongest association signal
  - another GWAS shows at least nominal association, default `< 0.05`
  - the paired GWASs point in the same effect direction, checked through the
    matched group-level Z-score signs in the SAS plotting path and matched
    group-level beta signs in the independent Perl verifier

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --spec your_spec.json \
  --step plot_local_manhattan \
  --get-common-associations=5e-8
```

Example manual local-Manhattan rerun for label tuning:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --spec your_spec.json \
  --step plot_local_manhattan \
  --local-max-hits-per-fig 30 \
  --local-manhattan-xgrp-y-pos -2.5 \
  --local-manhattan-yoffset-top 14 \
  --local-manhattan-yoffset-bottom 0.5 \
  --local-manhattan-angle4xaxis-label 60 \
  --local-manhattan-fontsize 3.0
```

Performance note for many common top hits:

- the local top-hit SAS runners now choose the first usable threshold with a
  cheap candidate count and only then call `get_top_signal_within_dist`
- `DiffGWASDeps/get_top_signal_within_dist.sas` now uses a grouped greedy
  interval-selection pass instead of the older broad self-join pattern, which
  materially reduces time when many loci pass the initial association filter
- for common-association local Manhattan / local GTF runs, explicitly upload
  and `%include "~/get_top_signal_within_dist.sas";` from `DiffGWASDeps`
  inside the submitted SAS script. Otherwise SAS ODA can silently reuse an
  older home-directory macro and fall back into the old `MAP_DUPS2TOPS` /
  insufficient-WORK-space path
- inside that local macro, use the plain data-step overlap check
  `(&pos_var >= _sel_dis_st) and (&pos_var <= _sel_dis_end)`;
  `between` is SQL syntax and will break the data-step greedy selector

and relaxes automatically until at least one eligible top hit is found.

Performance note:

- local top-hit runners no longer call the expensive distance-dedup macro once
  per threshold in the ladder
- they first choose the threshold cheaply by counting qualifying candidates, and
  then call `get_top_signal_within_dist` once with the chosen threshold

- local GTF-backed plots:
  - include standardized differential P for all pairwise comparisons
  - include raw association P for each individual GWAS
  - exclude raw differential P tracks
  - pair those association tracks with matching Z-score variables:
    standardized differential Z for pairwise tracks and raw GWAS Z for
    individual-GWAS tracks

For ancestry-style comparisons where the same GWAS appears in multiple pairs,
prefer aliasing one representative association track per unique GWAS, for
example `ASN_P`, `AFR_P`, `EUR_P`, rather than plotting duplicated pair-local
copies of the same association track.

## Example commands

Wide subset extraction from a config:

```bash
perl DiffGWASDeps/extract_significant_diff_gwas.pl \
  --config ./configs/preset_pgc_scz_sex_diff.json
```

Override just the threshold while keeping the preset:

```bash
perl DiffGWASDeps/extract_significant_diff_gwas.pl \
  --config ./configs/preset_pgc_scz_sex_diff.json \
  --threshold 5e-8
```

Single-locus extraction from a config:

```bash
perl DiffGWASDeps/extract_single_snp_wide_diff_gwas.pl \
  --config ./configs/preset_pgc_scz_sex_diff.json \
  --target-snp rs17425819
```

An ancestry-oriented example:

```bash
perl DiffGWASDeps/extract_significant_diff_gwas.pl \
  --config ./configs/preset_generic_ancestry_diff.json \
  --input /mnt/e/path/to/SCZ_ancestry_diff.stdized.tsv.gz \
  --output /mnt/e/path/to/SCZ_ancestry_diff.wide_subset.final.tsv.gz \
  --manifest /mnt/e/path/to/SCZ_ancestry_diff.wide_subset.final.manifest.tsv
```

## Design rule

Keep the reusable science logic in:

- config files
- versioned Perl/SAS/Bash scripts
- workflow notes/skills

Keep the Perl MCP server limited to stable execution primitives such as:

- run Perl or bash
- upload/download/delete/list SAS ODA files
- run SAS ODA code

That split keeps future ancestry, cohort, platform, or study-version comparisons much easier to support.

## Runner config layer

The SAS ODA shell wrappers now also support uppercase runner presets through:

- [configs/runner_pgc_scz_sex_diff.json](/G:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/configs/runner_pgc_scz_sex_diff.json)
- [configs/runner_generic_ancestry_diff.json](/G:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/configs/runner_generic_ancestry_diff.json)
- [emit_diff_gwas_runner_env.pl](/G:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/DiffGWASDeps/emit_diff_gwas_runner_env.pl)
- [generate_sas_wide_import_include.pl](/G:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/DiffGWASDeps/generate_sas_wide_import_include.pl)
- [render_sas_template.pl](/G:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/DiffGWASDeps/render_sas_template.pl)

Use them like this:

```bash
RUNNER_CONFIG_JSON=./configs/runner_pgc_scz_sex_diff.json ./run_sas_oda_manhattan4diffgwas_download_png.sh
RUNNER_CONFIG_JSON=./configs/runner_pgc_scz_sex_diff.json ./run_sas_oda_local_top_hits_manhattan_download_png.sh
RUNNER_CONFIG_JSON=./configs/runner_pgc_scz_sex_diff.json ./run_sas_oda_local_top_hits_with_gtf_download_html.sh
TARGET_SNP=rs17425819 RUNNER_CONFIG_JSON=./configs/runner_pgc_scz_sex_diff.json ./run_sas_oda_single_snp_with_gtf_download_html.sh
```

These runner presets can set shell-facing variables such as:

- `PROJECT_TAG`
- `DATA_GZ`
- `SOURCE_LONG_GZ`
- `EXTRACTOR_CONFIG_JSON`
- `MANHATTAN_P_VAR`
- `MANHATTAN_OTHER_P_VARS`
- `MANHATTAN_GWAS_LABEL_NAMES`
- `TOP_HIT_FOCUS_PVAR`
- `TOP_HIT_FILTER_EXPR`
- `TOP_HIT_MAF_THRESHOLD`
- `TOP_HIT_GNOMAD_FREQ_FILE`
- `TOP_HIT_GNOMAD_POP_MAP`
- `GTF_ASSOC_PVARS`
- `GTF_ZSCORE_VARS`
- `GTF_LABELS`

For the single-SNP GTF runner, `EXTRACTOR_CONFIG_JSON` lets the shell wrapper
call `extract_single_snp_wide_diff_gwas.pl --config ...` rather than relying on
a project-specific hard-coded long GWAS path.

## Dependency layout

The main automation entry point stays at the repo root:

- `auto_prepare_and_run_diff_gwas.pl`

Its helper layer now lives under:

- `DiffGWASDeps/`

That folder contains the reusable Perl/Bash/SAS dependency scripts such as:

- `merge_pgc_vcf_sumstats_long.pl`
- `diff_pairwise_gwas.pl`
- `standardize_diff_gwas_zscore.pl`
- `extract_significant_diff_gwas.pl`
- `extract_single_snp_wide_diff_gwas.pl`
- `emit_diff_gwas_runner_env.pl`
- `generate_sas_wide_import_include.pl`
- `render_sas_template.pl`
- `sort_long_gwas_by_coord.sh`
- `run_sas_oda_manhattan4diffgwas_download_png.sh`
- `run_sas_oda_local_top_hits_manhattan_download_png.sh`
- `run_sas_oda_local_top_hits_with_gtf_download_html.sh`
- `run_sas_oda_single_snp_with_gtf_download_html.sh`
- `run_sas_oda_manhattan4diffgwas.sas`
- `run_sas_oda_local_top_hits_manhattan.sas`
- `run_sas_oda_local_top_hits_with_gtf.sas`
- `run_sas_oda_single_snp_with_gtf.sas`
- `Manhattan4DiffGWASs_png.sas`
- `Lattice_gscatter_over_bed_track.sas`

The shell runners and MCP wrapper were updated to resolve those helpers from
`DiffGWASDeps/` automatically. The project-root shell runner names are kept as
small compatibility wrappers that forward into the real implementations in
`DiffGWASDeps/`.

## SAS ODA macro-loading behavior

The Perl SAS ODA runner now clearly follows this rule:

- when a new SAS ODA session is created, it auto-loads macros from `~/Macros`
  once via `importallmacros_ue(MacroDir=%sysfunc(pathname(HOME))/Macros,fileRgx=.,verbose=0)`
- when a persistent session is reused, it does not rerun `%importallmacros_ue`
- if a dead session must be recreated, the macro load happens again for that new session

Because of that, submitted SAS scripts in this pipeline should not call
`%importallmacros_ue` themselves. The redundant calls were removed from the
GTF-backed local runners, and the shell wrappers no longer upload
`importallmacros_ue.sas` for those runs.

## SAS schema layer

The SAS import blocks are now generated from config too.

The new helper:

- `generate_sas_wide_import_include.pl`

builds the `data scz_mh; infile ...; input ...;` block from:

- `base_cols`
- `value_fields`
- `pair_map`
- `prefix_order`
- `char_lengths`
- `alias_map`

It also auto-derives Z-score variables for any `<name>_BETA` plus matching
`<name>_SE` pair. The shell runners render that generated block into the SAS
templates through:

- `render_sas_template.pl`

That means the extractor layer, the runner layer, and the SAS wide-import layer
now all follow the same config contract.

## Current boundary

The remaining project-specific part is upstream of the wide schema: how a raw
long differential GWAS is produced from the source cohort files and what pair
tags it uses. Once that long differential GWAS exists, the downstream wide
extraction and SAS plotting layers are now config-driven.

## PGC ancestry notes

The raw files currently present in:

- `E:\LongCOVID_HGI_GWAS\PGC_Large_GWASs\PGC_SCZ_Ancestry_Stratified_GWASs`

show these available ancestry/group strata in their filenames:

- `afram`
- `asian`
- `european`
- `latino`
- `core`
- `primary`

So future ancestry-difference presets should be built around the actual pairwise
contrasts you generate from those strata, for example `ASN_vs_EUR`,
`AFR_vs_EUR`, or `LAT_vs_EUR`, rather than reusing sex-specific field names.

## Tested ancestry run

The generalized pipeline was tested on:

- `PGC3_SCZ_wave3.afram.autosome.public.v3.vcf.tsv.gz`
- `PGC3_SCZ_wave3.asian.autosome.public.v3.vcf.tsv.gz`
- `PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz`

with pairwise contrasts:

- `SCZ_W3_ASN_vs_EUR -> ASN_EUR`
- `SCZ_W3_AFR_vs_EUR -> AFR_EUR`
- `SCZ_W3_ASN_vs_AFR -> ASN_AFR`

Validated generated artifacts:

- merged long table:
  `PGC_SCZ_ancestry_diff_effects_merged_long.tsv.gz`
- sorted/indexed long table:
  `PGC_SCZ_ancestry_diff_effects_merged_long.sorted.coord.tsv.gz`
- differential table:
  `PGC_SCZ_ancestry_diff_effects.tsv.gz`
- standardized differential table:
  `PGC_SCZ_ancestry_diff_effects.stdized.tsv.gz`
- wide subset:
  `PGC_SCZ_ancestry_diff_effects.stdized.wide_beta_se_p_p_lt_0p05.final.tsv.gz`
- auto-generated configs:
  `configs/auto_PGC_SCZ_ancestry_diff_effects_{merge,diff,preset,runner}.json`
- SAS outputs:
  `PGC_SCZ_ancestry_SAS_manhattan.png`
  `PGC_SCZ_ancestry_SAS_local_top_hits_manhattan.png`
  `PGC_SCZ_ancestry_SAS_local_top_hits_with_gtf.html`

Updated plotting-spec test results:

- refreshed ancestry wide subset now includes `STD_DIFF_Z` columns in addition
  to `STD_DIFF_P`
- updated ancestry runner config uses:
  - genome-wide/local Manhattan:
    `ASN_EUR_STD_P`, `AFR_EUR_STD_P`, `ASN_AFR_STD_P`, `AFR_P`, `ASN_P`, `EUR_P`
  - local GTF:
    `ASN_EUR_STD_P AFR_EUR_STD_P ASN_AFR_STD_P AFR_P ASN_P EUR_P`
    with matching
    `ASN_EUR_STD_Z AFR_EUR_STD_Z ASN_AFR_STD_Z AFR_Z ASN_Z EUR_Z`
- ancestry rerun outputs were regenerated successfully:
  - `PGC_SCZ_ancestry_SAS_manhattan.png`
  - `PGC_SCZ_ancestry_SAS_local_top_hits_manhattan.png`
  - `PGC_SCZ_ancestry_SAS_local_top_hits_with_gtf.html`
- the updated local GTF HTML contains `4` embedded images

Regression check on the original sex-stratified schizophrenia workflow:

- the updated genome-wide runner spec also completed successfully for
  `configs/runner_pgc_scz_sex_diff.json`
- that run produced a fresh `PGC_SCZ_SAS_manhattan.png`

Key ancestry issues found and resolved:

1. `afram` header schema is slimmer than `asian/european`.
   Its file lacked `NCAS`, `NCON`, and `NEFF`, so `merge_pgc_vcf_sumstats_long.pl`
   was updated to treat those fields as optional.

2. The sorted long table is `bgzip`, not plain single-member gzip.
   `diff_pairwise_gwas.pl` was updated to stream `/mnt/...` gzip inputs through
   `gzip -dc`, so it can read full multi-member `bgzip` files instead of only the
   first header block.

3. Automation shelling from Cygwin Perl must use Cygwin bash explicitly.
   `auto_prepare_and_run_diff_gwas.pl` now defaults to `/bin/bash` instead of
   plain `bash`, which had resolved to WSL on this machine.

4. Local-top-hit Manhattan output must not reuse the genome-wide output prefix.
   The local runner now honors `LOCAL_OUTPUT_PREFIX` and `LOCAL_HTML_TITLE`,
   preventing overwrites of the genome-wide Manhattan files.

5. The refreshed local GTF plotting spec required `STD_DIFF_Z` to survive the
   wide-subset extraction, not only `STD_DIFF_P`. The config defaults and wide
   import generator were updated so `STD_DIFF_Z` is carried through and exposed
   to SAS through `post_alias_map` where needed.

Remaining caveat from the tested ancestry local/GTF path:

- The local-top-hit and GTF SAS logs still include macro-level warnings and some
  `ERROR:` lines, even though the output files were generated successfully.
  Future runs should verify both file existence and rendered-image counts, not
  assume that a noisy SAS log means total failure.

- In the updated ancestry test, the local Manhattan and local GTF logs still
  report duplicate-key and recursive `CREATE TABLE` warnings, and the GTF log
  still reports `WORK.FINAL_LABEL_SORTED` errors inside the shared macro stack.
  Those messages did not stop HTML/PNG generation, but they remain good cleanup
  targets for future macro debugging.

## Adaptive local-label spacing

For split local top-hit Manhattan panels, do not treat one fixed pair like
`xgrp_y_pos=-1.8` and `offset=(20,0.5)` as a universal solution.

The bottom `SNP:gene` labels behave better when spacing is guessed from the
actual panel density:

- how many GWAS tracks are stacked
- how many top-hit columns are in the current batch

The current implementation now does this inside:

- `DiffGWASDeps/run_sas_oda_local_top_hits_manhattan.sas`

Current rule of thumb:

- sparse panels keep milder bottom spacing
- dense panels, especially about 6 tracks x 6 columns, get:
  - more negative `xgrp_y_pos`
  - smaller bottom `offset=(...)`
  - slightly smaller `gwas_label_y_frac` so the bottom track title sits higher
  - wider SNP-vs-gene separation in the rotated two-line label split

For the most crowded local panels in this project, a better rule of thumb was
closer to:

- `xgrp_y_pos=-2.8` or lower
- `yoffset_setting=offset=(15,0.5)` or larger when the bottom box needs more height

Also, do not rely on spacing tweaks alone when the panel is still too crowded.
For local Manhattan figures with many stacked GWAS tracks, reduce the number of
top-hit columns per panel as well:

- about 6 stacked GWAS tracks -> cap local panels at 4 columns
- about 4-5 stacked GWAS tracks -> cap local panels at 5 columns

The current runner now does that automatically before calculating
`n_hit_batches`.

One more label-placement lesson from debugging:

- for rotated `SNP:gene` labels in `Manhattan4DiffGWASs_png.sas`, the original
  annotate `position` logic matters
  - first line: `position='B'`
  - second line: `position='E'`
- treat `xgrp_y_pos` and `yoffset_setting` as ways to move the whole label box
  up or down
- do not replace the `B/E` split with ad hoc x/y overrides unless there is a
  very specific reason, because that can make the two vertical lines harder to
  distinguish

Trial-and-error rule that worked best here:

1. if bottom `SNP:gene` labels are too close to the scatter band, first extend
   the bottom box by increasing the first `offset=(...)` value
2. then move the label block farther down with a more negative `xgrp_y_pos`
3. only after that consider reducing `local_max_hits_per_fig` further

The automation entrypoint exposes the requested panel width directly:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --spec your_spec.json \
  --step plot_local_manhattan \
  --local-max-hits-per-fig 4
```

That value is treated as an upper bound. The CLI/spec layer accepts up to
10 columns per panel, and the SAS runner may still tighten it further for dense
multi-track figures.

This is specifically to reduce overlap between:

- the bottom `SNP:gene` labels
- the lowest scatter band
- the lowest stacked-track title

## Recent bug fix (2026-05-04)

- Issue: the common-association mode sometimes produced an over-restrictive
  SAS filter expression containing clauses like `(..._STD_P>=0.5)` which
  excluded valid single-GWAS top signals that also had nominal support in a
  partner GWAS.
- Fix: `auto_prepare_and_run_diff_gwas.pl` now strips accidental `_STD_P>=0.5`
  (and similarly formatted) clauses when `--get-common-associations` is used.
  This prevents valid common-association loci from being dropped.
- Verification: the automation will call
  `DiffGWASDeps/verify_common_association_loci.pl` in common-association
  mode and write two outputs near the artifact directory:
  - `<artifact_stem>.common_assoc_candidates.tsv` (full candidate list)
  - `<artifact_stem>.common_assoc_verify.tsv` (distance-pruned loci)
- Repro command (example):

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --spec configs/test_PGC_GWAS4testing_all3_manual.json \
  --force --get-common-associations
```

- Observed result for `PGC_GWAS4testing`: verifier found 6,535 candidates
  at `5e-8` and the pruned list contained 14 loci (see generated
  `.common_assoc_verify.tsv`).
- Observed result for the bundled PGC schizophrenia sex-stratified wide file
  during the 2026-06-14 review: verifier found `11,703` candidates at
  `5e-8`; with the corrected `1e6` total pruning span and the default
  `MAF > 0.01` QC, the pruned list contained `103` loci
  (`tmp_common_verify_postfix_1e6.tsv`). The earlier undercount came from an
  overly broad `1e8` pruning span and a wrapper-level `15`-locus cap, not
  from a flaw in the grouped greedy SAS macro itself.

If you prefer, I can wire the SAS runner to use the verifier-derived CSV as
the plotting input to guarantee the same set of loci is plotted as verified.
