---
name: gwas-sas-oda-workflow
description: "Reusable scientific workflow for GWAS summary statistics: merge sex/population/dataset-stratified GWAS files, sort/bgzip/tabix-index combined GWAS tables, compute/subset data for multiple Manhattan plots, run SAS ODA plotting through a Perl MCP server, download PNG/HTML outputs, and clean temporary SAS ODA files. Use when Codex is asked to orchestrate GWAS data processing across local Cygwin/Perl scripts, htslib/tabix, and SAS OnDemand for Academics."
---

# GWAS SAS ODA Workflow

## Strategy

Use this skill as the workflow brain. Keep MCP tools small and reliable.

- Put scientific orchestration, file naming, validation, and decision rules in this skill and project scripts.
- Put atomic capabilities in the Perl MCP server: run Cygwin bash/Perl, run/upload/download/delete/list SAS ODA files, inspect directories, and optionally wrap bgzip/tabix.
- Do not put one-off cohort-specific GWAS logic into `server.pl`; keep that logic in versioned Perl/SAS/Bash scripts.

## Required Context

Before running a workflow, identify:

- Input GWAS directory and compression format.
- Whether the inputs are separate raw GWAS tables or one merged-wide GWAS
  table that already carries shared locus columns plus cohort-level
  `BETA_*` / `SE_*` / `P_*` blocks.
- Column mapping for chromosome, position, alleles, SNP ID, beta or OR, SE, P, and dataset tag.
- Whether chrX files must be vertically combined with sex/population-matched autosomal files.
- Whether the GWAS tables already contain allele-frequency columns that can be
  used for top-hit MAF filtering, and if not, whether a local gnomAD lookup
  file is available. If neither is available, the current safeguard will
  conservatively drop those candidates instead of keeping them without an MAF
  check.
- Desired final products: long merged table, tabix-indexed coordinate table,
  differential effects, Manhattan subset, local GTF plot, forest plot, and
  SAS ODA or gunplot PNG/HTML outputs.
- Cygwin path equivalents for Windows drives, usually `/mnt/e/...` and `/mnt/g/...`.
- Whether the user wants a fast install/render validation or a full
  genome-wide rerun. In Ubuntu Docker, those are very different time costs.
- On Windows, whether the user is installing through the portable Cygwin
  bootstrap or an existing Cygwin shell. The portable bootstrap defaults to
  `%USERPROFILE%\CygwinPortablePipeline`; pass `-PortableRoot` only when a
  different isolated root is needed.
- Whether the local network intercepts TLS. If Cygwin `curl` fails with a
  self-signed certificate chain during bootstrap, rerun the Windows portable
  installer with `-AllowInsecureDownloads` so phase-2 repo-local downloads use
  `PIPELINE_CURL_INSECURE=1`.
- SAS ODA size limits; avoid uploading large raw GWAS files.

## Workflow

1. Inspect inputs with Cygwin commands.
   Prefer `zcat file.gz | head`, `zcat file.gz | wc -l`, and header checks through the Perl MCP `run_perl_or_bash_cmd` tool.

   For merged-wide inputs in this repository:

   - first preview the inferred spec with:
     - `perl auto_prepare_and_run_diff_gwas.pl --gwas-dir AOA_GWAS_Data/ --preview-spec --generate-spec-only`
   - confirm the detector chose `source_mode=merged_gwas_table`
   - confirm the selected source file is the original merged study table, not a
     generated `*.merged_plotwide.tsv.gz` or prior plot artifact
   - confirm the inferred wide-column blocks match the real header, for example:
     - `BETA_DS_ALL`, `SE_DS_ALL`, `P_DS_ALL`
     - `BETA_MP2PRT`, `SE_MP2PRT`, `P_MP2PRT`
   - if the merged table uses nonstandard headers, add a
     `--raw-column-alias-config` JSON before attempting a full run

   For Ubuntu Docker validation in this repository:

   - keep `PIPELINE_WORKDIR=/opt/MultiGWAS-Explorer` so the wrapper uses the
     Linux-installed repo inside the image
   - do not mount a host checkout over `/opt/MultiGWAS-Explorer`, because that can
     hide the image's Linux `.venv-pipeline/` and repo-local Perl modules
   - expect the first uncached image build to take several minutes because the
     container installs Ubuntu packages, Java, gnuplot, ImageMagick, Python
     headers, and repo-local Perl/Python dependencies
   - after the image exists, run `bash install/check_pipeline_install.sh`
     first; that should be quick compared with a full plot rerun
   - for a quick renderer validation, prefer one inquiry SNP with:
     - `--plots local_manhattan,local_gtf`
     - `--target-snps rs185665940`
   - only add `manhattan` after the local plots succeed, because the
     genome-wide gunplot stage is the main runtime hotspot
   - in the 2026-06-08 validation of this repository, the Ubuntu Docker
     genome-wide gunplot step scanned 2,417,954 wide rows and wrote an
     intermediate `.plot.tsv` of about 110 MB; that single stage took about
     9 minutes, whereas the one-locus local Manhattan and local GTF panels each
     finished in about 10-12 seconds

   For Windows portable Cygwin installation in this repository:

   - run the PowerShell bootstrap from the repository root:
     - `powershell -NoProfile -ExecutionPolicy Bypass -File .\install\install_windows_portable_cygwin.ps1`
   - expect the installer to refresh Cygwin packages including the htslib build
     headers: `libbz2-devel`, `libcurl-devel`, `liblzma-devel`,
     `openssl-devel`, and `zlib-devel`
   - if `bgzip` and `tabix` are absent, the installer downloads htslib 1.20
     into `tools/` and builds repo-local tools under `local/bin/`
   - if TLS verification fails while downloading `cpanm` or htslib, rerun with:
     - `powershell -NoProfile -ExecutionPolicy Bypass -File .\install\install_windows_portable_cygwin.ps1 -AllowInsecureDownloads`
   - if the installer warns that it could not resolve `java.exe`, treat the
     dependency smoke test as valid but set `SASPY_JAVA_WIN` before launching
     SAS ODA sessions

2. Merge GWASs locally.
   Write or reuse a Perl script that streams gz files and outputs a normalized long-format table. Use stable keys such as `CHR`, `BP`, `A1`, `A2` only when those columns are explicitly available. Do not relabel A1/A2 as REF/ALT unless the source documentation says so.

3. Sort and index.
   Sort by chromosome and position, remove or quarantine non-coordinate rows, compress with `bgzip`, and index with `tabix -s <chr_col> -b <pos_col> -e <pos_col>`. For custom genome-wide Manhattan subsets, do not treat the original `CHR` text as already plot-safe; first normalize any `chr` prefix and sex-chromosome labels, then sort on the normalized coordinate fields that SAS will actually read.

4. Prepare plot subset locally.
   Transform long data to the wide format needed by the plotting macro. Keep only required columns and filter to a reasonable P threshold such as any selected P `< 0.05`, unless the user requests otherwise. Normalize chromosome labels before the final export that SAS will read:
   - strip an optional `chr` prefix
   - map `X -> 23`
   - map `Y -> 24`
   - drop rows whose chromosome or base-pair position still cannot be parsed
   - sort by the normalized chromosome/position fields, not the raw text labels
   A practical failure mode from this repository: if SAS imports `CHR` as numeric while the compact subset still contains a literal `X`, SAS converts that row to missing `.`, and the genome-wide Manhattan plot can show a false extra block before `chr1`.
   When top-hit reliability matters, preserve GWAS frequency/info fields in the
   wide subset too:
   - `GROUP1_FRQ_A`
   - `GROUP1_FRQ_U`
   - `GROUP2_FRQ_A`
   - `GROUP2_FRQ_U`
   - `GROUP1_INFO`
   - `GROUP2_INFO`
   Then apply the top-hit safeguard:
   - use GWAS-derived MAF first
   - fall back to a local gnomAD lookup only when GWAS frequency columns are
     absent
   - if neither source is available, expect `maf_source=UNKNOWN` and
     conservative filtering at the default threshold
   - keep hits only when `MAF > 0.01`, unless the user requests a different
     threshold through `top_hit_maf_threshold`
   The requested top-hit CSVs should preserve:
   - `selected_maf`
   - `maf_source`
   - `gwas_group1_maf`
   - `gwas_group2_maf`
   - `gwas_pair_maf_min`
   - `gnomad_maf`
   - `gnomad_pops`
   - `maf_filter_decision`
   - `maf_filter_reason`
   For merged-wide studies, prefer the in-repo converter path rather than
   hand-rolling a temporary wide table:
   - `DiffGWASDeps/convert_merged_gwas_to_plotwide.pl`
   - `DiffGWASDeps/generate_sas_wide_import_include.pl`

5. Run SAS ODA only on small subsets.
   Upload the gz subset and macro files. In SAS, read gz files using:

   ```sas
   filename mhgz zip "~/file.tsv.gz" gzip;
   ```

6. Generate PNG, not inline SVG/HTML.
   Close SASPy's internal inline destination with `ods _all_ close;`, use `ods listing` for SAS/GRAPH PNG output, and create a tiny HTML wrapper that references the PNG.

7. Download and verify outputs.
   Download the PNG and HTML explicitly after the plot run. Check both local files exist and are non-empty. Open the HTML only after verification. When a local GTF run produces a PNG, prefer a small figure-first final HTML that embeds that PNG and keep any raw SAS HTML as a sidecar instead of opening the sparse helper HTML directly.
   For forest plots, also verify that:
   - a single-SNP run produced one manifest row with `track_id=single_snp`
   - a multi-SNP run produced one manifest row per displayed GWAS / cohort
   - the manifest panel labels are not blank
   - the PNG y-axis labels match the intended mode:
     - single SNP: cohorts on the y-axis
     - multi SNP: SNP IDs on the left, nearby genes on the right

8. Clean SAS ODA temporary data.
   Delete uploaded gz subsets and transient macro files after successful download, unless the user asks to keep them.
   Prefer the low-level SAS ODA wrapper's bulk operations when several remote
   files must be managed together. In this project,
   `run_sas_codes_or_script_in_ODA.pl` now accepts repeated upload, download,
   file-info, and delete arguments in one invocation, plus regex-based delete
   patterns. Regex deletes match both basenames like `plot.png` and resolved
   remote paths like `~/plot.png`, so either `.*\.png$` or `~\/.*\.png` is
   acceptable depending on what is clearer in context. When calling from
   PowerShell or another shell layer, quote remote home paths like
   `'~/plot.png'` so the literal `~` reaches the helper unchanged. The helper
   now also normalizes both `~/...` and absolute SAS home paths such as
   `/home/...` consistently for `download-file`, `delete-file`, and
   `file-info`, and it verifies deletes by checking that the target no longer
   resolves afterward.

9. Preserve local GTF gene-track intent explicitly.
   When building local GTF plots, prefer a region-limited Gencode subset over a
   large prebuilt library table so the workflow can control whether
   non-protein-coding genes are included. In this project, local GTF plots are
   protein-coding-only by default; enable non-coding genes only by setting
   `include_non_protein_coding_genes_in_local_gtf: 1` in the spec JSON.
   For large local windows, prefer the project path that pre-extracts the GTF
   subset locally and uploads that compact subset to SAS ODA, rather than
   asking SAS `WORK` to materialize a very large GTF region on demand.
   In this project, `--local-gtf-window-bp` now controls both the extracted GTF
   half-window and the displayed local GTF plot half-window.
   The subset upload path now gzip-compresses that local GTF table before
   transfer, and the SAS import block reads it through
   `filename ... zip ... gzip`. For a first stress-test rerun, prefer a window
   such as `1e8` before trying `1e9`.
   When many genes overlap one local locus, keep the overall figure size fixed
   and instead let the pipeline slightly increase the lower gene-track share by
   auto-tuning the SAS `pct4neg_y` parameter.
   In this repository, both the batched and the single-SNP SAS local-GTF
   wrappers now start from the same more readable defaults before that
   auto-tuning path runs:
   - `GTF_PCT4NEG_Y=1.4`
   - `GTF_DESIGN_HEIGHT=1000`
   Prefer raising `GTF_PCT4NEG_Y` first when the bottom gene/exon track still
   looks cramped, instead of immediately editing the SAS macro by hand.
   Prefer the in-repo `DiffGWASDeps/SNP_Local_Manhattan_With_GTF.sas` over an
   untracked remote copy in SAS ODA so local GTF reruns use the repository's
   patched plotting logic.
   When a large local GTF window causes the underlying macro to expand the gene
   search across chromosome start, preserve the final displayed x-axis by
   forcing it back to the min/max association-signal positions from the GWAS
   subset instead of accepting a left boundary of `0`.
   If SAS already completed remotely but the wrapper had trouble downloading the
   final HTML, prefer recovering the helper-saved `sas_res_*.html` artifact and
   the final PNG path mentioned in the SAS log before declaring the rerun
   failed. If the PNG was recovered successfully, rebuild the user-facing final
   HTML around that PNG and keep the raw SAS HTML as a sidecar.
   For the direct single-SNP SAS local-GTF wrapper, remember that
   `TARGET_SNP` still defines the centered local window while
   `GTF_LABEL_SNPS` can supply additional comma-separated rsIDs to label inside
   that same locus. Prefer this path when several inquiry SNPs fall in one
   shared window and you want one centered local-GTF figure instead of one
   locus per target. When the label list contains three or fewer rsIDs, the
   default `auto_rotate2zero=1` path now keeps those top labels horizontal in
   the local-GTF headroom.
   The single-SNP extractor manifest now also reports whether the target row
   survived and which prefix blocks were present or missing. Use that manifest
   before blaming a sparse locus: a target row can still be valid for local GTF
   centering even when one ancestry / pair block is blank in the emitted wide
   row.
   If SAS reports that the target SNP was not found in the uploaded GWAS subset
   but the local manifest already says `target_row_found_in_window=1`, treat
   that first as a remote upload / helper integrity problem. In this project,
   the wrapper now verifies remote uploaded file sizes and can keep recovering
   the remote HTML when SAS finished but the helper failed late.

10. Treat local Manhattan layout tuning as a first-class rerun option.
   If the top-hit SNP-gene labels are crowded or visually unbalanced, rerun the
   local Manhattan stage with the pipeline overrides instead of patching SAS by
   hand. Prefer:
   - `--local-max-hits-per-fig` up to `30`
   - `--local-manhattan-angle4xaxis-label`
   - `--local-manhattan-xgrp-y-pos`
   - `--local-manhattan-yoffset-top`
   - `--local-manhattan-yoffset-bottom`
   - `--local-manhattan-fontsize`
   - `--local-manhattan-y-axis-label-size`
   - `--local-manhattan-y-axis-value-size`

11. Prefer GTF fallback gene assignment for local top-hit plots.
   When HaploReg does not return a valid nearby gene for a top hit, use the
   pipeline's region-limited Gencode fallback rather than leaving the label as
   `NA`. The same principle now applies to local Manhattan top-hit labeling and
   local GTF bottom tracks. If the top-hit CSV is missing locally during a
   local-GTF rerun, use the verified common-association or differential top-hit
   table as the fallback source of locus coordinates before building the local
   GTF subset.
   The same fallback principle now applies to multi-SNP forest plots when the
   right-side gene labels would otherwise be missing or `NA`.

12. Treat forest plots as first-class outputs in both backends.
   In this repository:
   - SAS ODA forest plots are driven by:
     - `DiffGWASDeps/beta2OR_forest_plot.sas`
     - `DiffGWASDeps/run_sas_oda_top_hits_forest_plot.sas`
   - gunplot forest plots are driven by:
     - `DiffGWASDeps/gunplot/pdl_gunplot_forest.pl`
     - `auto_prepare_and_run_diff_gwas_with_gunplot.pl --plots forest`
   Prefer these conventions:
   - single-SNP forest plot:
     - cohorts on the y-axis
     - `OR and 95% CI` on the x-axis
   - multi-SNP forest plot:
     - one panel per cohort / displayed GWAS track
     - SNP IDs on the left y-axis
     - nearby gene labels on the right y-axis
     - horizontal dashed divider between differential and common hits
     - significance star for `P < 5e-8`
   The gunplot renderer is intentionally styled to stay close to the SAS ODA
   forest output, so small remaining differences should mostly be backend
   rasterization details rather than different ordering, labels, or grouping.

13. Regenerate manuscript tables from the same validated top-hit sources.
   When the user needs manuscript-ready supplementary tables, prefer:
   - `DiffGWASDeps/regenerate_manuscript_hit_tables.pl`
   Current default behavior in this repository:
   - `Table_S1_all_common_association_loci.csv` is now a full-strata common-hit
     table, not a pooled-only summary
   - `Table_S2_differential_loci.csv` is now a full-strata differential-hit
     table, not a pooled-only summary
   - both tables retain pooled and ancestry-specific association `P`, `BETA`,
     and `SE`, plus pairwise differential `P`, `BETA`, and `SE` for
     `ALL`, `EUR`, and `ASN`
   - both tables also preserve the MAF QC provenance needed for later review
   The helper also writes mirrored convenience copies:
   - `Table_S1_all_common_association_loci_full_strata.csv`
   - `Table_S2_differential_loci_full_strata.csv`
   Keep the narrower manuscript main-text outputs separate:
   - `Table_2_representative_common_loci.csv`
   - `Table_1_top_differential_locus.csv`
   If the user only has an older pooled-only `Table_S1` export and does not
   want to rerun the full regeneration step, use:
   - `DiffGWASDeps/augment_common_hits_table_s1.pl`
   If a Windows-editable workbook is needed afterward, use:
   - `DiffGWASDeps/export_augmented_table_s1_excel.ps1`

## Project Scripts

When working in this user's GWAS/SAS ODA project, prefer these existing scripts if present:

- `merge_scz_sex_stratified_long.pl`
- `bgzip_tabix_diff_gwas.pl`
- `standardize_diff_gwas_zscore.pl`
- `extract_significant_diff_gwas.pl`
- `prepare_sas_manhattan_subset.pl`
- `convert_merged_gwas_to_plotwide.pl`
- `run_sas_oda_manhattan4diffgwas.sas`
- `Manhattan4DiffGWASs_png.sas`
- `run_sas_oda_manhattan4diffgwas_download_png.sh`
- `run_sas_oda_local_top_hits_manhattan.sas`
- `run_sas_oda_local_top_hits_manhattan_download_png.sh`
- `run_sas_oda_local_top_hits_with_gtf_download_html.sh`
- `run_sas_oda_single_snp_with_gtf_download_html.sh`
- `regenerate_manuscript_hit_tables.pl`
- `augment_common_hits_table_s1.pl`
- `auto_prepare_and_run_diff_gwas.pl`

For a new project, copy/adapt the script templates in this skill's `scripts/` folder only after checking the target column names and file sizes.

## Validation Checklist

- Report row counts at each major stage.
- Validate headers and column mappings before a full run.
- For merged-wide auto-detection, verify that `--gwas-dir` selected the
  original merged table instead of a generated `merged_plotwide`, `gunplot`,
  `png`, or `html` artifact.
- For merged-wide spec generation, verify that the inferred `wide_columns`
  blocks match the cohort names you expect to display and compare.
- If top-hit MAF filtering is enabled, report whether the final hit list used
  GWAS frequencies or gnomAD fallback, and whether any candidates were removed
  because `maf_source=UNKNOWN`.
- Use `perl DiffGWASDeps/test_top_hit_maf_filter.pl --no-real --keep-workdir`
  for a quick synthetic regression of the MAF safeguard branches.
- For the bundled PGC schizophrenia example, expect the real differential and
  common-association validations to take minutes rather than seconds; run them
  separately if you want clearer timing/logs.
- For tabix output, test at least one autosomal interval and chrX/23 interval.
- For SAS ODA output, verify local PNG and HTML sizes. A tiny HTML wrapper is expected; a hundreds-of-MB HTML file means inline ODS output leaked into the result.
- For custom genome-wide Manhattan subsets, visually verify that the x-axis
  starts with `chr1` and does not contain a phantom pre-`chr1` block caused by
  unsorted rows or `CHR='X'` being imported as numeric missing `.` in SAS.
- For local GTF reruns, also verify the SAS log reports a non-zero exon count.
- For forest plots, verify the requested top-hit CSV, manifest TSV, PNG panel
  count, and left/right y-axis labeling mode.
- For manuscript-table regeneration, verify that:
  - `Table_S1_all_common_association_loci.csv` carries pooled plus ancestry
    strata columns instead of only pooled `ALL_*` fields
  - `Table_S2_differential_loci.csv` carries pooled plus ancestry strata
    columns instead of only pooled `ALL_*` fields
  - the full-strata mirrored CSVs were emitted
  - the row counts still match the selected common and differential loci
- Confirm SAS ODA cleanup by listing or checking the deleted paths.
- For `run_sas_codes_or_script_in_ODA.pl` incidents, prefer this quick AI
  debug ladder before moving into `%include` or plotting macros:
  - `perl DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl --check-sas-oda-login-only`
  - `perl DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl --output-prefix oda_file_probe --dir4listing '~' --file-info '~/importallmacros_ue.sas'`
  - `perl DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl --output-prefix oda_transfer_probe --upload-file ./small_test.txt --download-file '~/small_test.txt' --download-local-path ./small_test.roundtrip.txt`
  - `perl DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl --output-prefix oda_code_probe --code "%put HELLO_FROM_ODA;"`
  - a minimal file-based `proc print` smoke script before any real plotting submit
  - `SAS_ODA_AUTOLOAD_MACROS=1 ... --file codex_put_smoke.sas --run-timeout-seconds 30` to isolate the default `~/Macros` bootstrap path
  - `SAS_ODA_AUTOLOAD_MACROS=0 ... --file codex_manual_importallmacros_smoke.sas` to time manual in-SAS `~/Macros` loading separately from the helper bootstrap path
- For the newer bootstrap instrumentation, inspect these artifacts together:
  - `<output-prefix>/output.run.status.json`
  - `<output-prefix>/output.macro_bootstrap.log.txt`
  - `<output-prefix>/output.html.info.txt`
- Interpret the newer debug lines carefully:
  - `Upload step: macro bootstrap helper: importallmacros_ue.sas ...`
    means the helper file upload finished; it is not the actual bootstrap
    submit
  - `SAS ODA macro bootstrap started at ...`
    means the wrapper has entered the real bootstrap submit
  - if the manual in-SAS macro probe is fast but the forced helper bootstrap
    probe stalls, treat that as evidence that the helper bootstrap path is the
    bottleneck rather than the `~/Macros` library contents
- When the user is debugging through PowerShell into Cygwin on Windows, prefer
  `--file` over inline `--code` for anything larger than the trivial
  `%put HELLO_FROM_ODA;` probe because shell quoting can silently truncate the
  SAS payload before submit.

## References

Read `references/tooling-strategy.md` when deciding whether a workflow step belongs in a Codex skill, a project script, or a new Perl MCP server tool.
