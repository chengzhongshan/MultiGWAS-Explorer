# PGC Sex Pipeline Benchmark

Date: `2026-04-28`

Dataset:
- `E:\LongCOVID_HGI_GWAS\PGC_Large_GWASs\PGC_SCZ_Sex_Stratified_GWASs\PGC_SCZ_female_vs_male_diff_effects.stdized.wide_beta_se_p_p_lt_0p05.final.tsv.gz`

Scope:
- benchmark the operational plotting pipeline on the validated sex-specific wide subset
- avoid re-extracting from the full `2.05 GB` long standardized file during timing

Config updates applied for the sex preset before benchmarking:
- `TOP_HIT_FOCUS_PVAR=ASN_STD_DIFF_P`
- `TOP_HIT_FILTER_EXPR=((ASN_STD_DIFF_P>0) and (ASN_STD_DIFF_P<1E-4)) or ((EUR_STD_DIFF_P>0) and (EUR_STD_DIFF_P<1E-4))`
- `TOP_HIT_SIGNAL_THRSHD=1e-4`

## Cold timings

| Step | Seconds | Status | Notes |
|---|---:|---|---|
| Delete remote shared gz | 4.485 | OK | cleanup before cold upload |
| Upload shared sex wide gz | 28.055 | OK | persistent ODA session, fresh upload |
| Genome-wide Manhattan | 29.175 | OK | reused uploaded wide gz |
| Single-SNP GTF (`rs17425819`) | 185.055 | OK | includes helper macro uploads, wide gz upload, SAS run, HTML download, cleanup |

## Local top-hit Manhattan status

The sex local top-hit Manhattan workflow did not produce a stable benchmarkable PNG in this session.

What was learned:
- the original sex preset was too strict at `1e-6`
- the runner config was corrected to use `ASN_STD_DIFF_P` / `EUR_STD_DIFF_P` with `1e-4`
- after correction, remaining failures were dominated by SAS ODA session-helper instability and rerun sequencing around remote input reuse

Practical status:
- this stage is **not yet reliable enough** to treat as a final benchmark number

## Warm-reuse benchmark status

Warm-reuse reruns were attempted with:
- `SKIP_DATA_UPLOAD=1`
- `KEEP_REMOTE_PLOT_DATA=1`
- persistent SAS ODA session reuse

Observed issue:
- the persistent session helper intermittently failed with:
  - `Upload failed: PYTHON ERROR: failed to read response header from session server`

Interpretation:
- the plotting pipeline itself is not the only moving part
- the long-lived session-server transport remains unstable enough that warm-reuse timings are not yet trustworthy as final benchmark numbers

## Output artifacts confirmed

- `benchmark_pgc_sex_plot_manhattan.png`
- `benchmark_pgc_sex_plot_manhattan_png.html`
- `benchmark_pgc_sex_plot_local_top_hits_with_gtf.html`

## Recommendation

Use these as the current reliable benchmark numbers for the sex-stratified pipeline:
- cold shared upload: about `28 s`
- genome-wide Manhattan: about `29 s`
- single-SNP local GTF: about `185 s`

Before publishing a warm-reuse benchmark, first stabilize the persistent session server used by `run_sas_codes_or_script_in_ODA.pl`.
