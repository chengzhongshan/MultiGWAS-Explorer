# Common-Association Top-Hit Review

Date: `2026-06-14`

## Question

Why did the bundled PGC schizophrenia sex-stratified common-association path
return far fewer loci than expected from the published schizophrenia GWAS
literature?

## Findings

- The SAS macro `DiffGWASDeps/get_top_signal_within_dist.sas` is not the main
  problem. Its grouped greedy pruning logic is consistent with the older macro
  family and uses a total exclusion span of `pos_dist_thrshd`, applied as
  `BP +/- 0.5 * pos_dist_thrshd`.
- The undercount mainly came from pipeline configuration and wrapper behavior:
  - the bundled common-association spec used `top_hit_dist_bp=1e8`
  - the local-top-hit wrappers hard-coded `--max-hits 15`
- The `1e8` span is too coarse for the bundled PGC schizophrenia
  common-association selector and collapses many distinct loci into one region.
- The hidden `15`-locus cap further truncated wrapper-driven outputs even when
  more loci were available.

## Published Reference

The 2022 Nature schizophrenia GWAS paper reports `287` distinct genomic loci in
the combined primary meta-analysis:

- Trubetskoy et al., Nature 2022
- https://www.nature.com/articles/s41586-022-04434-5

That published count is not expected to match this pipeline exactly because the
pipeline's common-association mode is stricter: it requires a strongest signal
in one displayed GWAS, nominal support in another displayed GWAS, and
same-direction effects across the qualifying pair.

## Real Validation on the Bundled PGC Sex-Stratified Wide File

Input:

- `E:\LongCOVID_HGI_GWAS\PGC_Large_GWASs\PGC_SCZ_Sex_Stratified_GWASs\PGC_SCZ_female_vs_male_diff_effects.stdized.wide_beta_se_p_p_lt_0p05.final.tsv.gz`

Verifier command:

```bash
perl DiffGWASDeps/verify_common_association_loci.pl \
  --spec configs/spec_pgc_scz_sex_common_automation.json \
  --output tmp_common_verify_postfix_1e6.tsv \
  --candidates-out tmp_common_verify_postfix_1e6.candidates.tsv
```

Observed metrics:

- raw qualifying rows across the ladder: `1,161,042`
- candidates with `COMMON_ASSOC_P < 5e-8`: `11,703`
- chosen threshold: `5e-8`
- candidates filtered by MAF: `0`
- selected loci after distance pruning with total span `1e6`: `103`

Saved artifacts:

- `tmp_common_verify_postfix_1e6.tsv`
- `tmp_common_verify_postfix_1e6.candidates.tsv`
- `tmp_common_verify_postfix_1e6.log`

## Distance Sensitivity on the Same 11,703 Genome-Wide Significant Candidates

Using the saved candidate table from the same run:

- total span `1e8` -> `39` loci
- total span `1e7` -> `75` loci
- total span `1e6` -> `103` loci
- total span `5e5` -> `135` loci

Interpretation:

- `1e8` is too broad for this use case and materially undercounts loci
- `1e6` restores a publication-scale count while remaining far below the raw
  candidate total

## Code Changes Applied

- changed the common-association default pruning span to `1e6`
- added `top_hit_max_loci` / `TOP_HIT_MAX_LOCI`
- removed the wrapper hard-cap of `15` loci
- clarified in the SAS macro comment that `pos_dist_thrshd` is the total span

## Practical Recommendation

For bundled PGC schizophrenia common-association reruns:

- keep `COMMON_ASSOC_P` at genome-wide significance (`5e-8`) when possible
- use `top_hit_dist_bp=1e6` as the default pruning span
- leave `top_hit_max_loci=0` unless you intentionally want a hard cap
