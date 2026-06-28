# Evaluation: Same-Direction Common-Association Mode

Date: 2026-05-07

## Goal

Reevaluate the `--get-common-associations` feature so that local top-hit
selection requires:

1. a strong single-GWAS association in one GWAS
2. at least one nominal association in the paired comparison GWAS
3. the same effect direction across that paired comparison

This evaluation used PGC schizophrenia datasets already used in the pipeline:

- sex-stratified PGC schizophrenia GWASs
- the mixed 3-GWAS testing set in `PGC_GWAS4testing`

## Implementation Changes

### 1. Pipeline logic

Updated:

- [auto_prepare_and_run_diff_gwas.pl](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/auto_prepare_and_run_diff_gwas.pl)

Common-association mode now builds `TOP_HIT_FILTER_EXPR` using:

- strongest single-GWAS association P (`COMMON_ASSOC_P`)
- partner nominal association P `< 0.05`
- concordant direction based on matched group-level Z-score signs

Example pattern now emitted into runner configs:

```text
((G1_P>0) and (G1_P=COMMON_ASSOC_P))
and ((G2_P>0) and (G2_P<0.05))
and (((G1_Z>0) and (G2_Z>0)) or ((G1_Z<0) and (G2_Z<0)))
```

### 2. Independent verifier

Added:

- [DiffGWASDeps/verify_common_association_loci.pl](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/DiffGWASDeps/verify_common_association_loci.pl)
- [verify_common_association_loci.pl](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/verify_common_association_loci.pl)

The verifier reads the wide GWAS file directly and independently checks:

- common-association thresholds
- paired nominal replication
- same-direction requirement
- distance pruning behavior matching `get_top_signal_within_dist.sas`

It uses:

- group-level beta signs by default
- optional `--direction-metric z`

and supports both wide-table layouts seen in this project:

- label-style columns such as `ALL_FEMALE_P`
- pairwise columns such as `CORE_FEMALE_GROUP1_P`

## Tests

### A. Sex-stratified PGC schizophrenia

Spec:

- [spec_pgc_scz_sex_common_automation.json](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/configs/spec_pgc_scz_sex_common_automation.json)

Regenerated runner config:

- [auto_PGC_SCZ_female_vs_male_diff_effects_runner.json](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/configs/auto_PGC_SCZ_female_vs_male_diff_effects_runner.json)

Independent verification outputs:

- [PGC_SCZ_common_assoc_same_dir_5e8.tsv](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/PGC_SCZ_common_assoc_same_dir_5e8.tsv)
- [PGC_SCZ_common_assoc_same_dir_ladder.tsv](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/PGC_SCZ_common_assoc_same_dir_ladder.tsv)

Observed result:

- `0` candidates at `5e-8`
- `0` candidates at `1e-6`
- `0` candidates at `1e-5`

Interpretation:

- under the stricter "nominal replication plus same direction" definition,
  this sex-stratified wide dataset did not retain any common-association loci
  at the current threshold ladder

This is a scientifically useful result, because it shows that the earlier
single-locus common hit was driven by the weaker earlier definition and does
not survive the new concordance rule.

### B. PGC_GWAS4testing 3-GWAS example

Spec:

- [test_PGC_GWAS4testing_all3_common_manual.json](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/configs/test_PGC_GWAS4testing_all3_common_manual.json)

Independent verification outputs:

- [PGC_GWAS4testing_all3_common_same_dir_ladder.tsv](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/PGC_GWAS4testing_all3_common_same_dir_ladder.tsv)
- [PGC_GWAS4testing_all3_common_same_dir_ladder.candidates.tsv](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/PGC_GWAS4testing_all3_common_same_dir_ladder.candidates.tsv)

Observed result:

- raw qualifying candidates across thresholds: `694,217`
- candidates with `COMMON_ASSOC_P < 5e-8`: `23,846`
- selected loci after `1e8` total distance pruning: `49`

Interpretation:

- the new rule is not over-restrictive in general
- it can recover many concordant common loci when the dataset genuinely
  contains them
- the current `top_hit_dist_bp = 1e8` still performs very aggressive locus
  collapsing

## Performance Notes

### What improved

- the verifier runs locally and avoids SAS ODA entirely
- it gives a quick sanity check before expensive local-Manhattan or local-GTF
  plotting

### What remains expensive

- `PGC_GWAS4testing_all3` is still heavy enough that full common-mode local
  Manhattan plotting may exceed the outer shell timeout even when the logic is
  correct

## Issues Encountered

### 1. Schema mismatch across wide datasets

Problem:

- the sex-stratified wide file uses alias columns like `ALL_FEMALE_P`
- the `PGC_GWAS4testing_all3` wide file exposes pairwise columns like
  `CORE_FEMALE_GROUP1_P`

Resolution:

- the verifier now resolves the first matching column from candidate sets
  derived from spec labels and pair prefixes

### 2. Earlier common-mode semantics were too permissive

Problem:

- earlier logic only required nominal association in the paired GWAS

Resolution:

- common-mode now also requires concordant direction in the paired GWAS

## Practical Recommendation

Use the verifier first:

```bash
perl verify_common_association_loci.pl --spec ./configs/spec_pgc_scz_sex_common_automation.json
```

Then adjust thresholds or distance pruning before invoking SAS ODA if needed.

For datasets with many true common loci, consider reducing `top_hit_dist_bp`
from `1e8` to something closer to `1e6` or `5e6` so biologically distinct
loci are not collapsed too aggressively.
