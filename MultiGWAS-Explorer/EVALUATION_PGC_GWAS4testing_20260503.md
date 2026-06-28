# Evaluation: PGC_GWAS4testing Pipeline and MCP Integration

Date: 2026-05-03

Dataset under evaluation:

- `/mnt/e/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_GWAS4testing/`

Scope:

- `auto_prepare_and_run_diff_gwas.pl`
- local SAS ODA plot runners in `DiffGWASDeps/`
- Perl MCP server exposure through `server.pl`
- documentation and workflow skill coverage

## Summary

The refactored pipeline is in a usable and mostly robust state for the
`PGC_GWAS4testing` schizophrenia test set.

What is solid:

- differential-analysis pipeline generation and rerun controls
- wide-subset generation and output correctness
- differential local Manhattan rendering and split-panel output
- top-hit CSV export with stacked single-GWAS and differential signals
- MCP wrapper support for targeted reruns and common-association mode arguments

What remains externally fragile:

- common-association local Manhattan rendering depends on SAS ODA availability
- recent common-mode failures were caused by remote SAS ODA login errors, not by
  local script syntax

## Inputs Evaluated

The 3-GWAS manual test spec:

- [test_PGC_GWAS4testing_all3_manual.json](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/configs/test_PGC_GWAS4testing_all3_manual.json)

The common-association variant:

- [test_PGC_GWAS4testing_all3_common_manual.json](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/configs/test_PGC_GWAS4testing_all3_common_manual.json)

## Performance

### Differential path

The pipeline produces a large but manageable wide subset for plotting:

- wide manifest:
  [PGC_GWAS4testing_all3_diff_effects.stdized.wide_beta_se_p_p_lt_0p05.final.manifest.tsv](</E:/LongCOVID_HGI_GWAS/PGC_Large_GWASs/PGC_GWAS4testing/PGC_GWAS4testing_all3_diff_effects.stdized.wide_beta_se_p_p_lt_0p05.final.manifest.tsv>)

Key values from that manifest:

- `rows_read`: `21823411`
- `rows_written`: `2016998`
- `groups_seen`: `7320339`
- `pair_prefixes`: `CORE_FEMALE,CORE_MALE,FEMALE_MALE`

The MCP-driven standardization step also completed successfully on this test
set and logged:

- `Numeric N`: `22734149`
- `Rows written`: `22734149`

Source:

- [server_8083.err.log](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/server_8083.err.log)

### Local top-hit selection performance

The main performance improvement implemented and evaluated here is:

1. do not call `get_top_signal_within_dist` once per threshold in the
   threshold ladder
2. count cheap candidates first
3. call the expensive top-hit deduplication macro only once with the chosen
   threshold

The macro itself was also refactored:

- [get_top_signal_within_dist.sas](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/DiffGWASDeps/get_top_signal_within_dist.sas)

Instead of the older broad SQL overlap pattern, it now:

- filters candidates first
- sorts once
- performs a grouped greedy interval-selection pass

This is the right direction for cases where many loci pass the initial
association filter.

### Common-association path

The common-association path is still the most sensitive performance area,
because it depends on:

- many more candidate loci passing the single-GWAS association threshold
- SAS ODA availability during the actual plot render

During evaluation, one representative common-mode run took a long time and
ultimately failed because SAS ODA itself rejected the session:

- [tmp1777830819/output.html.info.txt](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/tmp1777830819/output.html.info.txt)

The key failure was:

- `The application could not log on to the server "odaws01-usw2.oda.sas.com:8591". The server configuration is invalid.`

So for common-mode performance, the local code is improved, but full-end
timings are still confounded by remote ODA instability.

## Usability

### Strengths

The current user-facing workflow is substantially better than the earlier
state:

- step-aware reruns are available
- common-association mode is exposed through both CLI and MCP
- local panel width is configurable
- common-mode differential non-significance is configurable
- docs now explain that `--get-common-associations=THR` applies to
  single-GWAS association P, not differential P

Relevant files:

- [auto_prepare_and_run_diff_gwas.pl](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/auto_prepare_and_run_diff_gwas.pl)
- [server.pl](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/server.pl)
- [PerlMCP_Server_Instruction4CODEX.txt](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/PerlMCP_Server_Instruction4CODEX.txt)
- [GENERALIZE_DIFF_GWAS_PIPELINE.md](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/GENERALIZE_DIFF_GWAS_PIPELINE.md)
- [SKILL.md](</C:/Users/cheng/.codex/skills/gwas-sas-oda-workflow/SKILL.md>)

### MCP integration status

The MCP surface is present and documented.

The server schema now exposes:

- `get_common_associations`
- `common_association_top_hit_threshold`
- `common_association_diff_nonsig_threshold`
- `local_max_hits_per_fig`
- step-aware rerun arguments such as `list_steps`, `step`, `from_step`,
  `to_step`

During evaluation, the MCP config-generation path worked and produced the
expected config outputs for the PGC schizophrenia test specs.

### Usability issue found and fixed

When SAS ODA failed upstream, the external wrapper could still emit a misleading
"SAS job is completed!" line and uninitialized-value warnings.

That was cleaned up in:

- [run_sas_codes_or_script_in_ODA.pl](</g:/NGS_lib/Linux_codes_SAM/run_sas_codes_or_script_in_ODA.pl>)

The wrapper now:

- avoids undefined-value warnings around missing `htmlfilename`
- reports a clearer message when the job finishes without a downloadable HTML
  artifact

This makes failure logs much easier to interpret.

## Correctness

### Differential outputs

The differential local Manhattan outputs are present and coherent:

- [PGC_GWAS4testing_all3_diff_SAS_local_top_hits_manhattan.html](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/PGC_GWAS4testing_all3_diff_SAS_local_top_hits_manhattan.html)
- [PGC_GWAS4testing_all3_diff_SAS_local_top_hits_manhattan.png](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/PGC_GWAS4testing_all3_diff_SAS_local_top_hits_manhattan.png)
- parts 2 through 10 were also generated

The CSV export is also present:

- [PGC_GWAS4testing_all3_diff_SAS_local_top_hits_manhattan_top_hits.csv](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/PGC_GWAS4testing_all3_diff_SAS_local_top_hits_manhattan_top_hits.csv)

It contains:

- hit ordering and panel assignment
- SNP / chromosome / position / allele context
- per-GWAS association signals
- pairwise differential signals
- standardized differential signals
- matching Z-score fields

This is consistent with the intended "plot plus shareable summary table"
behavior.

### Common-association logic

The common-association logic was refined during this evaluation.

The generated runner config now correctly shows:

- `TOP_HIT_FOCUS_PVAR = COMMON_ASSOC_P`
- `COMMON_ASSOC_P = min(single-GWAS association P variables)`
- one GWAS must provide that strongest signal
- another GWAS must be nominally associated, default `< 0.05`
- at least one pairwise standardized differential P must be relaxed
  non-significant, default `>= 0.5`

Verified in:

- [auto_PGC_GWAS4testing_all3_diff_effects_runner.json](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/configs/auto_PGC_GWAS4testing_all3_diff_effects_runner.json)

### Remaining correctness caveat

Because the recent common-mode render failed at the SAS ODA login layer, we do
not yet have a final successful common-mode PNG from this evaluation pass.

That means:

- common-mode config generation and runner wiring are validated
- common-mode final render correctness is not fully re-confirmed in this pass
  because of the remote ODA outage

## Changes Made During This Evaluation

1. Common-association semantics corrected in
   [auto_prepare_and_run_diff_gwas.pl](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/auto_prepare_and_run_diff_gwas.pl)
   so that:
   - the optional threshold is clearly a single-GWAS association threshold
   - the relaxed non-differential default is `0.5`

2. MCP schema updated in
   [server.pl](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/server.pl)
   to expose:
   - `common_association_diff_nonsig_threshold`

3. Failure-path usability improved in
   [run_sas_codes_or_script_in_ODA.pl](</g:/NGS_lib/Linux_codes_SAM/run_sas_codes_or_script_in_ODA.pl>)
   by removing misleading completion output when HTML artifacts are missing

4. Documentation refreshed in:
   - [PerlMCP_Server_Instruction4CODEX.txt](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/PerlMCP_Server_Instruction4CODEX.txt)
   - [GENERALIZE_DIFF_GWAS_PIPELINE.md](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/GENERALIZE_DIFF_GWAS_PIPELINE.md)
   - [Add_perlmcp_server2Codex.sh](/g:/NGS_lib/Linux_codes_SAM/Conda_and_Docker_Related_Scripts/perlMCP4Gemini_Paper/Add_perlmcp_server2Codex.sh)
   - [SKILL.md](</C:/Users/cheng/.codex/skills/gwas-sas-oda-workflow/SKILL.md>)

## Overall Assessment

Performance:

- good for differential-mode preparation and plotting
- improved for local top-hit candidate selection
- still dependent on SAS ODA health for end-to-end common-mode runs

Usability:

- strong improvement from earlier revisions
- MCP exposure is now real and documented
- targeted reruns and common-mode controls make the tool much easier to use

Correctness:

- differential-mode outputs look correct and complete for this test set
- common-mode selection logic is now better aligned with the intended biology
- final common-mode PNG generation still needs one clean SAS ODA pass after the
  remote service stabilizes

## Recommended Next Steps

1. Re-run the common-association local Manhattan step once SAS ODA logins are
   stable again, using:
   - `--get-common-associations=5e-8`
   - `--common-association-diff-nonsig-threshold 0.5`

2. Add a lightweight dry-run summary mode that reports:
   - candidate counts per common-mode threshold
   - number of loci retained after the relaxed non-differential filter
   before invoking SAS ODA

3. Restart any long-lived MCP daemon after schema changes so the newest
   arguments are actually available to clients.
