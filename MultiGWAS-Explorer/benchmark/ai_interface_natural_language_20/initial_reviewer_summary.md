# Reviewer-response benchmark summary

Repository commit: `1eb07b8633ab17b30a920a66ddc1f1de1b07f863`
Evaluated prompts: 20 across 5 task families
Replicates represented: 1

## Ready-to-adapt response to the reviewer

We conducted a focused AI-interface evaluation using a frozen MultiGWAS-Explorer revision (`1eb07b8633ab17b30a920a66ddc1f1de1b07f863`). We evaluated 20 natural-language prompts representing 5 task families with fresh Codex sessions. The command-line workflow was treated as the scientific reference. We recorded whether Codex selected the prespecified repository action, whether execution reached a terminal outcome, source-edit violations, and end-to-end wall time. For families configured with explicit scientific artifact paths, task-concordant outputs were additionally canonicalized and compared with the corresponding CLI reference hashes.

**Codex.** Overall benchmark success was 85.0% (64.0-94.8%); positive-task completion was 83.3% (55.2-95.3%); prespecified action/invocation concordance was 85.0% (64.0-94.8%); paraphrase consistency was 60.0% (23.1-88.2%); safe-failure detection was 75.0% (30.1-95.4%); clarification handling was 100.0% (51.0-100.0%); median end-to-end latency was 29.088 s (IQR 23.309-45.484 s; p95 59.223 s). Canonical artifact parity was observed in 2/2 task-concordant assessable runs.

### MCP subset

The `preferred_interface=mcp` task families are reported separately to characterize Codex orchestration through the local MCP server.
**Codex MCP subset.** Positive-task completion was 50.0% (15.0-85.0%); prespecified invocation concordance was 62.5% (30.6-86.3%); safe-failure detection was 75.0% (30.1-95.4%); median end-to-end latency was 33.027 s (IQR 24.412-42.830 s).

## Failed-prompt diagnostics

- `F02_P01`: the executed action did not match the preregistered action. The first observed action was `sas_pipeline` with spec `E:/LongCOVID_HGI_GWAS/PGC_Large_GWASs/Manuscript4AI_with_PGC_GWASs/FR_Reviewers_Comments/perlMCP4Gemini_Paper/configs/spec_pgc_scz_ancestry_diff_automation.json`.
- `F02_P02`: the executed action did not match the preregistered action. The first observed action was `sas_pipeline` with spec `configs/spec_pgc_scz_ancestry_diff_automation.json`.
- `F21_P03`: the executed action did not match the preregistered action. The first observed action was `sas_pipeline` with spec `configs/THIS_FILE_SHOULD_NOT_EXIST.json`.

Repeated-run consistency was not assessed because this focused benchmark used one repetition per prompt.

These results evaluate the reliability of the AI orchestration layer in invoking the versioned MultiGWAS-Explorer workflow; they do not imply that the language model performs or improves the underlying GWAS statistical calculations.

## Suggested manuscript Methods sentence

The AI interface was evaluated with 20 preregistered natural-language prompts spanning 5 task families. Codex was run in fresh headless sessions against the same frozen repository and local MCP endpoint. The direct command-line workflow served as the reference. The experiment quantified prespecified action concordance, task completion, source-edit violations, and end-to-end execution time; task-concordant outputs were canonicalized and compared with the CLI reference by SHA-256 hash. Repeated-run consistency was not assessed because each prompt was executed once. Deliberately invalid prompts were scored for safe-failure behavior. Contradictory prompts were required to request clarification without execution.
