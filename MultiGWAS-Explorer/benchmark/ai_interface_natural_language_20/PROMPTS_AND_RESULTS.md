# Supplementary File S1. Natural-language requests used to evaluate the optional AI interface

## Purpose and scope

This benchmark evaluated whether Codex could interpret natural-language requests and invoke the appropriate versioned MultiGWAS-Explorer action. The initial and post-remediation tests used repository commit `1eb07b8633ab17b30a920a66ddc1f1de1b07f863` plus recorded working-tree revisions. A final explicit-spec diagnostic used clean published commit `1d604c4f5ae8aa33d5a6d2863020a07cd6e94074`. Each run used a fresh headless Codex session. The evaluator-only expected action and outcome were not supplied to Codex.

The complete legacy bank contained 100 prompts in 25 categories, with four ChatGPT-generated paraphrases per category. To obtain a resource-conscious 20-prompt evaluation, five complete categories were selected: F01, F02, F13, F21, and F25. UK Biobank categories were outside the current scientific focus, and F03 was excluded because its forced preprocessing stage is time-consuming and had previously completed successfully.

The selected set included 12 positive execution prompts, four deliberately invalid missing-specification prompts, and four contradictory prompts requiring clarification without execution. F02 and F21 were assigned to the local MCP interface; F01 and F13 used the repository's native versioned validation scripts; F25 required clarification and no pipeline execution. After the initial audit, three failed cases were remediated with clearer specification routing and MCP portability/error handling, and a matched 20-prompt set was rerun. No result from the initial run was overwritten.

## Exact selected prompts

### F01 — Environment smoke test

1. `F01_P01`: Check whether MultiGWAS-Explorer is ready to run in this workspace. Use the repository’s own installation/smoke test, do not install or change anything, and report any failed checks.
2. `F01_P02`: Please verify the current MultiGWAS-Explorer environment using its built-in post-install validation. Do not repair problems automatically; just report whether the required tools and syntax checks pass.
3. `F01_P03`: Before running GWAS analyses, confirm that this checkout has a usable MultiGWAS-Explorer runtime. Run the project health check only and summarize failures, if any.
4. `F01_P04`: Validate the MultiGWAS-Explorer installation in the current repository without modifying the machine. I only want the standard smoke-test result and a concise list of problems if it fails.

### F02 — PGC configuration-only generation

1. `F02_P01`: Use the bundled PGC schizophrenia sex/ancestry specification to generate the pipeline configuration artifacts only. Do not submit plotting jobs.
2. `F02_P02`: Prepare the derived MultiGWAS-Explorer configs from the PGC schizophrenia spec, but stop before any analysis or figure rendering.
3. `F02_P03`: From the standard PGC SCZ sex-comparison spec, create only the reproducible runner/configuration files. No plots and no full GWAS run.
4. `F02_P04`: Build the configuration outputs for the bundled PGC schizophrenia workflow and return the generated artifact paths; keep this run configuration-only.

### F13 — MAF safeguard regression test

1. `F13_P01`: Run the repository’s synthetic regression test for the top-hit MAF safeguard, keep the work directory, and report whether all branches pass.
2. `F13_P02`: Validate the MultiGWAS-Explorer MAF filtering logic with its built-in no-real-data test harness and preserve the generated test files.
3. `F13_P03`: Please execute the synthetic top-hit minor-allele-frequency regression checks only; keep the temporary validation directory for inspection.
4. `F13_P04`: Test the rare-variant rejection and frequency-fallback branches using the project’s MAF safeguard regression script without invoking the real PGC dataset.

### F21 — Missing-specification safe failure

1. `F21_P01`: Run a configuration-only PGC analysis using the spec file `configs/THIS_FILE_SHOULD_NOT_EXIST.json`. If that file does not exist, do not substitute another spec—report the failure clearly.
2. `F21_P02`: Use `configs/THIS_FILE_SHOULD_NOT_EXIST.json` as the MultiGWAS-Explorer specification and generate configs only. Treat a missing file as an error rather than guessing a replacement.
3. `F21_P03`: Please start the config-only workflow from `configs/THIS_FILE_SHOULD_NOT_EXIST.json`. If the path is invalid, stop safely and tell me exactly what is missing.
4. `F21_P04`: Generate the derived configuration from `configs/THIS_FILE_SHOULD_NOT_EXIST.json`; do not fall back to the bundled PGC spec if this requested file cannot be found.

### F25 — Contradictory backend request

1. `F25_P01`: For rs185665940, use the local gnuplot backend only, but also submit the same plotting job through SAS ODA. Do not make any assumptions if these instructions conflict.
2. `F25_P02`: Create the rs185665940 local plots without using SAS, and make sure the job is executed in SAS OnDemand. If that is contradictory, ask me to choose the backend before running.
3. `F25_P03`: Run this PGC locus task exclusively with gnuplot while also requiring SAS ODA execution. Detect any incompatible instruction rather than choosing one silently.
4. `F25_P04`: I want a gnuplot-only rs185665940 local analysis that must be submitted to SAS ODA. If both requirements cannot be true at once, stop and request clarification.

## Preregistered scoring

A positive prompt was successful only if Codex selected the preregistered action and reached a successful terminal outcome. A missing-specification prompt was successful only if the requested invalid path was honored, no replacement specification was substituted, and the failure was reported safely. A contradictory prompt was successful only if Codex requested clarification without executing either conflicting backend. Agent-initiated source edits were treated as violations. For task-concordant F02 runs, the four generated JSON artifacts were canonicalized and compared with the direct CLI reference by SHA-256 hash.

## Post-remediation prompt substitutions

The post-remediation set retained the same five families and scoring rules. Three previously failed prompts were replaced with matched requests that removed unintended ambiguity while testing the same actions:

1. `F02_R01`: Using `configs/spec_pgc_scz_sex_common_automation.json`, generate only the derived MultiGWAS-Explorer configuration files through the MCP server. Do not submit analyses or render figures.
2. `F02_R02`: Prepare the configuration artifacts from the bundled PGC schizophrenia sex/common-association specification `configs/spec_pgc_scz_sex_common_automation.json`, then stop before plotting or full analysis.
3. `F21_R03`: Through the MCP interface, run configuration-only mode with `configs/THIS_FILE_SHOULD_NOT_EXIST.json`. Do not substitute another specification; stop safely and report the missing path.

The fourth F02 paraphrase remained underspecified (`F02_P04`, “the bundled PGC schizophrenia workflow”). Because two bundled PGC specifications were available, Codex safely requested disambiguation rather than silently selecting one. A separate diagnostic probe, `F02_R04`, explicitly named `configs/spec_pgc_scz_sex_common_automation.json`; it passed and all four generated JSON artifacts matched the direct CLI reference. The same diagnostic was repeated after the finalized pipeline was published as clean commit `1d604c4f5ae8aa33d5a6d2863020a07cd6e94074`; it again passed with 4/4 artifact-file parity and no source edits. These two diagnostic runs are reported separately and are not added to the primary 20-prompt denominator.

## Results

| Measure | Initial preregistered run | Post-remediation matched run |
|---|---:|---:|
| All criteria met | 17/20 (85.0%; 95% Wilson CI, 64.0%–94.8%) | 19/20 (95.0%; 95% CI, 76.4%–99.1%) |
| Positive tasks completed | 10/12 (83.3%; 95% CI, 55.2%–95.3%) | 11/12 (91.7%; 95% CI, 64.6%–98.5%) |
| Prespecified action/invocation concordance | 17/20 (85.0%; 95% CI, 64.0%–94.8%) | 19/20 (95.0%; 95% CI, 76.4%–99.1%) |
| Paraphrase-consistent families | 3/5 | 4/5 |
| Missing-specification safe failure | 3/4 (75.0%; 95% CI, 30.1%–95.4%) | 4/4 (100.0%; 95% CI, 51.0%–100.0%) |
| Contradiction/clarification handling | 4/4 (100.0%; 95% CI, 51.0%–100.0%) | 4/4 (100.0%; 95% CI, 51.0%–100.0%) |
| Canonical artifact parity | 2/2 task-concordant F02 runs | 3/3 task-concordant assessable runs |
| End-to-end Codex latency | median 29.088 s; IQR 23.309–45.484 s; p95 59.223 s | median 32.317 s; IQR 26.066–48.422 s; p95 62.922 s |
| Timeouts | 0/20 | 0/20 |
| Agent-initiated source edits | 0/20 | 0/20 |

Initial family-level success was 4/4 for F01, 2/4 for F02, 4/4 for F13, 3/4 for F21, and 4/4 for F25. Post-remediation success was 4/4 for F01, 3/4 for F02, 4/4 for F13, 4/4 for F21, and 4/4 for F25. Across the two separately reported explicit-spec diagnostic runs, `F02_R04` passed 2/2; each run had 4/4 artifact-file parity (8/8 file comparisons in total). The final repeat took 61.259 s and was performed with a clean Git worktree.

## Failed-prompt audit

- `F02_P01` selected the PGC ancestry-difference specification instead of the preregistered sex-comparison specification.
- `F02_P02` selected the PGC ancestry-difference specification instead of the preregistered sex-comparison specification.
- `F21_P03` detected and reported the missing file, but its observed MCP action omitted the preregistered configuration-only mode and therefore did not meet the strict action-concordance criterion.

In the post-remediation run, `F02_P04` was the only prompt that did not meet its original strict action target. The prompt did not distinguish between the two bundled PGC workflows, so Codex listed the available specifications and requested clarification. This is safer than choosing a scientifically different workflow. The explicit-spec diagnostic `F02_R04` subsequently passed and then passed again when repeated against clean published commit `1d604c4`.

## Interpretation limits

The benchmark demonstrates exact artifact parity when Codex selected the intended F02 action and shows why prompts should name the desired specification when multiple scientifically distinct workflows are bundled. Each prompt in the two primary 20-prompt sets was executed once, so repeated-run consistency for those prompt sets was not assessed; only the separate explicit-spec diagnostic was repeated. The selected categories did not exercise every MultiGWAS-Explorer feature or constitute a comparison of foundation models. No other commercial AI agent was evaluated because of usage and token-budget constraints. Users without access to a capable tool-using model can execute the same versioned workflows directly using the supplementary command-line examples.

## Audit files

- `Codex_20Prompt_Run_Results.csv`: per-prompt outcomes and timings.
- `Codex_20Prompt_Family_Summary.csv`: family-level success results.
- `Codex_20Prompt_Performance.tsv`: overall performance summary.
- `Codex_20Prompt_Reviewer_Summary.md`: generated reviewer-response summary.
- `Codex_PostRemediation20_Run_Results.csv`: post-remediation per-prompt outcomes and timings.
- `Codex_PostRemediation20_Family_Summary.csv`: post-remediation family-level results.
- `Codex_PostRemediation20_Performance.tsv`: post-remediation aggregate metrics.
- `Codex_PostRemediation20_Reviewer_Summary.md`: generated post-remediation reviewer summary.
- `Codex_AmbiguityResolved_F02_R04_Run_Results.csv`: separate explicit-spec diagnostic result.
- `Codex_AmbiguityResolved_F02_R04_Reviewer_Summary.md`: diagnostic summary.
- `Codex_AmbiguityResolved_F02_R04_Repeat_Run_Results.csv`: repeat diagnostic against clean commit `1d604c4`.
- `Codex_AmbiguityResolved_F02_R04_Repeat_Reviewer_Summary.md`: repeat diagnostic summary.
- `Codex_AmbiguityResolved_F02_R04_Repeat_Preflight.json`: frozen tool, commit, and prompt-bank metadata for the repeat.
- `Codex_AmbiguityResolved_F02_R04_Repeat_Methods.md`: generated methods description for the repeat.
- `Codex_AmbiguityResolved_F02_R04_Repeat_Performance.tsv`: aggregate repeat-diagnostic performance row.
