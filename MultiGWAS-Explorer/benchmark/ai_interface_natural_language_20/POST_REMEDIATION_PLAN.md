# Post-remediation 20-prompt benchmark plan

Preregistered on 2026-08-29 before execution.

This validation preserves the original 20-prompt benchmark as an immutable initial result (17/20 successful). It reuses the 17 prompts that initially succeeded without changing their wording and replaces only the three failed prompts (`F02_P01`, `F02_P02`, and `F21_P03`) with matched prompts (`F02_R01`, `F02_R02`, and `F21_R03`).

The F02 replacements identify the intended PGC schizophrenia sex-stratified/common-association profile because the original phrase “sex/ancestry specification” was compatible with two different bundled specifications. The F21 replacement retains the intentionally absent path and configuration-only intent. All three use the same evaluator-only family gold actions as the original prompts.

The remediation consists of clearer MCP schema semantics, a read-only bundled-spec catalog, and validation that prevents `generate_spec_only` (GWAS-directory spec inference) from being used as a substitute for `mode=configs` (configuration-only execution of an existing spec).

Primary endpoint: all-criteria success among the 20 post-remediation prompts. Results must be labeled post-remediation and reported alongside, not in place of, the original 17/20 result.
