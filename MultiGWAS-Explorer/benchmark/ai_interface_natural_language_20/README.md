# Natural-language Codex benchmark

This directory contains the machine-readable evidence for the focused
natural-language evaluation of the optional MultiGWAS-Explorer AI interface.
It is separate from `benchmark/agent_interface/`, which measures deterministic
JSON-RPC dispatch without a language-model interpretation step.

## Design

- A legacy bank of 100 requests contained 25 task families with four
  ChatGPT-generated paraphrases per family.
- The resource-conscious evaluation selected 20 prompts from five complete
  families: environment validation (F01), PGC configuration-only generation
  (F02), the synthetic MAF-safeguard regression (F13), missing-specification
  safe failure (F21), and contradictory-backend clarification (F25).
- UK Biobank families were outside the current scientific focus. F03 was not
  rerun because its forced preprocessing is time-consuming and had previously
  completed successfully.
- Each primary prompt was executed once in a fresh headless Codex session.
  Evaluator-only gold actions were hidden from Codex. Positive requests had to
  select the preregistered action and reach a successful terminal outcome;
  invalid and contradictory requests were scored for safe failure and
  clarification without execution, respectively.
- Task-concordant F02 outputs were canonicalized and compared with four direct
  command-line reference JSON artifacts by SHA-256.

## Results

The immutable initial benchmark met all criteria for 17/20 prompts (85.0%;
95% Wilson CI, 64.0%-94.8%). Before its execution, the matched remediation set
was defined to retain 17 prompt wordings and replace the three failed requests
with clearer requests testing the same actions. That run met 19/20 criteria
(95.0%; 95% CI, 76.4%-99.1%), including:

- positive-task completion: 11/12;
- prespecified action concordance: 19/20;
- missing-specification safe failure: 4/4;
- contradictory-request clarification: 4/4;
- canonical artifact parity: 3/3 assessable task-concordant runs; and
- median latency: 32.317 s (IQR 26.066-48.422; p95 62.922), with no timeouts
  or agent-initiated source edits.

The only unmatched post-remediation prompt said only "the bundled PGC
schizophrenia workflow" even though two scientifically different bundled
specifications were available. Codex requested clarification instead of
silently selecting one. The explicit-spec diagnostic `F02_R04` was therefore
reported outside the 20-prompt denominator. It passed in the original probe
and again in a fresh session against clean published commit
`1d604c4f5ae8aa33d5a6d2863020a07cd6e94074` (2/2); each run matched all four
reference artifacts (8/8 artifact-file comparisons in total).

## Interpretation

These results evaluate orchestration: whether the agent maps a request to the
intended versioned repository action. They do not show that a language model
performs or improves the GWAS statistics. The primary prompts were each run
once, so repeated-run consistency for those sets was not assessed. The five
selected families do not exercise every pipeline feature, particularly the
long-running SAS ODA and genome-wide rendering paths, and this was not a
head-to-head comparison of foundation models. The same workflows remain
available directly through the documented command-line interface.

## Files

- `PROMPTS_AND_RESULTS.md`: exact selected prompts, scoring rules, results,
  failed-prompt audit, and limitations.
- `post_remediation_20_prompts.jsonl`: exact matched post-remediation prompt
  set.
- `explicit_spec_F02_R04.jsonl`: explicit-spec diagnostic prompt.
- `initial_*`: initial 20-prompt results and frozen preflight metadata.
- `post_remediation_*`: matched post-remediation results and metadata.
- `explicit_spec_repeat_*`: repeat diagnostic results from clean commit
  `1d604c4`.
- `gold_actions.json` and `run_codex_gemini_benchmark.pl`: evaluator action
  definitions and the frozen Perl runner snapshot used to score the runs.
- `SHA256SUMS.tsv`: byte sizes and SHA-256 hashes for the versioned audit
  bundle.

Absolute run paths in the CSV files are provenance records from the evaluation
machine; they are not required to execute the pipeline.
