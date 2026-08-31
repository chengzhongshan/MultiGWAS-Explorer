# Reviewer-response benchmark summary

Repository commit: `1eb07b8633ab17b30a920a66ddc1f1de1b07f863`
Tracked working-tree diff SHA-256: `c9f4ee1b890313577703ab88d755b940b2e97a24e70c2285575cc7f02cfe2cef`
Evaluated prompts: 20 across 5 task families
Replicates represented: 1

## Ready-to-adapt response to the reviewer

We conducted a focused AI-interface evaluation using repository commit `1eb07b8633ab17b30a920a66ddc1f1de1b07f863` plus recorded working-tree modifications (tracked-diff SHA-256 `c9f4ee1b890313577703ab88d755b940b2e97a24e70c2285575cc7f02cfe2cef`). We evaluated 20 natural-language prompts representing 5 task families with fresh Codex sessions. The command-line workflow was treated as the scientific reference. We recorded whether Codex selected the prespecified repository action, whether execution reached a terminal outcome, source-edit violations, and end-to-end wall time. For families configured with explicit scientific artifact paths, task-concordant outputs were additionally canonicalized and compared with the corresponding CLI reference hashes.

**Codex.** Overall benchmark success was 95.0% (76.4-99.1%); positive-task completion was 91.7% (64.6-98.5%); prespecified action/invocation concordance was 95.0% (76.4-99.1%); paraphrase consistency was 80.0% (37.6-96.4%); safe-failure detection was 100.0% (51.0-100.0%); clarification handling was 100.0% (51.0-100.0%); median end-to-end latency was 32.317 s (IQR 26.066-48.422 s; p95 62.922 s). Canonical artifact parity was observed in 3/3 task-concordant assessable runs.

### MCP subset

The `preferred_interface=mcp` task families are reported separately to characterize Codex orchestration through the local MCP server.
**Codex MCP subset.** Positive-task completion was 75.0% (30.1-95.4%); prespecified invocation concordance was 87.5% (52.9-97.8%); safe-failure detection was 100.0% (51.0-100.0%); median end-to-end latency was 27.818 s (IQR 26.905-32.609 s).

## Failed-prompt diagnostics

- `F02_P04`: the executed action did not match the preregistered action. The first observed action was `sas_pipeline`.

Repeated-run consistency was not assessed because this focused benchmark used one repetition per prompt.

These results evaluate the reliability of the AI orchestration layer in invoking the versioned MultiGWAS-Explorer workflow; they do not imply that the language model performs or improves the underlying GWAS statistical calculations.

## Suggested manuscript Methods sentence

The AI interface was evaluated with 20 preregistered natural-language prompts spanning 5 task families. Codex was run in fresh headless sessions against the same repository commit `1eb07b8633ab17b30a920a66ddc1f1de1b07f863` plus recorded working-tree modifications (tracked-diff SHA-256 `c9f4ee1b890313577703ab88d755b940b2e97a24e70c2285575cc7f02cfe2cef`) and local MCP endpoint. The direct command-line workflow served as the reference. The experiment quantified prespecified action concordance, task completion, source-edit violations, and end-to-end execution time; task-concordant outputs were canonicalized and compared with the CLI reference by SHA-256 hash. Repeated-run consistency was not assessed because each prompt was executed once. Deliberately invalid prompts were scored for safe-failure behavior. Contradictory prompts were required to request clarification without execution.
