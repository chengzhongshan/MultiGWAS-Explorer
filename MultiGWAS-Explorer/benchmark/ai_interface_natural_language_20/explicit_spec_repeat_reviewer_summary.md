# Reviewer-response benchmark summary

Repository commit: `1d604c4f5ae8aa33d5a6d2863020a07cd6e94074`
Evaluated prompts: 1 across 1 task family
Replicates represented: 1

## Ready-to-adapt response to the reviewer

We conducted a focused AI-interface evaluation using frozen repository commit `1d604c4f5ae8aa33d5a6d2863020a07cd6e94074`. We evaluated 1 natural-language prompts representing 1 task family with fresh Codex sessions. The command-line workflow was treated as the scientific reference. We recorded whether Codex selected the prespecified repository action, whether execution reached a terminal outcome, source-edit violations, and end-to-end wall time. For families configured with explicit scientific artifact paths, task-concordant outputs were additionally canonicalized and compared with the corresponding CLI reference hashes.

**Codex.** Overall benchmark success was 100.0% (20.7-100.0%); positive-task completion was 100.0% (20.7-100.0%); prespecified action/invocation concordance was 100.0% (20.7-100.0%); median end-to-end latency was 61.259 s (IQR 61.259-61.259 s; p95 61.259 s). Canonical artifact parity was observed in 1/1 task-concordant assessable runs.

### MCP subset

The `preferred_interface=mcp` task families are reported separately to characterize Codex orchestration through the local MCP server.
**Codex MCP subset.** Positive-task completion was 100.0% (20.7-100.0%); prespecified invocation concordance was 100.0% (20.7-100.0%); median end-to-end latency was 61.259 s (IQR 61.259-61.259 s).

Repeated-run consistency was not assessed because this focused benchmark used one repetition per prompt.
Safe-failure detection was not assessed because the selected prompts contained no deliberately invalid requests.

These results evaluate the reliability of the AI orchestration layer in invoking the versioned MultiGWAS-Explorer workflow; they do not imply that the language model performs or improves the underlying GWAS statistical calculations.

## Suggested manuscript Methods sentence

The AI interface was evaluated with 1 preregistered natural-language prompts spanning 1 task family. Codex was run in fresh headless sessions against the same frozen repository commit `1d604c4f5ae8aa33d5a6d2863020a07cd6e94074` and local MCP endpoint. The direct command-line workflow served as the reference. The experiment quantified prespecified action concordance, task completion, source-edit violations, and end-to-end execution time; task-concordant outputs were canonicalized and compared with the CLI reference by SHA-256 hash. Repeated-run consistency was not assessed because each prompt was executed once. Safe-failure behavior was not assessed because no deliberately invalid prompt was selected.
