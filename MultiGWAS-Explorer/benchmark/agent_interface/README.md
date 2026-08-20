# Optional MCP-interface evaluation

This benchmark compares the command-line path with deterministic dispatch
through the local Model Context Protocol (MCP) server. It does **not** submit a
free-form natural-language question to a large language model and does not
score generated answers.

The 10 successful timed MCP trials sent this JSON-RPC request. Only the
two-digit trial number in `output_file` changed:

```json
{"jsonrpc":"2.0","id":100,"method":"tools/call","params":{"name":"auto_prepare_and_run_diff_gwas","arguments":{"spec_file":"configs/spec_pgc_scz_sex_common_automation.json","mode":"configs","skip_plots":"true","output_file":"benchmark/agent_interface/agent_run_XX.log"}}}
```

For readability, the equivalent task was: “Using
`configs/spec_pgc_scz_sex_common_automation.json`, generate the four workflow
configuration files only; do not execute plots; save the run log as
`benchmark/agent_interface/agent_run_XX.log`.” That sentence is an equivalent
description, not an input used for timing.

Three error-handling trials changed `spec_file` to
`configs/INTENTIONALLY_MISSING_REVIEWER_BENCHMARK.json`. Corrected retries
restored the valid specification. Fully expanded records are in
`AI_Interface_Evaluation_Requests.jsonl`.

The timer begins immediately before the HTTP `tools/call` and stops after log
polling detects completion or error. MCP-server startup is excluded. The test
therefore evaluates deterministic tool dispatch, configuration generation,
logging, checksum parity, missing-specification detection, and corrected
retry—not natural-language comprehension, prompt robustness, hallucination,
visual interpretation, or autonomous scientific decisions.

Run the benchmark from the repository root:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass `
  -File benchmark/run_agent_interface_benchmark.ps1
```

Recorded result: CLI 10/10, MCP 10/10, exact four-artifact parity 10/10,
intentionally invalid requests detected 3/3, and corrected retries successful
3/3. See `SUMMARY.md` and the TSV files in this directory.
