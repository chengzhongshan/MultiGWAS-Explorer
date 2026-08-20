# Agent-interface benchmark

Predeclared task: generate the four deterministic configuration artifacts from `configs/spec_pgc_scz_sex_common_automation.json` in `configs` mode with plotting disabled.

The timed MCP path used an exact JSON-RPC `tools/call`; it did not submit or
score a free-form natural-language question. See `README.md` and
`AI_Interface_Evaluation_Requests.jsonl` for the payload and all expanded
trial records.

- CLI successes: 10/10
- MCP-agent successes: 10/10
- Exact four-artifact checksum parity: 10/10
- CLI elapsed seconds, median (IQR): 0.722 (0.720-0.743)
- MCP-agent elapsed seconds, median (IQR): 2.621 (2.619-2.630)
- Injected failures detected: 3/3
- Corrected retries recovered: 3/3

The MCP timing includes HTTP dispatch, background-process launch, and status polling; the CLI timing does not.

The repeated artifact-set SHA-256 was identical within and between paths.
Failure/recovery timing had median 2.620 s (IQR 2.620-2.670 s).
