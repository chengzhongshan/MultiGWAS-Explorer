# SAS ODA reliability and cleanup

The tested SAS ODA improvements from MCP4SAS have been ported to the explorer's
CLI runner and its `run_sas_codes_or_script_in_ODA` MCP tool.

## Changes

- Macro jump labels such as `%finished:` no longer trigger remote macro loaders.
- `ERROR:` inside an MPRINT source trace or ordinary data no longer produces a
  false failure. Actual SAS errors and unresolved-macro warnings still fail.
- Standalone CLI runs use owned `.sas_oda_work_*` scratch directories and remove
  them on exit, including handled failures and timeouts. Logs, HTML, status JSON,
  and downloads remain in their requested output directories.
- Timeout termination snapshots the worker's descendants and translates Cygwin
  PIDs to Windows PIDs before calling `taskkill`. This prevents surviving helpers
  from writing into a deleted scratch directory.
- Long-lived SASPy services use a stable temporary location instead of inheriting
  the short-lived CLI scratch path.
- The explorer's SAS MCP callback now uses `DiffGWASDeps/SAS_ODA_MCP.pm`.
  Completed jobs return an exit code and an MCP error flag. Missing completion
  records report `UNKNOWN` rather than unconfirmed success.
- MCP jobs retain logs and HTML in `SAS_Output_ODA_*`, remove their `tmp*` worker
  scratch on completion, and remove PID/completion markers after polling.
  Custom log destinations receive a copy of the diagnostic log.

The differential-GWAS pipeline, scientific defaults, other MCP tools, and existing
SAS macros remain part of the existing explorer implementation.

## Perl and Bash commands in Cygwin

From an updated explorer checkout:

```bash
export PATH=/usr/bin:$PATH
export SAS_ODA_AUTOLOAD_MACROS=0
export OPEN_RESULT=0

perl run_sas_codes_or_script_in_ODA.pl \
  --code 'proc means data=sashelp.class n mean; var age; run;' \
  --output-prefix SAS_Output_summary
```

Disable global macro autoload only for self-contained examples. Plotting scripts
that use the remote macro library can continue using their normal settings.

To retain surviving intermediates for debugging:

```bash
perl run_sas_codes_or_script_in_ODA.pl \
  --code 'proc print data=work.missing_table; run;' \
  --output-prefix SAS_Output_debug --keep-intermediate-files
```

`SAS_ODA_KEEP_INTERMEDIATES=1` also enables standalone debug retention.
Regular per-file cleanup inside helpers still occurs.

Preview and remove inactive, marked scratch directories older than 24 hours:

```bash
perl cleanup_sas_oda_intermediates.pl --directory .
perl cleanup_sas_oda_intermediates.pl --directory . --apply
```

Use `--older-than-hours 0` to include recent retained scratch. The cleanup tool
verifies ownership, resolved parent paths, age, and PID. Active, unmarked, and
symlink directories are preserved. This command targets `.sas_oda_work_*`, not
final outputs or unrelated legacy directories.

## Repeatable checks

### Compatibility filename and file-based smoke test

The canonical runner is `run_sas_codes_or_script_in_ODA.pl`. The compatibility
entry point `run_sas_codes_or_files_in_SAS_ODA.pl` now delegates to that runner,
preserving arguments, standard input, and the exit code. Both filenames are
available in the three active copies.

From your explorer or benchmark checkout in Cygwin:

```bash
export PATH=/usr/bin:$PATH
SAS_ODA_AUTOLOAD_MACROS=0 OPEN_RESULT=0 \
  perl run_sas_codes_or_files_in_SAS_ODA.pl \
  --file examples/sas_oda_runner_smoke.sas \
  --output-prefix SAS_Output_runner_file_verified --run-timeout-seconds 120

perl tests/run_live_cli_cleanup.pl --live \
  --runner run_sas_codes_or_files_in_SAS_ODA.pl
```

The compatibility filename passed all 14 original live cleanup assertions
(`SAS_Output_cli_cleanup_6f4t7T/`). The generated file-based smoke test also ran
in the active benchmark checkout and produced the expected **N=3, sum=60**,
with exit 0 and no SAS errors or warnings.

This test exposed a stale path in the final status JSON: the HTML had moved to
the results directory, but `htmlfilename` still referred to deleted scratch.
The runner now updates that field after saving the HTML. A targeted rerun passed
five assertions, including checking that the recorded HTML path exists
(`SAS_Output_cli_cleanup_fN6dhr/`). The live suite now includes this check for
both successful and debug-retained cases.

A final file-based run with the default macro-autoload setting also completed
with exit 0 and zero SAS errors/warnings. The self-contained script was detected
automatically, and its result again showed N=3 and sum=60. The benchmark
artifacts are in `SAS_Output_runner_default/`; the final status JSON points to
the retained HTML file. The ordinary command tested was:

```bash
OPEN_RESULT=0 perl run_sas_codes_or_files_in_SAS_ODA.pl \
  --file examples/sas_oda_runner_smoke.sas \
  --output-prefix SAS_Output_runner_default --run-timeout-seconds 180
```

Offline checks:

```bash
perl tests/test_cli_cleanup.pl
perl DiffGWASDeps/test_sas_log_and_labels.pl
perl DiffGWASDeps/test_sas_oda_mcp.pl
perl -c run_sas_codes_or_script_in_ODA.pl
perl -c server.pl
```

Real ODA checks using the existing configured account:

```bash
perl tests/run_live_cli_cleanup.pl --live
perl tests/run_live_mcp_oda.pl --live
```

The CLI suite generates a successful summary, a missing-table failure, a timeout,
and a debug-retained run. The MCP worker suite invokes the same callback used by
the explorer server and checks real ODA success/failure, logs, and HTML. It does
not exercise the server's complete MCP transport/handshake.

Offline validation passed 15 cleanup assertions, 16 log/label assertions, and
18 MCP worker assertions, including dispatch through the installed MCP SDK.
The explorer SDK passes both the tool and argument hash; its callback preserves
that two-parameter signature and forwards the argument hash to the shared worker.
The real CLI suite passed 14 assertions on
October 12, 2026; its artifacts are in `SAS_Output_cli_cleanup_ByNLhj/`.
The real MCP worker callback suite passed nine assertions. Successful and
intentionally failing ODA jobs retained their results under
`SAS_Output_ODA_4JEHIO/` and `SAS_Output_ODA_QcMdct/` respectively.

## Updated active copies

The port was applied to the explorer Git checkout, the working pipeline copy
beside MCP4SAS, and the active benchmark runner checkout shown in the IDE.
The relevant helper implementations are synchronized across those three copies.
Backup and packaging snapshots were not used as active targets. No remote push
was performed.
