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

## Real local-GTF comparison: association heights flattened

On October 12, 2026, the rs2070788/rs383510 SAS ODA local-GTF figure
showed twelve horizontal rows of points, while the gnuplot figure showed
varying association heights. The SAS runner correctly passed
`makeheatmapdotintooneline=0`, and the association mapper computed
`-log10(P)` separately from signed LD colors. However,
`Lattice_gscatter_over_bed_track.sas` subsequently changed the setting to
`1` when the number of groups exceeded its default threshold of ten.
That replaced the association heights with the constant `1.5`.

The automatic compression threshold now defaults to `0` (disabled).
The threshold condition also requires a positive threshold. This preserves
the requested association heights for plots with twelve or more tracks.
Explicit `makeheatmapdotintooneline=1` remains supported; callers can also
opt into automatic compression with a positive
`drawdotintooneline_if_totsc_gt` value. Both the main lattice macro and its
alternate local copy were corrected in all three active pipeline copies.

Reproduce the real validation in Cygwin from the explorer checkout:

```bash
perl auto_prepare_and_run_diff_gwas.pl \
  --spec configs/spec_pgc_scz_sex_common_automation.json \
  --target-snps rs2070788,rs383510 \
  --plots local_gtf \
  --ld-display-mode heatmap \
  --ld-population EUR \
  --ld-r2-threshold 0 \
  > sas_gtf_height_fix_validation.log 2>&1

perl DiffGWASDeps/test_selected_signed_ld_gtf.pl
```

The dispatcher regression passed all 39 assertions, covering locus grouping,
plot reuse, resume after failure, and status reporting. The dispatcher hashes
the plotting macro sources, so changing the lattice macro invalidates the
completed-plot cache without manually deleting the final HTML.

A Perl/Gunzip comparison of the real SAS and gnuplot input tables found 5,413
overlapping variants and 70,195 nonmissing P-value comparisons across fifteen
shared P columns, with no differences above 1e-10. This isolates the flat
association tracks to rendering rather than different association inputs.

The live rerun exited successfully. The final progress report has
`complete=1`, `remaining=0`, and an empty failure message. The regenerated
`PGC_SCZ_SAS_local_top_hits_with_gtf_common_targets_e25248ab57b2.png`
is 374,903 bytes; visual inspection confirmed varying association heights in
all twelve tracks and retained both SNP labels and the gene track. The final
HTML is the corresponding 509-byte figure wrapper. The pre-fix PNG was retained
locally as `*.before_height_fix.png`. The runner downloaded the explicit ODS
output and completed remote input/output cleanup. Validation artifacts remain
local and are not required to reproduce this correction.

## Shared SAS preflight corrections from real examples

On October 12, 2026, the shared source-analysis fixes validated through
MCP4SAS's 90 new example scenarios were applied to the explorer's canonical
ODA runner. `DiffGWASDeps/SAS_Source_Analysis.pm` now supplies the preflight
source views. The compatibility CLI continues to delegate to the canonical
runner, so both names use the corrected checks.

Preflight accepts macro definitions after other statements on the same
line, raw inline-data records containing apostrophes/comment delimiters,
quoted comment tokens, and statement comments after semicolons. Commented
PUT examples and macro-like text do not become compile blockers. Missing
inline-data terminators and genuine quote/comment/macro errors remain
detectable. Analysis preserves source locations and does not modify the
SAS program submitted to ODA.

The following Cygwin regressions passed 45 assertions (8, 24, and 13):

```bash
export PATH=/usr/bin:$PATH
perl tests/test_inline_macro_preflight.pl
perl tests/test_inline_data_preflight.pl
perl tests/test_comment_quote_preflight.pl
perl -c run_sas_codes_or_script_in_ODA.pl
```

The tracked smoke example combines the three quote/comment cases:

```bash
SAS_ODA_AUTOLOAD_MACROS=0 OPEN_RESULT=0 \
  perl run_sas_codes_or_files_in_SAS_ODA.pl \
  --file examples/sas_oda_preflight_smoke.sas \
  --output-prefix SAS_Output_preflight_smoke \
  --run-timeout-seconds 180
```

The real ODA run exited 0, emitted `PREFLIGHT_TOKEN_OK` and
`PREFLIGHT_VALUE_OK`, retained nonempty HTML, and reported zero SAS ERROR
and WARNING lines. Its local log is
`SAS_Output_preflight_smoke/output.html.info.txt`; the generated HTML is
`SAS_Output_preflight_smoke/sas_res_wsfdc.html`. Artifact names can vary on
rerun. An additional live check of the quoted-token case also passed.

The canonical runner, source-analysis module, regression tests, smoke
example, and this documentation were synchronized to the three active
pipeline copies. Generated PGC configuration overrides, LD caches, and
runtime results were retained locally and excluded from publication.

## SAS library validation, October 13, 2026

The shared runner was improved using 100 new library-derived live SAS ODA
cases in MCP4SAS (IDs 91–190, based on 80 distinct source files). The final
run passed all 100 validations: 98 successful calculations and two
intentional failures correctly reported. Source hashes, each calculation,
initial failures, corrections and Cygwin commands are recorded in
[the MCP4SAS report](https://github.com/chengzhongshan/MCP4SAS/blob/main/RealExamples4ImprovingMCP4SAS.md).

The companion runner now accepts percent-escaped STR delimiters and BQUOTE
apostrophes without false preflight failures. Macro dependency scanning
ignores NRSTR literals and ordinary single-quoted DATA strings while
retaining evaluated/deferred macro code. PROC FORMAT names ending in
digits receive advisory warnings; SAS still reports actual validity.
INPUT flowover and empty-dataset notes are exposed as nonfatal diagnostics,
including structured status, exit code, log path and diagnostics in the
MCP completion response.

Active include discovery preserves apostrophes inside double-quoted paths
and excludes masked literal include text. Filenames containing apostrophes
are uploaded via a safely named ZIP and extracted to their original names
with verified byte counts. Both persistent and one-shot upload routing
use this transfer path. The live include test submits a uniquely named
file, checks its calculation and HTML, and deletes only that remote fixture.

Run from this repository's script directory in Cygwin:

```bash
export PATH=/usr/bin:$PATH
prove tests/test_macro_quoting_preflight.pl tests/test_format_name_advisory.pl \
  tests/test_sas_log_diagnostics.pl tests/test_quoted_filename_upload.pl \
  DiffGWASDeps/test_sas_oda_mcp.pl \
  DiffGWASDeps/test_sas_include_preflight_lookup.pl \
  DiffGWASDeps/test_sas_log_and_labels.pl \
  DiffGWASDeps/test_sas_oda_upload_timeout.pl
SAS_ODA_AUTOLOAD_MACROS=0 OPEN_RESULT=0 \
  perl run_sas_codes_or_files_in_SAS_ODA.pl \
  --file examples/sas_library_macro_quote_smoke.sas \
  --output-prefix SAS_Output_library_macro_quotes \
  --run-timeout-seconds 180
perl tests/test_live_include_literal_paths.pl --live
```

The offline checks passed 105 assertions (97 across the first seven checked
files and eight upload-routing assertions). The live smoke run exited 0,
printed `LIBRARY_MACRO_QUOTE_OK` and retained nonempty HTML at
`SAS_Output_library_macro_quotes/sas_res_4foQ9.html`. Artifact names vary on
rerun. The live include check also passed all five assertions, retaining its
HTML and verifying cleanup of the unique apostrophe-containing remote fixture.
Runtime results, generated configuration overrides and LD caches remain
local. The shared fixes, tests, smoke example and documentation are
synchronized to all three active pipeline copies.

## Additional library batch: cases 191–290

A further 100 new examples in MCP4SAS use 49 additional library sources.
The case catalogue, source hashes, numerical/text oracles, discovery
failures and full Cygwin commands are documented in
[the MCP4SAS report, section 11](https://github.com/chengzhongshan/MCP4SAS/blob/main/RealExamples4ImprovingMCP4SAS.md).
Duplicate validation normalizes generated IDs/whitespace, and the new
batch rejects unexpected advisory diagnostics in successful jobs.
The complete final MCP4SAS rerun passed all 100 validations: 97 successful
jobs and three intentional SAS failures correctly reported.

Cases 191 and 192 exposed unnecessary targeted macro loaders for NRSTR
literals nested inside BQUOTE and inside double-quoted DATA strings.
The shared scanner now collects nested macro quote spans and recognizes
evaluated macro quoting inside double quotes. NRSTR literal calls remain
masked while genuine nested evaluated calls retain their dependencies.
Escaped macro quotes also preserve the surrounding SAS string boundary.

The expanded regression suite passed 113 assertions across eight Perl
test files. A live compatibility run through the public wrapper exited
0, emitted `LIBRARY_NESTED_QUOTE_OK`, injected no loaders for its literal
macro names and produced nonempty HTML at
`SAS_Output_library_nested_quotes/sas_res_8BVcF.html`.

```bash
export PATH=/usr/bin:$PATH
prove tests/test_macro_quoting_preflight.pl tests/test_format_name_advisory.pl \
  tests/test_sas_log_diagnostics.pl tests/test_quoted_filename_upload.pl \
  DiffGWASDeps/test_sas_oda_mcp.pl \
  DiffGWASDeps/test_sas_include_preflight_lookup.pl \
  DiffGWASDeps/test_sas_log_and_labels.pl \
  DiffGWASDeps/test_sas_oda_upload_timeout.pl
SAS_ODA_AUTOLOAD_MACROS=0 OPEN_RESULT=0 \
  perl run_sas_codes_or_files_in_SAS_ODA.pl \
  --file examples/sas_library_nested_quote_smoke.sas \
  --output-prefix SAS_Output_library_nested_quotes \
  --run-timeout-seconds 180
```

The scanner, regressions, smoke example and documentation are synchronized
to the three active pipeline copies. Runtime results and generated
configuration/cache files remain local.

## Additional library batch: cases 291–390

A further 100 examples use 42 additional SAS library files. The complete
catalogue, source hashes, discovery failures, corrections and Cygwin
commands are recorded in
[MCP4SAS report section 12](https://github.com/chengzhongshan/MCP4SAS/blob/main/RealExamples4ImprovingMCP4SAS.md).
The final full rerun passed all 100 validations: 98 successful jobs and two
intentional SAS failures correctly reported. All 42 source hashes and
normalized duplicate checks against the preceding 290 examples passed.

The shared runner now analyzes decoded CALL EXECUTE/DOSUBL code literals,
preserving source offsets and masking nested NRSTR/ordinary quoted DATA
text while retaining real evaluated dependencies. Explicit PROC PRINTTO
jobs reset stale SYSCC before execution, restore logging afterwards and
surface error condition codes of at least 8. This fixes a discovered
false-success response for a redirected missing-dataset error. The
redirected file itself is not automatically downloaded; the default log
contains the condition-code summary. Ordinary submissions are unchanged.

Long literal filerefs now receive a nonfatal preflight advisory, and actual
uninitialized-variable NOTES appear as structured advisory diagnostics.
Generated fixture issues involving FILEVAR bounds, implicit NOBS, floating
precision, IML orientation, SQL syntax and fileref length were corrected.
The live harness refreshes its next session after timeout/session loss and
retains the failed job without replaying it.

The companion offline checks passed 141 assertions across 11 Perl files.
The public-wrapper live checks passed all seven assertions, preserving
deferred calculations and HTML while correctly failing a real redirected
error with a SYSCC summary and no generated syntax error. Artifacts are
retained locally in `SAS_Output_deferred_printto_tdbPRQ`.

```bash
export PATH=/usr/bin:$PATH
prove tests/test_macro_quoting_preflight.pl tests/test_format_name_advisory.pl \
  tests/test_sas_log_diagnostics.pl tests/test_quoted_filename_upload.pl \
  DiffGWASDeps/test_sas_oda_mcp.pl \
  DiffGWASDeps/test_sas_include_preflight_lookup.pl \
  DiffGWASDeps/test_sas_log_and_labels.pl \
  DiffGWASDeps/test_sas_oda_upload_timeout.pl \
  tests/test_deferred_macro_literals.pl tests/test_printto_status_guard.pl \
  tests/test_fileref_name_advisory.pl
perl tests/test_live_deferred_and_printto.pl --live
```

The shared fixes, tests, standalone smoke example and this documentation
are synchronized to all three active pipeline copies. Generated
configuration overrides, caches and runtime artifacts remain local.

## Additional library batch: cases 391–490

Another 100 new examples draw on 41 additional library files. The full
case table, oracles, hashes, observed failures, repairs and Cygwin commands
are in [MCP4SAS report section 13](https://github.com/chengzhongshan/MCP4SAS/blob/main/RealExamples4ImprovingMCP4SAS.md).
The final full rerun passed all 100 validations: 98 successful jobs and two
intentional failures correctly reported. All 41 source hashes, normalized
duplicate checks against the preceding 390 cases and remote fixture cleanup
also passed.

Deferred `%include` directives with single-quoted paths were previously
missed because CALL EXECUTE/DOSUBL strings double the quote delimiters.
The shared extractor now decodes executable strings recursively and reads
SAS doubled delimiters within filenames. Dependency-path rewriting keeps
the correct escaping in direct and deferred source. Execution-context
detection follows parentheses, handles concatenation and stops treating
ordinary text after a completed DOSUBL call as executable code.

The new include fixtures use unique names. The harness retains local
programs and logs, deletes only remote fixtures submitted by the selected
cases and records cleanup success separately; cleanup failure returns a
nonzero harness exit. A name-literal fixture's VALIDVARNAME reset was moved
after its assertions and PRINT, preserving the caller's intended SAS policy.

Offline verification passed 152 assertions across 12 Perl test files.
The public-wrapper live test passed all seven assertions in
`SAS_Output_deferred_includes_g6G3xm`: both DOSUBL and CALL EXECUTE loaded a
single-quoted filename containing an apostrophe and space, produced the
correct value and nonempty HTML, and deleted only their unique remote fixture.

```bash
export PATH=/usr/bin:$PATH
prove tests/test_macro_quoting_preflight.pl tests/test_format_name_advisory.pl \
  tests/test_sas_log_diagnostics.pl tests/test_quoted_filename_upload.pl \
  DiffGWASDeps/test_sas_oda_mcp.pl \
  DiffGWASDeps/test_sas_include_preflight_lookup.pl \
  DiffGWASDeps/test_sas_log_and_labels.pl \
  DiffGWASDeps/test_sas_oda_upload_timeout.pl \
  tests/test_deferred_macro_literals.pl tests/test_printto_status_guard.pl \
  tests/test_fileref_name_advisory.pl tests/test_deferred_include_paths.pl
perl tests/test_live_deferred_includes.pl --live
```

The shared module fixes, unit/live tests and documentation are synchronized
to all three active pipeline copies. Generated configurations, caches and
runtime artifacts remain local.

## Validation with another 100 library-derived examples (491–590)

The MCP4SAS catalogue adds 100 distinct examples based on 40 further SAS
library files, with hashes and duplicate checks against the previous 490.
Discovery reproduced a dependency defect with multiple quoted local files
in one `%include` statement: only the first was uploaded, so SAS could not
open the second. The shared parser now collects each quoted operand,
including decoded DOSUBL and CALL EXECUTE literals. Comment separators,
escaped quotes and semicolons within paths are handled; options and later
statements are not mistaken for dependency paths. Ten focused regressions
cover these boundaries and inactive literal text.

The new Perl public-wrapper live test creates two uniquely named local SAS
files, runs their include list directly and through DOSUBL and CALL EXECUTE,
checks a computed total of 33 and nonempty HTML, and deletes only its two
remote fixtures. Its owned local output directory is ignored by Git.

```bash
export PATH=/usr/bin:$PATH
prove tests/test_macro_quoting_preflight.pl tests/test_format_name_advisory.pl \
  tests/test_sas_log_diagnostics.pl tests/test_quoted_filename_upload.pl \
  DiffGWASDeps/test_sas_oda_mcp.pl \
  DiffGWASDeps/test_sas_include_preflight_lookup.pl \
  DiffGWASDeps/test_sas_log_and_labels.pl \
  DiffGWASDeps/test_sas_oda_upload_timeout.pl \
  tests/test_deferred_macro_literals.pl tests/test_printto_status_guard.pl \
  tests/test_fileref_name_advisory.pl tests/test_deferred_include_paths.pl \
  tests/test_include_file_lists.pl
perl tests/test_live_include_file_lists.pl --live
```

The offline suite passed **162 assertions across 13 Perl files**. The shared
parser and Perl unit/live tests are synchronized to all three active
pipeline copies. The full case catalogue, Cygwin commands, original failures,
example corrections and final MCP validation evidence are recorded in
[RealExamples4ImprovingMCP4SAS.md, section 14](https://github.com/chengzhongshan/MCP4SAS/blob/main/RealExamples4ImprovingMCP4SAS.md).

The public-wrapper live test passed all ten assertions in
`SAS_Output_include_lists_Tgn3QJ`, including removal of its two remote
fixtures. This supplemental validation is separate from the 100 numbered
MCP examples.

The final complete MCP run passed **100/100 validations** in
`SAS_Output_real30_biOUtX`: 98 successful jobs and two expected failures.
Its four remote include fixtures were deleted successfully. The report
also retains the discovery run's 94/100 result and distinguishes runner
repairs from corrected example assertions and SAS syntax assumptions.
