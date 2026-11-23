# Portable macro queries and SAS runtime update

This incremental port uses MCP4SAS source revision
`58c2fd4378b27b65e9d8d2c684d3b0c208790a66`. It retains Explorer's
explicit-profile routing, target/configuration isolation, `close_session`
operation, result paths and Perl compatibility aliases. The original
research macros, GWAS inputs and plotting parameters are unchanged.

## Query available macros

From the inner `MultiGWAS-Explorer` directory, with SASPy configured:

```bash
perl run_sas_codes_or_files_in_ODA.pl \
  --code '%macroparas(macrorgx=GTEx,numlines2print=1000,print_outdsd=1);' \
  --persistent --session-id my_macros --open-html

python run_sas_codes_or_files_in_ODA.py \
  --code '%macroparas(macrorgx=GTEx,numlines2print=1000,print_outdsd=1);' \
  --persistent --session-id my_python_macros --open-html
```

The bundled portable `%macroparas` reads source files without executing the
macros or their demonstration programs. Its HTML uses one full-path heading per
macro, followed by left-aligned source text without printed line numbers.
Original dark-blue text and light-blue source cells are retained. Blank source
rows contain a TAB rather than a nonbreaking space, so copied code compiles.
ODS can trim leading indentation; this preview is not a byte-for-byte export.
The result dataset retains full paths programmatically. `print_outdsd=1` prints
the macro-name and parameter table, and a following bare `PROC PRINT;` refers
to that result table.

The query searches an explicit SAS-side `dir=`, then `MCP4SAS_MACRO_DIR`,
configured SASAUTOS locations and the SAS HOME/Macros directory. A client-side
`--macro-dir` controls dependency discovery and uploading; it is not a way to
make an inaccessible network path readable by SAS. No hard-coded institutional
network share is searched by default. An explicit user macro definition,
include or local macro file takes precedence over the bundled fallback.

`DiffGWASDeps/importallmacros_ue.sas` includes the same portable definition.
The importer records its selected `MacroDir` for subsequent queries and skips
legacy `macroparas.sas`/`macroparas_portable.sas` during bulk loading, so an old
definition cannot replace the portable one. Run
`python tools/sync_macroparas_bundle.py --check` to detect bundle drift.

**Bulk-import prerequisite:** the original Linux/ODA branch of
`%importallmacros_ue` still calls its library's enumeration/deletion helpers.
A complete macro library must contain those helpers. Bundling the portable
query does not turn an arbitrary folder of unrelated `.sas` files into a
complete import library. An initial test against a deliberately minimal
two-file fixture exposed these missing helpers; it failed visibly, rather
than being counted as successful. The subsequent import checks use the
existing complete local/ODA library. Plain portable macro queries need none
of those legacy helper macros.

## Cygwin and Python

Explorer now includes the [native Python CLI subset](../MCP4SAS_Use_Python/README.md).
It can hand off from Cygwin to an explicitly configured native Windows Python
3.10+ with SASPy. It converts client file/configuration/output paths, leaving
SAS code and remote SAS paths unchanged. Invalid explicit interpreters remain
actionable errors; they are not silently replaced. No virtual environment,
credentials or SAS JAR files are distributed in this port.

## Submission performance and state

Persistent Perl sessions keep a stable local HOME per profile, configuration
and session ID even when output folders change. Existing compatible helpers
are not restarted just to install performance changes: use an unused helper
port/new session when applying a new helper version, preserving needed WORK
tables before deliberately closing an older helper. Normal transport-failure
recovery retains its existing behavior; this is not a guarantee against WORK
loss after a real connection failure.

Internal HOME/template setup and metadata queries use text results. File
metadata is returned in one SAS response instead of four extra `symget`
round trips. Large source strings (default threshold 32 KiB) use one staged
file and a short `%INCLUDE` in the same SAS session; small programs remain
direct submissions. Failed staging does not submit or replay user code.
`submitted.sas` retains the processed source, and the SAS log records the
staged byte count and SHA-256. Windows HOME paths are accepted for include
debugging instead of being replaced by `/home`.

Local SAS writes new GTL templates into WORK and retains the configured read
paths, avoiding a writable SASUSER requirement. Template/HOME setup is cached.
These changes reduce avoidable operations; no controlled performance benchmark
or token/cost comparison is claimed by this port.

## Repeat the checks

```bash
bash install/check_sas_runtime.sh

# Explicitly opt in to real SAS; keep private credentials in SASPy config.
python tests/run_live_macro_query.py --live --runner perl --saspy-cfgname oda
python tests/run_live_macro_query.py --live --runner python --saspy-cfgname winlocal \
  --saspy-cfgfile /path/to/private/sascfg.py
```

The live test creates uniquely named sessions, verifies source-preview HTML,
recompiles a deliberately broken legacy query without clearing WORK, extracts
and executes the actual HTML source cells, checks numeric/string values, and
closes only its own resources. Add `--macro-library` with a complete existing
SAS-side library to test the bulk importer too. Raw logs and submitted programs
are retained in its `SAS_Output_RandID_macroquery_live_*` output folder.

Exact outcomes and scope are recorded in
[SAS_RUNTIME_MACRO_QUERY_VALIDATION.json](SAS_RUNTIME_MACRO_QUERY_VALIDATION.json).
The local runtime check passed 620 Perl assertions in 36 test plans, 50 Python
unit tests, the source/plan checks and both manifest checks. Eight handoff tests
include an actual Cygwin-to-native-Windows subprocess, without submitting SAS.
Live SAS checks passed for Perl and Python against both local Windows SAS and
ODA (five submissions per combination). The committed opt-in harness was also
executed successfully against local Windows SAS with Python and ODA with Perl.
Linux SAS routing was checked offline; a Linux SAS installation was not launched.

The [incremental manifest](SAS_RUNTIME_MACRO_QUERY_MANIFEST.json) records source
hashes and the pinned upstream revision. The earlier Nov21 validation remains a
historical record, not evidence that its full GWAS pipeline was rerun here.
