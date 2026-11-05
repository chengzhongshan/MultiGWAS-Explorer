# SAS ODA macro and open-code improvements — 2026-11-02

The canonical `run_sas_codes_or_script_in_ODA.pl` and its compatibility entry
point `run_sas_codes_or_files_in_SAS_ODA.pl` now use the shared macro preflight
analyzer. The changes were ported from the verified MCP4SAS evaluation work while
preserving this project's temporary-directory handling, HTML paths and runner
configuration.

The analyzer catches open-code LOCAL, indexed/WHILE/UNTIL/standalone DO, orphan
END/ELSE/THEN, nested open IF, an IF action without DO, and missing-semicolon
recursion. Existing GOTO/RETURN/orphan MEND detection remains. It respects macro
definition scope, comments, quoted and macro-quoted literals, and inline records.
Supported top-level IF/THEN/ELSE DO groups remain accepted on SAS 9.4M5 or later.

The shared loader classifier also recognizes quoting/scope builtins and the
abbreviated `%INC` statement. Include discovery preserves source order and handles
literal deferred includes in single- and double-quoted execution APIs. RESOLVE
dependencies and directly unquoted NRSTR calls are recognized. The canonical
runner keeps all these rules in one shared source-analysis module.

Macro job reports contain specific categories and practical hints for open-code
errors, missing terminators, unresolved references and expression/argument issues.
The existing strict unresolved-macro-invocation policy remains: a SAS WARNING can
still deliberately fail the wrapper. The evaluation harness tests that policy
separately from an actual SAS ERROR. Timeout terminal messages invalidate an
expected-error check even if the wrapper exits 1.

## Validation

The bundled `tests/macro_open_code50` corpus contains 50 distinct synthetic
boundary examples with book-source/line hashes and numerical/error/warning
oracles. The source books are read-only inputs; no private analyses or book
datasets are copied into the programs.

The initial 50-case evidence was captured in MCP4SAS: 50 final outcomes verified
across four runs, with four original stale attempts excluded and preserved.
That reference audit is clearly named
`tests/macro_open_code50/reference_results_mcp4sas_2026-11-02.json`.
The original comparison against 4,000 previous generated programs was performed
there. This checkout validates uniqueness within the bundled 50 and also compares
older suites if they are present.

Project-specific checks on this port passed:

- **411 offline assertions across 22 test files**.
- **30 normal-helper preflight assertions**, covering both entry points and five
  confirmed structural failure types. Each returns exit 1 with complete failed
  metadata and an explicit message that SAS submission was skipped.
- **4/4 isolated SAS ODA outcomes** in this checkout: supported open conditional,
  invalid LOCAL, invalid IF action and strict unresolved-invocation policy. The
  preserved run directory is `SAS_Output_library500_live_0iEH3T`.

From this project's checkout in Cygwin:

```bash
export PATH=/usr/bin:$PATH
prove tests/test_macro_open_code50.pl tests/test_open_code_macro_diagnostics.pl \
  tests/test_include_source_order.pl tests/test_include_abbreviation.pl \
  tests/test_resolve_macro_dependencies.pl tests/test_unquote_macro_dependencies.pl \
  tests/test_double_quoted_deferred_includes.pl
prove -v tests/test_live_open_code_preflight.pl :: --live
perl tests/run_live_macro_open_code50.pl --live --raw \
  --only M4001,M4026,M4035,M4047
```

The raw probe deliberately bypasses static preflight and uses a fresh one-shot
connection per case, with a 120-second timeout. Normal user commands keep their
preflight protection. Raw case IDs above cover a supported conditional, invalid
LOCAL, an invalid IF action and strict unresolved-invocation policy respectively.
The reference audit does not claim that all 50 raw jobs were rerun in this port.

The corpus README includes full generation/collection commands. Current source,
submitted source, evaluation-contract hashes, actual logs, process exits and
terminal status must all agree before an outcome is credited. Explicit
`--skip-stale` excludes an entire changed batch and preserves its old metadata;
it cannot promote stale evidence to a passing result.

Runtime logs, credentials, automatic GWAS configurations and generated benchmark
outputs are excluded from this source update. No Python evaluation scripts are
required.

The later [2026-11-03 adaptation record](SAS_ODA_ADAPTATION_20261103.md) covers
the Windows CPAN downloader fix, cache-only GTEx alternative and Explorer's own
500-case round10 validation.
