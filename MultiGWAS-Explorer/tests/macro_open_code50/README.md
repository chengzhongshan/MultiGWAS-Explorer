# Open-code macro evaluation: M4001–M4050

These are **50 distinct boundary examples**, not ten parameter variations of five
templates. Their normalized executable bodies differ from all 4,000 earlier
programs. They use synthetic data and actual source/line evidence from the macro
books in the requested SAS book sample directory. Original book files are unchanged.

Forty-nine source anchors come from *Carpenter's Complete Guide to the SAS Macro
Language, Second Edition* and one from *SAS Macro Programming Made Easy, Second
Edition*. Three `derived_extension` entries explicitly anchor modern boundary
experiments to a real macro definition instead of claiming that the exact newer
statement occurs in the older book. The manifest records titles, filenames,
lines, source/line hashes, executable hashes, numerical/error/warning oracles
and evaluation-contract hashes in the run evidence.

| IDs | Examples |
| --- | --- |
| M4001–M4003 | Supported open-code IF/DO, ELSE and sequential conditionals |
| M4004–M4006 | Indexed, WHILE and UNTIL loops correctly wrapped in macros |
| M4007–M4011 | LOCAL scope, RETURN, GOTO, inert definitions and nested IF inside macros |
| M4012–M4018 | STR/NRSTR, quoted controls, comments and DATALINES4 records |
| M4019–M4025 | Option restoration, EVAL/SYSEVALF, SUPERQ, delimiters, indirect references and SYMPUTX timing |
| M4026–M4030 | Invalid open LOCAL and indexed/WHILE/UNTIL/standalone DO |
| M4031–M4033 | Orphan END, ELSE and THEN |
| M4034–M4036 | Nested open IF, missing DO action and missing-semicolon recursion |
| M4037–M4039 | Invalid open GOTO/RETURN and orphan MEND |
| M4040–M4045 | EVAL division/decimal errors, unknown SYSFUNC, duplicate parameters, argument ordering and empty IF |
| M4046–M4050 | Unresolved variables/invocations, missing symbol deletions and recursive defaults |

The final outcomes are **29 successful MCP jobs, 20 intended SAS-error jobs and
one intended MCP policy failure**. Five cases have required actual warnings.
M4047 is special: SAS produces a warning and the expected numeric marker, while
MCP4SAS deliberately fails unresolved macro invocations. This existing protective
policy is retained and verified separately from SAS ERROR severity.

The shared source analyzer now catches 14 structural problems before submission,
compared with three in the preserved baseline. Thirty positive/warning examples
remain unblocked. Six expression/parameter errors still require actual SAS
evaluation; the static scanner is not a complete macro interpreter. Supported
open-code conditional groups assume SAS 9.4M5 or later.

MCP job summaries now include categorized macro diagnostics and hints. Source
echoes and submitted strings cannot fabricate these diagnostics. Warning-only
oracles require actual WARNING lines plus a runtime marker and correct status;
the strict unresolved-invocation case additionally requires the failed wrapper
exit. Timeout messages invalidate expected-error checks even when the wrapper
uses exit 1.

From the project checkout in Cygwin:

```bash
export PATH=/usr/bin:$PATH
perl tools/build_macro_open_code50.pl \
  '/cygdrive/h/F_Queens/360yunpan/SASCodesLibrary/SAS PDF/sas online samples'
prove tests/test_macro_open_code50.pl tests/test_open_code_macro_diagnostics.pl
perl tests/test_live_open_code_preflight.pl --live
perl tests/run_live_macro_open_code50.pl --live --raw
```

The `--raw` runner deliberately bypasses static preflight to measure SAS runtime
behavior in fresh isolated connections (120-second timeout per case). Normal
helper commands retain preflight. Use this option only for diagnostic evaluation.

`reference_results_mcp4sas_2026-11-02.json` is the preserved 50/50 evidence from
the MCP4SAS implementation of these shared components. It records four runs and
four excluded stale attempts; it is not a claim that all 50 raw jobs were rerun
in this checkout. The four-thousand-case novelty comparison was performed in
MCP4SAS; this checkout checks uniqueness within the bundled 50 cases.

The project-specific validation and deployment notes are in
`../../docs/SAS_MACRO_PREFLIGHT_IMPROVEMENTS.md`. To collect a new raw run here:

```bash
perl tools/collect_library_evaluation_results.pl --round 9 \
  --run SAS_Output_library500_live_REPLACE_WITH_YOUR_RUN \
  --output tests/macro_open_code50/verified_results_DATE.json
```

Explicit `--skip-stale` excludes entire changed batches and retains their original
metadata in `stale_attempts`; excluded attempts never count as passes. Without
that option, changed source/contract evidence fails collection. Raw logs and
credentials remain outside commits. No Python evaluation scripts are required.