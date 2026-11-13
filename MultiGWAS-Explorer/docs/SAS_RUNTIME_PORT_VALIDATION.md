# MCP4SAS runtime port: validation and remaining issues

Date: 12 November 2026. Ported from MCP4SAS
`351e7905908dc026f65ddc2da712467d6ccc7812` onto MultiGWAS-Explorer `d10985a`.
Testing used Cygwin Perl and native Windows Python 3.13.11 / SASPy 5.104.0.
Credentials, account paths and raw protocol captures remain local.

## Passing checks

- `bash install/check_sas_runtime.sh`: all 20 Perl files passed, with 278
  top-level assertions, plus all 10 Python archive tests. The check includes
  one-shot failure reporting, actual Cygwin/native file I/O, supplied/custom
  macro classification, literal contexts, include ordering, alias delegation,
  persistent bulk RPC, local homes and Windows COM thread ownership.
- Initial seven-file port-specific regression run: 82 assertions passed.
- Source profiler: all nine operation-result hashes agree before and after on
  the same frozen clustergram, regional gene/GTF and forest programs. Three
  repetitions after the port produced stable hashes. Local timings and source
  hashes are in [SAS_RUNTIME_SOURCE_PROFILE.json](SAS_RUNTIME_SOURCE_PROFILE.json).
- Generated Python helper matches the checked-in helper exactly. The archive
  implementation and source analyzer match the upstream bytes recorded in
  [SAS_RUNTIME_PORT_MANIFEST.json](SAS_RUNTIME_PORT_MANIFEST.json).
- Perl source, shell syntax and Python compilation checks passed for the
  changed runtime and helper files. Git whitespace checks passed.

## Actual SAS ODA smoke test

The canonical MultiGWAS-Explorer CLI created a fresh persistent ODA session with
syntax preflight enabled and global macro autoload disabled for a SAS-only
fixture. It uploaded two small TSVs through an actual bulk request, submitted a
program that required two rows and a sum of 33, generated PRINT output, and
downloaded both files through a second bulk request. Both returned SHA-256
hashes matched the inputs. The saved SAS log contains the completion marker
and zero SAS error lines. Archive verification messages and the two actual
bulk RPCs were checked separately. Both temporary remote files were deleted
and verified; only the owned test helper was shut down afterward.

The initial test launcher stopped on ordinary native STDERR diagnostics; after
correcting its capture, the SAS job completed. Its first completion-marker
assertion incorrectly inspected console output instead of the saved SAS log.
That assertion was corrected and the completed job was verified from retained
artifacts without rerunning it. These harness mistakes were not SAS failures.
The record does not report a complete-workflow time from the interrupted
harness. No full GWAS pipeline or new agent-token/cost benchmark was run.

## Broader offline suite: four baseline failures remain

The broader non-live Test::More selection ran 61 files and 1,440 assertions.
Fifty-seven files passed; four files failed. All four failures were reproduced
from a separately exported **unmodified `d10985a`** snapshot under the same
local environment (34 assertions; the same nine assertion failures and two
incomplete plans). They therefore remain baseline issues outside this runtime
port, rather than being reported as a clean full-suite pass.

| Existing test | Observed failure and follow-up |
|---|---|
| `install/test_github_sync.pl` | Native Git could not create the Cygwin `/tmp/...` bare repository. Review the test's local path boundary for native Git before counting this scenario as validated on Windows. |
| `install/test_gnuplot_chromosome_bounds.pl` | Four assertions expect an end-clipped start of 236,321,905; current fixed-span behavior returns 228,956,422 for hg38 and 229,250,621 for hg19. Reconcile the intended span contract with the existing test; this port changes neither implementation nor expectations. |
| `install/test_merged_meta_import.pl` | The Cygwin fixture fails its direct-import assertion after `cygpath failed`, then cannot read its expected spec. Reproduce and correct that fixture/import boundary separately. |
| `install/test_sas_manhattan_compact_input.pl` | Four assertions expect a chromosome-X fixture row that is absent from current output. Determine whether current filtering or the expected contract should change before modifying scientific behavior. |

An extra invocation with a broad `tests/test_*.pl` glob also selected live-only
programs without `--live`; they refused to run. That was an invalid test
selection, not additional executed SAS evidence, and is not counted above.

The focused runtime check is wired into the existing Installation preflight
and dependency smoke checks. Platform/container outcomes must be read from
the actual GitHub Actions run for the published commit; this report does not
infer those results from the workflow definition.

## Scope

The original working checkout's edited generated configurations and retained
results were preserved. No scientific SAS macro, genomic dataset, LD policy,
plot parameter or figure was changed. Existing warnings were not hidden.
The upstream 72 numerical/pixel equivalence checks are retained MCP4SAS evidence;
they are not described as a fresh 72-check MultiGWAS execution. This port has
an actual transfer/SAS smoke test and source-analysis evidence, not full
scientific or publication validation.
