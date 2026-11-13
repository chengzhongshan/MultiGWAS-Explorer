# SAS runtime performance: MCP4SAS port

The November 2026 MCP4SAS runtime fixes are integrated into MultiGWAS-Explorer's
SAS execution path. The upstream source is
[MCP4SAS 351e790](https://github.com/chengzhongshan/MCP4SAS/commit/351e7905908dc026f65ddc2da712467d6ccc7812);
the MultiGWAS-Explorer base is `d10985a`. The
[port manifest](SAS_RUNTIME_PORT_MANIFEST.json) records exact source revisions,
file hashes and which implementations were adapted.

## Fixes and integration boundaries

| Boundary | Problem | Ported correction |
|---|---|---|
| MCP argument dispatch | File-existence probes were applied to multiline or very large SAS source | `DiffGWASDeps/SAS_ODA_MCP.pm` now checks the shared filename-candidate predicate first. The MCP callback in server.pl continues to delegate to this module. |
| Include preflight | Unresolved SAS macro-variable paths triggered local searches and remote lookup | The canonical CLI defers targets containing `&` or `%` to SAS execution. Literal includes remain inspected; normal SAS failures remain visible. |
| Persistent transfers | The archive flag was ignored while files were transferred individually | `transfer_many` sends one bulk RPC to the existing SASPy session. The Python archive helper batches transfers and metadata, validates member sizes/CRCs, stages downloads and cleans temporary files. Explicit no-archive requests retain individual transfers. |
| Macro discovery | SAS-supplied LEFT/TRIM were treated as custom missing macros | Shared classification uses normal SAS autocall resolution for these helpers and retains custom dependencies. Standard SASAUTOS configuration is assumed. |
| Source inspection | The preceding source was rescanned for each quoted literal | One ordered context scan per analysis preserves CALL EXECUTE, DOSUBL, RESOLVE, nesting, escaped quotes and include order. |
| Cygwin/native Python | Nested transfer paths and one-shot status paths crossed the boundary without conversion | Convert local paths for the selected Python; preserve remote SAS paths and caller payloads. |

This project keeps its own architecture: `run_sas_codes_or_script_in_ODA.pl`
is the canonical CLI. The `...file...`, `...files...` and `...files_in_SAS...`
spellings are thin aliases, not independently maintained runner copies. Its
`SAS_ODA_*` temporary-directory configuration, retained HTML paths, COM thread
ownership, local SAS home support and GWAS-specific server tools are preserved.
The port does not replace server.pl with the smaller MCP4SAS server or change
GWAS calculations, LD policy, scientific macros or figure parameters.

The persistent protocol advances to `2026-11-12-persistent-archive`; the generated
and checked-in Python session helpers are tested for exact equality. An older
helper must be restarted through the normal lifecycle rather than silently
accepting a bulk request it does not support. Use distinct session ports when
running independent checkouts or clients with different runtime implementations.

## Local measurements in this checkout

The same frozen source bytes were profiled before and after this port, without
connecting to SAS. Each pair below sums one macro dependency scan and one
include-target scan. Before values are single observations; after values sum
the medians of three samples per function.

| Source | Bytes | Before, seconds | After, seconds |
|---|---:|---:|---:|
| Clustergram | 57,624 | 1.24 | 0.29 |
| Regional gene/GTF | 380,785 | 24.51 | 2.31 |
| Forest | 375,348 | 21.96 | 2.26 |

All nine function-result hashes, including the separate lexical-analysis
operation, match before and after. These are local source-inspection timings,
not complete GWAS runtimes or estimates of SAS CPU time. The sanitized
[profile record](SAS_RUNTIME_SOURCE_PROFILE.json) contains samples and hashes.

The upstream MCP4SAS evaluation observed a complete workflow reduction from
816.52 to 149.49 seconds; its retained direct SASPy baseline was 164.27 seconds.
Those were single runs at different times with different transfer policies.
They are **upstream evidence**, not new end-to-end measurements of this
MultiGWAS-Explorer port. See the upstream
[performance investigation](https://github.com/chengzhongshan/MCP4SAS/blob/351e7905908dc026f65ddc2da712467d6ccc7812/docs/PERFORMANCE_DEBUGGING.md)
for the full diagnostic history and the 72 retained output-equivalence checks.
Model tokens, billed cost and human effort were not measured.

## Diagnosing a slow request

1. Retain source/input hashes, source revisions, interpreter versions, session
   reuse, macro autoload, syntax-preflight and archive settings. Keep credentials
   outside shared records. Compare identical workloads and transfer policies.
2. Timestamp local preparation, session-server request receipt/return, actual
   SAS submission and result collection. A wrapper-submit status does not prove
   SAS has started. Profile a long gap before the server receives the request.
3. Verify the actual helper, protocol and `bulk_transfer` RPC. An enabled flag
   alone did not prove batching in the earlier implementation. Record fallback
   behavior rather than assuming every successful job used the intended path.
4. Run the local profiler from the inner pipeline directory:

   ```bash
   perl tools/profile_sas_source.pl --repeat 3 \
     --output source-profile.json /path/to/analysis.sas
   ```

   Hashing is outside the measured intervals. Compare result hashes as well as
   time. Preserve dependency semantics; do not disable safeguards for speed.
5. Reproduce the failing boundary with the offline checks, then verify a small
   live job before a long scientific run. Distinguish launch, network and login
   failures from SAS errors. On Windows use Cygwin-relative or converted script
   paths and the Python selected by the installation.
6. Compare numerical tables independently of rendered appearance after any
   live scientific rerun. Keep original warnings and failed attempts visible.

## Regression checks

After selecting the installed Perl/Python environment, run:

```bash
bash install/check_sas_runtime.sh
bash install/check_pipeline_install.sh
```

The focused check is offline and needs no SAS account or Java session. It tests
source boundaries, the actual persistent bulk dispatch, archive failures,
metadata batching, canonical aliases, executable-string semantics, Windows path
I/O, one-shot failure reporting, local homes and COM thread ownership. The full
installation check also runs it, extending the existing platform/container CI.
See [SAS_RUNTIME_PORT_VALIDATION.md](SAS_RUNTIME_PORT_VALIDATION.md) for actual
results and limits.

Saved programs still run without an AI agent:

```bash
perl run_sas_codes_or_script_in_ODA.pl --file analysis.sas \
  --persistent --session-id myproject
```

The configured SASPy/Java/SAS backend is still required. Executing saved code
this way makes no model request. Archive staging validates files before
promotion, but multiple final renames are not a transaction across destinations.
