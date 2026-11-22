# SAS execution validation — 22 November 2026

## Changes and verified results

The plotting workflow and its MCP tool now accept `saspy_cfgname`,
`saspy_cfgfile` and `dir4local_sas_home`; the command-line equivalents are
documented in [SELECT_SAS_RUNTIME.md](SELECT_SAS_RUNTIME.md). ODA stays the
default. Explicit local profiles are not retried through ODA.

Persistent sessions are separated by profile and configuration-file path.
Local plotting steps share a dedicated staging home. Figure-path recovery
accepts Windows paths and spaces. Local session initialization prepends a
writable WORK template store instead of relying on a writable SASUSER store.
Original scientific macros and their plotting parameters were preserved.

The follow-up MCP4SAS discovery correction is pinned to commit
`79ada136979f0e91694d3191e02191a421235abe`. Windows CRLF endings are removed
from interpreter records and package-path probes. This prevents an installed
SASPy environment from being rejected and an unrelated Python being selected.
The source manifest records original and adapted hashes separately.

| Check | Actual outcome |
|---|---|
| Complete offline runtime suite | 600 passing TAP assertions in 34 Perl scripts; 15 archive-transfer and 11 browser Python tests passed |
| MCP SDK argument forwarding | 34 assertions passed, including profile/config/home paths with spaces |
| Real configured ODA PRINT | Exit 0; 19 SASHELP.CLASS rows verified; SAS platform LIN; full SAS log and HTML retained |
| Real configured ODA listing | Exit 0; no Python activation or `PIPELINE_PYTHON_BIN` override required |
| Earlier isolated Windows IOM integration | 18 assertions passed: transfers, PRINT/SGPLOT, persistent WORK/process identity, recursive macro staging, exact-file cleanup |
| Later Windows IOM/COM and native batch probes | Did not complete reliably; excluded from successful validation claims |
| Local Linux SAS | Profile, factory, request routing and failure contracts passed offline; no licensed Linux SAS executable was available |

GitHub runtime CI passed on Ubuntu 24.04, macOS 15 and Windows Server 2022
for the local-target implementation (runs `44990827663`, `44991896887`).
The discovery follow-up is checked separately by the same matrix.

## Unresolved local Windows checks

**LOCAL-WIN-STARTUP:** later IOM/COM requests and a native SAS batch probe
timed out or terminated during initialization/execution on this host. A prior
18-check Windows run succeeded; that is not evidence that every subsequent
configuration or graphical workflow succeeded. The native batch probe created
an HTML file but did not provide verified completion, so it is not counted as
a passing test. No assumption about licensing, authentication or memory is
presented as an established cause.

**LOCAL-WIN-GTF:** the original saved multi-GWAS GTF example initially could
not update `BedGraph` in the local template store. SAS then found an older
SASUSER template; these outputs are not accepted as scientific validation.
The WORK template-store safeguard was added. A rerender in the owned session
subsequently exceeded its 600-second limit. Its numerical and graphical
outputs therefore remain unvalidated on local Windows SAS. The original
macro definitions, failed log and rerender program were retained. The
previous successful ODA validation of the same frozen example remains a
separate result; a full GWAS pipeline was not rerun here.

**HOST-PATH:** unavailable mapped-drive/tool entries in the inherited Windows
PATH could stall even a trivial Cygwin Perl command. Tests used a scoped local
tool PATH; the user's system PATH was not changed. This is distinct from the
CRLF discovery defect. Failed startup attempts and partial/interrupted tests
are retained and excluded from the passing totals.

Use a short dedicated staging path on Windows. A deeply nested staging path
failed remote file-size verification; repeating the earlier integration test
with a shorter staging path passed all 18 assertions. Original input folders
must not be used as staging folders because cleanup removes staged copies.

## Evidence

Local raw evidence is retained under
`SAS_Output_RandID_local_selection_20261121/`: submitted programs, complete SAS
logs, execution/status JSON, HTML, round-trip files, offline TAP logs and
failed attempts. Private runtime configuration and credentials are not added
to Git. These results validate execution plumbing and limited examples;
they are not full scientific validation or a performance benchmark.
