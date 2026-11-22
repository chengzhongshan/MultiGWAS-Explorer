# SAS ODA runtime update from MCP4SAS

The current MCP4SAS Perl runtime has been ported to MultiGWAS-Explorer. The
source revision is `cd0caa928bd37da3a7720a136d6949c261cbb545`; the target started
at `3a7a5fa1235dfbf80b11e2fbc01550e85291dc3d`. The
[manifest](SAS_RUNTIME_PORT_MANIFEST.json) records file hashes and adaptations.
The upstream [MIT notice](../DiffGWASDeps/MCP4SAS-MIT-License.txt) is retained.

## What changed

- Local macro directories are normalized before use. An inaccessible directory
  falls back to the working directory with a visible warning. Recursive local
  `%include` files are staged together with their dependencies, including paths
  containing spaces. Quoted upload filenames retain their intended spelling.
- Persistent uploads/downloads retain archive batching. New-file metadata is
  queried in batches; existing-file timestamp checks remain. Failed upload
  verification can retry once from a private staged copy. Integrity failures
  remain failures rather than being treated as completed transfers.
- Large SAS source preflight no longer copies the entire remaining source at
  every character when finding macro calls. Nested calls, quoted/comment text,
  missing parentheses and diagnostic line numbers retain their meaning.
  SAS-supplied constructs, including `%sysexec`, do not trigger custom macro
  lookup.
- Java startup messages are preserved and classified. Missing encryption
  classes, missing Java, SSL, DNS, heap and Windows IOM authentication failures
  receive different guidance. Java errors are not automatically relabeled as
  wrong passwords or assumed SAS WORK exhaustion. Windows network-path errors
  remain non-retryable infrastructure failures.
- `--open-html` explicitly requests the system default browser;
  `--no-open-html` suppresses it. macOS uses `open`, Windows uses its registered
  file association, and Linux uses available desktop launchers. Firefox is not
  required. Failed launchers are not reported as successful merely because a
  process was forked. Opening a window still depends on the host's desktop
  session and browser configuration.
- The Explorer SAS MCP tool waits up to 20 seconds for a verified terminal
  result by default. Longer jobs return a structured `RUNNING` result and PID;
  subsequent calls poll that accepted job without resubmitting SAS. Set
  `wait_seconds=0` for immediate background delivery, or a value up to 30.
  `max_log_chars` defaults to 16000; zero returns the full log. The complete
  saved log remains available regardless of response excerpting. Macro
  directory, local SAS home and browser options reach the actual CLI.

The Perl entry points remain `run_sas_codes_or_script_in_ODA.pl` (canonical),
`run_sas_codes_or_files_in_ODA.pl`, `run_sas_codes_or_file_in_ODA.pl` and
`run_sas_codes_or_files_in_SAS_ODA.pl`. Thin aliases are preserved, so fixes
apply to existing Explorer plotting wrappers without duplicating a large CLI.
The shared Python SASPy bridge and transfer/browser helpers were updated as
required by the Perl runner; the independent MCP4SAS Python runtime and VS Code
extension were not copied into Explorer.

## Running saved SAS

From the inner `MultiGWAS-Explorer` directory, after configuring an ODA profile:

```bash
perl run_sas_codes_or_files_in_ODA.pl --saspy-cfgname oda \
  --code 'proc print data=sashelp.class(obs=5); run;' \
  --output-prefix SAS_Output_RandID_oda_example --open-html
```

If the profile is in a separate private file, add
`--saspy-cfgfile /path/to/private/sascfg_personal.py`. The selected Python/SASPy
environment must be able to read it. Do not put passwords in committed scripts
or command examples. A named SASPy profile selects its target explicitly; do
not rely on a failed local profile falling back to ODA.

Use `--persistent --session-id analysis_name` consistently to retain WORK
datasets between calls. Existing helper processes do not automatically reload
new source code. For a fresh test, select an unused `SAS_ODA_SESSION_PORT` and
a new session name; preserve or finish other work before intentionally closing
its sessions. Do not use a global process-kill option to refresh unrelated jobs.

## Java and platform setup

For ODA, use the Java executable and IOM classpath configured for the selected
SASPy profile. If encryption files are missing, obtain them through the
[official SAS download](https://support.sas.com/downloads/package.htm?pid=2494)
using a SAS Profile, extract the package and locate `sas.rutil.jar`,
`sas.rutil.nls.jar` and `sastpj.rutil.jar`. Add their absolute paths alongside the
existing bridge/IOM dependencies in that private profile's `classpath`.
Paths must be readable by the Python/Java process; native Windows Python needs
Windows paths, which Cygwin can produce with `cygpath -am`.

An encryption class error is different from an SSPI authentication failure.
For local Windows IOM, confirm the intended local profile, its localhost target
and the Windows identity running SAS. Do not enter ODA credentials into a
local-SAS profile as a workaround. A configured Windows COM profile must keep
COM operations on their owning thread. Linux local SAS requires installed,
licensed SAS and its own STDIO profile; a unit contract is not a live Linux SAS
execution. This update retains these target-selection boundaries.

Check a configured connection with:

```bash
perl run_sas_codes_or_files_in_ODA.pl --saspy-cfgname oda --check-saspy-connection-only
```

If a company firewall blocks outbound TCP 8591 to the configured SAS ODA host,
ask IT/IS to allow that connection. The runner cannot bypass a blocked port.

## Reproducible checks

The offline suite runs without ODA credentials, Java sessions, an AI agent or
scientific downloads:

```bash
bash install/check_sas_runtime.sh
```

The script verifies that Perl actually executes a probe before counting tests.
On Windows, use Bash and Perl from the same Cygwin installation. A stale PATH
pointing to another installation can otherwise invalidate a test run. The new
GitHub `SAS ODA runtime` workflow uses Cygwin Perl with native Python on Windows,
and the platform's Perl/Python on Ubuntu and macOS.

An opt-in live regression exercises two-file transfer, persistent WORK reuse
across aliases, recursive macro staging, retained PRINT HTML, exact remote
cleanup and actual persistent bulk requests:

```bash
perl tests/test_live_runtime_updates.pl --live \
  --saspy-cfgfile /path/to/private/sascfg_personal.py
```

It creates uniquely named files and its own SAS session/server, retains logs
under `SAS_Output_RandID_runtime_live_*`, and closes only its own resources.
The [validation record](SAS_RUNTIME_VALIDATION_20261121.json) contains the
sanitized outcomes and source/output hashes. Raw execution artifacts and private
configuration remain local.

## Validation scope

The local offline check passed **561 assertions in 32 Perl test files**, plus
**15 archive tests and 11 browser tests** in Python. The existing actual MCP
SDK dispatch fixture passed 31 assertions. Server syntax passed after its
existing optional prerequisites were supplied in an isolated test directory;
no packages were installed globally. Java error fixtures run copied unchanged
CLI/modules from owned temporary directories, so they also work in read-only
images and cannot accidentally import an installed SASPy in place of the stub.

The live synthetic regression passed **17/17 checks** against SAS ODA. The
separate frozen multi-GWAS plotting payload passed **11,295 checks**, using
1,880 variants across three panels. Its four numerical CSV tables matched the
retained independent MCP4SAS direct-Python outputs at tolerance `1e-10`; its
decoded PNG pixels were identical. The plotting payload includes its own macro
definitions: this validates the Explorer runner, rather than replacing or
rerunning Explorer's entire production scientific pipeline. All five original
SAS plotting warnings remain recorded.

One full synthetic test took 130.084 seconds; the plotting/transfer/checking
example took 147.760 seconds. These are single execution observations, not a
paired speed comparison. Local Linux SAS and local Windows SAS were not newly
executed; their configuration, HOME and COM-threading contracts were checked
offline.

The initial runtime branch passed the targeted GitHub workflow on
[all three platforms](https://github.com/chengzhongshan/MultiGWAS-Explorer/actions/runs/44987080093).
Read subsequent workflow runs for the exact final publishing commit; this
record does not infer new platform results from the workflow definition.

Cleanup used only the test-owned session/server and exact uploaded synthetic
filenames. The intentional `endsas;` command returns exit 74; upstream recovery
can prewarm an owned replacement session before the helper is explicitly shut
down. That cleanup behavior is retained as a limitation, not counted as a
failed plotting execution or hidden as a warning-free close. Shared prepared
PoC fixtures were retained instead of deleting potentially pre-existing inputs.
Raw logs preserve the original execution paths, including the isolated
worktree used during testing.

This port changes execution, staging, diagnostics and result delivery. It does
not change original scientific SAS macros, GWAS inputs, allele harmonization,
LD selection or plotting parameters. Numerical and appearance comparisons are
reported separately. Successful rendering does not independently validate
biology, assembly labels or a complete GWAS pipeline. No AI token use or cost
savings were measured.
