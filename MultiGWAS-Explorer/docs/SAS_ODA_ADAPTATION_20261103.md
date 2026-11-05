# SAS ODA and installation improvements adapted on 2026-11-03

The Windows dependency downloader fix, cache-only GTEx plotting alternative and
500-case round10 evaluation corpus have been adapted from MCP4SAS to Explorer.
The tests below executed Explorer's own helpers and modules. The 500 cases are
the same round10 cases already generated for MCP4SAS, not an additional corpus.

## Installation

`install/common.sh` now routes all four CPAN installation paths through
`pipeline_cpanm`: ordinary modules, Cygwin compatibility repair, Devel::CheckLib
and PDL. Cygwin explicitly selects curl with `--no-lwp --no-wget --curl --verbose`
and the configured HTTPS mirror. This addresses the downloader selection that
previously failed on Windows even after switching to an HTTPS mirror.

The wrapper invokes the repo-local cpanminus script with the active Perl. It
retains Explorer's Perl ABI checks, local library directories, compatibility
repairs, PDL build timeouts and MAKEFLAGS. Its options array is always nonempty
so Bash 3's nounset handling does not fail on macOS.

`install/test_cpan_downloader.sh` passed 29 assertions covering Cygwin arguments,
Linux and Darwin argument selection, paths containing spaces, a custom mirror,
PDL timeout arguments and failure propagation. These mock argument checks do not
claim a real Linux/macOS installation. The existing first-install environment,
ABI rotation and Java path checks also passed under local Cygwin. The downloader
test now runs in the existing Unix environment check and before the Windows
portable installer in GitHub Actions. No new Explorer Actions run was triggered
by these local checks.

## GTEx plotting

The new `examples/sas_macros/Boxplots4GenesInGTExV8_Cached.sas` requires a prepared
long expression table and has no download path. It preserves the cache and
expression scale, validates requested genes and schema, and refuses rebuild
requests. The original local and online macros have not been replaced.

`tools/subset_gtex_expression.pl` streams an existing plain/gzip matrix to a small
gene subset before SAS import. Twelve assertions checked exact values, gzip,
missing genes, refusal to overwrite, malformed rows and truncated input.
Three fresh SAS ODA checks passed through Explorer's helper: correct two-gene
output without changing the source, missing-cache refusal and rebuild refusal.
See [GTEx recovery instructions](GTEX_LARGE_DOWNLOAD_RECOVERY.md) for usage and
the limits of the original resource-failure diagnosis.

## Evaluation evidence

All **500/500 cases** passed in Explorer run
`SAS_Output_library500_live_yqnfuI`: 450 expected successes and 50 expected SAS
errors across 50 parameterized families. The error cases require the expected
diagnostic and ordinary exit 1; a timeout or remote session termination cannot
pass. Distribution tolerances were `1e-10*max(1,abs(expected))` from the start.

Independent collection verified all 500 current source bodies, submitted batches,
evaluation contracts, actual SAS logs, exits and terminal status files, with no
unsuccessful attempts. The audit is
`tests/library500_round10/verified_results_2026-11-03.json`. The separate
`reference_results_mcp4sas_2026-11-03.json` preserves original MCP4SAS evidence;
it is not counted as an Explorer run. See the
[corpus README](../tests/library500_round10/README.md) for source provenance and
the distinction between source-derived examples and distribution API extensions.

The complete offline Perl suite passed **928 top-level checks in 24 files**.
It includes source/manifest hashes and executable novelty checks against earlier
suites present in this checkout. MCP4SAS's larger historical novelty comparison
is recorded as reference history, rather than claimed as rerun here.

## Cygwin commands

Run from the `MultiGWAS-Explorer` project directory:

```bash
export PATH=/usr/bin:$PATH
bash install/test_install_environment.sh

# Offline suite: live tests require explicit --live and are run separately.
offline=()
for test_file in tests/test_*.pl; do
  case "$test_file" in tests/test_live_*) continue ;; esac
  offline+=("$test_file")
done
prove "${offline[@]}"

prove -v tests/test_live_gtex_cache_only.pl :: --live
perl tests/run_live_library500.pl --live --round 10

# Use the run directory printed by the preceding command for a fresh rerun.
perl tools/collect_library_evaluation_results.pl --round 10 \
  --run SAS_Output_library500_live_yqnfuI \
  --output tests/library500_round10/verified_results_2026-11-03.json
```

Raw test artifacts are local ignored outputs under
`SAS_Output_adaptation_20261103`, `SAS_Output_library500_live_yqnfuI` and
`SAS_Output_gtex_cache_JW2JzC`. Only synthetic data and owned WORK resources were
used by the live checks. The original book library is read-only input. New tools
and evaluation scripts use Perl and Bash.

Persistent-session startup also uses the selected interpreter's path convention.
A native Windows Python session created through the corrected launcher was
reused by this runner with the same SAS PID 37624 and the same five-row WORK
dataset. HTML: `SAS_Output_RandID_persistent_path_shared_20261104/sas_res_FIVdu.html`.
The canonical runner's help now includes the persistent `mc1` create/read example.

## Installation workflow repair: Bash 3.2 compatibility (2026-11-04)

The previous GitHub Installation run 41922257401 failed in all six installation
jobs. Every job reached the existing LD heatmap contract check and stopped with
`run_sas_oda_local_top_hits_manhattan_download_png.sh uses Bash 4 case conversion
and will fail in macOS Bash 3.2`. Dependency installation and SASPy imports had
already passed; this failure was separate from local Windows SAS authentication.

The new automatic PLINK2 lead-selection branches used Bash 4 uppercase expansion.
Both branches now use the existing portable `printf | tr` pattern, preserving
case-insensitive LD/PLINK2_1KG matching. The existing contract check reproduced
the failure before the change and passed afterwards. Shell syntax, installer
environment checks (including all 29 CPAN downloader assertions), and the 12
requested-hit gene checks passed locally.

The workflow now runs the existing plotting contract and CPAN argument checks
in a short preflight job before any of the six installation jobs starts. Full
installation and synthetic plotting checks remain required on all platforms.
The workflow tests installed dependencies and synthetic plots without SAS ODA
credentials or a licensed local SAS instance. The previously published Windows
persistent-session runtime fix remains in commit adb2b7e.
