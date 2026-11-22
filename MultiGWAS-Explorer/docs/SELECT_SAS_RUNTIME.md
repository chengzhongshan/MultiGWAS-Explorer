# Choose SAS ODA or licensed local SAS

SAS OnDemand for Academics remains the default. To execute a plotting workflow
on an installed copy of SAS, explicitly select its SASPy profile. The same
selection reaches genome-wide Manhattan, regional Manhattan, GTF and forest
plotting steps, including their upload, download and cleanup operations.

## Windows

Run these commands from the installed MultiGWAS-Explorer directory in Cygwin.
Use native Windows Python with SASPy and a working local SAS configuration.
`winlocal` is the Windows IOM profile; see
[Windows IOM/SSPI setup](WINDOWS_LOCAL_IOM_SSPI.md). `winlocal_com` is an
alternative for an installation with a registered SAS COM provider and pywin32.
Neither profile supplies a SAS license or installs SAS Foundation.

```bash
export PIPELINE_PYTHON_BIN='C:/path/to/python.exe'

# First verify the connection without running an analysis.
perl run_sas_codes_or_files_in_ODA.pl \
  --saspy-cfgname winlocal \
  --saspy-cfgfile 'C:/path/to/sascfg_winlocal_iom.py' \
  --check-saspy-connection-only

# Run your existing analysis specification using local Windows SAS.
perl auto_prepare_and_run_diff_gwas.pl --spec configs/my_analysis.json \
  --saspy-cfgname winlocal \
  --saspy-cfgfile 'C:/path/to/sascfg_winlocal_iom.py' \
  --plots manhattan,local_manhattan,local_gtf,forest \
  --no-gnuplot-fallback-on-sas-failure --no-gnuplot-fallback-on-sas-space
```

The last two flags require SAS rendering to succeed instead of allowing the
existing gnuplot recovery path. Gnuplot recovery never changes the requested
SAS profile to ODA. Existing reference-data, scientific selection and plot
options in your specification retain their meanings.

## Linux

Use Linux Python with SASPy and a licensed Linux SAS executable. The supplied
STDIO profile avoids Windows Java and COM configuration.

```bash
export PIPELINE_PYTHON_BIN="$(command -v python3)"
export SASPY_LOCAL_SAS_PATH='/path/to/SASFoundation/9.4/sas'

perl run_sas_codes_or_files_in_ODA.pl \
  --saspy-cfgname linuxlocal --saspy-cfgfile "$PWD/install/sascfg_linux_local.py" \
  --check-saspy-connection-only

perl auto_prepare_and_run_diff_gwas.pl --spec configs/my_analysis.json \
  --saspy-cfgname linuxlocal --saspy-cfgfile "$PWD/install/sascfg_linux_local.py" \
  --plots local_manhattan,local_gtf,forest \
  --no-gnuplot-fallback-on-sas-failure --no-gnuplot-fallback-on-sas-space
```

Local Linux execution requires SAS installed on that Linux host. This is not
native SAS execution on macOS; a Mac can use ODA or a separately configured
remote SAS connection. Live Linux SAS was unavailable for this change;
Linux profile/routing/error behavior is covered by offline tests.

## Save the selection in your specification

The workflow accepts `saspy_cfgname`, `saspy_cfgfile` and
`dir4local_sas_home` in the analysis JSON. Command-line options take precedence,
then specification values, then `SASPY_CFGNAME` (or `SASPY_CONFIG_NAME`),
`SASPY_CFGFILE` and `SAS_LOCAL_HOME_DIR`. Without a configured profile, it uses
`oda`. A configuration file must actually contain the selected profile; a
failed local connection is reported without retrying ODA.

```json
{
  "saspy_cfgname": "winlocal",
  "saspy_cfgfile": "C:/path/to/sascfg_winlocal_iom.py",
  "dir4local_sas_home": "E:/SAS staging/my_analysis"
}
```

Add these keys to your complete analysis specification. This fragment alone
does not describe GWAS inputs or study groups. Keep credentials out of shared
specifications and generated configuration files.

## How files and sessions are kept consistent

For known local profiles, the workflow defaults to a dedicated child directory
`SAS_Local_Home_<artifact_stem>` under `output_dir`. Every plotting helper uses
this same home even when its output prefix changes. Alternatively set
`--local-sas-home DIR` (alias `--dir4local_sas_home DIR`). Choose a writable,
dedicated staging directory, not your original data directory: requested
cleanup removes staged copies of inputs. For Windows SAS, use a short path
such as `C:/SAS_staging/my_analysis`; deeply nested paths can exceed SAS/Windows
file limits. Custom local profile names require an
explicit home. Remote IOM/SSH sessions reject local staging homes.

Generated runner JSON pins the profile, configuration path and home. It can be
reused through `RUNNER_CONFIG_JSON` by the existing shell plotting helpers.
Local templates use the configured SAS HOME rather than an ODA user's shared
library. Figure-path recovery accepts Windows drive paths, Linux paths and
paths containing spaces.

Persistent connections are isolated by both profile and configuration-file
path, even when a caller reuses the same `--session-id`. Switching between
ODA and Windows/Linux SAS therefore cannot reuse the other target's WORK
tables or connection. Changing the contents of a configuration file at the
same path still requires a fresh session ID. The new helper protocol restarts
older helper versions; save work before upgrading and recreate needed WORK
tables afterward.

To return to ODA, select `--saspy-cfgname oda` and your ODA configuration file.
An inherited local staging home is cleared. Explicitly supplying a local home
and `oda` together is an error.

`--local-sas-only` keeps its previous meaning: emit desktop SAS debug programs
without executing them. It does **not** select a local execution backend.

## Regression tests

`tests/test_explicit_sas_target.pl` exercises the workflow's actual config
emitter, shell exports, framed persistent handler, target separation, failed
local connections without ODA fallback, and output paths with spaces.
`tests/test_local_sas_home.pl` covers Linux STDIO factories and HOME setup.

The opt-in live test can be repeated on a configured Windows installation:

```bash
perl tests/test_live_runtime_updates.pl --live --saspy-cfgname winlocal \
  --saspy-cfgfile 'C:/path/to/sascfg_winlocal_iom.py' --with-graph \
  --outdir SAS_Output_RandID_local_check
```

It uses an isolated helper port and retains submitted programs, complete SAS
logs, PRINT/SGPLOT HTML, round-trip files and execution metadata. It checks
two rows summing to 33, byte-identical transfers, persistent WORK reuse,
recursive macro staging and exact-file cleanup. This is execution plumbing
validation, not a rerun of a scientific GWAS pipeline.
