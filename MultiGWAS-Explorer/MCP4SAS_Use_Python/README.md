# Python SAS runner for MultiGWAS-Explorer

This is the native CLI dependency subset of MCP4SAS revision
`58c2fd4378b27b65e9d8d2c684d3b0c208790a66`, under the retained
[MIT license](../DiffGWASDeps/MCP4SAS-MIT-License.txt). The full MCP4SAS Python
MCP server, extension, installation tools and developer corpora are not bundled.
Explorer's production plotting commands continue to use its Perl runner.

Run from the inner `MultiGWAS-Explorer` directory:

```bash
python run_sas_codes_or_files_in_ODA.py --saspy-cfgname oda \
  --code 'proc print data=sashelp.class;run;' --open-html
python run_sas_codes_or_files_in_ODA.py \
  --code '%macroparas(macrorgx=GTEx,numlines2print=1000,print_outdsd=1);' \
  --persistent --session-id my_macros --open-html
```

Both commands require an installed, configured SASPy environment. If launched
from Cygwin, the CLI can hand off to a native Windows Python 3.10+ with SASPy.
It discovers existing environments and accepts `PIPELINE_PYTHON_BIN` explicitly;
it does not silently replace an invalid explicit interpreter setting. Client
paths are converted during handoff; SAS source and remote SAS paths are retained.

Use `--saspy-cfgname winlocal --saspy-cfgfile /path/to/private/sascfg.py` for
configured Windows SAS, or `--saspy-cfgname linuxlocal` for configured Linux SAS.
A requested local profile does not fall back to ODA on connection failure.
Perl and Python have separate persistent-session implementations: use the same
runner, profile, configuration and session ID to reuse WORK tables.

`--help`, `--dry-run`, and log classification do not require a SAS connection.
The Python CLI implements the corresponding submission, transfer, macro,
profile, session and browser controls; it is not a drop-in implementation of
every internal Perl-only worker option. See `--help` for its supported options.

See [the port and validation notes](../docs/SAS_RUNTIME_MACRO_QUERY_20261122.md)
for macro-query formatting, bulk-import prerequisites and test coverage.
