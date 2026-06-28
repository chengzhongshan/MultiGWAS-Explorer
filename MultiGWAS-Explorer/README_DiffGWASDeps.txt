DiffGWASDeps Layout
===================

This project keeps the user-facing differential-GWAS entrypoint at the repo root:

  perl auto_prepare_and_run_diff_gwas.pl ...

Its required Perl, shell, and SAS dependency files now live in:

  DiffGWASDeps/

Keep these together when moving or packaging the workflow:

1. auto_prepare_and_run_diff_gwas.pl
2. DiffGWASDeps/
3. configs/
4. server.pl                  (if MCP access is needed)
5. Vendored SAS ODA helper stack inside `DiffGWASDeps/`:
   - `run_sas_codes_or_script_in_ODA.pl`
   - `SAS_ODA_Runner.pm`
   - `sas_oda_session_server.py`
   - `GetPIDs.pl` (used by the MCP server helper path)

Minimum portable deployment checklist
-------------------------------------

To keep the workflow runnable after copying it to a new place or exposing
`auto_prepare_and_run_diff_gwas.pl` on PATH, verify all of the following:

[ ] `auto_prepare_and_run_diff_gwas.pl` is present
[ ] `DiffGWASDeps/` is a sibling of the main entry script
[ ] `configs/` is a sibling of the main entry script
[ ] `DiffGWASDeps/` contains the helper Perl scripts, shell runners, and SAS
    templates/macros
[ ] `DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl` is present
[ ] `DiffGWASDeps/SAS_ODA_Runner.pm` is present
[ ] `DiffGWASDeps/sas_oda_session_server.py` is present if persistent SAS
    ODA sessions are used
[ ] the target machine can run Perl, Bash/Cygwin, Python+saspy, and the SAS
    ODA helper path

If you want old direct helper commands to keep working, keep the root-level
Perl wrappers too.

What lives in DiffGWASDeps
--------------------------

- Perl helper scripts for merge/diff/standardize/extract
- shell runners for SAS ODA plotting
- SAS templates/macros used by those runners

Compatibility wrappers
----------------------

Small root-level wrapper scripts are kept for older manual commands such as:

  perl extract_significant_diff_gwas.pl ...
  perl extract_single_snp_wide_diff_gwas.pl ...

Those wrappers forward into the real implementations in DiffGWASDeps/.

Portability goal
----------------

The intended steady-state layout is:

- users run auto_prepare_and_run_diff_gwas.pl from PATH or from the repo root
- the script resolves its helpers from the sibling DiffGWASDeps/ directory
- shell runners inside DiffGWASDeps resolve WORKDIR relative to themselves
- if historical global `PERL5LIB` helper copies are absent, the vendored
  `DiffGWASDeps/` helper stack is used first

If the startup dependency check fails, verify that DiffGWASDeps contains the
required files and sits beside auto_prepare_and_run_diff_gwas.pl.

Quick SAS ODA helper debug order
--------------------------------

When `run_sas_codes_or_script_in_ODA.pl` misbehaves, keep the first probes
small and layered:

1. login only:
   - `perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl --check-sas-oda-login-only`
2. file-only listing / file-info:
   - `perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl --output-prefix oda_file_probe --dir4listing '~' --file-info '~/importallmacros_ue.sas'`
3. small upload/download round trip:
   - `perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl --output-prefix oda_transfer_probe --upload-file ./small_test.txt --download-file '~/small_test.txt' --download-local-path ./small_test.roundtrip.txt`
4. smallest SAS submit:
   - `perl ./DiffGWASDeps/run_sas_codes_or_script_in_ODA.pl --output-prefix oda_code_probe --code "%put HELLO_FROM_ODA;"`
5. then only after that move into `%include`, plotting macros, or forced
   macro-bootstrap tests

Helpful debug artifacts now include:

- `<output-prefix>/output.run.status.json`
- `<output-prefix>/output.macro_bootstrap.log.txt`
- `<output-prefix>/output.html.info.txt`

Interpretation tips:

- `Upload step: macro bootstrap helper: importallmacros_ue.sas ...`
  means the tiny helper upload finished; it is not the real macro bootstrap
- `SAS ODA macro bootstrap started at ...`
  means the wrapper has entered the default `~/Macros` autoload submit
- if a manual in-SAS `%importallmacros_ue(...)` probe is fast but the forced
  helper bootstrap still stalls, the bottleneck is the helper bootstrap path,
  not the `~/Macros` directory contents

Crowding control for Manhattan plots
------------------------------------

The generalized runner config now supports:

- `manhattan_fig_height`
- `local_manhattan_fig_height`
- `local_max_hits_per_fig`

When these are not supplied, the pipeline derives reasonable defaults from the
number of plotted tracks. Local top-hit Manhattan plots also support batching:
if the number of selected top loci exceeds `local_max_hits_per_fig`, the first
image keeps the base PNG name and later batches are written as `_part2`,
`_part3`, and so on. The downloaded HTML wrapper includes all generated panels.

Step-aware reruns
-----------------

`auto_prepare_and_run_diff_gwas.pl` now supports targeted reruns so you do not
have to rerun the whole workflow when only one stage matters.

Useful CLI patterns:

- `perl auto_prepare_and_run_diff_gwas.pl --spec your_spec.json --list-steps`
- `perl auto_prepare_and_run_diff_gwas.pl --spec your_spec.json --step plot_local_gtf`
- `perl auto_prepare_and_run_diff_gwas.pl --spec your_spec.json --from-step extract_wide_subset`
- `perl auto_prepare_and_run_diff_gwas.pl --spec your_spec.json --plot-manhattan --force`

Supported named steps currently include:

- `merge_raw`
- `sort_long`
- `diff_pairs`
- `standardize_diff`
- `extract_wide_subset`
- `plot_manhattan`
- `plot_local_manhattan`
- `plot_local_gtf`
- `cleanup_shared_plot_data`

The Perl MCP tool `auto_prepare_and_run_diff_gwas` exposes the same ideas with
arguments like `list_steps`, `step`, `from_step`, `to_step`, and convenience
booleans such as `plot_local_gtf=true`.
