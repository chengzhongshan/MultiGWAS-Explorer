# Debugging the AOA GWAS comparison pipeline

This is a reproducible record of the 2026-10-09 AOA rerun and the checks to
make before trusting its figures. Run the commands from `AOA_GWAS_Data` unless
another directory is shown. Keep failed `run_*` and `upload_*` directories
until their logs have been examined.

## Input and expected outputs

The source is `DS_ALL_and_MP2PRT_plus_meta.gz`, a merged-wide GWAS table with
DS_ALL, MP2PRT, and precomputed meta-analysis P values. The inferred spec is
`MultiGWAS-Explorer/MultiGWAS-Explorer/configs/auto_AOA_GWAS_Data_diff_merged_from_dir.spec.json`.
The wide plot table contains 13,465,870 rows. The genome-wide Manhattan input
retains 1,125,786 rows having **at least one displayed GWAS P < 0.05**. This
filter applies only to the genome-wide Manhattan plot; local windows keep all
available SNPs so their shape and LD context are preserved.

```bash
../MultiGWAS-Explorer/MultiGWAS-Explorer/auto_prepare_and_run_diff_gwas.pl --gwas-dir .
../MultiGWAS-Explorer/MultiGWAS-Explorer/auto_prepare_and_run_diff_gwas_with_gnuplot.pl --gwas-dir .
```

The first command uses SAS ODA and can take hours when many local GTF loci
remain. The second uses local gnuplot. An overall exit status of zero does not
prove every SAS stage succeeded: the wrapper can continue with gnuplot after a
SAS failure. Inspect the named SAS output and the stage log.

## What failed, and what the rerun established

| Observation | Diagnosis and action |
| --- | --- |
| `merged_plotwide.manifest.tsv` was missing or incomplete | The generated wide file alone does not prove it matches the current spec. Let the extractor regenerate the manifest and wide table, then verify both exist and the header has the requested pair and meta columns. The 2026-10-09 regeneration read and wrote 13,465,870 rows. |
| Genome-wide SAS and gnuplot plots succeeded | Verify the row count in the compact-subset log and open the PNG. The AOA SAS plot has four tracks: Meta, MP2PRT, DS_ALL, and their differential P. |
| The original combined SAS local Manhattan submit lost its SAS session | Its diagnostic file reported `SAS_ODA_REMOTE_SESSION_TERMINATED` (exit 74), followed by `SASIOConnectionTerminated: No SAS process attached`. A SAS WORK, quota, or memory problem is plausible, but the saved diagnostics do **not** prove which limit was reached. Preserve the run directory and log. |
| Automatic local GTF plotting completed | `AOA_GWAS_Data_diff_SAS_local_top_hits_with_gtf.progress.json` reported 47 of 47 loci complete. Check the final HTML and each linked PNG; do not infer completion from the top-level wrapper alone. |
| A one-locus compact SAS local Manhattan test succeeded | For `rs12028518` and a ±1 Mb window, tabix extracted 12,329 rows into a 696,229-byte gzip instead of uploading the 721 MB whole-genome wide table. SAS returned no `ERROR` lines and the PNG and HTML downloaded. |
| An explicit-target download asked for a nonexistent LD audit | Explicit targets bypass LD clumping, so SAS does not create the audit file. The runner now omits that optional download for explicit targets. This was a download-manifest bug, not a plot failure. |
| The first 47-locus SAS run said `0 ERROR lines`, yet its log contained `CHR=.` and `_ERROR_=1` for chromosome-X rows | SAS imported literal `X` into a numeric `CHR` field. The combined SAS upload now converts `X` to `23` (and `Y` to `24`); gnuplot locus files retain their original chromosome strings. The corrected 47-locus input contains 16,191 chr23 rows and no literal `X` rows. Check for `Invalid data for CHR` and `_ERROR_=1`, not only lines beginning `ERROR:`. |
| The chrX meta track was blank although `META_P` existed | The SAS Manhattan macro kept only rows with nonmissing differential P before building all tracks. Near `rs145027345`, 2,984 rows have meta P; near `rs145598300`, 2,581 do. None of those rows shares a differential P value. Retain a row when **any displayed P track** is present, then plot each track from its own P column. Do not infer absent meta-analysis signals from a blank track without counting nonmissing source values. |
| Two chromosome-X rsIDs appeared as duplicate lead entries | The merged-wide table has more than one allele row for each ID. For explicitly requested leads, select one representative allele row per request using the smallest focus P; leave all rows in the surrounding scatter data. |
| The combined SAS local Manhattan plot lacked chromosome names, then the first added labels overlapped gene symbols | Build each ordered group label from `SNP:gene:chrN` (using `chrX` for numeric 23). Render the three pieces vertically and offset the chromosome label to the right by an amount derived from font and figure size. Inspect a one-locus plot and a crowded panel after changing annotation positions. |

For the 47 GTF lead loci, tabix returned 426,369 unique rows across the ±1 Mb
windows. Their combined gzip was 23,778,835 bytes, about 3.3% of the original
721,502,597-byte wide gzip. This is the local Manhattan input used for the
multi-locus verification run.

The gnuplot full comparison completed, including genome-wide Manhattan, local
Manhattan, and local GTF figures. Its local extraction reused the indexed BGZF
cache. A tabix cache is valid only for the matching source size and modification
time; let the indexer rebuild it after source changes.

The chrX correction does not merge meta and differential P values across
allele rows. Each plotted point retains the P value supplied on its own row.

The final 2026-10-09 SAS rerun used the 47 completed GTF leads and the ±1 Mb
compact local Manhattan input. It returned 47 distinct lead rows and four
nonempty, readable PNG panels, with meta points visible at the chrX loci and
vertical `SNP`, gene, and `chrN` labels under every panel. The SAS submit
reported zero `ERROR` lines. The runner downloaded all four panels and cleaned
its successful remote temporary files. The gnuplot full rerun and its local
plots also completed. The original automatic SAS local Manhattan submission
still has the large-upload failure described above; use the completed-lead
compact rerun below for that stage.

## Inspect a failed SAS stage

From the pipeline directory, locate the latest stage and inspect its status:

```bash
cd ../MultiGWAS-Explorer/MultiGWAS-Explorer
ls -dt run_local_hits_manhattan_png_* | head
cat run_local_hits_manhattan_png_*/output.non_retryable_remote_termination.txt
cat run_local_hits_manhattan_png_*/output.run.status.json
```

Use the specific latest directory rather than a wildcard when multiple runs
exist. Search its SAS log for `ERROR:`, `WARNING:`, and `SAS_ODA_REMOTE_SESSION_TERMINATED`.
If SAS never returned a complete log, record the remote termination as an
unknown cause; do not label it a proven WORK exhaustion. For the GTF sequence:

```bash
cat AOA_GWAS_Data_diff_SAS_local_top_hits_with_gtf.progress.json
test -s AOA_GWAS_Data_diff_SAS_local_top_hits_with_gtf.html
```

The progress file lists each locus and whether it completed. Resume the
pipeline without deleting its progress file or completed per-locus HTML/PNG
files. Upload and plot one locus at a time to limit SAS ODA storage and memory.

## Keep local SAS input small

The local Manhattan runner now extracts requested loci from the tabix-indexed
merged-wide table before SAS upload when `TARGET_SNP_LIST` is set. This path
retains all SNPs in the selected windows; the P < 0.05 genome-wide filter is
not used. A ±1 Mb window is a practical first rerun for AOA. A very large
`--local-gtf-window-bp` (especially over 5 Mb) increases upload, SAS WORK,
memory, and rendering costs. Try a smaller window and record the actual window
used when comparing figures.

For one target GTF plot, run:

```bash
../MultiGWAS-Explorer/MultiGWAS-Explorer/auto_prepare_and_run_diff_gwas.pl \
  --spec ../MultiGWAS-Explorer/MultiGWAS-Explorer/configs/auto_AOA_GWAS_Data_diff_merged_from_dir.spec.json \
  --step plot_local_gtf --target-snps rs12028518 \
  --local-gtf-window-bp 1e6 --no-gnuplot-fallback-on-sas-failure
```

For a combined local Manhattan rerun of the 47 completed GTF leads, run from
the pipeline directory:

```bash
python3 - <<'PY'
import json
from pathlib import Path

progress = json.loads(Path('AOA_GWAS_Data_diff_SAS_local_top_hits_with_gtf.progress.json').read_text())
config = json.loads(Path('configs/auto_AOA_GWAS_Data_diff_merged_runner.json').read_text())
config['TARGET_SNP_LIST'] = ','.join(locus['snp'] for locus in progress['loci'] if locus['status'] == 'complete')
config['LOCAL_WINDOW_BP'] = '1e6'
config['LOCAL_OUTPUT_PREFIX'] = 'AOA_GWAS_Data_diff_SAS_local_top_hits_manhattan'
out = Path('cache/local_manhattan_reuse/aoa_completed_gtf_leads.json')
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text(json.dumps(config, indent=2))
print(f'{len(config["TARGET_SNP_LIST"].split(","))} loci -> {out}')
PY
RUNNER_CONFIG_JSON="$PWD/cache/local_manhattan_reuse/aoa_completed_gtf_leads.json" \
  OPEN_RESULT=0 KEEP_REMOTE_PLOT_DATA=0 \
  bash DiffGWASDeps/run_sas_oda_local_top_hits_manhattan_download_png.sh
```

Check that the progress JSON describes the current AOA spec before using its
leads. Target SNPs and their window width are part of the compact-cache key;
a changed target set or width creates a new subset. The runner downloads plots
and removes its successful temporary SAS ODA files afterward.

## Verification before publishing figures

1. Confirm the wide table and manifest are nonempty and the manifest matches
   the current source/spec. Confirm the compact Manhattan count and that Meta
   P appears as its own track.
2. Check each SAS stage log independently for `ERROR:` or remote termination.
   Search for `Invalid data for CHR` and `_ERROR_=1` as well: SAS can report
   bad imported rows while the helper still summarizes `0 ERROR lines`. A
   gnuplot fallback does not repair a failed SAS artifact.
3. Check GTF progress reaches `complete == total` and `remaining == 0`.
   Verify the final HTML links resolve to nonempty files; open representative
   PNGs to inspect SNP/gene/chromosome labels, genes, signed-LD colorbar, and
   window boundaries.
4. Keep the failed SAS diagnostic folders until the failure is explained. Clean
   only successful temporary folders after all output files have been checked.

The scientific rule is to preserve every SNP inside a local window, but to use
the nominal-significance subset for the genome-wide Manhattan plot. The
operational rule is to verify each plot stage separately and minimize the data
uploaded to SAS ODA.
