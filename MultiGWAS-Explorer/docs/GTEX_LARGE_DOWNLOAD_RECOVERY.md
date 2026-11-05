# GTEx boxplots after an abruptly terminated SAS ODA session

The original MCP4SAS investigation's saved
`SAS_Output_RandID979354/output.html.info.txt` contains the exact
MAP3K19/R3HDM1 request with `UseGeneratedDsd=0`. Its terminal status records
`sas_oda_remote_session_termination`, exit 74, and no definitive SAS runtime log.
The macro bootstrap completed. A resource failure is plausible, but the log
does not prove whether disk, memory, a service limit or another cause ended SAS.
The failed large job was not resubmitted.

## Why the original macro takes the expensive path

The original macro in the user's library:

- Treats `UseGeneratedDsd=0` as a request to build the dataset. `PreviousDsd=tgt`
  is only used for cached plotting when `UseGeneratedDsd=1`.
- Checks for `GTEx.exp_all`; if it is absent, invokes the UCSC importer on
  `https://cells.ucsc.edu/gtex8/exprMatrix.tsv.gz`.
- Imports the whole matrix before selecting the two genes. The importer defaults
  to up to one million selected columns and an extremely large record length.
  Putting the file in WORK does not remove all resource limits.
- Uses an outer condition that can rebuild from `exp_all` even when cached
  plotting was requested, modifies the selected previous dataset in place, and
  contains `%sysfunc(exists(...))` where SAS's dataset function is `exist`.

A LIBNAME assignment made in SAS Studio is not automatically an assignment in a
new SASPy/ODA connection. Explicitly assign `GTEx` to the folder containing the
prepared table in the same submission. Do not assume an old Studio libref exists.

## Cache-only alternative

### Run a self-contained learning example

From the project directory in Cygwin:

```bash
bash examples/run_gtex_boxplot_demo.sh
```

The script prints the SAS learning example and exact submission command. It
combines the local macro definition with `examples/sas/gtex_boxplot_demo.sas`
before submission; an online SAS `%include` cannot read a local Windows file.
The example generates synthetic data for three genes, three tissues, two sexes
and eight observations per group, then plots only MAP3K19 and R3HDM1. It also
prints 12 observations and a summary table into the HTML result. The command
prints the generated HTML and log paths at the end.

On 2026-11-04 this command passed through SAS ODA with no SAS ERROR or WARNING
lines and the actual runtime message:

```text
GTEX_DEMO_PASS: source=144 selected=96 genes=2
```

Artifacts are retained locally in `SAS_Output_gtex_cache_demo_Qtc9xe`: the HTML
is `plot/sas_res__rLfA.html`, the log is `plot/output.html.info.txt`, and the
embedded plot was extracted unchanged as `gtex_boxplot.png` for visual review.
The reviewed plot has two gene panels, three tissue categories, and female/male
boxplots. These expression values are synthetic learning data, not GTEx findings.

### Use your prepared GTEx table

`examples/sas_macros/Boxplots4GenesInGTExV8_Cached.sas` defines
`%Boxplots4GenesInGTExV8_Cached` without modifying the user's original local or
online macro. It requires a prepared long table with:

- `rownames`: gene symbol, or an identifier ending in `|SYMBOL`;
- numeric `exp`: expression already transformed as intended for plotting;
- `cluster`; and
- `sex` for `bysex=1`, or `Tissue` otherwise.

It filters the requested genes, checks each is present after optional WHERE
filtering, writes a separate output table and draws the panels. It preserves the
expression scale and the cache. Missing input, invalid columns and requests with
`UseGeneratedDsd=0` fail with an actionable error. There is no HTTP or full-matrix
download path. Existing outputs require explicit `replace=1`, and the output
member name must differ from the source cache member.

From the `MultiGWAS-Explorer` project directory in Cygwin, set the real
prepared-data folder below:

```bash
{
  cat examples/sas_macros/Boxplots4GenesInGTExV8_Cached.sas
  cat <<'SAS'
libname GTEx '/home/YOUR_ACCOUNT/PREPARED_DATA_FOLDER';
%Boxplots4GenesInGTExV8_Cached(
  genes=MAP3K19 R3HDM1,
  dsdout=exp,
  bysex=1,
  UseGeneratedDsd=1,
  PreviousDsd=tgt,
  Lib4PreviousDsd=GTEx,
  WhereFilters4Boxplot=%str(),
  boxplot_width=1200,
  boxplot_height=1000
);
SAS
} | perl run_sas_codes_or_script_in_ODA.pl --code -
```

### Tested real-data learning script

On 2026-11-04, `examples/run_real_gtex_boxplots.sh` successfully plotted real
MAP3K19 and R3HDM1 expression through `run_sas_codes_or_file_in_ODA.pl`. That
spelling is a compatibility alias for the main runner and uses its configuration.

Run from the project directory in Cygwin:

```bash
# Default: read the user's existing full local expression matrix.
bash examples/run_real_gtex_boxplots.sh

# Faster reuse: read the small two-gene subset saved by the successful run.
bash examples/run_real_gtex_boxplots.sh \
  SAS_Output_gtex_cache_real_ATeUu4/gtex_two_genes.tsv \
  '~/data/GTEx_V8'
```

The first optional argument is a local matrix/subset TSV; the second is the
online directory containing `umap.sas7bdat`. The script prints its exact runner
command and creates a fresh output directory on each invocation. The submitted
SAS program, input manifest, small expression subset, upload table, HTML, terminal
status and SAS log remain together for review. The original book/library files
and existing GTEx tables are read-only inputs.

The full local source used for this run was
`E:/LongCOVID_HGI_GWAS/GTEx_Exp_SAS_Data/exp.matrix/exprMatrix_ensgSymbolCombined.tsv`
(4,261,582,197 bytes). Perl streamed its 56,200 gene rows and retained just:

- `ENSG00000176601.12|MAP3K19`;
- `ENSG00000048991.16|R3HDM1`.

The two-gene subset is 661,175 bytes, with SHA-256
`fa2e5ced0edbba1b1e23d8247c748641641c2f86cba43ca33083f002ef5afa59`.
The saved learning script was executed against that subset. It uploaded a
2,020,843-byte long table, rather than transferring the full matrix, and joined
sample IDs directly to `~/data/GTEx_V8/umap.sas7bdat` using `sample_id=seq_id`.
It checks for duplicate phenotype IDs, unmatched samples and invalid expression
values before plotting. The uniquely named upload is deleted after SAS reads it.

The preparation step applies `log2(raw_exp+1)` once, matching the original macro's
transformation. The cached plotting macro then preserves those transformed
values. The original rule excluding clusters containing `Outlier` is retained;
none of this run's matched rows were removed by that rule.

The actual runtime result was:

```text
GTEX_REAL_PASS: input=34764 samples=17382 joined=34764 plotted=34764 genes=2
```

Each gene has 17,382 observations across 54 tissue categories. The SAS log had
zero ERROR and WARNING lines, and terminal JSON recorded completed success.
The HTML includes female/male boxplots and tissue-specific descriptive summaries.
These are descriptive expression plots; the script does not perform a sex
differential-expression significance test.

Artifacts are retained locally in `SAS_Output_gtex_cache_real_ATeUu4`:

- `real_gtex_boxplot.sas`: complete submitted SAS code, including the macro;
- `input_manifest.json`: input checksum, counts and transformation/join rules;
- `plot/sas_res_g7BwS.html`: real boxplots and descriptive tables;
- `plot/output.html.info.txt`: actual SAS log;
- `real_gtex_boxplot.png`: the embedded PNG extracted unchanged for visual review.

The plot was visually checked for two gene panels, readable tissue labels and
the female/male legend. New learning scripts use Bash, Perl and SAS.

## If a prepared cache is unavailable

Prepare a gene subset on the local computer before importing into SAS. The Perl
utility streams an existing local UCSC-style tab-delimited matrix or gzip file
and retains its sample header and every matching gene row:

```bash
perl tools/subset_gtex_expression.pl \
  --input /path/to/exprMatrix.tsv.gz \
  --genes 'MAP3K19 R3HDM1' \
  --output gtex_two_genes.tsv
```

It never overwrites the input or an existing output. Missing genes, ragged selected
rows and corrupt/truncated gzip files fail before the output is published. Only
one matrix row is held at a time. This utility does not download the full file
or construct sex/cluster metadata. Import the small subset and join matching
sample annotations to create the long cache before calling the cache-only plot.
It deliberately does not guess missing phenotype information.

For the original macro, a minimal repair is to prioritize `UseGeneratedDsd=1`,
validate the prepared dataset before starting work, replace `exists` with `exist`,
filter into a separate WORK table, and require explicit opt-in before any full
matrix download. Increasing MCP's timeout alone does not bound matrix resource
use or recover an abruptly terminated SAS session.

## Verification

`tests/test_gtex_expression_subset.pl` passed 12 offline assertions for plain/gzip
input, exact expression preservation, missing genes, existing output, ragged rows
and truncation. `prove -v tests/test_live_gtex_cache_only.pl :: --live` passed three isolated
SAS ODA cases: a 24-row synthetic source produced the correct 16-row two-gene
output while preserving the source, missing-cache refusal, and rebuild refusal.
Explorer's passing live run is retained in the ignored
`SAS_Output_gtex_cache_JW2JzC` folder, with the test transcript in
`SAS_Output_adaptation_20261103/gtex_live.tap`.

Development exposed a FINDW detail: its `T` modifier trims the delimiter argument
as well as the search strings, so a space-only delimiter becomes empty. The final
macro strips the gene token and uses only `I` for case-insensitive matching. See
the [SAS FINDW documentation](https://support.sas.com/documentation/cdl/en/lefunctionsref/63354/HTML/default/p16rdsa30vmm43n1ej4936nwa01t.htm).
Earlier failed prototypes from the original investigation are not counted as
Explorer evidence. No large expression matrix or original user dataset was used
in these tests. See the [adaptation record](SAS_ODA_ADAPTATION_20261103.md) for
the commands used to validate this port.
