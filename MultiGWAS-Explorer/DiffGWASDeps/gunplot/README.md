PDL + gnuplot alternative plotting helpers

This folder contains initial proof-of-concept scripts to generate Manhattan
and local locus plots using gnuplot. They are intentionally lightweight and
meant to serve as a drop-in alternative to SAS ODA plotting for the
`auto_prepare_and_run_diff_gwas` pipeline.

Files
- `pdl_gunplot_manhattan.pl` - generate genome-wide multi-track Manhattan PNG
- `pdl_gunplot_local_gtf.pl` - generate a local locus PNG with optional GTF gene track

Notes
- Scripts read a wide-format TSV.gz with columns `CHR`, `BP`, `SNP`, and one or
  more P-value columns (detected automatically or passed via `--pcols`).
- `gnuplot` must be available in PATH. These scripts generate `.gp` files and
  run gnuplot to produce PNG output.
- The local GTF plot requires a local GTF file passed with `--gtf` to render
  gene rectangles; otherwise the gene track is skipped.
- The repository wrapper now tunes the genome-wide Manhattan palette and
  top-of-panel GWAS labels to follow the SAS ODA multi-track style more
  closely, while still prioritizing reliability and portability.
- Small visual differences can still remain because gnuplot and SAS ODA do not
  rasterize points identically.

Usage examples

Generate a Manhattan from a prepared wide subset:

  perl pdl_gunplot_manhattan.pl --data /path/to/wide_subset.tsv.gz --outdir ./gunplot_out

Generate a local locus plot for SNP `rs12345` with a GTF:

  perl pdl_gunplot_local_gtf.pl --data /path/to/wide_subset.tsv.gz --snp rs12345 --gtf /path/to/gencode.gtf --outdir ./gunplot_out

Wrapper

`auto_prepare_and_run_diff_gwas_with_gunplot.pl` (in repository root) wraps
the existing data-preparation logic (optionally calling the original
`auto_prepare_and_run_diff_gwas.pl` to produce a wide subset) and then calls
these plotting helpers. Use `--data-gz` to supply an existing wide subset.
