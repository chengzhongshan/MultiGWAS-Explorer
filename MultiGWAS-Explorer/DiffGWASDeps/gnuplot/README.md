# PDL + gnuplot alternative plotting helpers

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
- Nearby query SNPs whose SNP-centered windows overlap are rendered in one
  shared local-Manhattan/local-GTF locus. Each query is labeled once with a
  separate leader line, and the stable combined PNG points to that single
  coordinate-aware plot instead of a stale multi-panel composite.
- High-LD highlighting is disabled by default to keep wide locus plots
  readable. Enable it with `--highlight-high-ld-snps`, or supply an explicit
  list with `--ld-snps` (which also enables highlighting).
  High-LD variants can then be supplied directly or resolved automatically.
  The resolver checks a normalized local HaploReg TSV/SQLite cache first and
  queries HaploReg only for cache misses, saving successful web results for
  reuse. The default population is EUR and the default display threshold is
  r2 >= 0.8; use `--ld-population`, `--ld-r2-threshold`, and `--ld-cache`
  to override them. Use `--ld-marker-symbol` and `--ld-marker-color` to choose
  star/plus/cross/circle/square/triangle/diamond and a named or `#RRGGBB`
  color. Markers are overlaid on every association track in both
  local-Manhattan and local-GTF plots.
- A reusable tabix-indexed GTF cache is queried for each locus; it is not
  rebuilt for every SNP when the cached BGZF file and index remain valid.
- The repository wrapper now tunes the genome-wide Manhattan palette and
  top-of-panel GWAS labels to follow the SAS ODA multi-track style more
  closely, while still prioritizing reliability and portability.
- Small visual differences can still remain because gnuplot and SAS ODA do not
  rasterize points identically.

Usage examples

Generate a Manhattan from a prepared wide subset:

  perl pdl_gunplot_manhattan.pl --data /path/to/wide_subset.tsv.gz --outdir ./gnuplot_out

Generate a local locus plot for SNP `rs12345` with a GTF:

  perl pdl_gunplot_local_gtf.pl --data /path/to/wide_subset.tsv.gz --snp rs12345 --gtf /path/to/gencode.gtf --outdir ./gnuplot_out

Wrapper

`auto_prepare_and_run_diff_gwas_with_gunplot.pl` (in repository root) wraps
the existing data-preparation logic (optionally calling the original
`auto_prepare_and_run_diff_gwas.pl` to produce a wide subset) and then calls
these plotting helpers. Use `--data-gz` to supply an existing wide subset.
The wrapper name is retained for backward compatibility, but the canonical
implementation directory is `DiffGWASDeps/gnuplot/` and the default reusable
GTF cache is `.gnuplot_gtf_cache`.

For example, render two nearby query SNPs in one locus and highlight a known
LD-linked SNP:

  perl auto_prepare_and_run_diff_gwas_with_gunplot.pl \
    --spec configs/spec_pgc_scz_sex_common_automation.json \
    --target-snps rs2070788,rs383510 \
    --plots local_gtf \
    --ld-snps rs2837646 \
    --ld-marker-symbol diamond --ld-marker-color '#6A0DAD'

The SAS automation entry point also defaults to a one-time gnuplot fallback
when a plotting job is specifically classified as SAS ODA WORK/storage
exhaustion (exit 73), or when SASPy terminates with native exit 134/SIGABRT.
Use `--no-gnuplot-fallback-on-sas-space` or
`--no-gnuplot-fallback-on-sas-failure` to disable these policies.
Authentication, syntax, data, and ordinary SAS program errors do not trigger
fallback.

For a local-first high-LD workflow, prepare a small candidate file containing
the query SNPs and extract all of their proxies from previously downloaded
`LD_EUR.tsv.gz` (or another population archive):

  HAPLOREG_LD_ARCHIVE_DIR=/path/to/haploreg/data \\
    bash DiffGWASDeps/build_haploreg_ld_candidate_cache.sh \\
      --candidates query_snps.csv --populations EUR --min-r2 0.8 \\
      --retain-all-proxies --output cache/query_high_ld.tsv

Pass the result with `--ld-cache cache/query_high_ld.tsv`. If a query is not
present, the default web fallback requests only that SNP from HaploReg and
adds the returned rows to the reusable local web cache. Use
`--no-ld-web-fallback` (gnuplot) or `--no-local-ld-web-fallback` (SAS entry
point) for strictly offline execution.
