# PLINK2 / 1000 Genomes LD validation

The LD-pruning and local-plot paths use direct phased PLINK2 calculations from
1000 Genomes Phase 3 as the primary method. HaploReg4 is fallback-only when the
configured genotype reference is unavailable. Missing or monomorphic variants
are reported as `NOT_ESTIMABLE`; absence is never converted to `r2=0`.

## Common-association reanalysis

Using PGC schizophrenia common-association candidates at `P < 5e-8`, phased
1000 Genomes EUR and EAS panels, `r2 >= 0.10`, a 1-Mb window, and the `ANY`
population rule:

- candidates evaluated: 11,703
- retained leads: 265
- candidates pruned by directly calculated LD: 11,438
- retained leads absent from the 1KG variant set and labeled not estimable: 145
- distance-fallback decisions: 0

The versioned result files are:

- `pgc_common_plink2_1kg_leads.tsv`
- `pgc_common_plink2_1kg_audit.tsv`
- `pgc_common_plink2_1kg_cache.tsv`
- `pgc_common_plink2_1kg_reference_status.tsv`

For the nominated pair, phased EUR `r2=0.00122378` for rs185665940 and
rs10166057. The EAS value is not estimable because rs185665940 is monomorphic
(alternate-allele count 0 of 1,008). rs753634 and rs7755143 are directly linked
at `r2=0.826822` in EUR and `r2=0.736677` in EAS, so rs7755143 is pruned.

## Three-SNP local-GTF heatmap validation

The reference SNP was rs10166057 and the requested targets were rs185665940,
rs10166057, and rs4852780. With `MAJOR4` (`EUR+AFR+AMR+EAS`), a 1-Mb PLINK2
window, and minimum `r2=0`, the cache contained 47,870 estimable proxy records
plus the reference row. The pooled-reference values were `r2=0.00173277` for
rs185665940 and `r2=0.134473` for rs4852780.

The Gnuplot local-GTF HTML/PNG render completed successfully. The SAS ODA path
completed numeric `LD_R2` augmentation and emitted its runnable local SAS script
in dry-run mode. Both renderers map `sign(Z) * r2` to one continuous diverging
scale, while `--ld-display-mode none` retains the original Z-score colormap.

The compact nominated-pair results are in `three_snp_local_ld_results.tsv`.
