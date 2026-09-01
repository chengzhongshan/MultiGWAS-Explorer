# HaploReg LD analysis for nominated variants

Date: 19 August 2026

The requested SAS macro
`H:/F_Queens/360yunpan/SASCodesLibrary/SAS-Useful-Codes/Macros/QueryMulti_LD_SNPs_at_Haploreg4.sas`
was run in SAS OnDemand for Academics for `rs185665940 rs10166057`, separately
for HaploReg populations `EUR` and `ASN`.

## High-LD proxy extraction (`r2 >= 0.5`)

| Query | EUR proxies excluding query | ASN proxies excluding query | Strongest observed proxy |
|---|---:|---:|---|
| rs185665940 | 1 | 0 | rs139205793, EUR r2=0.59, D'=0.83 |
| rs10166057 | 20 | 23 | EUR: rs3768641, r2=0.88; ASN: rs3768638 and rs3768639, r2=1.00 |

Complete annotation tables are in
`haploreg_ld_eur_r2ge0p5.tsv` and
`haploreg_ld_asn_r2ge0p5.tsv`. They include population frequencies,
HaploReg hg38 positions, r2, D', and functional annotations returned by the
macro.

## Pairwise relationship between the two nominated variants

Each variant was also queried separately with `ldThresh=0.01`, and the result
was filtered to the two nominated rsIDs. HaploReg returned only the query SNP
itself in EUR and ASN; the other nominated SNP was absent in all four queries.
Therefore, under this HaploReg/1000 Genomes query the supported statement is
`r2 < 0.01` in both EUR and ASN. The exact value below 0.01 is not reported.

This demonstrates very low pairwise LD but does not replace GCTA-COJO or
individual-level conditional regression. No ancestry-matched genotype
reference was available in the revision workspace, so conditional effect
estimates were not fabricated.

## Integration into top-hit selection

The local Manhattan/GTF and forest-plot workflow no longer uses a physical
distance window as its default definition of separate top hits. MAF-passing
candidates are ranked by raw differential P or the minimum within-stratum
association P, as appropriate, and greedily clumped through the
repository-local HaploReg query macro. The revision default prunes r2>=0.10 in
any configured population panel. Distance is retained only as an explicit
legacy mode or a labeled fallback when a lead has no usable LD response.

The real differential PGC benchmark contained three MAF-passing candidates at
raw P<1e-6. EUR/ASN clumping retained rs185665940 and rs753634 and pruned
rs7755143 under rs753634 at ASN r2=0.82. Both retained leads had complete EUR
and ASN responses. A separate four-SNP fixture demonstrated why distance is
insufficient: a 1-Mb window retained one lead, whereas LD clumping retained two
and pruned only documented proxy pairs. The audit files are
`ld_fixture_audit.tsv` and `pgc_diff_ld_audit.tsv`.

These outputs define LD-clumped lead SNPs; they do not establish conditional
association.

## Archive-backed acceleration for common hits

The official HaploReg v4 EUR and ASN LD archives were streamed in parallel and
filtered to the 11,665 unique rsIDs represented by 11,701 MAF-passing
common-association candidates. The 402-s scan retained 7,023,934
candidate-to-candidate edges at r2>=0.20. Building a reusable SQLite index and
running the first complete local clump took 118.832 s; exact repeat runs took
1.742-1.751 s. Archive-only processing retained 175 leads.

The SAS ODA hybrid validation imported the same cache, used live HaploReg only
when an archive query was absent, and completed in 1,126 s after upload and
bootstrap. It retained 206 leads, with 8,947 LD-pruned candidates and 2,548
labeled distance-fallback prunes. A live-query-only run did not complete in a
3,458-s monitored execution. Because these v4 archive files are based on the
older downloadable resource and begin at r2=0.20, the live path remains
available for lower thresholds and is the default fallback for archive misses.

## Reproducibility artifacts

- `run_haploreg_ld_benchmark.sas`: multi-query high-LD extraction.
- `run_haploreg_pairwise_ld_benchmark.sas`: low-threshold pairwise check.
- `haploreg_pairwise_ld.tsv`: filtered pairwise result.
- `ld_fixture_audit.tsv`: four-SNP LD-selector regression audit.
- `pgc_diff_ld_audit.tsv`: real differential top-hit clumping audit.
- `pgc_common_haploreg_archive_coverage.tsv`: concise archive coverage audit.
- `pgc_common_ld_cached_r2ge0p2_leads.tsv`: SAS hybrid common-hit leads.
- `pgc_common_ld_cached_r2ge0p2_audit.tsv`: full SAS hybrid audit.
- `haploreg_ld_benchmark_submit/` and `haploreg_pairwise_ld_submit2/`: SAS
  logs and HTML artifacts.
