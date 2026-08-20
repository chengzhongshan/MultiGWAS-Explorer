# LD-clumped top-hit benchmark

Date: 19 August 2026

The revised top-hit stage ranks MAF-passing variants by the applicable raw P
value and greedily removes candidates in HaploReg LD with each selected lead.
Differential mode uses raw `DIFF_P`; common-association mode uses the minimum
within-stratum association P. The default is r2>=0.10 in any configured
population panel. Distance is retained only as an explicit legacy method or a
labeled fallback after a failed LD query.

## Four-SNP regression fixture

The fixture contains two documented EUR LD pairs on chromosome 2. The former
1-Mb distance selector retained only rs10166057. At r2>=0.50, the LD selector
retained rs10166057 and rs185665940 and pruned rs3768644 under rs10166057
(r2=0.52) and rs139205793 under rs185665940 (r2=0.59). Both lead queries
completed successfully.

## Real differential PGC candidates

The MAF safeguard and raw-P<1e-6 criterion produced three candidates. EUR/ASN
clumping at r2>=0.10 under the ANY rule retained rs185665940 and rs753634 and
pruned rs7755143 under rs753634 at ASN r2=0.82. Both selected leads returned
complete EUR and ASN responses; no distance fallback was used.

The results are in `ld_top_hit_benchmark_summary.tsv`,
`reviewer_fixture_ld_independent_leads.tsv`, `reviewer_fixture_ld_audit.tsv`,
`reviewer_pgc_diff_ld_independent_leads.tsv`, and
`reviewer_pgc_diff_ld_audit.tsv`. “LD-clumped lead” is an operational pruning
label and is not evidence that a conditional effect survives.

## Downloadable-archive acceleration

HaploReg publishes population-specific LD archives. The EUR and ASN downloads
are 7.7 GB and 6.7 GB, respectively, and contain one query rsID per row followed
by semicolon-delimited proxy,r2,D-prime records. The new streaming cache builder
downloads/decompresses each population once and writes only edges for rsIDs in
the MAF-passing candidate table. SAS imports the combined cache, creates a
query-rsID/population index, and falls back to the web endpoint only for cache
misses.

The official archive and HaploReg documentation support r2>=0.2. Consequently,
the accelerated path rejects a requested threshold below its declared cache
minimum; the existing live-query path remains available when a lower cutoff is
required.

For the real common-association candidate table, the parallel EUR/ASN archive
scan completed in 402 seconds. It matched 11,235 of 11,665 candidate rsIDs in
EUR and 10,363 in ASN and retained 7,023,934 candidate edges. Building the
SQLite index and clumping all 11,701 candidate rows took another 118.832
seconds. The archive-only result contained 175 leads: 116 had complete cache
coverage, 20 had one population represented, and 39 had no archive row. It
pruned 8,839 candidates through LD and 2,687 through the explicitly labeled
distance fallback under cache-absent leads. Because the downloadable files are
the 2015 v4.0/1000 Genomes Phase 1 resource, current web queries remain the
default fallback for rsIDs absent from those archives.

Reusing the completed SQLite index reduced repeated full 11,701-candidate
clumps to 1.742-1.751 seconds. The repeated lead CSV and audit TSV had exact
SHA-256 parity (`326b481a...` and `f9db4fa3...`, respectively).

The same combined cache was then imported into SAS ODA and indexed by query
rsID and population. The hybrid run used current web queries only for archive
misses and completed the clumping stage in 1,126 seconds after dependency
upload and macro bootstrap. It retained 206 leads, pruned 8,947 candidates by
LD, and pruned 2,548 through the labeled distance fallback. Among retained
leads, 123 had complete cache responses, 20 combined cache and web responses,
39 had complete web responses, and 24 had no usable response and therefore
used the fallback. The 206-lead result differs from the archive-only 175-lead
result because current web data resolve some archive misses as low-LD rather
than applying the more aggressive 1-Mb fallback.

For comparison, the exhaustive web-only 11,701-candidate SAS job did not
finish during its 3,458-second monitored execution and wrote no lead/audit
result. The cache therefore converts the archive-covered majority of network
calls into indexed lookups while preserving an auditable web fallback.
