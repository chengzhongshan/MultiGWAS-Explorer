# Harmonization-filtered unfiltered schizophrenia QC

Input: `PGC_SCZ_female_vs_male_diff_effects.stdized.tsv.gz`

All 22,292,345 direct paired rows were streamed. Before genomic-inflation and candidate calculations, 3,386,012 A/T or C/G strand-ambiguous rows were excluded and the aligned mean effect-allele-frequency threshold of 0.20 was applied. No row exceeded the frequency threshold. This left 18,906,333 eligible rows. Lambda GC used raw `GROUP1_Z`, `GROUP2_Z`, and `DIFF_Z` on eligible autosomes; histogram bin width was 0.0001.

| Comparison | Female lambda GC | Male lambda GC | Difference lambda GC | Eligible autosomal rows | Eligible chromosome-X rows | Strand-ambiguous rows excluded |
|---|---:|---:|---:|---:|---:|---:|
| ALL | 1.3148 | 1.3662 | 1.0200 | 6,379,921 | 0 | 1,138,546 |
| EUR | 1.2900 | 1.3574 | 1.0211 | 6,456,348 | 194,288 | 1,187,692 |
| ASN | 1.0831 | 1.0651 | 0.9974 | 5,715,993 | 159,783 | 1,059,774 |

No raw differential P value reached 5e-8. Eligible row counts at raw `DIFF_P < 1e-5` were 41 ALL, 39 EUR, and 59 ASN. These are correlated row-level candidates, not independent loci. The minimum P variants remained rs185665940 in ALL and EUR (P=4.321342e-7) and rs42067 in ASN (P=4.146884e-7).

The earlier `benchmark/unfiltered_qc/` output is retained as a pre-exclusion sensitivity analysis. Its difference lambda values were almost identical, while its P<1e-5 counts were higher because strand-ambiguous rows were included.
