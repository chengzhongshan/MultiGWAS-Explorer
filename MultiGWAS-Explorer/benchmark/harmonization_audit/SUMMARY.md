# Full-source allele harmonization audit

The audit scanned 44,930,352 long-format rows and evaluated 22,638,007 pair-specific variant keys. Both female and male strata were present for 22,292,345 keys. Every paired record was a direct A1/A2 match; there were no swapped, complemented, incompatible, duplicate-tag, or effect-allele-frequency difference greater than 0.20 pairs.

| Pair | Keys evaluated | Both present/direct | Female only | Male only | A/T or C/G ambiguous | EAF difference >0.10 | Revised eligible |
|---|---:|---:|---:|---:|---:|---:|---:|
| ALL | 7,606,736 | 7,518,467 | 25,796 | 62,473 | 1,138,546 | 0 | 6,379,921 |
| EUR | 7,919,740 | 7,838,328 | 39,643 | 41,769 | 1,187,692 | 0 | 6,650,636 |
| ASN | 7,111,531 | 6,935,550 | 86,701 | 89,280 | 1,059,774 | 1 | 5,875,776 |
| Total | 22,638,007 | 22,292,345 | 152,140 | 193,522 | 3,386,012 | 1 | 18,906,333 |

The revised default excludes the 3,386,012 strand-ambiguous direct matches. The single ASN pair with effect-allele-frequency difference greater than 0.10 remained below the prespecified 0.20 exclusion threshold.

Machine-readable counts are in `harmonization_audit.tsv`; up to 100 examples per exception class and pair are in `harmonization_exceptions_sample.tsv`.
