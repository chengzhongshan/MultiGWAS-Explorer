# Nominated-variant detectable-effect analysis

Date: 19 August 2026

The calculation uses the normal-Wald approximation, two-sided alpha, 80%
target power, independent female and male strata, and the observed standard
errors at each variant. Effects are log odds ratios. For a test with standard
error `SE`, the minimum detectable absolute effect is
`(qnorm(1-alpha/2) + qnorm(0.80)) * SE`. This preserves imputation and
variant-availability effects reflected in the supplied summary statistics.

## Pooled-stratum results at alpha = 5e-8

| SNP | Test | MAF | Variant-specific cases/controls | 80%-power detectable effect |
|---|---|---:|---|---|
| rs185665940 | Female log(OR) | 0.0135 | 13,456 / 28,297 | OR <=0.618 or >=1.617 |
| rs185665940 | Male log(OR) | 0.0130 | 25,906 / 26,707 | OR <=0.649 or >=1.541 |
| rs185665940 | Female-minus-male log(OR) | 0.0130 | both rows above | OR ratio <=0.524 or >=1.909; observed-effect power 34.6% |
| rs10166057 | Female log(OR) | 0.0720 | 23,031 / 45,626 | OR <=0.855 or >=1.169 |
| rs10166057 | Male log(OR) | 0.0755 | 39,967 / 41,159 | OR <=0.876 or >=1.141 |
| rs10166057 | Female-minus-male log(OR) | 0.0720 | both rows above | OR ratio <=0.815 or >=1.227; observed-effect power 0.46% |

The result explains why unequal stratum-specific P values cannot establish a
sex-specific effect and why the nominated contrasts require substantially
larger effective sample sizes for genome-wide detection. Complete pooled,
EUR, ASN, genome-wide-alpha, and exploratory-alpha results are in
`top_variant_detectable_effects.tsv`; the extracted inputs and executable R
script are retained beside it.
