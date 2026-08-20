# Reviewer-driven differential-GWAS revision

Date: 19 August 2026

## Scientific scope

The schizophrenia analysis is an exploratory demonstration of MultiGWAS-Explorer, not a locus-discovery or genotype-by-sex interaction study. Raw differential P values drive inference and top-hit selection. A P value between 5e-8 and 1e-5 is suggestive and is not described as genome-wide significant, robust, female-specific, novel, causal, or conditionally independent.

## Implemented pipeline safeguards

- `diff_pairwise_gwas.pl` reports beta difference, SE, Z, two-sided raw P, allele relation, ambiguity, aligned mean-frequency difference, and harmonization status.
- The standalone script and automated configuration now default to exclusion of A/T and C/G pairs and aligned mean effect-allele-frequency differences greater than 0.20.
- The workflow preserves A1/A2 terminology and does not silently reinterpret the source alleles as REF/ALT.
- Raw `DIFF_P` is used for Manhattan tracks and top-hit thresholds. `STD_DIFF_Z` can be used for color mapping; the legacy standardized P field is not inferential.
- Top-hit selection defaults to MAF greater than 0.01, records the frequency source and group-specific MAF/INFO, and rejects unknown MAF unless explicitly overridden.
- Top-hit selection now defaults to greedy HaploReg LD clumping instead of physical-distance windows. Candidates are ranked by the applicable raw P value; after each lead is selected, candidates with r-squared at least 0.10 in any requested ancestry panel are pruned. Common-association mode ranks by the minimum within-stratum association P value, whereas differential mode ranks by raw `DIFF_P`.
- The ancestry panels, r-squared threshold, ANY/ALL population rule, and query-failure policy are explicit configuration fields. A physical-distance window is retained only as a labeled fallback after a failed HaploReg query or as an explicitly selected legacy mode.
- Every retained lead and removed proxy is written to an LD audit table with lead rank, ancestry panel, r-squared, query status, and selection method. These are LD-clumped lead variants, not proof of conditional association.
- Large candidate sets can use candidate-restricted caches built by streaming the official HaploReg v4 EUR/ASN LD archives. The archives support r-squared cutoffs of at least 0.20; cache misses fall back to the current web query, and lower cutoffs continue to use the live-query path.
- The unfiltered QC path rejects P-selected filenames unless a debugging override is supplied and now applies the same ambiguity/frequency eligibility policy before lambda GC and candidate enumeration.
- Chromosome X is counted separately and excluded from autosomal lambda GC.
- Manifests record parameters, row counts, exclusions, and output paths.

## Full-source allele audit

The audit scanned 44,930,352 long-format rows and evaluated 22,638,007 pair-specific keys. Both sex strata were present for 22,292,345 keys, all with direct A1/A2 agreement. There were no swapped, complemented, incompatible, duplicate-tag, or frequency-discordant (>0.20) pairs. Female-only and male-only counts were 152,140 and 193,522. Exclusion of 3,386,012 strand-ambiguous direct pairs left 18,906,333 eligible pairs.

## Harmonization-filtered unfiltered QC

| Pair | Female lambda GC | Male lambda GC | Difference lambda GC | Eligible autosomal rows | Eligible raw-P<1e-5 rows |
|---|---:|---:|---:|---:|---:|
| ALL | 1.3148 | 1.3662 | 1.0200 | 6,379,921 | 41 |
| EUR | 1.2900 | 1.3574 | 1.0211 | 6,456,348 | 39 |
| ASN | 1.0831 | 1.0651 | 0.9974 | 5,715,993 | 59 |

No eligible raw differential result reached 5e-8. The leading variants remained rs185665940 in ALL/EUR and rs42067 in ASN. The older pre-exclusion QC output is retained only as a sensitivity analysis.

## Official comparator benchmark

GWAMA 2.2.2 and EasyStrata 8.6 were executed in the pinned `multigwas-comparators:20260819` Ubuntu 24.04 image on the same 100,000 direct-allele rows. Each tool retained all rows. Maximum absolute P differences were 5.829e-7 for GWAMA `q_p-value` and 8.301e-8 for EasyStrata `CALCPDIFF`, with zero disagreements at P thresholds 0.05, 1e-5, and 5e-8. The same conclusion held separately for 86,473 non-ambiguous and 13,527 strand-ambiguous sensitivity rows.

Measured wall time/peak RSS were 2.73 s/17,484 KB for the direct MultiGWAS formula audit, 91.55 s/79,624 KB for GWAMA, and 4.78 s/97,520 KB for EasyStrata. These values describe one fixture and environment, not a general ranking.

## Agent-interface benchmark

All 10 CLI and 10 MCP-agent configuration-generation runs succeeded. Every agent run reproduced the same four-artifact SHA-256 set as the CLI path. Median elapsed time was 0.722 s (IQR 0.720-0.743) for CLI and 2.621 s (IQR 2.619-2.630) for the MCP path. Three of three injected missing-specification failures were detected and all corrected retries succeeded. The result is scoped to this deterministic task.

The timed MCP path did not submit or score a free-form natural-language question. It sent an exact JSON-RPC `tools/call` request to `auto_prepare_and_run_diff_gwas`; successful trials differed only in the output-log name. Thus, the benchmark evaluates deterministic dispatch, configuration generation, logging, checksum parity, error detection, and corrected retry behavior—not language-model comprehension or general AI performance. The exact payload and all expanded records are published under `benchmark/agent_interface/`.

## Nominated-variant validation

- rs185665940: pooled beta difference -0.519300, SE 0.102746, raw Z -5.05423, raw P 4.321e-7; female/male MAF 0.0135/0.0130 and INFO 0.788/0.782. At alpha 5e-8, its pooled sex contrast had 34.6% observed-effect power and an 80%-power detectable OR-ratio interval outside 0.524-1.909.
- rs10166057: pooled beta difference 0.0923935, SE 0.0324968, raw Z 2.84316, raw P 0.004467; female/male MAF 0.0720/0.0755 and INFO 0.951/0.961. Its corresponding observed-effect power was 0.46%, with an 80%-power detectable OR-ratio interval outside 0.815-1.227.
- The requested HaploReg SAS macro found pairwise r-squared below 0.01 between the two leads in EUR and ASN. High-LD proxy counts at r-squared at least 0.5 were 1 EUR/0 ASN for rs185665940 and 20 EUR/23 ASN for rs10166057.
- In the real differential top-hit benchmark, three MAF-passing candidates at raw P<1e-6 were reduced to two LD-clumped leads. Rs7755143 was pruned under rs753634 because HaploReg reported ASN r-squared=0.82; rs185665940 and rs753634 were retained with successful EUR/ASN queries.
- In a four-variant regression fixture, the former 1-Mb distance rule retained one lead, whereas EUR LD clumping at r-squared >=0.5 correctly retained two leads and pruned only the two documented high-LD proxies.
- The real common-association analysis contained 11,701 MAF-passing candidates at within-stratum P<5e-8. The one-time parallel EUR/ASN archive scan took 402 s and retained 7,023,934 candidate edges. SQLite indexing plus the first full local clump took 118.832 s, while indexed reruns took 1.742-1.751 s. Archive-only processing retained 175 leads. The SAS ODA hybrid cache-plus-web-fallback run retained 206 leads, pruned 8,947 candidates by LD and 2,548 by labeled distance fallback, and took 1,126 s after upload/bootstrap. A web-only run did not finish in its 3,458-s monitored execution and produced no result artifacts.

## Validation that could not be supported

Controlled individual-level genotypes and a suitable genotype reference were unavailable, so genotype-by-sex interaction, GCTA-COJO, and genotype-level conditional regression were not performed. Individual-level demographics/exclusions, complete cis-eQTL summary statistics, and adjusted GTEx expression/covariate data were also unavailable. Corresponding interaction, conditional-independence, causal-gene, colocalization, and sex-differential regulatory claims were removed.

## Reproducibility locations

- Full-source allele counts: `benchmark/harmonization_audit/`
- Harmonization-filtered unfiltered QC: `benchmark/unfiltered_qc_harmonized/`
- Official comparators: `benchmark/fixture_all_100k/` and `benchmark/official_tool_runs.tsv`
- Agent evaluation: `benchmark/agent_interface/`
- HaploReg outputs: `benchmark/reviewer_haploreg_*.tsv`
- LD top-hit benchmarks: `benchmark/reviewer_fixture_*leads.tsv`, `benchmark/reviewer_fixture_ld_audit.tsv`, `benchmark/reviewer_pgc_*_ld_*.tsv`, `benchmark/pgc_common_haploreg_archive_coverage.tsv`, and `benchmark/ld_top_hit_benchmark_summary.tsv`
- Power and nominated variants: `benchmark/top_variant_detectable_effects.tsv` and `benchmark/nominated_variant_audit.tsv`
- Reviewer-facing documents: `../Resubmission_Package/`
