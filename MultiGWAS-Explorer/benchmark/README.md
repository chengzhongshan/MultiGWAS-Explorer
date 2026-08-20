# GWAMA and EasyStrata comparison benchmark

## Reproduction driver

All commands used for the comparator benchmark and the other measured sections
of the revision report are collected in
`run_differential_gwas_revision_benchmarks.sh`. The script requires an explicit
target and writes reruns to `benchmark/reproduction/` by default, preserving the
submission artifacts in this directory.

```bash
# Inspect the resolved GWAS inputs before doing any work.
bash benchmark/run_differential_gwas_revision_benchmarks.sh show-inputs

# Focused reviewer-requested comparison with official GWAMA and EasyStrata.
bash benchmark/run_differential_gwas_revision_benchmarks.sh comparators

# Complete local/report stages (full scans and archive download are expensive).
bash benchmark/run_differential_gwas_revision_benchmarks.sh local-report

# Add the successful SAS ODA HaploReg/LD stages.
bash benchmark/run_differential_gwas_revision_benchmarks.sh all
```

The available targets separately expose raw-data preparation, the unfiltered
QC source, power calculations, MCP-interface trials, downloadable-archive LD
clumping, SAS ODA LD checks, the hybrid common-hit job, report rendering, and
compact result validation. Run the script with `--help` for inputs and
environment overrides.
The documented 3,458-second web-only noncompletion is an explicit
`web-only-control` target and is deliberately excluded from `all`.

This benchmark separates three questions that should not be conflated:

1. **Numerical agreement.** MultiGWAS-Explorer's beta-difference statistic is
   compared with EasyStrata `CALCPDIFF` option 1. With two independent strata,
   its squared Z statistic is also algebraically equal to the two-stratum
   Cochran Q statistic reported by GWAMA.
2. **Input/QC behavior.** Official GWAMA and EasyStrata runs must use the same
   allele-harmonized fixture and report retained/excluded markers separately.
3. **Operational scope.** Runtime, memory, plots, local gene tracks, manifests,
   CLI automation, and agent execution are compared as workflow features; they
   are not evidence that one statistical test is novel.

Run the in-repository numerical audit:

```bash
perl DiffGWASDeps/benchmark_diff_methods.pl \
  --input /path/to/unfiltered.wide.tsv.gz \
  --prefix ALL --rho 0 \
  --output benchmark/multigwas_numerical.tsv
```

The official comparator runs were completed in the pinned
`multigwas-comparators:20260819` Ubuntu 24.04 image. The image contains GWAMA
2.2.2 from Ubuntu and EasyStrata 8.6 from the official Regensburg distribution.
All three paths used the same 100,000 matched variants exported from the
unfiltered merged-long source. Commands, versions, elapsed time, and peak RSS
are recorded in `official_tool_runs.tsv`.

Results are in `fixture_all_100k/official_comparison.tsv`:

- Both official tools retained and matched all 100,000 variants.
- GWAMA `q_p-value` had maximum absolute P difference `5.83e-7`; this is
  expected from its six-decimal output formatting. There were zero threshold
  disagreements at 0.05, 1e-5, or 5e-8.
- EasyStrata `CALCPDIFF` had maximum absolute P difference `8.30e-8` and
  maximum absolute `-log10(P)` difference `4.07e-8`, with zero threshold
  disagreements.
- The fixture contained 86,473 non-ambiguous and 13,527 strand-ambiguous
  direct matches. Both tools retained every row and had zero threshold
  disagreements within either class. The ambiguous class is a formula
  sensitivity analysis only and is excluded by the revised biological default.
- Wall time / peak RSS were 2.73 s / 17,484 KB for the direct workflow audit,
  91.55 s / 79,624 KB for GWAMA, and 4.78 s / 97,520 KB for EasyStrata. These
  are environment-specific operational measurements, not a general ranking.

Acceptance criteria applied:

- absolute beta and SE differences below `1e-10` on unrounded fixtures;
- P-value agreement assessed on `-log10(P)` as well as absolute scale;
- variant counts reconciled before numerical comparison;
- ambiguous A/T and C/G formula results analyzed separately, with the
  full-source biological audit reported under `harmonization_audit/`;
- official-tool runtime and peak RSS measured with `/usr/bin/time -v`;
- chromosome X results reported separately from autosomes.

## Full-source revision results

The category audit under `harmonization_audit/` scanned 44,930,352 long rows.
It found 22,292,345 direct paired rows, 3,386,012 strand-ambiguous pairs, and
18,906,333 rows eligible under the revised default. No swapped, complemented,
incompatible, duplicate-tag, or EAF-difference-greater-than-0.20 pair occurred.

The final `unfiltered_qc_harmonized/` scan produced difference lambda GC values
of 1.0200 ALL, 1.0211 EUR, and 0.9974 ASN. Eligible row counts at raw
`DIFF_P < 1e-5` were 41, 39, and 59, respectively, and none reached 5e-8.
The older `unfiltered_qc/` output is a documented pre-exclusion sensitivity
analysis.

The requested HaploReg SAS runs, power calculations, nominated-variant audit,
and MCP-interface results are documented in `HAPLOREG_LD_SUMMARY.md`,
`POWER_ANALYSIS_SUMMARY.md`, `nominated_variant_audit.tsv`, and
`agent_interface/README.md`. The interface timing used deterministic JSON-RPC
tool calls, not free-form questions submitted to an LLM. Exact expanded
requests are in `agent_interface/AI_Interface_Evaluation_Requests.jsonl`.

## LD-clumped top-hit selection

The revised SAS top-hit stage ranks MAF-passing candidates by their applicable
raw P value and performs greedy HaploReg clumping. Its default excludes a
candidate when its LD with an already selected lead is r-squared >=0.10 in any
configured ancestry panel. Population panels, the r-squared cutoff, ANY/ALL
logic, and query-failure behavior are configuration fields. The former
base-pair window is available only as an explicit legacy method or a labeled
query-failure fallback. Audit tables distinguish `SELECTED_LEAD`, `PRUNED_LD`,
and fallback decisions; this operational LD pruning does not establish
conditional association.

Regression and real-data artifacts use the `reviewer_fixture_*` and
`reviewer_pgc_*_ld_*` filename prefixes. The differential real-data benchmark
reduced three MAF-passing candidates at raw P<1e-6 to two leads; rs7755143 was
pruned under rs753634 at ASN r-squared=0.82.

For the 11,701-candidate common-association analysis, the official HaploReg v4
EUR/ASN downloads were streamed once and reduced to a 7,023,934-edge local
cache at r-squared>=0.20. Extraction took 402 s, SQLite indexing plus the first
full local clump took 118.832 s, and indexed reruns took 1.742-1.751 s. The SAS
ODA cache-plus-web-fallback validation retained 206 leads and completed in
1,126 s after upload/bootstrap; a web-only attempt did not finish in its
3,458-s monitored execution. See `LD_TOP_HIT_BENCHMARK.md` and
`pgc_common_haploreg_archive_coverage.tsv` for scope and coverage.
