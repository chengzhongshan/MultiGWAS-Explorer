# Public sex-stratified GWAS validation (September 2026)

## Follow-up correction: forest text and output discoverability

User inspection found that the original gnuplot forest PNGs decoded correctly
but had invisible text in portable Cygwin. Hard-coded Arial selection was the
cause; the renderer now uses generic Sans, includes a Variant/Cohort y-axis
title, and refreshes cached panels after renderer changes or missing panels.
The real-data forest panels were regenerated and visually inspected. A Perl
pixel-region regression checks visible label text and passed on Windows and
Ubuntu with the pipeline-installed gnuplot. The earlier PNG-decoding pass
alone did not establish forest-label correctness.

Seven SAS PNGs existed in the pipeline working directory. A new results
gallery collects both backends under the selected output directory and
verifies the copied images using SHA-256. The public SAS test also disables
gnuplot fallback and forces fresh rendering to prevent ambiguous test results.

Test source: [PGC schizophrenia 2022](https://figshare.com/articles/dataset/scz2022/19426775) European female and male summary statistics,
including autosomes and X (GRCh37/hg19). All four published MD5 checksums passed.
The new Perl downloader also downloaded and verified the male X file from the
public Figshare endpoint.

## Scope and evidence

| Check | Result |
| --- | --- |
| Full source records | 15,758,068 |
| Full differential comparisons | 6,650,636 |
| Full independent numerical validation | PASS: all 6,650,636 comparisons, including 194,288 on X |
| Full filtered plotting table | 1,251,499 rows; no invalid filter values or missing/duplicate pair prefixes |
| Full gnuplot | PASS: all four plot families; 11 decoded PNGs |
| Strand-ambiguous groups excluded | 1,187,692 |
| Groups without a pair excluded | 81,412 |
| Real integration subset numerical checks | PASS: all 69,926 comparisons, chromosomes 1–22 and X |
| Subset gnuplot | PASS: genome-wide Manhattan, local Manhattan, GTF tracks, forest plots |
| Subset SAS ODA | PASS: same four plot families; PNG decoding verified |
| Real LD service and overlay | PASS: HaploReg lookup; 147 LD points plotted |
| Windows and Ubuntu path/index regressions | PASS: spaces in paths, autosomal/X tabix queries, nested precomputed configuration |
| Full SAS ODA | PASS: all four plot families; 7 decoded PNGs |

The integration subset retained source rows with BP modulo 100 equal to zero
or source P below 1e-5. It is a real-data subset, distinct from the unthinned
full run. Numerical validation independently computes female-minus-male effects,
SE, Z, and two-sided P using POSIX::erfc, assuming zero cross-sex covariance.

Local execution used Windows 10, portable Cygwin 3.6.10, Perl 5.44 and gnuplot
6.0.5. Ubuntu regression tests ran in WSL Ubuntu 24.04. GitHub Actions separately
tests fresh installation and synthetic rendering across its platform matrix;
those checks do not constitute a full public-GWAS run on every platform.

All six installation jobs passed for commit `f8780f7`: Windows, Ubuntu 24.04,
macOS 15 Apple Silicon, macOS 15 Intel, Docker, and Apptainer.
[GitHub Actions run](https://github.com/chengzhongshan/MultiGWAS-Explorer/actions/runs/35095368132).

The subsequent code revision `a551a5a` changes only inquiry/MAF log messages
and documentation. Its [repeat installation run](https://github.com/chengzhongshan/MultiGWAS-Explorer/actions/runs/35128834251)
had passed five jobs when this report was published; Intel macOS was still
running. The completed six-platform result above belongs to `f8780f7`.

Full-data inquiry targets were rs2232429 (common association), rs185665940
(minimum autosomal differential P), and rs62604261 (minimum X differential P).
The full numerical scan took 401 seconds; gnuplot genome-wide rendering took
142 seconds on this host. These are observed timings, not performance promises.

## Bugs corrected during testing

- Convert Cygwin data paths for native Windows tabix; require genuine indexes.
- Use platform-specific Perl dependencies and install Perl HTTPS dependencies.
- Correct gnuplot forest freshness-helper lookup.
- Create nested configuration directories and apply precomputed input paths
  before generating runner configuration.
- Pass the selected GENCODE cache to the SAS runner.
- Select the repository's SAS helper instead of an unrelated version on PATH.
- Replace macOS-incompatible `zcat` calls in standardization with list-form
  `gzip -dc` calls. The new CI regression caught this platform-specific failure.

The first full standardized-index attempt failed before the path correction;
the existing standardized table was subsequently indexed successfully without
recomputing its values. The original failed-run manifest is retained locally.

Inquiry plots use selected common, differential and X targets. These checks
do not establish biological significance or test genome-wide automatic LD
clumping. Test scripts are Perl; the SAS backend retains its SASPy dependency.

Explicit inquiry targets bypass automatic MAF filtering. Separate synthetic
MAF and allele-harmonization regressions passed; the inquiry plots are not
evidence of automatic MAF-filter execution.

SAS gene-track logs reported 346, 110 and 141 unique exons for the three loci.
Representative genome-wide and chromosome-X images were also inspected
visually. SAS verified deletion of temporary batch/forest artifacts and the
shared full plotting table after download.

The full run was resumed after an interruption during table extraction and
after a running shell was interrupted by a log-message edit. Completed data
stages were reused. The final resumed run and image validation passed; the
earlier failed logs remain in the local test directory. This is a completed
stage-by-stage validation, not a claim that the first invocation was error-free.

See [the Perl test instructions](PUBLIC_GWAS_TEST.md) to reproduce the test and
[the full numerical result](public_gwas_numeric_results.json) for chromosome
counts and target coordinates. Source files, credentials, and generated GWAS
tables are not included in Git.
