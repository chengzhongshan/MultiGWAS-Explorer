# Public schizophrenia-by-sex validation (2026-09-26)

The updated pipeline was rerun with the public female and male `ALL` GWAS
summary statistics from the PGC [scz2022 release](https://figshare.com/articles/dataset/scz2022/19426775/5).
The input files were `daner_PGC_SCZ_w3_81_0618a_all_female.gz` and
`daner_PGC_SCZ_w3_81_0618a_all_male.gz`. Both passed gzip integrity checks and
matched the release's published `cksum` values (female `2281346142`, male
`1867671634`). This run covers autosomes; the separate sex-stratified chrX
files were not included.

The local spec used `raw_pgc_vcf_sumstats`, one `ALL_FEMALE` versus `ALL_MALE`
pair, GRCh37/hg19, and a 500 kb half-window for the targeted GTF plot. Inputs,
temporary files, and gnuplot images remain in the ignored local
`cache/pgc_scz_sex_public/` directory. SAS images were downloaded to the
pipeline directory. The source files and genotype reference are not committed.

| Stage | Observed result |
| --- | ---: |
| Female source rows | 7,544,263 |
| Male source rows | 7,580,940 |
| Differential rows standardized | 6,379,921 |
| Standard wide rows (any configured P < 0.05) | 1,236,471 |
| Genome-wide displayed-track P < 0.05 rows | 1,221,905 |
| Full gnuplot Manhattan render | Passed; 1m 6s |
| Full SAS ODA Manhattan render | Passed; 3m 39s; 179,574-byte PNG |

The local chr2 test used rs185665940 at hg19 chr2:72,269,028. The SNP appears
in both public source files. The SAS single-target path read the tabix-indexed
standardized GWAS window, retained nonsignificant local SNPs, and added phased
EUR LD from the official
[PLINK2 1000 Genomes Phase 3 reference](https://www.cog-genomics.org/plink/2.0/resources).
The gnuplot local GTF manifest reports `ld_display_mode=heatmap`,
`signed_r2_coloring=1`, and five plotted LD points. The SAS ODA local GTF run
completed in 3m 46s and downloaded a 246,959-byte PNG with a gene track and a
fixed signed-LD colorbar from -1 to 1.

The explicit `--manhattan-all-snps` escape hatch was checked on a 100,000-row
slice from each public source. After pairing, its separate unfiltered wide
table contained 85,812 coordinate-valid SNPs; the default genome-wide plot
input contained 16,419 SNPs meeting the displayed-track P < 0.05 rule.
Returning to the default option restored the filtered cache. The full 6.38M-row
all-SNP render was not requested because it is much larger and slower than the
default scientific view.

This rerun exposed and fixed four portability/correctness issues: raw sorting
and standardized/locus extraction now discover `bgzip` and `tabix` in an
ancestor `local/bin` directory; explicit all-SNP mode now bypasses the earlier
wide-table P filter as well as gnuplot point thinning; GTF downloads use a
streamed curl transfer when available; and a transient Windows lock on a
temporary SAS file no longer changes a successful preparation into a failure.
