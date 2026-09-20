#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Temp qw(tempdir);
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw($GunzipError);
my $dir=tempdir('multigwas sort XXXXX',TMPDIR=>1,CLEANUP=>1);
my $input="CHR\tBP\tSNP\n23\t200\trs3\n1\t20\trs2\n1\t10\trs1\n1\tNA\trs4\n";
gzip(\$input=>"$dir/input.gz") or die $GzipError;
local $ENV{INPUT_GZ}="$dir/input.gz";
local $ENV{OUTPUT_GZ}="$dir/sorted.gz";
local $ENV{EXCLUDED_GZ}="$dir/excluded.gz";
local $ENV{TMPDIR_SORT}="$dir/tmp";
system('bash',"$Bin/../DiffGWASDeps/sort_long_gwas_by_coord.sh")==0 or die "Sorter failed\n";
my $fh=IO::Uncompress::Gunzip->new("$dir/sorted.gz",MultiStream=>1) or die $GunzipError;
my @lines=<$fh>;
$fh->close or die "Cannot close sorted gzip\n";
die "Wrong sorted content\n" unless join('',@lines) eq "#CHR\tBP\tSNP\n1\t10\trs1\n1\t20\trs2\n23\t200\trs3\n";
my $indexed="$dir/sorted.gz";
for my $region ('1:1-30','23:1-300') {
 open my $tab,'-|','tabix',$indexed,$region or die $!;
 my @rows=<$tab>;
 close $tab or die "tabix query failed\n";
 die "Empty region $region\n" unless @rows;
 die "Wrong query rows for $region\n" unless @rows==($region=~/^1:/?2:1);
}
print "PASS: coordinate sorting, excluded row, real tabix index, autosome/X queries, paths containing spaces\n";
my $diff="CHR\tBP\tSNP\tDIFF_Z\tDIFF_P\n1\t10\trs1\t1\t0.3\n1\t20\trs2\t-1\t0.3\n23\t200\trs3\t2\t0.05\n";
gzip(\$diff=>"$dir/diff.gz") or die $GzipError;
system($^X,"$Bin/../DiffGWASDeps/standardize_diff_gwas_zscore.pl",
 '--input',"$dir/diff.gz",'--output',"$dir/std.gz",'--manifest',"$dir/std.manifest.tsv")==0
 or die "Standardization failed\n";
open my $mf,'<',"$dir/std.manifest.tsv" or die $!;
my $metrics=do {local $/;<$mf>};close $mf;
die "Standardized output was not indexed\n" unless $metrics=~/^index_status\tcreated$/m;
print "PASS: standardized differential output has a tabix index\n";

# Force the target beyond the first BGZF member. The locus extractor must scan
# every member to discover its coordinates before issuing the indexed query.
my $bgzf_diff="CHR\tBP\tSNP\tPAIR_TAG\tDIFF_Z\tDIFF_P\n";
for my $i (1..12000) {
 $bgzf_diff .= "1\t$i\trs$i\tPAIR\t0.1\t0.9\n";
}
$bgzf_diff .= "6\t28359632\trsTargetAfterFirstBlock\tPAIR\t2\t0.01\n";
gzip(\$bgzf_diff=>"$dir/bgzf-diff.gz") or die $GzipError;
system($^X,"$Bin/../DiffGWASDeps/standardize_diff_gwas_zscore.pl",
 '--input',"$dir/bgzf-diff.gz",'--output',"$dir/bgzf-std.gz",'--manifest',"$dir/bgzf-std.manifest.tsv")==0
 or die "BGZF standardization failed\n";
system($^X,"$Bin/../DiffGWASDeps/extract_single_snp_wide_diff_gwas.pl",
 '--input',"$dir/bgzf-std.gz",'--target-snp','rsTargetAfterFirstBlock','--window-bp','100',
 '--output',"$dir/bgzf-target.gz",'--manifest',"$dir/bgzf-target.manifest.tsv",
 '--output-dir',$dir,'--base-cols','CHR,BP,SNP','--pair-col','PAIR_TAG',
 '--value-fields','DIFF_P','--pair-map','PAIR=PAIR','--prefix-order','PAIR')==0
 or die "BGZF target lookup failed\n";
open my $bgzf_manifest,'<',"$dir/bgzf-target.manifest.tsv" or die $!;
my $bgzf_metrics=do {local $/;<$bgzf_manifest>};close $bgzf_manifest;
die "Target after the first BGZF member was not preserved\n"
 unless $bgzf_metrics=~/^target_row_found_in_window\t1$/m;
print "PASS: target lookup scans all BGZF members before the tabix query\n";
system($^X,"$Bin/../DiffGWASDeps/extract_significant_diff_gwas.pl",
 '--input',"$dir/bgzf-std.gz",'--output',"$dir/bgzf-wide.gz",
 '--manifest',"$dir/bgzf-wide.manifest.tsv",'--threshold','0.05',
 '--base-cols','CHR,BP,SNP','--pair-col','PAIR_TAG','--value-fields','DIFF_P',
 '--filter-fields','DIFF_P','--pair-map','PAIR=PAIR','--prefix-order','PAIR')==0
 or die "BGZF wide-subset extraction failed\n";
open my $wide_manifest,'<',"$dir/bgzf-wide.manifest.tsv" or die $!;
my $wide_metrics=do {local $/;<$wide_manifest>};close $wide_manifest;
die "Wide-subset extraction stopped before the last BGZF member\n"
 unless $wide_metrics=~/^rows_read\t12001$/m && $wide_metrics=~/^rows_written\t1$/m;
print "PASS: wide-subset extraction scans all BGZF members\n";

my $gtf=<<'GTF';
##gtf-version 3
chr1	test	gene	100	200	.	+	.	gene_id "ENSG000001"; gene_name "GENE1"; gene_type "protein_coding";
chr2	test	gene	100	200	.	+	.	gene_id "ENSG000002"; gene_name "GENE2"; gene_type "protein_coding";
GTF
gzip(\$gtf=>"$dir/fixture.gtf.gz") or die $GzipError;
system($^X,"$Bin/../DiffGWASDeps/extract_gencode_gtf_subset.pl",
 '--gtf-gz',"$dir/fixture.gtf.gz",'--reference-build','hg38',
 '--region','1:50:250','--output',"$dir/gtf-subset.tsv")==0
 or die "Indexed GTF extraction failed\n";
open my $gtf_subset,'<',"$dir/gtf-subset.tsv" or die $!;
my $gtf_text=do {local $/;<$gtf_subset>};close $gtf_subset;
die "Expected GTF gene was not extracted\n" unless $gtf_text=~/\bGENE1\b/;
die "Out-of-region GTF gene was extracted\n" if $gtf_text=~/\bGENE2\b/;
die "GTF tabix index was not created\n" unless -s "$dir/fixture.gtf.bgz.tbi";
print "PASS: native bgzip/tabix build and query the indexed GTF with POSIX paths\n";
