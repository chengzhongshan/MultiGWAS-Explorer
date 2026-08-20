#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename qw(dirname);
use File::Spec;
use File::Temp qw(tempdir);
use IO::Compress::Gzip qw(gzip $GzipError);

my $tmp = tempdir(CLEANUP => 1);
my $input = File::Spec->catfile($tmp, 'synthetic.tsv.gz');
my $outdir = File::Spec->catdir($tmp, 'out');
my @header = qw(
    CHR BP A1 A2 SNP PAIR_TAG GROUP1_GWAS_TAG GROUP2_GWAS_TAG
    GROUP1_Z GROUP2_Z DIFF_Z DIFF_P IS_CHRX
    GROUP1_FRQ_A GROUP1_FRQ_U GROUP2_FRQ_A GROUP2_FRQ_U
);
my @rows = (
    [1, 100, 'A', 'G', 'rs_direct', 'TEST', 'TEST_FEMALE', 'TEST_MALE', 1, 0, 1, 0.000001, 0, 0.10, 0.11, 0.12, 0.13],
    [1, 200, 'A', 'T', 'rs_ambiguous', 'TEST', 'TEST_FEMALE', 'TEST_MALE', 2, 0, 2, 0.0001, 0, 0.10, 0.10, 0.10, 0.10],
    [1, 300, 'A', 'C', 'rs_frequency', 'TEST', 'TEST_FEMALE', 'TEST_MALE', 3, 0, 3, 0.0001, 0, 0.10, 0.10, 0.50, 0.50],
    ['X', 400, 'C', 'T', 'rs_x', 'TEST', 'TEST_FEMALE', 'TEST_MALE', 1, 0, 1, 0.000002, 1, 0.20, 0.20, 0.21, 0.21],
);
my $content = join("\t", @header) . "\n";
$content .= join("\t", @$_) . "\n" for @rows;
gzip \$content => $input or die "Cannot create synthetic gzip: $GzipError\n";

my $script = File::Spec->catfile(dirname(__FILE__), 'reviewer_qc_long_report.pl');
system($^X, $script,
    '--input', $input,
    '--outdir', $outdir,
    '--candidate-p', '1e-5',
    '--exclude-strand-ambiguous',
    '--max-eaf-abs-diff', '0.2',
) == 0 or die "reviewer_qc_long_report.pl failed\n";

open my $qfh, '<', File::Spec->catfile($outdir, 'unfiltered_genomic_inflation.tsv')
    or die "Missing QC result\n";
my $qheader = <$qfh>;
chomp $qheader;
my @qnames = split /\t/, $qheader, -1;
my %qidx;
@qidx{@qnames} = (0 .. $#qnames);
my $line = <$qfh> // die "Missing QC data row\n";
chomp $line;
my @q = split /\t/, $line, -1;
die "Wrong source count\n" unless $q[$qidx{SOURCE_PAIR_ROWS}] == 4;
die "Wrong eligible count\n" unless $q[$qidx{ELIGIBLE_VARIANT_ROWS}] == 2;
die "Wrong autosomal count\n" unless $q[$qidx{AUTOSOMAL_VARIANT_ROWS}] == 1;
die "Wrong chrX count\n" unless $q[$qidx{X_VARIANT_ROWS}] == 1;
die "Ambiguous exclusion failed\n" unless $q[$qidx{EXCLUDED_STRAND_AMBIGUOUS}] == 1;
die "Frequency exclusion failed\n" unless $q[$qidx{EXCLUDED_EAF_DIFF_GT_MAX}] == 1;
close $qfh;

open my $cfh, '<', File::Spec->catfile($outdir, 'unfiltered_difference_candidates_p_lt_1e_5.tsv')
    or die "Missing candidate result\n";
my @candidate_lines = <$cfh>;
close $cfh;
die "Excluded rows leaked into candidate output\n" unless @candidate_lines == 3;
die "Expected direct candidate missing\n" unless grep { /rs_direct/ } @candidate_lines;
die "Expected chrX candidate missing\n" unless grep { /rs_x/ } @candidate_lines;

print "reviewer_qc_long_report harmonization test: PASS\n";
