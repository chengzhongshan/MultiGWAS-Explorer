#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Temp qw(tempdir);
use JSON::PP qw(decode_json);

my $dir = tempdir(CLEANUP => 1);
open my $d, '>', "$dir/in.tsv" or die $!;
print {$d} join("\t", qw(CHR BP A1 A2 SNP ALL_GROUP1_BETA ALL_GROUP1_SE ALL_GROUP1_P ALL_GROUP1_FRQ_A ALL_GROUP1_INFO ALL_GROUP2_BETA ALL_GROUP2_SE ALL_GROUP2_P ALL_GROUP2_FRQ_A ALL_GROUP2_INFO ALL_DIFF_BETA ALL_DIFF_SE ALL_DIFF_P ALL_STD_DIFF_Z ALL_STD_DIFF_P)), "\n";
for my $i (1..2000) {
    my $snp = $i == 7 ? 'rsTest' : "rs$i";
    print {$d} join("\t", 1, $i, qw(A G), $snp, .1, .1, .3, ($i == 7 ? .01 : .2), ($i == 7 ? .7 : .95), -.1, .1, .3, .2, .95, .2, .141421, .1573, 1.41421, .1573), "\n";
}
close $d;
open my $c, '>', "$dir/spec.json" or die $!;
print {$c} '{"pairs":[{"prefix":"ALL"}]}'; close $c;
my $script = "$Bin/qc_report.pl";
system($^X, $script, '--input', "$dir/in.tsv", '--config', "$dir/spec.json", '--output-dir', "$dir/out", '--nominated-snps', 'rsTest', '--reservoir-size', 1000) == 0 or die "QC command failed\n";
open my $j, '<', "$dir/out/qc_summary.json" or die $!; local $/; my $x = decode_json(<$j>); close $j;
die "row count mismatch\n" unless $x->{rows_scanned} == 2000;
open my $h, '<', "$dir/out/nominated_variants.tsv" or die $!; my $txt = do { local $/; <$h> }; close $h;
die "expected MAF flag absent in:\n$txt\n" unless $txt =~ /MAF_AT_OR_BELOW_0\.01/;
die "expected INFO flag absent\n" unless $txt =~ /INFO_BELOW_0\.8/;
print "qc_report synthetic test: PASS\n";
