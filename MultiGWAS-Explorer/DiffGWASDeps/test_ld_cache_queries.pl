#!/usr/bin/env perl
use strict;
use warnings;

use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);

my $tmp = tempdir('multigwas_ld_queries_XXXXXX', TMPDIR => 1, CLEANUP => 1);
my $gwas = File::Spec->catfile($tmp, 'gwas.tsv');
my $cache = File::Spec->catfile($tmp, 'ld.tsv');
my $out = File::Spec->catfile($tmp, 'augmented.tsv');
open my $gf, '>', $gwas or die $!;
print {$gf} "CHR\tBP\tSNP\tP\n1\t100\trsA\t0.1\n1\t110\trsProxy\t0.2\n1\t120\trsB\t0.3\n";
close $gf;
open my $cf, '>', $cache or die $!;
print {$cf} "query_snp\tproxy_snp\tld_population\tproxy_r2\tsource\treference_panel\treference_build\tld_method\n";
print {$cf} "rsA\trsProxy\tEUR\t0.25\tPLINK2_1KG_DIRECT\t1000_GENOMES_PHASE_3\tGRCh37_hg19\tPLINK2_R2_PHASED\n";
print {$cf} "rsB\trsProxy\tEUR\t0.90\tPLINK2_1KG_DIRECT\t1000_GENOMES_PHASE_3\tGRCh37_hg19\tPLINK2_R2_PHASED\n";
close $cf;

my $augment = File::Spec->catfile($Bin, 'augment_gwas_with_ld_r2.pl');
system($^X, $augment, '--input', $gwas, '--ld-cache', $cache,
    '--reference-snp', 'rsA', '--output', $out) == 0 or die "LD augmentation failed\n";
open my $fh, '<', $out or die $!;
my $text = do { local $/; <$fh> };
close $fh;
die "The selected query did not control the proxy r2 value\n"
    unless $text =~ /^1\t110\trsProxy\t0\.2\t0\.25$/m;
die "The selected query SNP was not assigned r2=1\n"
    unless $text =~ /^1\t100\trsA\t0\.1\t1$/m;
die "An unrelated query SNP was incorrectly assigned r2=1\n"
    if $text =~ /^1\t120\trsB\t0\.3\t1$/m;

print "Multi-query LD cache filtering: PASS\n";
