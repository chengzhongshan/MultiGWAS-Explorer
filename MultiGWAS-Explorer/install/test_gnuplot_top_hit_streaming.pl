#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Spec;
use File::Temp qw(tempdir);
use IO::Compress::Gzip qw(gzip $GzipError);
use Test::More;

my $tmp = tempdir(CLEANUP => 1);
my $input = File::Spec->catfile($tmp, 'wide.tsv.gz');
my $source = join("\n",
    "CHR\tBP\tSNP\tDIFF_P\tOTHER_P",
    "1\t100\trsA\t1e-7\t2e-6",
    "1\t2000100\trsB\t2e-6\t3e-6",
    "2\t500\trsC\t0.8\t0.9",
) . "\n";
gzip \$source => $input or die "Cannot make fixture: $GzipError\n";
my $selector = File::Spec->catfile($Bin, '..', 'DiffGWASDeps', 'gnuplot',
    'select_top_hits_from_wide.pl');

sub selected_snps {
    my ($output) = @_;
    open my $fh, '<', $output or die "Cannot read $output: $!\n";
    <$fh>;
    my @snps = map { (split /\t/)[3] } <$fh>;
    close $fh;
    return \@snps;
}

my $primary = File::Spec->catfile($tmp, 'primary.tsv');
is(system($^X, $selector, '--input', $input, '--output', $primary,
    '--focus-pvar', 'DIFF_P', '--thresholds', '1e-6,1e-5',
    '--maf-threshold', '0'), 0, 'top-hit selector runs');
is_deeply(selected_snps($primary), ['rsA'],
    'first threshold with a candidate controls selection');

my $fallback = File::Spec->catfile($tmp, 'fallback.tsv');
is(system($^X, $selector, '--input', $input, '--output', $fallback,
    '--focus-pvar', 'OTHER_P', '--thresholds', '1e-6,1e-5',
    '--maf-threshold', '0'), 0, 'threshold fallback runs');
is_deeply(selected_snps($fallback), ['rsA', 'rsB'],
    'fallback threshold retains both distant eligible loci');

my $targeted = File::Spec->catfile($tmp, 'targeted.tsv');
is(system($^X, $selector, '--input', $input, '--output', $targeted,
    '--target-snps', 'rsC', '--maf-threshold', '0'), 0,
    'targeted selection runs');
is_deeply(selected_snps($targeted), ['rsC'],
    'targeted selection retains a nonsignificant requested SNP');

done_testing();
