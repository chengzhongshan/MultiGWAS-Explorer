#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Temp qw(tempdir);
use IO::Compress::Gzip qw(gzip $GzipError);
use Test::More;
use lib "$Bin/../DiffGWASDeps";
use ChromosomeBounds qw(chromosome_length locus_window);

is(chromosome_length('hg38', 'chr1'), 248956422, 'hg38 chromosome 1 length');
is(chromosome_length('GRCh37', '23'), 155270560, 'GRCh37 X alias');
is(chromosome_length('t2t', 'chr1'), 248387328, 'T2T/hs1 chromosome 1 length');
is_deeply([locus_window('hg38', '1', 246321905, 10000000)],
    [236321905, 248956422], 'right edge stops at chromosome end');
is_deeply([locus_window('hg38', '1', 120000000, 10000000)],
    [110000000, 130000000], 'interior window retains full span');
is_deeply([locus_window('hg19', '1', 246321905, 10000000)],
    [236321905, 249250621], 'hg19 uses its own chromosome end');

my $dir = tempdir('chromosome bounds XXXXX', TMPDIR => 1, CLEANUP => 1);
my $data = "$dir/wide.tsv.gz";
my $table = "CHR\tBP\tSNP\tP\n"
    . "1\t246321905\trsEnd\t1e-7\n"
    . "1\t248000000\trsNearby\t0.02\n"
    . "1\t120000000\trsInterior\t0.04\n";
gzip \$table => $data or die $GzipError;

my $extract = "$Bin/../DiffGWASDeps/gnuplot/extract_merged_locus_wide_batch.pl";
my $target = '{"snp":"rsEnd","chr":"1","bp":246321905}';
is(system($^X, $extract, '--input', $data, '--output-dir', $dir,
    '--window-bp', '10000000', '--reference-build', 'hg38',
    '--target', $target), 0, 'locus extraction succeeds at chromosome end');
my $locus_manifest = "$dir/gnuplot_locus_rsEnd_window_10000000.wide.manifest.tsv";
my $locus_metrics = metrics($locus_manifest);
is($locus_metrics->{window_end}, 248956422, 'tabix extraction window is clipped');

SKIP: {
    skip 'gnuplot is unavailable', 6 if system('gnuplot', '--version') != 0;
    my $render = "$Bin/../DiffGWASDeps/gnuplot/pdl_gnuplot_local_locus.pl";
    my $prefix = "$dir/end_plot";
    is(system($^X, $render, '--data', $data, '--snp', 'rsEnd',
        '--out-prefix', $prefix, '--window-bp', '10000000',
        '--reference-build', 'hg38', '--pcols', 'P'), 0,
        'gnuplot renders clipped chromosome-end window');
    my $plot_metrics = metrics("$prefix.manifest.tsv");
    is($plot_metrics->{window_start}, 236321905, 'plot records window start');
    is($plot_metrics->{window_end}, 248956422, 'plot records clipped end');
    open my $gp, '<', "$prefix.gp" or die $!;
    my $script = do { local $/; <$gp> };
    close $gp;
    like($script, qr/set xrange \[236321905:248956422\]/,
        'gnuplot x-axis stops at chromosome end');
    is(system($^X, $render, '--data', $data, '--snp', 'rsInterior',
        '--out-prefix', "$dir/interior_plot", '--window-bp', '10000000',
        '--reference-build', 'hg38', '--pcols', 'P'), 0,
        'gnuplot renders interior window');
    is(metrics("$dir/interior_plot.manifest.tsv")->{window_end}, 130000000,
        'interior plot still includes any real genomic gaps');
}

done_testing();

sub metrics {
    my ($path) = @_;
    open my $fh, '<', $path or die "Cannot read $path: $!";
    my %values = map { chomp; split /\t/, $_, 2 } <$fh>;
    close $fh;
    return \%values;
}
