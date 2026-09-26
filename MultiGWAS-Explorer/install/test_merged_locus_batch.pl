#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Temp qw(tempdir);
use JSON::PP qw(encode_json);
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Test::More;

my $dir = tempdir('merged locus XXXXX', TMPDIR => 1, CLEANUP => 1);
my $input = "$dir/wide.tsv.gz";
my $table = join("\n",
    "CHR\tBP\tSNP\tP\tBETA",
    "chr1\t100\trsA\t0.8\t0.1",
    "1\t900\trsFar\t0.001\t0.3",
    "1\t120\trsNearby\t0.9\t0.2",
    "2\t200\trsB\t0.7\t0.4",
) . "\n";
gzip \$table => $input or die $GzipError;
my $status = system($^X,
    "$Bin/../DiffGWASDeps/gnuplot/extract_merged_locus_wide_batch.pl",
    '--input', $input, '--output-dir', $dir, '--window-bp', '50',
    '--target', encode_json({snp => 'rsA', chr => '1', bp => 100}),
    '--target', encode_json({snp => 'rsB', chr => '2', bp => 200}));
is($status, 0, 'multiple merged-wide loci are extracted in one pass');
my $a = '';
gunzip "$dir/gunplot_locus_rsA_window_50.wide.tsv.gz" => \$a or die $GunzipError;
is($a, "CHR\tBP\tSNP\tP\tBETA\nchr1\t100\trsA\t0.8\t0.1\n1\t120\trsNearby\t0.9\t0.2\n",
    'local window keeps all rows and columns regardless of P');
my $b = '';
gunzip "$dir/gunplot_locus_rsB_window_50.wide.tsv.gz" => \$b or die $GunzipError;
is($b, "CHR\tBP\tSNP\tP\tBETA\n2\t200\trsB\t0.7\t0.4\n",
    'other chromosome receives only its own locus');
open my $mf, '<', "$dir/gunplot_locus_rsA_window_50.wide.manifest.tsv" or die $!;
my %manifest = map { chomp; split /\t/, $_, 2 } <$mf>;
close $mf;
is($manifest{target_snp}, 'rsA', 'manifest identifies the target SNP');
is($manifest{window_bp}, 50, 'manifest records the requested window');

my $indexed = "$dir/wide.bgz";
is(system($^X, "$Bin/../DiffGWASDeps/gnuplot/index_merged_wide_tabix.pl",
    '--input', $input, '--output', $indexed), 0,
    'merged-wide table is sorted and indexed with bgzip and tabix');
ok(-s "$indexed.tbi", 'tabix index was created');
is(system($^X, "$Bin/../DiffGWASDeps/gnuplot/index_merged_wide_tabix.pl",
    '--input', $input, '--output', $indexed), 0,
    'the matching tabix index can be reused');
is(system($^X,
    "$Bin/../DiffGWASDeps/gnuplot/extract_merged_locus_wide_batch.pl",
    '--input', $input, '--indexed-input', $indexed,
    '--output-dir', $dir, '--window-bp', '50',
    '--target', encode_json({snp => 'rsA', chr => '1', bp => 100}),
    '--target', encode_json({snp => 'rsB', chr => '2', bp => 200})), 0,
    'local windows are fetched through tabix');
my $indexed_a = '';
gunzip "$dir/gunplot_locus_rsA_window_50.wide.tsv.gz" => \$indexed_a or die $GunzipError;
my @indexed_lines = split /\n/, $indexed_a;
my @streamed_lines = split /\n/, $a;
my @indexed_rows = sort @indexed_lines[1, 2];
my @streamed_rows = sort @streamed_lines[1, 2];
is_deeply(\@indexed_rows, \@streamed_rows,
    'tabix extraction retains aliases and nonsignificant SNPs');

SKIP: {
    skip 'gnuplot is unavailable', 3 if system('gnuplot', '--version') != 0;
    my $prefix = "$dir/local_rsA";
    my $plot_status = system($^X,
        "$Bin/../DiffGWASDeps/gnuplot/pdl_gunplot_local_locus.pl",
        '--data', "$dir/gunplot_locus_rsA_window_50.wide.tsv.gz",
        '--snp', 'rsA', '--out-prefix', $prefix,
        '--window-bp', '50', '--pcols', 'P',
        '--width', '600', '--height', '400');
    is($plot_status, 0, 'local Manhattan renders from the batch-extracted window');
    ok(-s "$prefix.png", 'local Manhattan PNG is nonempty');
    open my $pmf, '<', "$prefix.manifest.tsv" or die $!;
    my %plot_manifest = map { chomp; split /\t/, $_, 2 } <$pmf>;
    close $pmf;
    is($plot_manifest{rows_in_window}, 2,
        'local Manhattan includes both nonsignificant SNPs');
}
done_testing();
