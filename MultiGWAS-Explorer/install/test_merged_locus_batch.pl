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
    "X\t300\trsX\t0.6\t0.5",
) . "\n";
gzip \$table => $input or die $GzipError;
my $status = system($^X,
    "$Bin/../DiffGWASDeps/gnuplot/extract_merged_locus_wide_batch.pl",
    '--input', $input, '--output-dir', $dir, '--window-bp', '50',
    '--target', encode_json({snp => 'rsA', chr => '1', bp => 100}),
    '--target', encode_json({snp => 'rsB', chr => '2', bp => 200}));
is($status, 0, 'multiple merged-wide loci are extracted in one pass');
my $a = '';
gunzip "$dir/gnuplot_locus_rsA_window_50.wide.tsv.gz" => \$a or die $GunzipError;
is($a, "CHR\tBP\tSNP\tP\tBETA\nchr1\t100\trsA\t0.8\t0.1\n1\t120\trsNearby\t0.9\t0.2\n",
    'local window keeps all rows and columns regardless of P');
my $b = '';
gunzip "$dir/gnuplot_locus_rsB_window_50.wide.tsv.gz" => \$b or die $GunzipError;
is($b, "CHR\tBP\tSNP\tP\tBETA\n2\t200\trsB\t0.7\t0.4\n",
    'other chromosome receives only its own locus');
open my $mf, '<', "$dir/gnuplot_locus_rsA_window_50.wide.manifest.tsv" or die $!;
my %manifest = map { chomp; split /\t/, $_, 2 } <$mf>;
close $mf;
is($manifest{target_snp}, 'rsA', 'manifest identifies the target SNP');
is($manifest{window_bp}, 50, 'manifest records the requested window');
is($manifest{columns}, 'CHR,BP,SNP,P,BETA',
    'manifest records columns needed to reuse a locus without extraction');

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
gunzip "$dir/gnuplot_locus_rsA_window_50.wide.tsv.gz" => \$indexed_a or die $GunzipError;
my @indexed_lines = split /\n/, $indexed_a;
my @streamed_lines = split /\n/, $a;
my @indexed_rows = sort @indexed_lines[1, 2];
my @streamed_rows = sort @streamed_lines[1, 2];
is_deeply(\@indexed_rows, \@streamed_rows,
    'tabix extraction retains aliases and nonsignificant SNPs');

my $target_csv = "$dir/common_top_hits.csv";
open my $target_fh, '>', $target_csv or die $!;
print {$target_fh} "SNP,CHR,BP\nrsA,1,100\nrsNearby,1,120\nrsB,2,200\n";
close $target_fh or die $!;
my $combined = "$dir/common_local_windows.tsv.gz";
is(system($^X,
    "$Bin/../DiffGWASDeps/gnuplot/extract_merged_locus_wide_batch.pl",
    '--input', $input, '--indexed-input', $indexed,
    '--output-dir', $dir, '--window-bp', '50',
    '--targets-csv', $target_csv, '--combined-output', $combined), 0,
    'common-hit CSV drives indexed local-window extraction');
my $combined_text = '';
gunzip $combined => \$combined_text or die $GunzipError;
my @combined_rows = sort grep { length } (split /\n/, $combined_text)[1 .. 3];
is_deeply(\@combined_rows, [sort(
    "chr1\t100\trsA\t0.8\t0.1",
    "1\t120\trsNearby\t0.9\t0.2",
    "2\t200\trsB\t0.7\t0.4",
)], 'combined upload keeps full locus rows and de-duplicates overlapping windows');

my $sas_combined = "$dir/sas_local_windows.tsv.gz";
is(system($^X,
    "$Bin/../DiffGWASDeps/gnuplot/extract_merged_locus_wide_batch.pl",
    '--input', $input, '--indexed-input', $indexed,
    '--output-dir', $dir, '--window-bp', '50',
    '--target', encode_json({snp => 'rsX', chr => '23', bp => 300}),
    '--combined-output', $sas_combined, '--combined-numeric-chr'), 0,
    'SAS combined input is extracted for a chromosome X lead');
my $sas_combined_text = '';
gunzip $sas_combined => \$sas_combined_text or die $GunzipError;
like($sas_combined_text, qr/^23\t300\trsX\t0\.6\t0\.5$/m,
    'SAS combined input encodes X as numeric chromosome 23');
my $x_locus = '';
gunzip "$dir/gnuplot_locus_rsX_window_50.wide.tsv.gz" => \$x_locus or die $GunzipError;
like($x_locus, qr/^X\t300\trsX\t0\.6\t0\.5$/m,
    'gnuplot locus preserves the source chromosome label');

my $sas_locus_dir = "$dir/sas locus";
my $sas_cache_dir = "$dir/sas index";
my @sas_locus_cmd = ($^X,
    "$Bin/../DiffGWASDeps/gnuplot/prepare_indexed_merged_locus.pl",
    '--input', $input, '--target-snp', 'rsA', '--window-bp', '50',
    '--output-dir', $sas_locus_dir, '--cache-dir', $sas_cache_dir);
is(system(@sas_locus_cmd), 0,
    'SAS GTF preparation builds a sorted tabix index and extracts the target window');
my $sas_locus = '';
gunzip "$sas_locus_dir/gnuplot_locus_rsA_window_50.wide.tsv.gz" => \$sas_locus
    or die $GunzipError;
my @sas_locus_lines = split /\n/, $sas_locus;
my @sas_locus_rows = sort @sas_locus_lines[1, 2];
is_deeply(\@sas_locus_rows, \@streamed_rows,
    'SAS GTF tabix locus retains nonsignificant SNPs');
open my $smf, '<', "$sas_locus_dir/gnuplot_locus_rsA_window_50.wide.manifest.tsv" or die $!;
my %sas_manifest = map { chomp; split /\t/, $_, 2 } <$smf>;
close $smf;
is($sas_manifest{access_mode}, 'TABIX', 'SAS GTF locus records indexed access');
is(system(@sas_locus_cmd), 0, 'SAS GTF indexed locus is reusable');

SKIP: {
    skip 'gnuplot is unavailable', 8 if system('gnuplot', '--version') != 0;
    my $prefix = "$dir/local_rsA";
    my $plot_status = system($^X,
        "$Bin/../DiffGWASDeps/gnuplot/pdl_gnuplot_local_locus.pl",
        '--data', "$dir/gnuplot_locus_rsA_window_50.wide.tsv.gz",
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

    my $gtf = "$dir/chrX_genes.tsv";
    open my $gf, '>', $gtf or die $!;
    print {$gf} "chr\tgenesymbol\tgene\tst\ten\ttype\nX\tXGENE\tXGENE\t295\t305\tgene\n";
    close $gf;
    my $x_prefix = "$dir/local_rsX_with_genes";
    is(system($^X, "$Bin/../DiffGWASDeps/gnuplot/pdl_gnuplot_local_locus.pl",
        '--data', "$dir/gnuplot_locus_rsX_window_50.wide.tsv.gz",
        '--snp', 'rsX', '--out-prefix', $x_prefix,
        '--window-bp', '50', '--pcols', 'P', '--gtf', $gtf,
        '--width', '600', '--height', '400'), 0,
        'chrX gene-track Manhattan renders in chromosome-color mode');
    open my $xgp, '<', "$x_prefix.gp" or die $!;
    my $xgp_text = do { local $/; <$xgp> };
    close $xgp;
    like($xgp_text, qr/using 1:2 with points pt 7 ps 0\.72 lc rgb '#d95319'/,
        'chrX local points use the genome-wide chrX orange');
    ok(-s "$x_prefix.png", 'chrX local gene-track PNG is nonempty');

    my $genome_prefix = "$dir/genome_colors";
    is(system($^X, "$Bin/../DiffGWASDeps/gnuplot/pdl_gnuplot_manhattan.pl",
        '--data', $input, '--out-prefix', $genome_prefix,
        '--pcols', 'P', '--min-logp', '0', '--thin-mod', '1',
        '--no-remove-x-chr', '--width', '600', '--height', '400'), 0,
        'genome-wide Manhattan renders chromosomes including X');
    open my $gpt, '<', "$genome_prefix.plot.tsv" or die $!;
    my $gpt_text = do { local $/; <$gpt> };
    close $gpt;
    like($gpt_text, qr/^\S+\t\S+\t\d+\tX\t300\t\S+\t2$/m,
        'genome-wide chrX data uses the same orange palette index');
}
done_testing();
