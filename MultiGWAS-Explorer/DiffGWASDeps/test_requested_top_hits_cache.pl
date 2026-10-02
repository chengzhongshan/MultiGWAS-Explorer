#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Spec;
use File::Temp qw(tempdir);
use IO::Compress::Gzip qw(gzip $GzipError);
use Test::More;

my $dir = tempdir('requested_hits_cache_XXXXXX', TMPDIR => 1, CLEANUP => 1);
my $input = File::Spec->catfile($dir, 'wide.tsv.gz');
my $output = File::Spec->catfile($dir, 'requested.csv');
my $wide = "CHR\tBP\tSNP\tA1\tA2\tTEST_P\n1\t100\trs111\tA\tG\t0.001\n";
gzip(\$wide => $input) or die "Cannot create fixture: $GzipError\n";
my $script = File::Spec->catfile($Bin, 'generate_requested_top_hits_csv.pl');

sub run_generator {
    my ($gene) = @_;
    open my $pipe, '-|', $^X, $script,
        '--input', $input, '--output', $output,
        '--target-snps', 'rs111', '--target-snp-genes', "rs111:$gene",
        '--top-hit-focus-pvar', 'TEST_P', '--reuse-cache'
        or die "Cannot run generator: $!\n";
    local $/;
    my $log = <$pipe> // '';
    close $pipe or die "Generator failed\n";
    return $log;
}

my $first = run_generator('GENE1');
ok(-s $output && -s "$output.request.md5", 'CSV and request key were written');
unlike($first, qr/\[skip\]/, 'first request generates the CSV');
my $second = run_generator('GENE1');
like($second, qr/\[skip\] Reusing requested top-hit CSV/, 'identical request skips selection');
my $third = run_generator('GENE2');
unlike($third, qr/\[skip\]/, 'changed gene request invalidates the cache');
open my $fh, '<', $output or die $!;
local $/;
my $csv = <$fh>;
close $fh;
like($csv, qr/GENE2/, 'regenerated output contains the changed gene');
done_testing();
