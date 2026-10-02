#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Spec;
use File::Temp qw(tempdir);
use Test::More;

my $dir = tempdir('selected_signed_ld_gtf_XXXXXX', TMPDIR => 1, CLEANUP => 1);
my $csv = File::Spec->catfile($dir, 'hits.csv');
my $config = File::Spec->catfile($dir, 'runner.json');
my $runner = File::Spec->catfile($dir, 'fake_single_runner.sh');
my $index = File::Spec->catfile($dir, 'plots.html');
my $calls = File::Spec->catfile($dir, 'calls.tsv');
open my $fh, '>', $csv or die $!;
print {$fh} "SNP,CHR,BP\nrs111,1,100\nrs222,2,200\n";
close $fh;
open $fh, '>', $config or die $!;
print {$fh} "{}\n";
close $fh;
open $fh, '>', $runner or die $!;
print {$fh} <<'SH';
#!/usr/bin/env bash
set -euo pipefail
printf '%s\t%s\t%s\t%s\n' "$TARGET_SNP" "$LOCAL_WINDOW_BP" "$KEEP_REMOTE_PLOT_DATA" "$CLEAN_ODA_INPUT" >> "$TEST_CALLS"
printf '<html>%s</html>\n' "$TARGET_SNP" > "$TEST_DIR/$OUTPUT_HTML_BASENAME"
printf 'PNG:%s\n' "$TARGET_SNP" > "$TEST_DIR/${OUTPUT_HTML_BASENAME%.html}.png"
SH
close $fh;

local $ENV{TEST_DIR} = $dir;
local $ENV{TEST_CALLS} = $calls;
local $ENV{KEEP_REMOTE_PLOT_DATA} = 1;
local $ENV{CLEAN_ODA_INPUT} = 0;
my @cmd = ($^X, File::Spec->catfile($Bin, 'run_selected_signed_ld_gtf.pl'),
    '--targets-csv', $csv, '--runner-config', $config,
    '--output-html', $index, '--single-runner', $runner,
    '--window-bp', '5000000', '--population', 'EUR', '--mode', 'heatmap');
is(system(@cmd), 0, 'both loci render successfully');
open $fh, '<', $calls or die $!;
my @calls = <$fh>;
close $fh;
is(scalar(@calls), 2, 'one runner call per locus');
like($calls[0], qr/^rs111\t5000000\t0\t1\s*$/, 'first locus cleans its remote inputs');
like($calls[1], qr/^rs222\t5000000\t0\t1\s*$/, 'second locus cleans its remote inputs');
ok(-s File::Spec->catfile($dir, 'plots_rs111.html'), 'first final plot exists');
ok(-s File::Spec->catfile($dir, 'plots_rs222.html'), 'second final plot exists');
ok(-s $index, 'plot index exists');
is(system(@cmd), 0, 'repeat request succeeds');
open $fh, '<', $calls or die $!;
@calls = <$fh>;
close $fh;
is(scalar(@calls), 2, 'repeat request reuses both completed loci');
open $fh, '>', $config or die $!;
print {$fh} "{\"FOREST_DOTSIZE\":20}\n";
close $fh;
is(system(@cmd), 0, 'unrelated forest setting is accepted');
open $fh, '<', $calls or die $!;
@calls = <$fh>;
close $fh;
is(scalar(@calls), 2, 'unrelated setting leaves per-locus plots cached');
open $fh, '>', $config or die $!;
print {$fh} "{\"GTF_ASSOC_PVARS\":\"TEST_P\"}\n";
close $fh;
is(system(@cmd), 0, 'changed GTF setting is accepted');
open $fh, '<', $calls or die $!;
@calls = <$fh>;
close $fh;
is(scalar(@calls), 4, 'changed GTF setting rerenders both loci');
done_testing();
