#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Spec;
use File::Temp qw(tempdir);
use JSON::PP qw(decode_json);
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
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$TARGET_SNP" "$LOCAL_WINDOW_BP" "$KEEP_REMOTE_PLOT_DATA" "$CLEAN_ODA_INPUT" "$SINGLE_SNP_REUSE_SHARED_MACROS" "$INCLUDE_PREFLIGHT_ENABLED" "$SAS_ODA_REUSE_VERIFIED_MACRO_BOOTSTRAP_HELPER" >> "$TEST_CALLS"
if [[ "${FAIL_ONCE_SNP:-}" == "$TARGET_SNP" && ! -e "$TEST_DIR/failed.once" ]]; then
  touch "$TEST_DIR/failed.once"
  exit 74
fi
printf '<html>%s</html>\n' "$TARGET_SNP" > "$TEST_DIR/$OUTPUT_HTML_BASENAME"
printf 'PNG:%s\n' "$TARGET_SNP" > "$TEST_DIR/${OUTPUT_HTML_BASENAME%.html}.png"
SH
close $fh;

local $ENV{TEST_DIR} = $dir;
local $ENV{TEST_CALLS} = $calls;
local $ENV{KEEP_REMOTE_PLOT_DATA} = 1;
local $ENV{CLEAN_ODA_INPUT} = 0;
local $ENV{INCLUDE_PREFLIGHT_ENABLED} = 1;
my @cmd = ($^X, File::Spec->catfile($Bin, 'run_selected_signed_ld_gtf.pl'),
    '--targets-csv', $csv, '--runner-config', $config,
    '--output-html', $index, '--single-runner', $runner,
    '--window-bp', '5000000', '--population', 'EUR', '--mode', 'heatmap');
is(system(@cmd), 0, 'both loci render successfully');
open $fh, '<', $calls or die $!;
my @calls = <$fh>;
close $fh;
is(scalar(@calls), 2, 'one runner call per locus');
like($calls[0], qr/^rs111\t5000000\t0\t1\t0\t1\t0\s*$/, 'first locus uploads and checks shared macros');
like($calls[1], qr/^rs222\t5000000\t0\t1\t1\t0\t1\s*$/, 'second locus reuses shared macros and skips repeated checks');
ok(-s File::Spec->catfile($dir, 'plots_rs111.html'), 'first final plot exists');
ok(-s File::Spec->catfile($dir, 'plots_rs222.html'), 'second final plot exists');
ok(-s $index, 'plot index exists');
is(system(@cmd), 0, 'repeat request succeeds');
open $fh, '<', $calls or die $!;
@calls = <$fh>;
close $fh;
is(scalar(@calls), 2, 'repeat request reuses both completed loci');
unlink File::Spec->catfile($dir, 'plots_rs222.png') or die $!;
is(system(@cmd), 0, 'missing second plot is regenerated');
open $fh, '<', $calls or die $!;
@calls = <$fh>;
close $fh;
is(scalar(@calls), 3, 'only the missing locus reruns');
like($calls[2], qr/^rs222\t5000000\t0\t1\t0\t1\t0\s*$/,
    'first fresh locus rechecks shared macros after earlier cached plots');
open $fh, '>', $config or die $!;
print {$fh} "{\"FOREST_DOTSIZE\":20}\n";
close $fh;
is(system(@cmd), 0, 'unrelated forest setting is accepted');
open $fh, '<', $calls or die $!;
@calls = <$fh>;
close $fh;
is(scalar(@calls), 3, 'unrelated setting leaves per-locus plots cached');
open $fh, '>', $config or die $!;
print {$fh} "{\"GTF_ASSOC_PVARS\":\"TEST_P\"}\n";
close $fh;
is(system(@cmd), 0, 'changed GTF setting is accepted');
open $fh, '<', $calls or die $!;
@calls = <$fh>;
close $fh;
is(scalar(@calls), 5, 'changed GTF setting rerenders both loci');
like($calls[3], qr/\t0\t1\t0\s*$/, 'first changed-GTF locus verifies macros again');
like($calls[4], qr/\t1\t0\t1\s*$/, 'second changed-GTF locus reuses them');

open $fh, '>', $csv or die $!;
print {$fh} "SNP,CHR,BP\nrs111,1,100\nrs222,2,200\nrs333,3,300\n";
close $fh;
open $fh, '>', $config or die $!;
print {$fh} "{\"GTF_ASSOC_PVARS\":\"TEST_P2\"}\n";
close $fh;
{
    local $ENV{FAIL_ONCE_SNP} = 'rs222';
    is(system(@cmd) >> 8, 74, 'ODA disconnect preserves its exit code');
}
open $fh, '<', $calls or die $!;
@calls = <$fh>;
close $fh;
is(scalar(@calls), 7, 'interrupted batch stops after the failing second locus');
ok(-s File::Spec->catfile($dir, 'plots_rs111.html.request.md5'), 'first success has durable checkpoint');
ok(!-e $index, 'incomplete batch has no full index');
ok(-s File::Spec->catfile($dir, 'plots.partial.html'), 'incomplete batch has partial index');
my $progress = File::Spec->catfile($dir, 'plots.progress.json');
open $fh, '<', $progress or die $!;
my $report = decode_json(do { local $/; <$fh> });
close $fh;
is($report->{complete}, 1, 'progress reports one completed hit');
is($report->{remaining}, 2, 'progress reports two remaining hits');
is($report->{loci}[1]{status}, 'failed', 'failed locus is identified');
is(system(@cmd), 0, 'rerun resumes incomplete batch');
open $fh, '<', $calls or die $!;
@calls = <$fh>;
close $fh;
is(scalar(@calls), 9, 'rerun invokes only second and third loci');
like($calls[7], qr/^rs222\t.*\t0\t1\t0\s*$/, 'first fresh locus after resume verifies macros');
like($calls[8], qr/^rs333\t.*\t1\t0\t1\s*$/, 'next locus reuses verified macros');
ok(-s $index, 'completed batch writes full index');
ok(!-e File::Spec->catfile($dir, 'plots.partial.html'), 'completed batch removes partial index');
utime(time() + 5, time() + 5, $runner) or die $!;
is(system(@cmd), 0, 'launcher timestamp change does not invalidate plots');
open $fh, '<', $calls or die $!;
@calls = <$fh>;
close $fh;
is(scalar(@calls), 9, 'no locus reruns after launcher timestamp change');

my $busy_index = File::Spec->catfile($dir, 'busy.html');
my $busy_partial = File::Spec->catfile($dir, 'busy.partial.html');
my $busy_progress = File::Spec->catfile($dir, 'busy.progress.json');
mkdir $busy_partial or die $!;
mkdir $busy_progress or die $!;
my @busy_cmd = ($^X, File::Spec->catfile($Bin, 'run_selected_signed_ld_gtf.pl'),
    '--target-snps', 'rs444,rs555', '--runner-config', $config,
    '--output-html', $busy_index, '--single-runner', $runner);
is(system(@busy_cmd), 0, 'unavailable optional progress pages do not stop plotting');
ok(-s $busy_index, 'completed plot index is still written');
ok(-s File::Spec->catfile($dir, 'busy_rs444.html'), 'first plot survives status-page failure');
ok(-s File::Spec->catfile($dir, 'busy_rs555.html'), 'second plot survives status-page failure');
my @failed_status_temps = glob(File::Spec->catfile($dir, 'busy.*.tmp.*'));
is(scalar(@failed_status_temps), 0, 'failed status replacements leave no temporary files');
done_testing();
