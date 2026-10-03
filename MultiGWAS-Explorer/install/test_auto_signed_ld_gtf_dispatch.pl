#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Temp qw(tempdir);
use File::Spec;
use Test::More;

my $dir = tempdir('signed_ld_dispatch_XXXXXX', TMPDIR => 1, CLEANUP => 1);
my $dispatcher = "$Bin/../DiffGWASDeps/run_selected_signed_ld_gtf.pl";
my $runner = "$dir/mock_single_runner.sh";
my $config = "$dir/runner.json";
open my $cfg, '>', $config or die $!;
print {$cfg} "{}\n";
close $cfg;
open my $mock, '>', $runner or die $!;
print {$mock} <<'MOCK';
#!/usr/bin/env bash
set -euo pipefail
[[ -z "${DATA_GZ+x}" && -z "${REMOTE_DATA_BASENAME+x}" ]]
[[ "${GTF_LD_DISPLAY_MODE}" == "heatmap" ]]
[[ "${GTF_LD_REFERENCE_SNP}" == "${TARGET_SNP}" ]]
printf '%s\t%s\t%s\n' "${TARGET_SNP}" "${LOCAL_WINDOW_BP}" "${GTF_LD_HEATMAP_LEGEND_TITLE}" >> "${WORKDIR}/calls.tsv"
printf '%s\t%s\n' "${TARGET_SNP}" "${GTF_LD_R2_CACHE}" >> "${WORKDIR}/cache_calls.tsv"
printf '<html><img src="%s.png"></html>\n' "${TARGET_SNP}" > "${WORKDIR}/${OUTPUT_HTML_BASENAME}"
printf 'PNG:%s\n' "${TARGET_SNP}" > "${WORKDIR}/${OUTPUT_HTML_BASENAME%.html}.png"
MOCK
close $mock;
local $ENV{WORKDIR} = $dir;
local $ENV{DATA_GZ} = 'whole_genome.tsv.gz';
local $ENV{REMOTE_DATA_BASENAME} = 'whole_genome.tsv.gz';
my $csv = "$dir/hits.csv";
open my $fh, '>', $csv or die $!;
print {$fh} "hit_order,panel_index,CHR,BP,SNP,gene\n";
print {$fh} "1,1,9,121741165,rs75453394,DAB2IP\n";
print {$fh} "2,2,1,100,rsOther,GENE\n";
close $fh;
my $index = "$dir/selected.html";
is(system($^X, $dispatcher, '--targets-csv', $csv, '--runner-config', $config,
    '--output-html', $index, '--single-runner', $runner,
    '--window-bp', 500000, '--population', 'EUR'), 0,
    'automatically selected loci use one signed-LD runner per reference SNP');
ok(-s $index, 'multiple loci produce an HTML plot index');
open my $index_fh, '<', $index or die $!;
my $index_html = do { local $/; <$index_fh> };
close $index_fh;
like($index_html, qr/selected_rs75453394\.html/, 'index links to first locus plot');
like($index_html, qr/selected_rsOther\.html/, 'index links to second locus plot');
open my $calls_fh, '<', "$dir/calls.tsv" or die $!;
my @calls = <$calls_fh>;
close $calls_fh;
is(scalar @calls, 2, 'each selected locus has its own LD calculation');
like($calls[0], qr/^rs75453394\t500000\tSigned LD r2 to rs75453394/, 'first locus uses its own signed-LD legend');

open $fh, '>', $csv or die $!;
print {$fh} "hit_order,panel_index,CHR,BP,SNP,gene\n1,1,9,121741165,rs75453394,DAB2IP\n";
close $fh;
my $single_html = "$dir/single.html";
is(system($^X, $dispatcher, '--targets-csv', $csv, '--runner-config', $config,
    '--output-html', $single_html, '--single-runner', $runner,
    '--window-bp', 500000), 0, 'one automatic hit keeps the normal local-GTF output filename');
ok(-s $single_html, 'single-locus HTML was generated');

my $cache = "$dir/first_ld.tsv";
open my $ld, '>', $cache or die $!;
print {$ld} "SNP\tLD_R2\nrs75453394\t1\n";
close $ld;
open $cfg, '>', $config or die $!;
print {$cfg} '{"GTF_LD_R2_CACHE_BY_SNP":{"rs75453394":"', $cache, '"}}', "\n";
close $cfg;
my $explicit_html = "$dir/explicit.html";
is(system($^X, $dispatcher, '--target-snps', 'rs75453394,rsOther',
    '--runner-config', $config, '--output-html', $explicit_html,
    '--single-runner', $runner, '--window-bp', 500000), 0,
    'explicit multiple-SNP request uses the resumable dispatcher');
ok(-s $explicit_html, 'explicit request writes a completed index');
open my $cache_fh, '<', "$dir/cache_calls.tsv" or die $!;
my @cache_calls = <$cache_fh>;
close $cache_fh;
is($cache_calls[-2], "rs75453394\t$cache\n", 'explicit first SNP receives its own LD cache');
is($cache_calls[-1], "rsOther\t\n", 'explicit second SNP does not inherit first LD cache');

done_testing();
