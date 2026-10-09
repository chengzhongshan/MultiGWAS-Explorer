#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Temp qw(tempdir);
use File::Spec;
use IO::Compress::Gzip qw(gzip $GzipError);
use JSON::PP qw(encode_json);
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
printf '%s\t%s\t%s\t%s\n' "${TARGET_SNP}" "${LOCAL_WINDOW_BP}" "${GTF_LD_HEATMAP_LEGEND_TITLE}" "${GTF_LABEL_SNPS}" >> "${WORKDIR}/calls.tsv"
printf '%s\t%s\t%s\t%s\n' "${TARGET_SNP}" "${GTF_LABEL_LAYOUT:-auto}" "${GTF_LABEL_POSITIONS:-}" "${GTF_LABEL_FONT_SIZE:-}" >> "${WORKDIR}/layout_calls.tsv"
printf '%s\t%s\t%s\n' "${TARGET_SNP}" "${GTF_LABEL_HEADROOM_FRAC:-}" "${GTF_LABEL_CENTER_OFFSET:-}" >> "${WORKDIR}/headroom_calls.tsv"
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
like($index_html, qr/chr9:121741165/, 'index identifies the first plot by genomic location');
open my $calls_fh, '<', "$dir/calls.tsv" or die $!;
my @calls = <$calls_fh>;
close $calls_fh;
is(scalar @calls, 2, 'each selected locus has its own LD calculation');
like($calls[0], qr/^rs75453394\t500000\tSigned LD r2 x sign\(Z\) to rs75453394/, 'first locus uses its own signed-LD legend');

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

my $wide = "$dir/nearby.tsv.gz";
my $wide_text = "CHR\tBP\tSNP\tP\n6\t27500000\trsNearA\t1e-9\n6\t27600000\trsNearB\t2e-9\n";
gzip(\$wide_text => $wide) or die "gzip $wide: $GzipError\n";
open $cfg, '>', $config or die $!;
print {$cfg} encode_json({ DATA_GZ => $wide }), "\n";
close $cfg;
my $nearby_html = "$dir/nearby.html";
open $calls_fh, '<', "$dir/calls.tsv" or die $!;
@calls = <$calls_fh>;
close $calls_fh;
my $calls_before = scalar @calls;
is(system($^X, $dispatcher, '--target-snps', 'rsNearA,rsNearB',
    '--runner-config', $config, '--output-html', $nearby_html,
    '--single-runner', $runner, '--window-bp', 100000), 0,
    'explicit nearby targets are resolved from the GWAS table and merged');
open $calls_fh, '<', "$dir/calls.tsv" or die $!;
@calls = <$calls_fh>;
close $calls_fh;
is(scalar(@calls), $calls_before + 1, 'overlapping target windows invoke the SAS runner once');
like($calls[-1], qr/^rsNearA\t200000\t.*\trsNearA,rsNearB\s*$/,
    'merged locus uses the first lead as LD reference, expands its window, and labels both SNPs');
ok(-s $nearby_html, 'merged locus keeps the requested final HTML name');
ok(!-e "$dir/nearby_rsNearB.html", 'no redundant second local-GTF plot is created');

my $auto_nearby_csv = "$dir/auto_nearby.tsv";
open $fh, '>', $auto_nearby_csv or die $!;
print {$fh} "SNP\tCHR\tBP\tcommon_assoc_p\n";
print {$fh} "rsNearA\t6\t27500000\t1e-6\n";
print {$fh} "rsNearB\t6\t27600000\t1e-10\n";
close $fh;
open $calls_fh, '<', "$dir/calls.tsv" or die $!;
@calls = <$calls_fh>;
close $calls_fh;
$calls_before = scalar @calls;
is(system($^X, $dispatcher, '--targets-csv', $auto_nearby_csv,
    '--runner-config', $config, '--output-html', "$dir/auto_nearby.html",
    '--single-runner', $runner, '--window-bp', 100000), 0,
    'automatic overlapping leads share one plot');
open $calls_fh, '<', "$dir/calls.tsv" or die $!;
@calls = <$calls_fh>;
close $calls_fh;
is(scalar(@calls), $calls_before + 1, 'automatic region is rendered once');
like($calls[-1], qr/^rsNearB\t200000\t.*\trsNearB\s*$/,
    'automatic region labels and uses only its smallest-P lead');

open $calls_fh, '<', "$dir/layout_calls.tsv" or die $!;
my @layout_calls = <$calls_fh>;
close $calls_fh;
like($layout_calls[-2], qr/^rsNearA\thorizontal\trsNearA=\d+ rsNearB=\d+\t10\s*$/,
    'explicit nearby targets use planned horizontal positions when they fit');
like($layout_calls[-1], qr/^rsNearB\thorizontal\trsNearB=\d+\t10\s*$/,
    'automatic single lead receives a centered horizontal label');
open $calls_fh, '<', "$dir/headroom_calls.tsv" or die $!;
my @headroom_calls = <$calls_fh>;
close $calls_fh;
like($headroom_calls[-2], qr/^rsNearA\t0\.\d+\t0\.\d+\s*$/,
    'explicit locus receives planned height and centering');

{
    local $ENV{GTF_LABEL_PLAN_BACKEND} = 'sas';
    is(system($^X, $dispatcher, '--target-snps', 'rsNearA,rsNearB',
        '--runner-config', $config, '--output-html', "$dir/nearby_sas_fallback.html",
        '--single-runner', $runner, '--window-bp', 100000), 0,
        'SAS fallback can be selected when Perl planning is unavailable');
}
open $calls_fh, '<', "$dir/layout_calls.tsv" or die $!;
@layout_calls = <$calls_fh>;
close $calls_fh;
like($layout_calls[-1], qr/^rsNearA\tauto\t\t10\s*$/,
    'fallback leaves positions empty for the SAS macro to plan');

open $cfg, '>', $config or die $!;
print {$cfg} encode_json({ DATA_GZ => $wide, GTF_LD_REFERENCE_SNP => 'rsNearB' }), "\n";
close $cfg;
is(system($^X, $dispatcher, '--target-snps', 'rsNearA,rsNearB',
    '--runner-config', $config, '--output-html', "$dir/nearby_override.html",
    '--single-runner', $runner, '--window-bp', 100000), 0,
    'explicit LD lead override is accepted');
open $calls_fh, '<', "$dir/calls.tsv" or die $!;
@calls = <$calls_fh>;
close $calls_fh;
like($calls[-1], qr/^rsNearB\t200000\t.*\trsNearA,rsNearB\s*$/,
    'explicit LD lead override changes the reference but keeps all labels');

my $prior_html = "$dir/prior.html";
for my $prior ([rsPriorA => 21, 42841988], [rsPriorB => 21, 42858367]) {
    my $prior_csv = "$dir/prior_$prior->[0]_top_hit.csv";
    open my $prior_fh, '>', $prior_csv or die $!;
    print {$prior_fh} "SNP,CHR,BP\n$prior->[0],$prior->[1],$prior->[2]\n";
    close $prior_fh;
}
open $cfg, '>', $config or die $!;
print {$cfg} "{}\n";
close $cfg;
open $calls_fh, '<', "$dir/calls.tsv" or die $!;
@calls = <$calls_fh>;
close $calls_fh;
$calls_before = scalar @calls;
is(system($^X, $dispatcher, '--target-snps', 'rsPriorA,rsPriorB',
    '--runner-config', $config, '--output-html', $prior_html,
    '--single-runner', $runner, '--window-bp', 650000), 0,
    'coordinates from completed per-target CSVs avoid a genome-wide rsID scan');
open $calls_fh, '<', "$dir/calls.tsv" or die $!;
@calls = <$calls_fh>;
close $calls_fh;
is(scalar(@calls), $calls_before + 1, 'prior nearby target CSVs merge into one resumed locus');
like($calls[-1], qr/^rsPriorA\t666379\t.*\trsPriorA,rsPriorB\s*$/,
    'prior coordinates produce the expected expanded shared window');

my $auto_path = "$Bin/../auto_prepare_and_run_diff_gwas.pl";
open my $auto_fh, '<', $auto_path or die $!;
my $auto_source = do { local $/; <$auto_fh> };
close $auto_fh;
like($auto_source, qr/--targets-csv "\$common_ld_artifacts\{leads\}"/,
    'automatic signed-LD plots consume the PLINK2 LD-pruned lead table');

done_testing();
