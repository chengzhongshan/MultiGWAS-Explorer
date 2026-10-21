#!/usr/bin/env perl
use strict;
use warnings;

use FindBin qw($Bin);
use File::Spec;
use File::Temp qw(tempdir);
use JSON::PP qw(decode_json encode_json);
use Test::More;

my $root = File::Spec->rel2abs(File::Spec->catdir($Bin, '..'));
my $driver = File::Spec->catfile($root, 'auto_prepare_and_run_diff_gwas.pl');
my $tmp = tempdir('track_regex_test_XXXXXX', TMPDIR => 1, CLEANUP => 1);
my $spec_path = File::Spec->catfile($tmp, 'spec.json');
my $runner_path = File::Spec->catfile($tmp, 'auto_track_regex_fixture_runner.json');

my @groups = map { +{ tag => $_, files => [lc($_) . '.gz'] } }
  qw(ALL_FEMALE ALL_MALE EUR_FEMALE EUR_MALE ASN_FEMALE ASN_MALE);
my @pairs = (
    { pair_tag => 'ALL_FEMALE_vs_MALE', group1 => 'ALL_FEMALE', group2 => 'ALL_MALE', prefix => 'ALL', label => 'All' },
    { pair_tag => 'EUR_FEMALE_vs_MALE', group1 => 'EUR_FEMALE', group2 => 'EUR_MALE', prefix => 'EUR', label => 'EUR' },
    { pair_tag => 'ASN_FEMALE_vs_MALE', group1 => 'ASN_FEMALE', group2 => 'ASN_MALE', prefix => 'ASN', label => 'ASN' },
);
my $spec = {
    source_mode => 'raw_pgc_vcf_sumstats',
    project_tag => 'track_regex_fixture',
    artifact_stem => 'track_regex_fixture',
    input_dir => $tmp,
    output_dir => $tmp,
    configs_dir => $tmp,
    workdir => $root,
    reference_build => 'hg19',
    groups => \@groups,
    pairs => \@pairs,
};
open my $spec_fh, '>', $spec_path or die "Cannot write $spec_path: $!\n";
print {$spec_fh} JSON::PP->new->canonical(1)->pretty(1)->encode($spec);
close $spec_fh or die "Cannot close $spec_path: $!\n";

local $ENV{PIPELINE_RESPECT_SPEC_WORKDIR} = 1;

my @matrix_args = (
    '--exclude-manhattan-tracks', '/_META$/i',
    '--manhattan-track-order', '/^EUR(?:_|$)/i',
    '--exclude-local-manhattan-tracks', '/^ASN_/i',
    '--local-manhattan-track-order', '/^ALL(?:_|$)/i',
    '--exclude-local-gtf-tracks', '/_META$/i',
    '--local-gtf-track-order', '/^ASN(?:_|$)/i',
);
is(system($^X, $driver, '--spec', $spec_path, '--mode', 'configs', '--plots', 'manhattan', @matrix_args), 0,
    'regex exclusion and ordering matrix generates configs');
my $runner = read_json($runner_path);
is($runner->{MANHATTAN_EXCLUDED_TRACKS}, 'ALL_META,EUR_META,ASN_META', 'genomewide regex excludes all meta tracks');
is($runner->{MANHATTAN_ORDER_PRIORITY}, 'EUR,EUR_FEMALE,EUR_MALE', 'genomewide regex orders EUR tracks');
is($runner->{LOCAL_MANHATTAN_EXCLUDED_TRACKS}, 'ASN,ASN_FEMALE,ASN_MALE,ASN_META', 'local Manhattan regex matches IDs and variables');
is($runner->{LOCAL_MANHATTAN_ORDER_PRIORITY}, 'ALL,ALL_FEMALE,ALL_MALE,ALL_META', 'local Manhattan regex orders ALL tracks');
is($runner->{GTF_EXCLUDED_TRACKS}, 'ALL_META,EUR_META,ASN_META', 'local GTF regex excludes all meta tracks');
is($runner->{GTF_ORDER_PRIORITY}, 'ASN,ASN_FEMALE,ASN_MALE', 'local GTF regex orders remaining ASN tracks');

is(system($^X, $driver, '--spec', $spec_path, '--mode', 'configs', '--plots', 'manhattan',
        '--exclude-manhattan-track-regex', '/_META$/i'), 0,
    'dedicated genomewide exclusion regex generates configs');
$runner = read_json($runner_path);
is($runner->{MANHATTAN_EXCLUDED_TRACKS}, 'ALL_META,EUR_META,ASN_META',
    'dedicated regex excludes all meta tracks');

my $resolution_output = qx{$^X "$driver" --spec "$spec_path" --mode configs --plots manhattan --exclude-manhattan-tracks META 2>&1};
is($?, 0, 'aggregate META exclusion generates a track-resolution summary');
like($resolution_output,
    qr/\[tracks\] genomewide Manhattan exclusion input condition\(s\): list='META'/,
    'track-resolution summary reports the supplied exclusion condition');
like($resolution_output,
    qr/\[tracks\] genomewide Manhattan excluded tracks: ALL_META, EUR_META, ASN_META/,
    'track-resolution summary reports every concretely excluded track');

is(system($^X, $driver, '--spec', $spec_path, '--mode', 'configs', '--plots', 'manhattan',
        '--display-gwas', '/^(?:ALL|EUR)(?:_FEMALE|_MALE)?$/i'), 0,
    'regex inclusion with alternation generates configs');
$runner = read_json($runner_path);
is($runner->{DISPLAY_GWAS}, 'ALL,EUR,ALL_FEMALE,ALL_MALE,EUR_FEMALE,EUR_MALE',
    'display regex selects the expected tracks in catalog order');

is(system($^X, $driver, '--spec', $spec_path, '--mode', 'configs', '--plots', 'manhattan',
        '--no-compute-meta'), 0,
    '--no-compute-meta generates configs');
$runner = read_json($runner_path);
unlike($runner->{DISPLAY_GWAS}, qr/(?:^|,)\w*_?META(?:,|$)/,
    '--no-compute-meta removes computed meta tracks from the catalog');

my $heartbeat_output = qx{$^X "$driver" --spec "$spec_path" --step extract_wide_subset --force --step-heartbeat-seconds 1 2>&1};
my $heartbeat_status = $?;
isnt($heartbeat_status, 0, 'missing-data fixture preserves extract failure status');
like($heartbeat_output, qr/\[progress\] extract_wide_subset is running/,
    'extract_wide_subset reports an immediate running heartbeat');

done_testing();

sub read_json {
    my ($path) = @_;
    open my $fh, '<', $path or die "Cannot read $path: $!\n";
    local $/;
    my $decoded = decode_json(<$fh>);
    close $fh;
    return $decoded;
}
