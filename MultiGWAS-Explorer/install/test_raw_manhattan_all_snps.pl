#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../DiffGWASDeps";
use FixedEffectMeta qw(fixed_effect_meta);
use File::Temp qw(tempdir);
use JSON::PP qw(encode_json decode_json);
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Test::More;

my $dir = tempdir('raw manhattan XXXXX', TMPDIR => 1, CLEANUP => 1);
my (undef, undef, undef, $correlated_p) = fixed_effect_meta(.2, .1, .2, .1, .5);
cmp_ok(abs($correlated_p - 0.02092), '<', 0.0001,
    'configured cohort correlation increases fixed-effect meta uncertainty');
my $source = "$dir/stdized.tsv.gz";
my $header = join "\t", qw(
    CHR BP A1 A2 SNP PAIR_TAG GROUP1_BETA GROUP2_BETA DIFF_BETA
    GROUP1_SE GROUP2_SE DIFF_SE GROUP1_P GROUP2_P DIFF_P STD_DIFF_Z STD_DIFF_P
);
my $table = "$header\n"
    . join("\t", 1, 100, 'A', 'G', 'rsStrong', 'F_vs_M',
        0.2, 0.1, 0.1, 0.1, 0.1, 0.14, 0.01, 0.8, 0.4, 0.2, 0.6) . "\n"
    . join("\t", 1, 150, 'G', 'T', 'rsMetaOnly', 'F_vs_M',
        0.16, 0.16, 0, 0.1, 0.1, 0.14, 0.11, 0.11, 0.9, 0, 0.9) . "\n"
    . join("\t", 1, 200, 'C', 'T', 'rsWeak', 'F_vs_M',
        0.1, 0.1, 0, 0.1, 0.1, 0.14, 0.6, 0.7, 0.8, 0.1, 0.9) . "\n";
gzip \$table => $source or die $GzipError;

my $spec_path = "$dir/spec.json";
my $spec = {
    source_mode => 'precomputed_diff_stdized',
    input_stdized => $source,
    artifact_stem => 'raw_all_test',
    project_tag => 'RAW_ALL_TEST',
    input_dir => $dir,
    output_dir => $dir,
    configs_dir => "$dir/configs",
    workdir => "$Bin/..",
    reference_build => 'hg19',
    threshold => 0.05,
    manhattan_fig_height => 1875,
    top_hit_focus_prefix => 'TEST',
    top_hit_selection_method => 'distance',
    top_hit_ld_source => 'HAPLOREG4',
    groups => [
        { tag => 'F', files => [] },
        { tag => 'M', files => [] },
    ],
    pairs => [
        { pair_tag => 'F_vs_M', group1 => 'F', group2 => 'M',
          prefix => 'TEST', label => 'Female vs male' },
    ],
};
open my $sfh, '>', $spec_path or die $!;
print {$sfh} encode_json($spec);
close $sfh;
my $driver = "$Bin/../auto_prepare_and_run_diff_gwas.pl";
my $filtered = "$dir/raw_all_test.stdized.wide_beta_se_p_p_lt_0p05.meta_ivw.final.manifest.tsv";
my $all = "$dir/raw_all_test.stdized.wide_beta_se_p_all_snps.meta_ivw.final.manifest.tsv";

is(system($^X, $driver, '--spec', $spec_path, '--step', 'extract_wide_subset'),
    0, 'default standardized-wide extraction succeeds');
is(manifest($filtered)->{rows_written}, 2,
    'default wide cache retains a SNP significant only in the computed meta-analysis');
my $runner = decode_json(read_file("$dir/configs/auto_raw_all_test_runner.json"));
like($runner->{DISPLAY_GWAS}, qr/(?:^|,)META(?:,|$)/,
    'individual-GWAS inputs display the computed meta P track by default');
is($runner->{MANHATTAN_FIG_WIDTH}, 4500,
    'tall SAS Manhattan plots keep a wide chromosome-friendly aspect ratio');
is($runner->{MANHATTAN_INCLUDE_X_CHR}, 0,
    'genomewide Manhattan excludes chromosome X by default');
my $preset = "$dir/configs/auto_raw_all_test_preset.json";
my $local = "$dir/meta_only_local.tsv.gz";
is(system($^X, "$Bin/../DiffGWASDeps/extract_single_snp_wide_diff_gwas.pl",
    '--config', $preset, '--target-snp', 'rsMetaOnly', '--window-bp', '50',
    '--output-dir', $dir, '--output', $local, '--manifest', "$dir/meta_only_local.manifest.tsv"),
    0, 'local SNP extraction also computes the meta track');
my $local_text = '';
gunzip($local => \$local_text) or die $GunzipError;
my @local_lines = split /\n/, $local_text;
my @local_head = split /\t/, $local_lines[0];
my @local_values = split /\t/, (grep { /\trsMetaOnly\t/ } @local_lines)[0];
my %local_row; @local_row{@local_head} = @local_values;
cmp_ok(abs($local_row{META_P} - 0.023725), '<', 0.0001,
    'local meta P uses inverse-variance fixed-effect combination');
open my $import_pipe, '-|', $^X,
    "$Bin/../DiffGWASDeps/generate_sas_wide_import_include.pl",
    '--config', $preset, '--input-file', $local, '--remote-basename', 'local.tsv.gz'
    or die $!;
my $sas_import = do { local $/; <$import_pipe> };
is(close($import_pipe), 1, 'SAS local import reads the actual meta-wide header');
like($sas_import, qr/^\s+META_P\s*$/m,
    'SAS imports computed meta P as a numeric column');
like($sas_import, qr/^\s+META_Z\s*$/m,
    'SAS imports computed meta Z as a numeric column');
unlike($sas_import, qr/\bMETA_Z\s*=/,
    'SAS retains the supplied or computed meta Z instead of recalculating it');
is(system($^X, $driver, '--spec', $spec_path, '--step', 'plot_manhattan',
    '--exclude-manhattan-tracks', 'META,DIFFERENTIAL', '--list-steps'),
    0, 'genomewide Manhattan can exclude meta and differential tracks');
$runner = decode_json(read_file("$dir/configs/auto_raw_all_test_runner.json"));
is($runner->{MANHATTAN_EXCLUDED_TRACKS}, 'TEST,META',
    'excluded track IDs are recorded in the runner');
unlike(join(',', $runner->{MANHATTAN_P_VAR}, @{$runner->{MANHATTAN_OTHER_P_VARS}}),
    qr/(?:META_P|TEST_DIFF_P)/, 'excluded P columns are absent from genomewide Manhattan');
like($runner->{GTF_ASSOC_PVARS}, qr/META_P/,
    'genomewide exclusion leaves local GTF tracks intact');
is(system($^X, $driver, '--spec', $spec_path, '--step', 'plot_manhattan',
    '--manhattan-track-order', 'META,M', '--list-steps'),
    0, 'genomewide Manhattan accepts a partial multi-track order');
$runner = decode_json(read_file("$dir/configs/auto_raw_all_test_runner.json"));
is($runner->{MANHATTAN_P_VAR}, 'META_P',
    'requested meta track becomes the first Manhattan panel');
is($runner->{MANHATTAN_TRACK_ORDER}, 'META,M,TEST,F',
    'unlisted genomewide tracks follow in their original order');
like($runner->{GTF_ASSOC_PVARS}, qr/^TEST_DIFF_P F_P M_P META_P$/,
    'custom genomewide order does not reorder local GTF tracks');
is(system($^X, $driver, '--spec', $spec_path, '--step', 'plot_local_manhattan',
    '--step', 'plot_local_gtf',
    '--exclude-local-manhattan-tracks', 'META,TEST',
    '--local-manhattan-track-order', 'M,F',
    '--exclude-local-gtf-tracks', 'TEST,F',
    '--local-gtf-track-order', 'META,M', '--list-steps'),
    0, 'local Manhattan and GTF accept independent multi-track exclusions and orders');
$runner = decode_json(read_file("$dir/configs/auto_raw_all_test_runner.json"));
is($runner->{LOCAL_MANHATTAN_TRACK_ORDER}, 'M,F',
    'local Manhattan scatter tracks follow the requested order');
is($runner->{LOCAL_MANHATTAN_P_VAR}, 'M_P',
    'local Manhattan first P column follows its order');
is($runner->{GTF_TRACK_ORDER}, 'META,M',
    'local GTF scatter tracks follow their independent order');
is($runner->{GTF_ASSOC_PVARS}, 'META_P M_P',
    'local GTF excludes two unwanted P columns');
is($runner->{MANHATTAN_TRACK_ORDER}, 'TEST,F,M,META',
    'local plot changes leave genomewide Manhattan order intact');
is(system($^X, $driver, '--spec', $spec_path, '--step', 'extract_wide_subset',
    '--manhattan-all-snps'), 0, 'explicit all-SNP wide extraction succeeds');
is(manifest($all)->{rows_written}, 3,
    'all-SNP wide cache restores the nonsignificant SNP upstream of plotting');
is(manifest($all)->{threshold}, 2,
    'all-SNP wide extraction does not apply the nominal P threshold');
ok(-s $filtered && -s $all, 'filtered and all-SNP manifests remain separate');
is(system($^X, $driver, '--spec', $spec_path, '--step', 'extract_wide_subset'),
    0, 'default mode is restorable after all-SNP mode');
is(manifest($filtered)->{rows_written}, 2,
    'restored default still uses the nominally filtered wide table');
my $x_input = "$dir/x_input.tsv.gz";
my $x_rows = "CHR\tBP\tP\n1\t100\t0.01\nX\t200\t0.01\n23\t300\t0.01\n";
gzip \$x_rows => $x_input or die $GzipError;
my $compact = "$Bin/../DiffGWASDeps/prepare_sas_manhattan_plot_input.pl";
my @compact_args = ($^X, $compact, '--input', $x_input, '--pvars', 'P',
    '--output', "$dir/x_compact.tsv.gz", '--schema-out', "$dir/x_schema.json",
    '--manifest', "$dir/x_manifest.json");
is(system(@compact_args), 0, 'autosomal compact Manhattan input succeeds');
my $x_manifest = decode_json(read_file("$dir/x_manifest.json"));
is($x_manifest->{rows_written}, 1, 'default compact input excludes X and chromosome 23');
is($x_manifest->{rows_x_removed}, 2, 'both X spellings are counted as excluded');
is(system(@compact_args, '--include-x-chr'), 0, 'explicit chromosome X request succeeds');
$x_manifest = decode_json(read_file("$dir/x_manifest.json"));
is($x_manifest->{rows_written}, 3, 'explicit request includes X and chromosome 23');
SKIP: {
    skip 'gnuplot is unavailable', 4 unless system('gnuplot', '--version') == 0;
    my $plotter = "$Bin/../DiffGWASDeps/gnuplot/pdl_gunplot_manhattan.pl";
    my @plot_args = ($^X, $plotter, '--data', "$dir/x_compact.tsv.gz",
        '--pcols', 'P', '--labels', 'Test', '--width', 500, '--height', 300,
        '--min-logp', 0, '--thin-mod', 1);
    is(system(@plot_args, '--out-prefix', "$dir/x_plot_default"), 0,
        'gnuplot default-X figure renders');
    is(manifest("$dir/x_plot_default.manifest.tsv")->{rows_kept}, 1,
        'gnuplot default figure excludes both X and chromosome 23');
    is(system(@plot_args, '--out-prefix', "$dir/x_plot_included",
        '--no-remove-x-chr'), 0, 'gnuplot explicit-X figure renders');
    is(manifest("$dir/x_plot_included.manifest.tsv")->{rows_kept}, 3,
        'gnuplot explicit figure includes both X spellings');
}
is(system($^X, $driver, '--spec', $spec_path, '--step', 'plot_manhattan',
    '--include-x-chr', '--list-steps'), 0, 'SAS runner accepts the X override');
$runner = decode_json(read_file("$dir/configs/auto_raw_all_test_runner.json"));
is($runner->{MANHATTAN_INCLUDE_X_CHR}, 1, 'X override reaches the SAS runner');
done_testing();

sub read_file {
    my ($path) = @_;
    open my $fh, '<', $path or die $!;
    local $/;
    return <$fh>;
}

sub manifest {
    my ($path) = @_;
    open my $fh, '<', $path or die "Cannot read $path: $!\n";
    my %values;
    while (<$fh>) {
        chomp;
        my ($key, $value) = split /\t/, $_, 2;
        $values{$key} = $value if defined $value;
    }
    close $fh;
    return \%values;
}
