#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Temp qw(tempdir);
use File::Path qw(make_path);
use File::Copy qw(copy);
use File::Basename qw(dirname basename);
use Cwd qw(abs_path getcwd);
use JSON::PP qw(encode_json decode_json);
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Test::More;

my $dir = tempdir('multigwas meta XXXXX', TMPDIR => 1, CLEANUP => 1);
my $driver = "$Bin/../auto_prepare_and_run_diff_gwas.pl";
# Deliberately use an arbitrary filename and a directory with spaces.
my $input = "$dir/results.gz";
my $table = join("\t", qw(CHR BP SNP A2 A1 BETA_MP2PRT SE_MP2PRT P_MP2PRT PR_meta BETAR_meta WEIGHTED_Z_meta PWZ_meta BETA_DS_ALL SE_DS_ALL P_DS_ALL)) . "\n"
    . join("\t", 1, 100, 'rsTest', 'G', 'A', .4, .1, .0001, .0123, -.234, -3.5, .0004, .1, .2, .6) . "\n"
    . join("\t", 1, 200, 'rsMissing', 'T', 'C', ('.') x 7, .2, .3, .5) . "\n"
    . join("\t", 2, 300, 'rsOther', 'G', 'T', .2, .1, .05, .02, .15, 2.3, .03, .1, .2, .6) . "\n";
my ($header_line, @data_lines) = split /(?<=\n)/, $table;
gzip \$header_line => $input or die $GzipError;
my $second_stream = '';
gzip \join('', @data_lines) => \$second_stream or die $GzipError;
open my $gz_append, '>>:raw', $input or die $!;
print {$gz_append} $second_stream;
close $gz_append or die $!;
my $relative_dir = "$dir/relative input";
make_path($relative_dir);
copy($input, "$relative_dir/combined_meta.gz") or die $!;
my $previous_cwd = getcwd();
chdir dirname($relative_dir) or die $!;
my $relative_spec_path = "$dir/relative.spec.json";
my $relative_status = system($^X, $driver, '--gwas-dir', basename($relative_dir),
    '--spec-out', $relative_spec_path, '--generate-spec-only');
chdir $previous_cwd or die $!;
is($relative_status, 0, 'relative --gwas-dir generates a spec');
my $relative_spec = read_json($relative_spec_path);
is(abs_path($relative_spec->{output_dir}), abs_path($relative_dir),
    'auto-detected output directory is absolute and independent of the script location');
$relative_spec->{configs_dir} = "$dir/relative_configs";
write_json($relative_spec_path, $relative_spec);
is(system($^X, $driver, '--spec', $relative_spec_path, '--skip-plots',
    '--step', 'extract_wide_subset'), 0, 'relative directory spec converts from another working directory');
ok(-s "$relative_dir/" . $relative_spec->{artifact_stem} . '.merged_plotwide.tsv.gz',
    'wide output is written alongside the input directory');
my $spec_path = "$dir/import.json";
is(system($^X, $driver, '--input-merged', $input, '--spec-out', $spec_path,
    '--generate-spec-only'), 0, 'direct import generates a spec without separate cohort files');
my $spec = read_json($spec_path);
is($spec->{source_mode}, 'merged_gwas_table', 'uses existing merged workflow');
is(scalar @{$spec->{pairs}}, 1, 'retains cohort comparison');
is($spec->{merged_extra_tracks}[0]{p_col}, 'PR_meta', 'preserves old AOA meta P selection');
is($spec->{extra_tracks}[0]{zvar}, 'META_Z', 'meta track is exposed for plotting');
$spec->{configs_dir} = "$dir/configs";
$spec->{output_dir} = "$dir/output";
$spec->{artifact_stem} = 'meta_test';
$spec->{top_hit_ld_source} = 'HAPLOREG4';
write_json($spec_path, $spec);
is(system($^X, $driver, '--spec', $spec_path, '--skip-plots', '--list-steps'), 0,
    'existing pipeline generates configs for direct import');
my $runner_path = "$dir/configs/auto_meta_test_runner.json";
is(system($^X, $driver, '--spec', $spec_path, '--step', 'plot_local_gtf',
    '--target-snps', 'rsTest,rsOther', '--list-steps'), 0,
    'multi-target SAS GTF configuration is generated');
is(read_json($runner_path)->{GTF_LD_DISPLAY_MODE}, 'heatmap',
    'explicit target GTF plots default to the signed-LD heatmap');
is(system($^X, $driver, '--spec', $spec_path, '--step', 'plot_manhattan',
    '--manhattan-all-snps', '--list-steps'), 0,
    'SAS Manhattan accepts the explicit all-SNP option');
is(read_json($runner_path)->{MANHATTAN_ALL_SNPS}, 1,
    'all-SNP override reaches the SAS runner');
my $preset_path = "$dir/configs/auto_meta_test_preset.json";
my $preset = read_json($preset_path);
is(system($^X, "$Bin/../DiffGWASDeps/convert_merged_gwas_to_plotwide.pl",
    '--config', $preset_path), 0, 'merged converter runs');
my $output = '';
gunzip $preset->{output} => \$output or die $GunzipError;
my @lines = split /\n/, $output;
my @header = split /\t/, shift @lines;
my @values = split /\t/, $lines[0], -1;
my %row; @row{@header} = @values;
cmp_ok(abs($row{META_P} - .0123), '<', 1e-12, 'supplied meta P retained rather than recalculated');
cmp_ok(abs($row{META_Z} + 3.5), '<', 1e-12, 'supplied signed weighted Z retained');
cmp_ok(abs($row{META_BETA} + .234), '<', 1e-12, 'supplied meta beta retained');
cmp_ok(abs($row{MP2PRT_DS_ALL_DIFF_BETA} - .3), '<', 1e-12, 'cohort differential calculation unchanged');
is(scalar @lines, 3, 'rows with missing cohort values retained');
@values = split /\t/, $lines[1], -1;
@row{@header} = @values;
ok($row{META_P} !~ /\d/, 'missing meta P is not fabricated');
open my $sas_import_pipe, '-|', $^X,
    "$Bin/../DiffGWASDeps/generate_sas_wide_import_include.pl",
    '--config', $preset_path, '--input-file', $preset->{output},
    '--remote-basename', 'merged.tsv.gz' or die $!;
my $sas_import = do { local $/; <$sas_import_pipe> };
ok(close($sas_import_pipe), 'merged SAS import reads the real wide header');
like($sas_import, qr/^\s+META_Z\s*$/m,
    'merged SAS import includes the supplied signed meta Z');
unlike($sas_import, qr/\bMETA_Z\s*=/,
    'merged SAS import does not overwrite the supplied signed meta Z');
SKIP: {
    skip 'gnuplot or PDL is unavailable', 14
        unless system('gnuplot', '--version') == 0 && eval { require PDL; 1 };
    my $gnu_spec = "$dir/gunplot.spec.json";
    is(system($^X, "$Bin/../auto_prepare_and_run_diff_gwas_with_gunplot.pl",
        '--input-merged', $input, '--spec-out', $gnu_spec,
        '--plots', 'manhattan', '--display-gwas', 'Meta'), 0,
        'gnuplot entry point imports and plots supplied meta results directly');
    my $gnu_runner = read_json("$dir/auto_" . (read_json($gnu_spec)->{artifact_stem}) . '_runner.json');
    is($gnu_runner->{MANHATTAN_P_VAR}, 'META_P', 'gnuplot uses the supplied meta P column');
    opendir my $plot_dh, $dir or die $!;
    my @png = map { "$dir/$_" }
        grep { /GUNPLOT_manhattan_display_META_.*\.png$/ } readdir $plot_dh;
    closedir $plot_dh;
    ok(@png && -s $png[0], 'gnuplot writes a nonempty meta Manhattan image');
    my $relative_gnu_spec = "$dir/relative_gunplot.spec.json";
    chdir dirname($relative_dir) or die $!;
    my $relative_gnu_status = system($^X,
        "$Bin/../auto_prepare_and_run_diff_gwas_with_gunplot.pl",
        '--gwas-dir', basename($relative_dir), '--spec-out', $relative_gnu_spec,
        '--plots', 'manhattan', '--display-gwas', 'Meta');
    chdir $previous_cwd or die $!;
    is($relative_gnu_status, 0, 'gnuplot accepts a relative GWAS directory from another working directory');
    my $relative_gnu_cfg = read_json($relative_gnu_spec);
    is(abs_path($relative_gnu_cfg->{output_dir}), abs_path($relative_dir),
        'gnuplot generates an absolute output directory');
    opendir my $relative_plot_dh, $relative_dir or die $!;
    my @relative_png = map { "$relative_dir/$_" }
        grep { /GUNPLOT_manhattan_display_META_.*\.png$/ } readdir $relative_plot_dh;
    closedir $relative_plot_dh;
    ok(@relative_png && -s $relative_png[0], 'gnuplot renders a meta plot from the relative directory');
    my $default_gnu_spec = "$dir/relative_gunplot_default.spec.json";
    chdir dirname($relative_dir) or die $!;
    my $default_gnu_status = system($^X,
        "$Bin/../auto_prepare_and_run_diff_gwas_with_gunplot.pl",
        '--gwas-dir', basename($relative_dir), '--spec-out', $default_gnu_spec,
        '--plots', 'manhattan');
    chdir $previous_cwd or die $!;
    is($default_gnu_status, 0, 'default merged gnuplot refreshes Meta-only runner and renders all tracks');
    my $default_runner = read_json("$dir/auto_" . $relative_spec->{artifact_stem} . '_gunplot_runner.json');
    is_deeply($default_runner->{MANHATTAN_OTHER_P_VARS},
        [qw(MP2PRT_DS_ALL_GROUP2_P MP2PRT_DS_ALL_GROUP1_P META_P)],
        'gnuplot translates cohort aliases to columns in the wide file');
    my $default_prefix = $default_runner->{OUTPUT_PREFIX};
    $default_prefix =~ s/_SAS_/_GUNPLOT_/g;
    ok(-s "$relative_dir/$default_prefix.png",
        'default merged gnuplot writes a nonempty Manhattan image');
    open my $default_manifest_fh, '<', "$relative_dir/$default_prefix.manifest.tsv" or die $!;
    my %default_manifest = map { chomp; split /\t/, $_, 2 } <$default_manifest_fh>;
    close $default_manifest_fh;
    is($default_manifest{rows_scanned}, 2,
        'genome-wide gnuplot scans only rows with a displayed P below 0.05');
    chdir dirname($relative_dir) or die $!;
    my $all_gnu_status = system($^X,
        "$Bin/../auto_prepare_and_run_diff_gwas_with_gunplot.pl",
        '--gwas-dir', basename($relative_dir), '--spec-out', $default_gnu_spec,
        '--plots', 'manhattan', '--manhattan-all-snps');
    chdir $previous_cwd or die $!;
    is($all_gnu_status, 0, 'gnuplot accepts the explicit all-SNP Manhattan option');
    open $default_manifest_fh, '<', "$relative_dir/$default_prefix.manifest.tsv" or die $!;
    %default_manifest = map { chomp; split /\t/, $_, 2 } <$default_manifest_fh>;
    close $default_manifest_fh;
    is($default_manifest{rows_scanned}, 3,
        'explicit all-SNP mode bypasses the nominal P filter');
    my $custom_spec = "$dir/custom_gunplot.spec.json";
    is(system($^X, "$Bin/../auto_prepare_and_run_diff_gwas_with_gunplot.pl",
        '--input-merged', $input, '--spec-out', $custom_spec, '--plots', 'manhattan',
        '--exclude-manhattan-tracks', 'MP2PRT_DS_ALL,MP2PRT',
        '--manhattan-track-order', 'META,DS_ALL',
        '--exclude-local-manhattan-tracks', 'META,MP2PRT',
        '--local-manhattan-track-order', 'DS_ALL,MP2PRT_DS_ALL',
        '--exclude-local-gtf-tracks', 'MP2PRT,MP2PRT_DS_ALL',
        '--local-gtf-track-order', 'META,DS_ALL',
        '--include-x-chr'), 0,
        'gnuplot wrapper forwards independent track exclusions and orders');
    my $custom_artifact = read_json($custom_spec)->{artifact_stem};
    my $custom_runner = read_json("$dir/auto_${custom_artifact}_runner.json");
    is($custom_runner->{MANHATTAN_TRACK_ORDER}, 'META,DS_ALL',
        'gnuplot genomewide panels follow the requested order');
    is($custom_runner->{LOCAL_MANHATTAN_TRACK_ORDER}, 'DS_ALL,MP2PRT_DS_ALL',
        'gnuplot local Manhattan panels retain their independent order');
    is($custom_runner->{GTF_TRACK_ORDER}, 'META,DS_ALL',
        'gnuplot local GTF scatter panels retain their independent order');
    is($custom_runner->{MANHATTAN_INCLUDE_X_CHR}, 1,
        'gnuplot X override is recorded in the shared runner');
    my $compact_manifest = "$Bin/../cache/sas_manhattan/"
        . $custom_runner->{OUTPUT_PREFIX} . '.p_lt_0_05.manifest.json';
    is(read_json($compact_manifest)->{include_x_chr}, 1,
        'gnuplot X override reaches compact genomewide input');
    is($default_manifest{rows_thinned}, 0,
        'explicit all-SNP gnuplot mode disables background point thinning');
}
ok(system($^X, $driver, '--input-merged', $input, '--spec', $spec_path,
    '--generate-spec-only') != 0, 'ambiguous input options rejected');
my $invalid = "$dir/invalid.tsv";
open my $fh, '>', $invalid or die $!;
print {$fh} "CHR\tBP\tSNP\tPR_meta\n1\t100\trsTest\t0.1\n";
close $fh;
ok(system($^X, $driver, '--input-merged', $invalid, '--generate-spec-only') != 0,
    'unsupported meta-only table fails clearly');
done_testing();

sub read_json {
    open my $fh, '<', $_[0] or die $!;
    return decode_json(do { local $/; <$fh> });
}
sub write_json {
    open my $fh, '>', $_[0] or die $!;
    print {$fh} encode_json($_[1]);
    close $fh or die $!;
}
