#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Temp qw(tempdir);
use JSON::PP qw(encode_json decode_json);
use IO::Compress::Gzip qw(gzip $GzipError);
use Test::More;

my $dir = tempdir('raw manhattan XXXXX', TMPDIR => 1, CLEANUP => 1);
my $source = "$dir/stdized.tsv.gz";
my $header = join "\t", qw(
    CHR BP A1 A2 SNP PAIR_TAG GROUP1_BETA GROUP2_BETA DIFF_BETA
    GROUP1_SE GROUP2_SE DIFF_SE GROUP1_P GROUP2_P DIFF_P STD_DIFF_Z STD_DIFF_P
);
my $table = "$header\n"
    . join("\t", 1, 100, 'A', 'G', 'rsStrong', 'F_vs_M',
        0.2, 0.1, 0.1, 0.1, 0.1, 0.14, 0.01, 0.8, 0.4, 0.2, 0.6) . "\n"
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
my $filtered = "$dir/raw_all_test.stdized.wide_beta_se_p_p_lt_0p05.final.manifest.tsv";
my $all = "$dir/raw_all_test.stdized.wide_beta_se_p_all_snps.final.manifest.tsv";

is(system($^X, $driver, '--spec', $spec_path, '--step', 'extract_wide_subset'),
    0, 'default standardized-wide extraction succeeds');
is(manifest($filtered)->{rows_written}, 1,
    'default wide cache retains only a nominally significant SNP');
is(system($^X, $driver, '--spec', $spec_path, '--step', 'extract_wide_subset',
    '--manhattan-all-snps'), 0, 'explicit all-SNP wide extraction succeeds');
is(manifest($all)->{rows_written}, 2,
    'all-SNP wide cache restores the nonsignificant SNP upstream of plotting');
is(manifest($all)->{threshold}, 2,
    'all-SNP wide extraction does not apply the nominal P threshold');
ok(-s $filtered && -s $all, 'filtered and all-SNP manifests remain separate');
is(system($^X, $driver, '--spec', $spec_path, '--step', 'extract_wide_subset'),
    0, 'default mode is restorable after all-SNP mode');
is(manifest($filtered)->{rows_written}, 1,
    'restored default still uses the nominally filtered wide table');
done_testing();

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
