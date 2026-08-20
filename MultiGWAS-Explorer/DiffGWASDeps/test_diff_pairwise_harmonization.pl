#!/usr/bin/env perl
use strict;
use warnings;
use Test::More;
use File::Temp qw(tempdir);
use File::Spec;
use FindBin qw($Bin);
use JSON::PP qw(encode_json);

my $tmp = tempdir(CLEANUP => 1);
my $input = File::Spec->catfile($tmp, 'merged.tsv');
my $output = File::Spec->catfile($tmp, 'diff.tsv');
my $manifest = File::Spec->catfile($tmp, 'manifest.tsv');
my $config = File::Spec->catfile($tmp, 'config.json');

open my $ih, '>', $input or die $!;
print {$ih} join("\t", qw(CHR BP A1 A2 SNP GWAS_TAG SOURCE_FILE CHR_ORIGINAL IS_CHRX FRQ_A FRQ_U INFO BETA SE P)), "\n";
for my $row (
    [1, 100, 'A', 'C', 'rs_direct', 'F', 'f', 1, 0, .20, .22, .9, .10, .05, .04],
    [1, 100, 'A', 'C', 'rs_direct', 'M', 'm', 1, 0, .21, .23, .9, .02, .04, .6],
    [1, 200, 'A', 'T', 'rs_ambiguous', 'F', 'f', 1, 0, .40, .41, .9, .10, .05, .04],
    [1, 200, 'A', 'T', 'rs_ambiguous', 'M', 'm', 1, 0, .42, .41, .9, .02, .04, .6],
    [1, 300, 'G', 'A', 'rs_freq', 'F', 'f', 1, 0, .10, .10, .9, .10, .05, .04],
    [1, 300, 'G', 'A', 'rs_freq', 'M', 'm', 1, 0, .70, .70, .9, .02, .04, .6],
) {
    print {$ih} join("\t", @$row), "\n";
}
close $ih;

open my $ch, '>', $config or die $!;
print {$ch} encode_json({
    input => $input, output => $output, manifest => $manifest,
    base_cols => [qw(CHR BP A1 A2 SNP)], pairs => { F_vs_M => ['F', 'M'] },
    exclude_strand_ambiguous => 1, max_eaf_abs_diff => 0.2, rho => 0,
});
close $ch;

my $script = File::Spec->catfile($Bin, 'diff_pairwise_gwas.pl');
is(system($^X, $script, '--config', $config), 0, 'harmonized differential run succeeds');

open my $oh, '<', $output or die $!;
my @lines = <$oh>;
close $oh;
is(scalar(@lines), 2, 'only header and one accepted direct non-ambiguous locus remain');
like($lines[0], qr/ALLELE_RELATION\tSTRAND_AMBIGUOUS\tEAF_ABS_DIFF\tHARMONIZATION_STATUS/, 'audit columns are emitted');
like($lines[1], qr/rs_direct.*\tDIRECT\t0\t0\.0?1\taccepted_exact_allele_match\t/, 'accepted row records harmonization status');

open my $mh, '<', $manifest or die $!;
my $manifest_text = do { local $/; <$mh> };
close $mh;
like($manifest_text, qr/skipped_strand_ambiguous\t1/, 'ambiguous pair exclusion counted');
like($manifest_text, qr/skipped_eaf_discordance\t1/, 'frequency discordance exclusion counted');
like($manifest_text, qr/accepted_exact_allele_match\t1/, 'accepted exact match counted');

done_testing();
