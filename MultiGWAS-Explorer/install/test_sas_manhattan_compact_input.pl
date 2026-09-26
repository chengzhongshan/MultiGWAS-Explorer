#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Temp qw(tempdir);
use JSON::PP qw(encode_json decode_json);
use IO::Compress::Gzip qw(gzip $GzipError);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Test::More;

my $dir = tempdir('sas manhattan XXXXX', TMPDIR => 1, CLEANUP => 1);
my $input = "$dir/wide.tsv.gz";
my $schema = "$dir/schema.json";
my $output = "$dir/compact.tsv.gz";
my $compact_schema = "$dir/compact.schema.json";
my $manifest = "$dir/compact.manifest.json";
my $source = join("\n",
    join("\t", qw(CHR BP SNP DIFF_P GROUP1_P GROUP2_P META_P BETA)),
    join("\t", qw(chrX 200 rsX 0.8 0.04 0.9 0.9 0.2)),
    join("\t", qw(2 100 rsBoundary 0.05 0.8 0.9 0.7 0.3)),
    join("\t", qw(1 50 rsMeta 0.8 0.8 0.9 0.049 0.4)),
    join("\t", qw(1 10 rsDiff 0.001 0.8 0.9 0.8 0.5)),
    join("\t", qw(chrY bad rsBad 0.001 0.8 0.9 0.8 0.6)),
) . "\n";
gzip \$source => $input or die $GzipError;
open my $cfg, '>', $schema or die $!;
print {$cfg} encode_json({alias_map => {
    DS_ALL_P => 'GROUP2_P', MP2PRT_P => 'GROUP1_P',
}});
close $cfg;
my @cmd = ($^X, "$Bin/../DiffGWASDeps/prepare_sas_manhattan_plot_input.pl",
    '--input', $input, '--schema-config', $schema,
    '--pvars', 'DIFF_P DS_ALL_P MP2PRT_P META_P',
    '--output', $output, '--schema-out', $compact_schema,
    '--manifest', $manifest);
is(system(@cmd), 0, 'compact Manhattan input is prepared');
my $plain = '';
gunzip $output => \$plain or die $GunzipError;
my @rows = split /\n/, $plain;
is(shift @rows, "CHR\tBP\tDIFF_P\tDS_ALL_P\tMP2PRT_P\tMETA_P",
    'only coordinates and displayed P columns are uploaded');
is_deeply(\@rows, [
    "1\t10\t0.001\t0.9\t0.8\t0.8",
    "1\t50\t0.8\t0.9\t0.8\t0.049",
    "23\t200\t0.8\t0.9\t0.04\t0.9",
], 'any displayed P below 0.05 is retained and normalized coordinates are sorted');
my $stats = read_json($manifest);
is($stats->{rows_read}, 5, 'all source rows were scanned');
is($stats->{rows_written}, 3, 'nominal filter reduces the data');
is($stats->{rows_bad_coordinate}, 1, 'unusable coordinates are excluded');
is_deeply(read_json($compact_schema)->{wide_columns},
    [qw(CHR BP DIFF_P DS_ALL_P MP2PRT_P META_P)],
    'SAS import schema matches the compact table');
is(system(@cmd), 0, 'identical compact input is reusable');
done_testing();

sub read_json {
    open my $fh, '<', $_[0] or die $!;
    return decode_json(do { local $/; <$fh> });
}
