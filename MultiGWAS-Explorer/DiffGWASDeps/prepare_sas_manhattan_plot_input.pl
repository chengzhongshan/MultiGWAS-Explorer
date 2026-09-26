#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use JSON::PP qw(decode_json encode_json);
use IO::Uncompress::Gunzip qw($GunzipError);
use IO::Compress::Gzip qw($GzipError);
use File::Basename qw(dirname);
use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);

my ($input, $schema_config, $pvars, $output, $schema_out, $manifest) = ('') x 6;
my $threshold = 0.05;
my $all_snps = 0;
GetOptions(
    'input=s'         => \$input,
    'schema-config=s' => \$schema_config,
    'pvars=s'         => \$pvars,
    'output=s'        => \$output,
    'schema-out=s'    => \$schema_out,
    'manifest=s'      => \$manifest,
    'threshold=f'     => \$threshold,
    'all-snps!'       => \$all_snps,
) or die "Invalid Manhattan subset options\n";
die "--input, --pvars, --output, --schema-out, and --manifest are required\n"
    unless length($input) && length($pvars) && length($output)
        && length($schema_out) && length($manifest);
die "Input not found: $input\n" unless -s $input;
die "Threshold must be between 0 and 1\n" unless $threshold > 0 && $threshold <= 1;

my @pvars = grep { length } split /[\s,]+/, $pvars;
die "At least one P column is required\n" unless @pvars;
for my $pvar (@pvars) {
    die "Invalid P column name: $pvar\n" unless $pvar =~ /^[A-Za-z_][A-Za-z0-9_]*$/;
}
my %seen;
@pvars = grep { !$seen{$_}++ } @pvars;

my $cfg = {};
if (length $schema_config) {
    open my $fh, '<', $schema_config or die "Cannot read $schema_config: $!\n";
    $cfg = decode_json(do { local $/; <$fh> });
    close $fh;
}
my %alias = (
    %{ $cfg->{alias_map} || {} },
    %{ $cfg->{post_alias_map} || {} },
);
my @source_stat = (stat($input))[7, 9];
my $cache_key = JSON::PP->new->canonical->encode({
    input => File::Spec->rel2abs($input),
    source_size => $source_stat[0],
    source_mtime => $source_stat[1],
    pvars => \@pvars,
    alias => \%alias,
    threshold => $threshold,
    all_snps => $all_snps ? 1 : 0,
});
if (-s $output && -s $schema_out && -s $manifest) {
    my $old = eval {
        open my $fh, '<', $manifest or die $!;
        my $value = decode_json(do { local $/; <$fh> });
        close $fh;
        $value;
    };
    if ($old && ($old->{cache_key} || '') eq $cache_key) {
        print "[skip] Reusing compact SAS Manhattan input: $output ($old->{rows_written} rows)\n";
        exit 0;
    }
}

my $in = IO::Uncompress::Gunzip->new($input, MultiStream => 1)
    or die "Cannot open $input: $GunzipError\n";
my $header = <$in>;
die "Input is empty: $input\n" unless defined $header;
$header =~ s/[\r\n]+$//;
$header =~ s/^#//;
my @header = split /\t/, $header, -1;
my %index;
for my $i (0 .. $#header) { $index{$header[$i]} = $i; }
for my $required (qw(CHR BP)) {
    die "Required column $required is missing from $input\n" unless exists $index{$required};
}
my @physical_pvars;
for my $pvar (@pvars) {
    my $physical = exists $index{$pvar} ? $pvar : ($alias{$pvar} || '');
    die "Displayed P column $pvar is absent from $input and its schema aliases\n"
        unless length($physical) && exists $index{$physical};
    push @physical_pvars, $physical;
}

for my $path ($output, $schema_out, $manifest) {
    make_path(dirname($path)) unless -d dirname($path);
}
my $tmpdir = tempdir('sas_manhattan_XXXX', DIR => dirname($output), CLEANUP => 1);
my $unsorted = File::Spec->catfile($tmpdir, 'unsorted.tsv');
my $sorted = File::Spec->catfile($tmpdir, 'sorted.tsv');
my $tmp_gz = File::Spec->catfile($tmpdir, 'subset.tsv.gz');
open my $raw, '>', $unsorted or die "Cannot write $unsorted: $!\n";
my ($rows_read, $rows_written, $bad_coord) = (0, 0, 0);
while (my $line = <$in>) {
    $rows_read++;
    $line =~ s/[\r\n]+$//;
    my @values = split /\t/, $line, -1;
    my $chr = $values[$index{CHR}] // '';
    $chr =~ s/^chr//i;
    $chr = 23 if uc($chr) eq 'X';
    $chr = 24 if uc($chr) eq 'Y';
    my $bp = $values[$index{BP}] // '';
    unless ($chr =~ /^\d+$/ && $chr >= 1 && $chr <= 24
        && $bp =~ /^\d+$/ && $bp > 0) {
        $bad_coord++;
        next;
    }
    my @p = map { $values[$index{$_}] // '' } @physical_pvars;
    my $keep = $all_snps ? 1 : 0;
    if (!$all_snps) {
        for my $value (@p) {
            next unless $value =~ /^(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;
            if ($value < $threshold) { $keep = 1; last; }
        }
    }
    next unless $keep;
    print {$raw} join("\t", 0 + $chr, 0 + $bp, @p), "\n";
    $rows_written++;
}
close $in or die "Cannot finish reading $input: $GunzipError\n";
close $raw or die "Cannot finish writing $unsorted: $!\n";
die $all_snps ? "No coordinate-valid rows in $input\n"
    : "No rows pass P < $threshold for the displayed tracks\n"
    unless $rows_written;

local $ENV{LC_ALL} = 'C';
my $sort_bin = $^O eq 'cygwin' ? '/usr/bin/sort' : 'sort';
system($sort_bin, '-s', '-t', "\t", '-k1,1n', '-k2,2n',
    '-o', $sorted, $unsorted) == 0 or die "Could not sort compact Manhattan rows\n";
my $out = IO::Compress::Gzip->new($tmp_gz)
    or die "Cannot write $tmp_gz: $GzipError\n";
print {$out} join("\t", 'CHR', 'BP', @pvars), "\n";
open my $sorted_fh, '<', $sorted or die "Cannot read $sorted: $!\n";
while (my $line = <$sorted_fh>) { print {$out} $line; }
close $sorted_fh;
close $out or die "Cannot finish $tmp_gz: $GzipError\n";

my $schema = { wide_columns => ['CHR', 'BP', @pvars], char_lengths => {} };
my $report = {
    cache_key => $cache_key,
    source => $input,
    output => $output,
    threshold => $threshold,
    all_snps => $all_snps ? 1 : 0,
    displayed_pvars => \@pvars,
    physical_pvars => \@physical_pvars,
    rows_read => $rows_read,
    rows_written => $rows_written,
    rows_removed => $rows_read - $rows_written,
    rows_bad_coordinate => $bad_coord,
};
write_json($schema_out, $schema);
replace_file($tmp_gz, $output);
write_json($manifest, $report);
print "Compact SAS Manhattan input: $output\n";
print $all_snps
    ? "Rows read: $rows_read; kept all coordinate-valid SNPs: $rows_written; bad coordinates: $bad_coord\n"
    : "Rows read: $rows_read; kept at P < $threshold: $rows_written; bad coordinates: $bad_coord\n";
print "Columns: CHR BP ", join(' ', @pvars), "\n";

sub write_json {
    my ($path, $value) = @_;
    my $tmp = "$path.tmp.$$";
    open my $fh, '>', $tmp or die "Cannot write $tmp: $!\n";
    print {$fh} JSON::PP->new->canonical->encode($value), "\n";
    close $fh or die "Cannot close $tmp: $!\n";
    replace_file($tmp, $path);
}

sub replace_file {
    my ($from, $to) = @_;
    unlink $to if -e $to;
    rename $from, $to or die "Cannot replace $to: $!\n";
}
