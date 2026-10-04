#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use IO::Uncompress::Gunzip qw($GunzipError);
use IO::Compress::Gzip qw($GzipError);

my ($input, $output) = ('', '');
GetOptions('input=s' => \$input, 'output=s' => \$output)
    or die "Usage: $0 --input wide.tsv.gz --output sas_wide.tsv.gz\n";
die "--input and --output must name different files\n"
    unless length($input) && length($output) && $input ne $output;

my $in = IO::Uncompress::Gunzip->new($input, MultiStream => 1)
    or die "Cannot read $input: $GunzipError\n";
my $out = IO::Compress::Gzip->new($output)
    or die "Cannot write $output: $GzipError\n";
my $header = <$in> // die "Missing wide-table header in $input\n";
print {$out} $header or die "Cannot write header to $output: $GzipError\n";
(my $header_text = $header) =~ s/[\r\n]+\z//;
$header_text =~ s/^#//;
my @columns = split /\t/, $header_text, -1;
my ($chr_index) = grep { uc($columns[$_]) eq 'CHR' } 0 .. $#columns;
die "Missing CHR column in $input\n" unless defined $chr_index;

my ($rows, $converted) = (0, 0);
while (my $line = <$in>) {
    $line =~ s/[\r\n]+\z//;
    my @fields = split /\t/, $line, -1;
    die "Wide row has no CHR field in $input at line " . ($rows + 2) . "\n"
        unless @fields > $chr_index;
    if ($fields[$chr_index] =~ /^(?:chr)?X\z/i) {
        $fields[$chr_index] = 23;
        ++$converted;
    }
    print {$out} join("\t", @fields), "\n"
        or die "Cannot write row to $output: $GzipError\n";
    ++$rows;
}
close $in or die "Cannot close $input: $GunzipError\n";
close $out or die "Cannot close $output: $GzipError\n";
print "Rows: $rows; X-to-23: $converted; SAS input: $output\n";
