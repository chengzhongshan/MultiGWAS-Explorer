#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my @inputs;
my $output;
GetOptions(
    'input=s@' => \@inputs,
    'output=s' => \$output,
) or die usage();
die usage() unless @inputs && defined $output;

open my $out, '>:raw', $output or die "Cannot write $output: $!\n";
my $expected_header;
my $rows = 0;
for my $path (@inputs) {
    open my $in, '<:raw', $path or die "Cannot read $path: $!\n";
    my $header = <$in>;
    die "Empty LD cache: $path\n" unless defined $header;
    $header =~ s/[\r\n]+\z//;
    if (!defined $expected_header) {
        $expected_header = $header;
        print {$out} "$header\n";
    }
    elsif ($header ne $expected_header) {
        die "LD cache header mismatch in $path\n";
    }
    while (my $line = <$in>) {
        print {$out} $line;
        ++$rows;
    }
    close $in;
}
close $out or die "Cannot close $output: $!\n";
print "OUTPUT\t$output\nROWS\t$rows\nINPUTS\t", scalar(@inputs), "\n";

sub usage {
    return "Usage: combine_haploreg_ld_caches.pl --input EUR.tsv --input ASN.tsv --output combined.tsv\n";
}
