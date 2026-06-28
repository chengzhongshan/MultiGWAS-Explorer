#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $template = '';
my $output = '';
my @replacements;
my @file_replacements;

GetOptions(
    'template=s'     => \$template,
    'output=s'       => \$output,
    'replace=s@'     => \@replacements,
    'replace-file=s@'=> \@file_replacements,
) or die usage();

die "--template is required\n" unless length $template;
die "--output is required\n" unless length $output;

my %map;
for my $entry (@replacements) {
    my ($key, $value) = split /=/, $entry, 2;
    die "Invalid --replace entry: $entry\n" unless defined $key;
    $map{$key} = defined $value ? $value : '';
}
for my $entry (@file_replacements) {
    my ($key, $path) = split /=/, $entry, 2;
    die "Invalid --replace-file entry: $entry\n" unless defined $key && defined $path;
    open my $fh, '<', $path or die "Cannot read replacement file $path: $!\n";
    local $/;
    $map{$key} = <$fh>;
    close $fh;
}

open my $in, '<', $template or die "Cannot read template $template: $!\n";
local $/;
my $text = <$in>;
close $in;

for my $key (keys %map) {
    my $token = "__${key}__";
    my $value = defined $map{$key} ? $map{$key} : '';
    $text =~ s/\Q$token\E/$value/g;
}

open my $out, '>', $output or die "Cannot write output $output: $!\n";
print {$out} $text;
close $out or die "Cannot close output $output: $!\n";

sub usage {
    return <<"USAGE";
Usage:
  perl render_sas_template.pl --template in.sas --output out.sas [--replace KEY=VALUE] [--replace-file KEY=file]
USAGE
}
