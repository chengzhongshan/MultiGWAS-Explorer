#!/usr/bin/env perl
use strict;
use warnings;
my $file = shift or die "Usage: perl extract_sas_figure_path.pl SAS-run.log\n";
open my $fh, '<', $file or die "Cannot read $file: $!\n";
my $remaining = 0;
while (my $line = <$fh>) {
    if ($line =~ s/^.*The final figure is put here:\s*//) { $remaining = 40; }
    next unless $remaining-- > 0;
    $line =~ s/^\s+|\s+$//g;
    $line =~ s/^['"](.*)['"]$/$1/;
    if ($line =~ m{^(?:[A-Za-z]:[\\/]|/|~/).+\.(?:png|jpg|jpeg|svg)$}i) {
        print $line;
        exit 0;
    }
}
exit 1;
