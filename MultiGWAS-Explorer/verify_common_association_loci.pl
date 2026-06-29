#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);

my $target = "$Bin/DiffGWASDeps/verify_common_association_loci.pl";
do $target or die "Unable to load $target: $@\n$!\n";
