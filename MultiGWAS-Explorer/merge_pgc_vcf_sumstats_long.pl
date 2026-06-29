#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Spec;

my $target = File::Spec->catfile($Bin, 'DiffGWASDeps', 'merge_pgc_vcf_sumstats_long.pl');
exec('perl', $target, @ARGV) or die "Failed to exec $target: $!";
