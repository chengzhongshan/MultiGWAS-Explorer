#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Spec;

my $target = File::Spec->catfile($Bin, 'DiffGWASDeps', 'standardize_diff_gwas_zscore.pl');
exec('perl', $target, @ARGV) or die "Failed to exec $target: $!";
