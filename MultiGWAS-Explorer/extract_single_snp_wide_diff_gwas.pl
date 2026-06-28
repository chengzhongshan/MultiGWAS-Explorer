#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Spec;

my $target = File::Spec->catfile($Bin, 'DiffGWASDeps', 'extract_single_snp_wide_diff_gwas.pl');
exec('perl', $target, @ARGV) or die "Failed to exec $target: $!";
