#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Spec;

my $target = File::Spec->catfile($Bin, 'DiffGWASDeps', 'emit_diff_gwas_runner_env.pl');
exec('perl', $target, @ARGV) or die "Failed to exec $target: $!";
