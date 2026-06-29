#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Spec;

my $target = File::Spec->catfile($Bin, 'DiffGWASDeps', 'generate_sas_wide_import_include.pl');
exec('perl', $target, @ARGV) or die "Failed to exec $target: $!";
