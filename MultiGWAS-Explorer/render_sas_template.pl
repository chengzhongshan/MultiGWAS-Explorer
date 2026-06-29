#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Spec;

my $target = File::Spec->catfile($Bin, 'DiffGWASDeps', 'render_sas_template.pl');
exec('perl', $target, @ARGV) or die "Failed to exec $target: $!";
