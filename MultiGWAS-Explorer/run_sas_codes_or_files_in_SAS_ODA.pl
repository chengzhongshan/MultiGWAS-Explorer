#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Spec;

# Compatibility entry point; keep all SAS ODA logic in the canonical runner.
my $runner = File::Spec->catfile($Bin, 'run_sas_codes_or_script_in_ODA.pl');
die "Cannot find SAS ODA runner: $runner\n" unless -f $runner;
exec { $^X } $^X, $runner, @ARGV;
die "Cannot start SAS ODA runner: $!\n";
