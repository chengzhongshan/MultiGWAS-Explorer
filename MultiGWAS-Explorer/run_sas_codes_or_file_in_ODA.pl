#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);

# Compatibility spelling; keep execution and configuration in the main runner.
exec $^X, "$Bin/run_sas_codes_or_script_in_ODA.pl", @ARGV;
die "Cannot start SAS ODA runner: $!\n";
