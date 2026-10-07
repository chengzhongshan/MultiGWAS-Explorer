#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Spec;
use Test::More;

my $helper = File::Spec->catfile($Bin, '..', 'run_sas_codes_or_script_in_ODA.pl');
open my $fh, '<', $helper or die "Cannot read $helper: $!\n";
local $/;
my $source = <$fh>;
close $fh;

like(
    $source,
    qr/File::Spec->catfile\(getcwd\(\), 'DiffGWASDeps', \$base\)/,
    'include preflight checks the project DiffGWASDeps directory directly',
);
like(
    $source,
    qr/version-controlled source lives in DiffGWASDeps.*before considering a recursive search/s,
    'direct macro lookup is documented as preceding recursive fallback',
);

done_testing();
