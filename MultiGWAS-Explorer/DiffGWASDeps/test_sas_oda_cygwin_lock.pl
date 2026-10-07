#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Spec;
use Test::More;

my $wrapper = File::Spec->catfile($Bin, 'run_sas_oda_local_top_hits_manhattan_download_png.sh');
open my $fh, '<', $wrapper or die "Cannot read $wrapper: $!\n";
local $/;
my $source = <$fh>;
close $fh;

like($source, qr/\[\[ "\$\{pipeline_uname\}" != CYGWIN\* \]\].*command -v flock/s,
    'Cygwin explicitly bypasses descriptor-based flock');
like($source, qr/atomic directory cannot be\s*# inherited by SASPy\/Java/s,
    'the lock implementation documents why Cygwin uses mkdir');
like($source, qr/kill -0 "\$\{lock_owner\}".*rmdir "\$\{SAS_ODA_PIPELINE_LOCK_DIR\}"/s,
    'directory lock reclaims a stale owner PID');

done_testing();
