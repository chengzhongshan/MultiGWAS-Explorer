#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use Test::More;
sub source {
    open my $fh, '<:raw', "$Bin/../$_[0]" or die $!;
    return do { local $/; <$fh> };
}
like(source('run_sas_codes_or_script_in_ODA.pl'),
    qr/if \(SAS_Source_Analysis::include_target_needs_runtime_resolution\(\$target\)\) \{.*?next;.*?my \$remote_target = normalize_remote_include_path/s,
    'canonical runner defers dynamic includes before remote or recursive lookup');
for my $alias (qw(run_sas_codes_or_file_in_ODA.pl run_sas_codes_or_files_in_ODA.pl run_sas_codes_or_files_in_SAS_ODA.pl)) {
    like(source($alias), qr/exec.*?\$\^X.*?run_sas_codes_or_script_in_ODA\.pl|run_sas_codes_or_script_in_ODA\.pl.*?exec/s,
        "$alias delegates to the canonical runner");
    unlike(source($alias), qr/sub run_include_preflight/,
        "$alias does not carry a stale second preflight implementation");
}
like(source('DiffGWASDeps/SAS_ODA_MCP.pm'),
    qr/is_local_file_candidate\(\$sas_input\).*?-f \$sas_input/s,
    'MCP dispatch guards its actual source-versus-file probe');
done_testing();
