#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../DiffGWASDeps";
use SAS_ODA_MCP;
use Cwd qw(abs_path);
use Test::More;
die "Usage: perl tests/run_live_mcp_oda.pl --live\n" unless @ARGV == 1 && $ARGV[0] eq '--live';
chdir abs_path("$Bin/..") or die $!;
local $ENV{OPEN_RESULT} = 0;
local $ENV{SAS_ODA_AUTOLOAD_MACROS} = 0;
for my $failure (0, 1) {
    my $result = SAS_ODA_MCP::run_sas_oda_tool({persistent => 0, run_timeout_seconds => 120,
        sas_codes_or_file => $failure ? 'proc print data=work.explorer_missing; run;'
            : 'proc means data=sashelp.class n mean; var age; run;'});
    my $text = $result->{content}[0]{text};
    my ($pid) = $text =~ /PID: (\d+)/;
    my ($dir) = $text =~ /Results directory: ([^\r\n]+)/;
    ok($pid && $dir, 'real ODA MCP worker started');
    waitpid($pid, 0);
    my $done = SAS_ODA_MCP::run_sas_oda_tool({pid => $pid});
    is(!!$done->{isError}, !!$failure, 'real ODA completion error flag');
    like($done->{content}[0]{text}, qr/Exit code: $failure\b/, 'real ODA exit code returned');
    ok(-s "$dir/output.html.info.txt", 'diagnostic log preserved');
    if (!$failure) {
        my @html = glob("$dir/sas_res_*.html");
        ok(@html && -s $html[0], 'real ODA HTML preserved');
    }
    diag("Real ODA MCP results: $dir");
}
done_testing();
