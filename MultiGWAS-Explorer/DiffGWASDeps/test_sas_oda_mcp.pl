#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use lib $Bin;
use SAS_ODA_MCP;
use File::Temp qw(tempdir);
use Cwd qw(getcwd);
use Test::More;
use MCP::Tool;

open my $server_file, '<', "$Bin/../server.pl" or die $!;
my $server_source = do { local $/; <$server_file> }; close $server_file;
my ($callback_source) = $server_source =~ /code => (sub \([^\n]+\) \{ return SAS_ODA_MCP::run_sas_oda_tool\(\$args\) \})/;
die 'Missing SAS ODA callback' unless $callback_source;
my $callback = eval "use feature 'signatures'; $callback_source";
die $@ if $@;
{
    no warnings 'redefine';
    my $received;
    local *SAS_ODA_MCP::run_sas_oda_tool = sub {
        $received = $_[0];
        return {content => [{type => 'text', text => 'SDK dispatch accepted'}]};
    };
    my $tool = MCP::Tool->new(code => $callback);
    my $args = {sas_codes_or_file => 'proc print; run;', persistent => 'false'};
    my $result = $tool->call($args, {});
    is($received, $args, 'installed MCP SDK forwards the actual argument hash');
    is($result->{content}[0]{text}, 'SDK dispatch accepted', 'MCP callback dispatch succeeds');
}

my $original = getcwd();
my $root = tempdir(CLEANUP => 1);
my $fake = "$root/fake_runner.pl";
open my $fh, '>', $fake or die $!;
print {$fh} <<'RUNNER';
use strict; use warnings;
my $dir;
while (@ARGV) { my $arg=shift; $dir=shift if $arg eq '--output-prefix' }
die 'Missing directory' unless $dir;
open my $log, '>', "$dir/output.html.info.txt" or die $!;
print {$log} "=== SAS Log ===\n";
print {$log} $ENV{ODA_TEST_FAIL} ? "ERROR: Intentional missing table\n" : "NOTE: Result ready\n";
close $log;
open my $html, '>', "$dir/sas_res_test.html" or die $!;
print {$html} '<html>Preserved result</html>'; close $html;
exit($ENV{ODA_TEST_FAIL} ? 1 : 0);
RUNNER
close $fh;
chdir $root or die $!;
{
    no warnings 'redefine';
    local *SAS_ODA_MCP::runner_path = sub { $fake };
    for my $failure (0, 1) {
        local $ENV{ODA_TEST_FAIL} = $failure;
        my $started = SAS_ODA_MCP::run_sas_oda_tool({
            sas_codes_or_file => 'proc print; run;', persistent => 0,
            output_file => "custom_$failure.txt",
        });
        my $text = $started->{content}[0]{text};
        my ($pid) = $text =~ /PID: (\d+)/;
        my ($dir) = $text =~ /Results directory: ([^\r\n]+)/;
        ok($pid && $dir, 'job starts with PID and retained result directory');
        waitpid($pid, 0);
        my $done = SAS_ODA_MCP::run_sas_oda_tool({pid => $pid});
        is(!!$done->{isError}, !!$failure, 'MCP error flag follows worker exit');
        like($done->{content}[0]{text}, qr/Exit code: $failure\b/, 'exit code returned');
        like($done->{content}[0]{text}, qr/=== SAS Log ===/, 'full log returned');
        ok(-f "custom_$failure.txt", 'custom log copied');
        ok(-f "$dir/sas_res_test.html", 'HTML artifact retained');
        ok(!-f "$dir/output.html.info.pid", 'completed PID marker removed');
        opendir my $dh, $root or die $!;
        my @scratch = grep { /^tmp[A-Za-z0-9_]+$/ } readdir $dh; closedir $dh;
        is(scalar @scratch, 0, 'worker scratch removed');
    }
}
chdir $original or die $!;
done_testing();
