#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../DiffGWASDeps";
use SAS_ODA_MCP;
use File::Temp qw(tempdir);
use File::Spec;
use Cwd qw(getcwd);
use JSON::PP qw(decode_json);
use Test::More;

my $original=getcwd();
my $root=tempdir('explorer-mcp-wait-XXXX',TMPDIR=>1,CLEANUP=>1);
my $fake=File::Spec->catfile($root,'runner.pl');
open my $program,'>',$fake or die $!;
print {$program} <<'WORKER';
use strict; use warnings; use JSON::PP qw(encode_json); use Time::HiRes qw(sleep);
my @args=@ARGV;my $dir;
while(@ARGV){my $arg=shift;$dir=shift if $arg eq '--output-prefix';}
die 'No output prefix' unless $dir;
sleep($ENV{MCP_TEST_DELAY} || 0);
if ($ENV{MCP_TEST_GATE}) {
    for (1..500) { last if -f $ENV{MCP_TEST_GATE}; sleep(0.01); }
}
open my $log,'>',"$dir/output.html.info.txt" or die $!;
print {$log} "=== SAS Log ===\nNOTE: start\n",('x' x 20000),"\nNOTE: end\n";
close $log;
open my $capture,'>',"$dir/arguments.json" or die $!;
print {$capture} encode_json(\@args);close $capture;
exit($ENV{MCP_TEST_EXIT} || 0);
WORKER
close $program;
chdir $root or die $!;
END { chdir $original if defined $original; }
{
    no warnings 'redefine';
    local *SAS_ODA_MCP::runner_path=sub{$fake};
    for my $exit (0,1) {
        local $ENV{MCP_TEST_EXIT}=$exit;
        my $result=SAS_ODA_MCP::run_sas_oda_tool({sas_codes_or_file=>'%put READY;',persistent=>0,
            macro_dir=>'macros with spaces',open_html=>JSON::PP::false,dir4local_sas_home=>'local home',max_log_chars=>1000});
        is($result->{structuredContent}{status},$exit?'FAILED':'COMPLETE','default bounded wait returns verified terminal status');
        is($result->{structuredContent}{exit_code},$exit,'worker exit preserved');
        is(!!$result->{isError},!!$exit,'failure flag agrees with exit');
        like($result->{content}[0]{text},qr/Middle output omitted/,'large reply is excerpted');
        my $path=$result->{structuredContent}{log_path};
        open my $log,'<',$path or die $!;my $full=do{local$/;<$log>};close $log;
        cmp_ok(length($full),'>',20000,'full saved log remains available');
        open my $capture,'<',File::Spec->catfile($result->{structuredContent}{results_directory},'arguments.json') or die $!;
        my $args=decode_json(do{local$/;<$capture>});close $capture;
        ok(grep($_ eq '--no-open-html',@$args),'explicit browser policy reaches CLI');
        ok(grep($_ eq 'macros with spaces',@$args),'macro directory is passed without shell splitting');
        ok(grep($_ eq 'local home',@$args),'local SAS home is passed intact');
    }
    local $ENV{MCP_TEST_EXIT}=0;
    local $ENV{MCP_TEST_GATE}=File::Spec->catfile($root,'release-worker');
    my $running=SAS_ODA_MCP::run_sas_oda_tool({sas_codes_or_file=>'%put LATER;',persistent=>0,wait_seconds=>0});
    ok(!-e $ENV{MCP_TEST_GATE},'explicit background mode returns before the blocked worker is released');
    is($running->{structuredContent}{status},'RUNNING','accepted background job is structured');
    open my $gate,'>',$ENV{MCP_TEST_GATE} or die $!;close $gate;
    my $done=SAS_ODA_MCP::run_sas_oda_tool({pid=>$running->{structuredContent}{pid},wait_seconds=>3,max_log_chars=>0});
    is($done->{structuredContent}{status},'COMPLETE','poll waits for the same accepted job');
    unlike($done->{content}[0]{text},qr/Middle output omitted/,'zero log limit requests full content');
    cmp_ok(length($done->{content}[0]{text}),'>',20000,'full response includes complete log');
    opendir my $folder,$root or die $!;
    my @before=sort readdir $folder;closedir $folder;
    for my $invalid ({wait_seconds=>31},{wait_seconds=>-1},{max_log_chars=>-1},{max_log_chars=>1.5}) {
        my $failure=SAS_ODA_MCP::run_sas_oda_tool($invalid);
        ok($failure->{isError} && !$failure->{structuredContent}{submitted},'invalid delivery option cannot submit a job');
    }
    opendir $folder,$root or die $!;my @after=sort readdir $folder;closedir $folder;
    is_deeply(\@after,\@before,'invalid options create no scratch or results directories');
}
done_testing;
