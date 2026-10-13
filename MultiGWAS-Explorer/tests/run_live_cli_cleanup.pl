#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use Cwd qw(abs_path);
use File::Temp qw(tempdir);
use File::Path qw(make_path);
use Test::More;
use Getopt::Long qw(GetOptions);
use JSON::PP qw(decode_json);

my ($live, $only, $timeout_seconds, $runner_name) = (0, undef, 3, 'run_sas_codes_or_script_in_ODA.pl');
GetOptions('live' => \$live, 'only=s' => \$only, 'timeout-seconds=i' => \$timeout_seconds,
    'runner=s' => \$runner_name) or die "Invalid options\n";
die "Usage: perl tests/run_live_cli_cleanup.pl --live [--only timeout]\n" unless $live;
die "Unknown case\n" if defined($only) && $only !~ /^(success|failure|timeout|retained)$/;
die "Timeout must be positive\n" unless $timeout_seconds > 0;
die "Unsupported runner name\n" unless $runner_name eq 'run_sas_codes_or_script_in_ODA.pl'
    || $runner_name eq 'run_sas_codes_or_files_in_SAS_ODA.pl';
my $root = abs_path("$Bin/..");
my $out = tempdir('SAS_Output_cli_cleanup_XXXXXX', DIR => $root, CLEANUP => 0);
local $ENV{SAS_ODA_AUTOLOAD_MACROS} = 0;
local $ENV{OPEN_RESULT} = 0;
local $ENV{SAS_ODA_WORKER_TEMP_DIR}; delete $ENV{SAS_ODA_WORKER_TEMP_DIR};
for my $case (qw(success failure timeout retained)) {
    next if defined($only) && $only ne $case;
    my $dir = "$out/$case"; make_path($dir);
    my $pid = fork(); die $! unless defined $pid;
    if (!$pid) {
        chdir $dir or die $!;
        open STDOUT, '>', 'runner.log' or die $!;
        open STDERR, '>&STDOUT' or die $!;
        my $code = $case eq 'failure' ? 'proc print data=work.cleanup_missing; run;'
            : $case eq 'timeout' ? 'data _null_; rc=sleep(30,1); run;' : 'proc print data=sashelp.class(obs=3); run;';
        my @args = ('--code', $code, '--output-prefix', 'result', '--run-timeout-seconds',
            $case eq 'timeout' ? $timeout_seconds : 120, '--run-timeout-grace-seconds', 1);
        push @args, '--keep-intermediate-files' if $case eq 'retained';
        exec $^X, "$root/$runner_name", @args;
        die $!;
    }
    waitpid($pid, 0); my $exit = $? >> 8;
    ok(($case eq 'failure' || $case eq 'timeout') ? $exit != 0 : $exit == 0, "$case expected exit");
    opendir my $dh, $dir or die $!;
    my @scratch = grep { /^\.sas_oda_work_/ } readdir $dh; closedir $dh;
    is(scalar @scratch, $case eq 'retained' ? 1 : 0, "$case scratch lifecycle");
    if ($case eq 'success' || $case eq 'retained') {
        my @html = glob("$dir/result/sas_res_*.html");
        ok(@html && -s $html[0], "$case final HTML preserved");
        ok(-s "$dir/result/output.html.info.txt", "$case diagnostic log preserved");
        open my $sf, '<', "$dir/result/output.run.status.json" or die $!;
        my $status = do { local $/; decode_json(<$sf>) }; close $sf;
        ok($status->{htmlfilename} && -s $status->{htmlfilename}, "$case status points to retained HTML");
    }
    if ($case eq 'retained') {
        is(system($^X, "$root/cleanup_sas_oda_intermediates.pl", '--directory', $dir,
            '--older-than-hours', 0, '--apply'), 0, 'retained test scratch cleaned by maintenance command');
        ok(!-d "$dir/$scratch[0]", 'maintenance removed only this test scratch');
    }
}
diag("CLI cleanup artifacts: $out");
done_testing();
