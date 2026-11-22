#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use File::Temp qw(tempdir);
use File::Spec;
use IPC::Open3;
use Symbol qw(gensym);
use JSON::PP qw(decode_json);
use Test::More;

my $root=tempdir(CLEANUP=>1);
my $log=File::Spec->catfile($root,'network.log');
open my $fh,'>',$log or die $!;
print {$fh} "ERROR: Windows error code: 67 in hx_disk_access\n";close $fh;
sub execute {
    my $error=gensym;
    my $pid=open3(undef,my $output,$error,$^X,@_);
    my $text=do{local$/;<$output>} // '';
    my $stderr=do{local$/;<$error>} // '';
    waitpid($pid,0);
    return ($?>>8,$text,$stderr);
}
for my $entry(qw(run_sas_codes_or_script_in_ODA.pl run_sas_codes_or_files_in_ODA.pl run_sas_codes_or_file_in_ODA.pl run_sas_codes_or_files_in_SAS_ODA.pl)) {
    my $path="$Bin/../$entry";
    my ($exit,$text,$error)=execute($path,'--help');
    is($exit,0,"$entry help succeeds");
    like($text,qr/--open-html/,'browser override documented');
    like($text,qr/--macro-dir/,'macro discovery documented');
    like($text,qr/--archive-transfers/,'batching control documented');
    ($exit,$text,$error)=execute($path,'--classify-sas-log',$log);
    is($exit,1,"$entry preserves failed-log exit");
    my $status=decode_json($text);
    is($status->{failure_class},'sas_windows_network_path_unavailable','alias preserves infrastructure classification');
    ok(!$status->{retryable},'classification cannot request a replay');
    ($exit,$text,$error)=execute($path,'--dry-run','--code','data _null_; run;',
        '--no-open-html','--no-archive-transfers','--macro-dir',File::Spec->catdir($root,'missing'),
        '--persistent','--session-id','alias_probe');
    is($exit,0,"$entry accepts new flags without SAS access");
    like($text,qr/DRY-RUN: persistent=1 session_id=alias_probe/,'session identity reaches canonical CLI');
    like($text,qr/using the current working directory/,'missing macro directory has visible fallback');
}
done_testing;
