#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../DiffGWASDeps";
use SAS_ODA_Runner;
use File::Spec;
use File::Path qw(make_path);
use File::Temp qw(tempdir);
use Cwd qw(abs_path);
use Getopt::Long qw(GetOptions);
use Digest::SHA qw(sha256_hex);
use IO::Socket::INET;
use JSON::PP qw(encode_json);
use Test::More;

my ($live,$cfgfile,$outdir);
GetOptions('live'=>\$live,'saspy-cfgfile=s'=>\$cfgfile,'outdir=s'=>\$outdir) or die 'Invalid test option';
plan skip_all=>'Use --live with configured SAS ODA credentials' unless $live;
my $root=abs_path("$Bin/..");
die "Choose a new --outdir; existing results must not be reused\n" if defined($outdir) && -e $outdir;
$outdir //= tempdir('SAS_Output_RandID_runtime_live_XXXXXX',DIR=>$root,CLEANUP=>0);
make_path($outdir) unless -d $outdir;
$outdir=abs_path($outdir);
sub read_bytes {open my $fh,'<:raw',$_[0] or die $!;local$/;return <$fh>;}
sub write_bytes {open my $fh,'>:raw',$_[0] or die $!;print {$fh} $_[1];close $fh or die $!;}
my $socket=IO::Socket::INET->new(LocalAddr=>'127.0.0.1',LocalPort=>0,Listen=>1,Proto=>'tcp') or die $!;
my $port=$socket->sockport;close $socket;
my $tag='mgx'.time.'_'. $$;
my $session=$tag.'_runtime';
my $used=0;
my $owner=$$;
local $ENV{SAS_ODA_SESSION_PORT}=$port;
local $ENV{SASPY_CFGNAME}='oda';
local $ENV{SASPY_CFGFILE}=$cfgfile if defined $cfgfile;
local $ENV{OPEN_RESULT}=0;
local $ENV{SAS_ODA_AUTOLOAD_MACROS}=0;
local $ENV{PERLMCP_SASPY_LOAD_MACROS}=0;
local $ENV{SAS_ODA_SESSION_DEBUG_LOG}=SAS_ODA_Runner::_normalize_local_path_for_python("$outdir/session.log",$ENV{PIPELINE_PYTHON_BIN} || 'python3');
my @common=('--persistent','--session-id',$session,'--saspy-cfgname','oda','--no-open-html','--run-timeout-seconds',180);
push @common,('--saspy-cfgfile',$cfgfile) if defined $cfgfile;
my @records;
sub action {
    my ($label,$entry,@args)=@_;
    my $prefix="$outdir/$label";
    my $child=fork();die $! unless defined $child;
    if(!$child){
        open STDOUT,'>',"$prefix.console.log" or die $!;
        open STDERR,'>&STDOUT' or die $!;
        exec {$^X} $^X,"$root/$entry",@common,'--output-prefix',$prefix,@args;
        die $!;
    }
    $used=1;waitpid($child,0);
    my $exit=$? & 127 ? 128+($? & 127) : $?>>8;
    push @records,{phase=>$label,entry=>$entry,exit=>$exit};
    return ($exit,(-f "$prefix/output.html.info.txt" ? read_bytes("$prefix/output.html.info.txt") : ''));
}
sub close_owned {
    return unless $used;
    local $ENV{SAS_ODA_SESSION_PORT}=$port;
    local $ENV{OPEN_RESULT}=0;
    local $ENV{SAS_ODA_AUTOLOAD_MACROS}=0;
    eval {action('close','run_sas_codes_or_files_in_ODA.pl','--code','endsas;');};
    eval {
        my $connection=IO::Socket::INET->new(PeerAddr=>'127.0.0.1',PeerPort=>$port,Proto=>'tcp',Timeout=>5) or die $!;
        my $body=encode_json({cmd=>'shutdown'});
        print {$connection} pack('N2',0,length($body)),$body;
        my $response;read($connection,$response,4096);close $connection;
    };
    $used=0;
}
END {close_owned() if defined($owner) && $$==$owner;}

my @remote=map {"${tag}_$_.".(/^[ab]$/?'tsv':'sas')} qw(a b parent child);
write_bytes("$outdir/$remote[0]","x\n11\n");
write_bytes("$outdir/$remote[1]","x\n22\n");
my $program=<<"SAS";
data work.port_a; infile "~/$remote[0]" firstobs=2; input x; run;
data work.port_b; infile "~/$remote[1]" firstobs=2; input x; run;
data work.port_values; set work.port_a work.port_b; run;
proc sql noprint; select count(*),sum(x) into :port_n trimmed,:port_sum trimmed from work.port_values; quit;
%macro check_port; %if &port_n ne 2 or &port_sum ne 33 %then %do; %put ERROR: Invalid port values; %abort cancel; %end; %mend; %check_port;
proc print data=work.port_values; run;
%put MGX_PORT_PASS=transfer;
%put MGX_PORT_PID=&sysjobid;
SAS
write_bytes("$outdir/transfer.sas",$program);
my ($exit,$log)=action('transfer','run_sas_codes_or_files_in_ODA.pl','--file',"$outdir/transfer.sas",
    '--upload-file',"$outdir/$remote[0]",'--upload-file',"$outdir/$remote[1]",
    '--download-file',"~/$remote[0]",'--download-local-path',"$outdir/returned-a.tsv",
    '--download-file',"~/$remote[1]",'--download-local-path',"$outdir/returned-b.tsv");
is($exit,0,'actual batched upload, SAS submit and download succeeds');
like($log,qr/^MGX_PORT_PASS=transfer\s*$/m,'runtime calculation verifies two rows and sum 33');
unlike($log,qr/^\s*ERROR(?:\s+\d+-\d+)?:/m,'no SAS error in retained full log');
my ($pid)=$log=~/^MGX_PORT_PID=(\S+)/m;
for my $i(0,1){
    my $returned="$outdir/returned-".($i?'b':'a').'.tsv';
    ok(-f $returned,'downloaded file exists');
    is(-f $returned?sha256_hex(read_bytes($returned)):'missing',sha256_hex(read_bytes("$outdir/$remote[$i]")),'round-trip bytes match');
}
my @html=glob("$outdir/transfer/sas_res_*.html");
ok(@html && -s $html[0],'actual PRINT HTML is retained');
($exit,$log)=action('reuse','run_sas_codes_or_file_in_ODA.pl','--code',q{proc sql noprint; select sum(x) into :again trimmed from work.port_values; quit; %macro check_reuse; %if &again ne 33 %then %do; %put ERROR: WORK was not retained; %abort cancel; %end; %mend; %check_reuse; %put MGX_PORT_PASS=reuse; %put MGX_PORT_PID=&sysjobid;});
is($exit,0,'another alias reuses the same persistent WORK');
like($log,qr/^MGX_PORT_PASS=reuse\s*$/m,'reused WORK has the expected values');
my ($again_pid)=$log=~/^MGX_PORT_PID=(\S+)/m;
ok(defined($pid) && defined($again_pid) && $pid eq $again_pid,'SAS process identity is unchanged');
write_bytes("$outdir/$remote[2]","%include \"$remote[3]\";\n");
write_bytes("$outdir/$remote[3]", "%macro ${tag}_calc; proc sql noprint; select sum(x) into :nested trimmed from work.port_values; quit; %if &nested ne 33 %then %do; %put ERROR: Nested calculation failed; %abort cancel; %end; %put MGX_PORT_PASS=nested; %mend;\n");
write_bytes("$outdir/nested.sas","%include \"$outdir/$remote[2]\";\n%${tag}_calc;\n");
($exit,$log)=action('nested','run_sas_codes_or_files_in_SAS_ODA.pl','--file',"$outdir/nested.sas",'--macro-dir',$outdir);
is($exit,0,'recursive local include/macro staging succeeds');
like($log,qr/^MGX_PORT_PASS=nested\s*$/m,'staged child macro executes with expected data');
unlike($log,qr/^\s*ERROR(?:\s+\d+-\d+)?:/m,'nested full SAS log contains no runtime errors');
my @cleanup;push @cleanup,('--delete-file',"~/$_") for @remote;
my @metadata;push @metadata,('--file-info',"~/$_") for @remote;
($exit,$log)=action('cleanup','run_sas_codes_or_script_in_ODA.pl',@cleanup,@metadata);
is($exit,0,'only this test’s exact uploaded files are deleted');
my $cleanup=read_bytes("$outdir/cleanup.console.log");
is(scalar(()=$cleanup=~/^EXISTS\t0\s*$/mg),4,'all four remote test files are verified absent');
if(-f "$outdir/session.log"){
    my $debug=read_bytes("$outdir/session.log");
    cmp_ok(scalar(()=$debug=~/request cmd=bulk_transfer session_id=\Q$session\E/mg),'>=',2,'actual persistent bulk RPCs are recorded');
}else{fail('dedicated session debug log is retained');}
close_owned();
write_bytes("$outdir/execution.json",JSON::PP->new->canonical->pretty->encode({session=>$session,records=>\@records,scope=>'Synthetic runner integration; not a scientific pipeline or performance benchmark.'}));
diag("Retained execution evidence: $outdir");
done_testing;
