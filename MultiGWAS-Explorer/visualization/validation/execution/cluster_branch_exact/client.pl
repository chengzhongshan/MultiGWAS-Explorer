#!/usr/bin/env perl
use strict;
use warnings;
use FindBin;
use Cwd qw(abs_path);
use File::Path qw(make_path);
use File::Basename qw(basename dirname);
use File::Copy qw(copy);
use JSON::PP;
use Digest::SHA qw(sha256_hex);
use LWP::UserAgent;
use POSIX qw(strftime);

# Reader client: no .audit-local, private configuration or author's paths.
chdir "$FindBin::Bin/.." or die $!;
sub bytes { open my $f,'<:raw',$_[0] or die "$_[0]: $!"; local $/; <$f> }
sub put { my($p,$s)=@_; make_path(dirname($p)); open my $f,'>:raw',$p or die "$p: $!"; print $f $s }
sub json { decode_json(bytes($_[0])) }
sub save { put($_[0],JSON::PP->new->canonical->pretty->encode($_[1])) }
sub stamp { strftime('%Y-%m-%dT%H:%M:%SZ',gmtime) }
sub own_file {
    my($p)=@_; my $root=abs_path('companion'); my $full=abs_path($p);
    die "Input must be an existing file inside companion/: $p\n"
        unless defined($full) && -f $full && index($full,"$root/")==0;
    return $full;
}
my($mode,$name,$variant)=@ARGV;
$mode //= 'help'; $variant //= 'baseline';
my $catalog_path=$ENV{BOOK_LESSON_CATALOG}//'companion/reader/lessons.json';
own_file($catalog_path);
my $catalog=json($catalog_path);
if($mode eq 'list') {
    print "$_\t$catalog->{$_}{title}\n" for sort keys %$catalog; exit;
}
if($mode eq 'help') {
    print "perl scripts/run_visualization.pl list\n",
          "perl scripts/run_visualization.pl discover\n",
          "perl scripts/run_visualization.pl run cluster_branch_exact [baseline|exercise]\n",
          "perl scripts/run_visualization.pl verify .runs/RUN_DIRECTORY\n",
          "BOOK_MCP_URL defaults to http://127.0.0.1:18789/mcp\n",
          "BOOK_LESSON_CATALOG optionally selects a catalog inside companion/\n";
    exit;
}
if($mode eq 'verify') {
    my $p=abs_path($name // ''); my $base=abs_path('.runs');
    die "Verify requires an owned .runs/ directory\n" unless $p && $base && index($p,"$base/")==0;
    my $m=json("$p/manifest.json"); my @bad;
    for my $f(@{$m->{artifacts}//[]}) {
        push @bad,$f->{name} unless -f "$p/artifacts/$f->{name}"
            && sha256_hex(bytes("$p/artifacts/$f->{name}")) eq $f->{sha256};
    }
    push @bad,'submitted.sas' if $m->{submitted_sha256}
        && sha256_hex(bytes("$p/submitted.sas")) ne $m->{submitted_sha256};
    for my $pair(['source.sas','source_sha256'],['client.pl','client_sha256']) {
        push @bad,$pair->[0] if $m->{$pair->[1]}
            && (!-f "$p/$pair->[0]" || sha256_hex(bytes("$p/$pair->[0]")) ne $m->{$pair->[1]});
    }
    push @bad,'terminal outcome' unless $m->{expected_outcome_observed};
    push @bad,'complete capture' unless $m->{capture_complete} || ($m->{scope}//'') =~ /^Tool discovery/;
    push @bad,'sas.log' if $m->{log_sha256}
        && (!-f "$p/sas.log" || sha256_hex(bytes("$p/sas.log")) ne $m->{log_sha256});
    print @bad ? "FAIL: @bad\n" : "PASS: recorded submission, terminal outcome and artifact hashes\n";
    print "Numerical meaning requires the lesson's separate expected-value checks.\n";
    exit(@bad ? 1 : 0);
}
die "Use help, list, discover, run or verify\n" unless $mode eq 'discover' || $mode eq 'run';
my $s;
if($mode eq 'run') {
    die "Unknown lesson; use list\n" unless defined($name) && ($s=$catalog->{$name});
    die "Use baseline or exercise\n" unless $variant =~ /^(baseline|exercise)$/;
    die "No executed exercise source for $name\n" if $variant eq 'exercise' && !$s->{exercise};
    own_file($s->{$variant eq 'exercise' ? 'exercise' : 'source'});
    own_file($_) for @{$s->{inputs}//[]},@{$s->{definitions}//[]};
}
my $url=$ENV{BOOK_MCP_URL}//'http://127.0.0.1:18789/mcp';
die "This reader client requires a loopback MCP endpoint\n"
    unless $url =~ m{^http://(?:127\.0\.0\.1|localhost):\d+/mcp$};
my $token='mgxviz_'.time.'_'.$$;
my $folder='.runs/'.$token.'_'.($name//'discovery').'_'.$variant;
die "Refusing to reuse an existing result directory\n" if -e $folder;
make_path($folder);
put("$folder/client.pl",bytes("$FindBin::Bin/run_visualization.pl"));
print "Run directory: $folder\n";
my $ua=LWP::UserAgent->new(timeout=>240); my($seq,@rpc)=(0);
sub rpc {
    my($method,$params)=@_; my $request={jsonrpc=>'2.0',id=>++$seq,method=>$method,params=>$params//{}};
    my $at=stamp(); my $h=$ua->post($url,'Content-Type'=>'application/json',Content=>encode_json($request));
    my $response=eval {decode_json($h->decoded_content)};
    save("$folder/rpc_$seq.json",{request=>$request,response=>$response,http_status=>$h->code,sent_utc=>$at,received_utc=>stamp()});
    push @rpc,"rpc_$seq.json";
    die "HTTP/JSON failure: ".$h->status_line."\n" unless $h->is_success && $response;
    die encode_json($response->{error})."\n" if $response->{error};
    return $response->{result};
}
sub result_text { join "\n",map {$_->{text}//''} @{$_[0]{content}//[]} }
sub job {
    my($args)=@_; my $r=rpc('tools/call',{name=>'run_sas_codes_or_files_in_ODA',arguments=>$args});
    my($pid)=result_text($r)=~/^PID:\s*(\d+)/m;
    die "No job identifier returned: ".result_text($r)."\n" unless $pid;
    for my $n(0..30) {
        sleep 30 if $n;
        $r=rpc('tools/call',{name=>'run_sas_codes_or_files_in_ODA',arguments=>{pid=>0+$pid}});
        my $t=result_text($r);
        next if $t =~ /STATUS:\s*RUNNING/;
        my($status)=$t=~/STATUS:\s*(COMPLETE|FAILED)/;
        my($exit)=$t=~/Exit code:\s*(\d+)/;
        die "Unrecognized final reply; inspect $folder\n" unless $status;
        return {status=>$status,exit_code=>defined($exit)?0+$exit:undef,
                is_error=>$r->{isError}?JSON::PP::true:JSON::PP::false,pid=>0+$pid,text=>$t};
    }
    die "Polling budget exhausted; preserve the run records before deciding whether to retry\n";
}
rpc('initialize',{protocolVersion=>'2024-11-05',capabilities=>{},clientInfo=>{name=>'SAS-book-reader',version=>'1.0'}});
my $tools=rpc('tools/list',{}); save("$folder/tools.json",$tools);
die "Required SASPy-backed tool absent\n" unless grep {$_->{name} eq 'run_sas_codes_or_files_in_ODA'} @{$tools->{tools}//[]};
if($mode eq 'discover') {
    print "$_->{name}\n" for @{$tools->{tools}};
    save("$folder/manifest.json",{endpoint=>$url,client=>'SAS-book-reader 1.0',rpc=>\@rpc,expected_outcome_observed=>JSON::PP::true,scope=>'Tool discovery; no SAS submitted'}); exit;
}
my $source=$s->{$variant eq 'exercise' ? 'exercise' : 'source'};
my $lesson=bytes($source); my $source_snapshot=$lesson; my @inputs; my @uploads;
for my $input(@{$s->{inputs}//[]}) {
    my $base=basename($input); my $uploaded="${token}_$base";
    my $p="$folder/input/$uploaded"; put($p,bytes($input)); push @uploads,abs_path($p);
    push @inputs,{source=>$input,sha256=>sha256_hex(bytes($input)),remote=>"~/$uploaded"};
    $lesson =~ s{"\Q$base\E"}{"~/$uploaded"}g;
    $lesson =~ s{zip=\Q$base\E,}{zip=~/$uploaded,}g;
}
my %target=(persistent=>'false',saspy_cfgname=>$ENV{BOOK_SAS_TARGET}//'oda',run_timeout_seconds=>$s->{timeout}//300);
if(@uploads) { my $u=job({upload_file=>\@uploads,skip_upload_if_same=>'false',%target}); die "Input upload failed; inspect $folder\n" unless $u->{status} eq 'COMPLETE' && !$u->{is_error} }
my $remote="~/$token";
my $prefix=qq{options dlcreatedir notes source nodate nonumber;\nlibname bookout "$remote";\nlibname bookout clear;\n%let BOOK_OUT=%sysget(HOME)/$token;\n%let P2_OUT=%sysget(HOME)/$token;\n%let P3_OUT=%sysget(HOME)/$token;\n%let AF_OUT=%sysget(HOME)/$token;\n%let AF_REMOTE=$token;\n};
my @dependencies;
for my $definition(@{$s->{definitions}//[]}) { $prefix.="\n".bytes($definition)."\n";push @dependencies,{file=>$definition,sha256=>sha256_hex(bytes($definition))} }
for my $support(@{$s->{support}//[]}) { own_file($support);$prefix.="\n".bytes($support)."\n";push @dependencies,{file=>$support,sha256=>sha256_hex(bytes($support))} }
if($s->{upload_macrovar}) { die 'No uploaded input' unless @inputs;$prefix.="%let $s->{upload_macrovar}=$inputs[0]{remote};\n" }
$prefix.=qq{\nods listing gpath="&BOOK_OUT";\n};
$lesson =~ s/gpath='\.'/gpath="&BOOK_OUT"/g;
my $suffix="\n/* Reader capture: export values without changing the analysis. */\n";
my @files;
for my $table(@{$s->{tables}//[]}) {
    die 'Unsafe export table' unless $table =~ /^(?:\w+\.)?\w+$/;
    my $short=$table; $short=~s/^.*\.//;
    $suffix.=qq{data work.reader_export;set $table;format _numeric_ best32.;run;\nproc export data=work.reader_export outfile="$remote/$short.csv" dbms=csv replace;run;\n};
    push @files,"$short.csv";
}
$suffix.=qq{data work.reader_environment;length sas_version \$64 platform \$32 process_id \$40; sas_version="&sysvlong";platform="&sysscp";process_id="&sysprocessid";run;\nproc export data=work.reader_environment outfile="$remote/environment.csv" dbms=csv replace;run;\n};
push @files,'environment.csv',@{$s->{figures}//[]};
if($s->{display}) {
    die 'Unsafe display table' unless $s->{display}=~/^(?:\w+\.)?\w+$/;
    $suffix.=qq{options papersize=("6.5in","2.2in") leftmargin=.1in rightmargin=.1in topmargin=.1in bottommargin=.1in;\nods noproctitle;title;footnote;\nods printer printer=png300 file="$remote/table.png" style=HTMLBlue dpi=300;\nproc print data=$s->{display} noobs;run;\nods printer close;\n};
    push @files,'table.png';
}
if($s->{collect_images}) {
    $suffix.=qq{filename bookdir "$remote";\ndata work.reader_images;length Name \$256;handle=dopen('bookdir');if handle>0 then do i=1 to dnum(handle);Name=dread(handle,i);if lowcase(scan(Name,-1,'.')) in ('png','svg') then output;end;rc=dclose(handle);keep Name;run;\nproc export data=work.reader_images outfile="$remote/reader_images.csv" dbms=csv replace;run;\n};
    push @files,'reader_images.csv';
}
my $submitted=$prefix."\n/* LESSON_BEGIN */\n$lesson\n/* LESSON_END */\n".$suffix;
put("$folder/source.sas",$source_snapshot);put("$folder/submitted.sas",$submitted);
my $r=job({sas_codes_or_file=>$submitted,%target});put("$folder/result.txt",$r->{text});
my($log)=$r->{text}=~/^SAS log for debugging saved to:\s*([^\r\n]+)/m;
copy($log,"$folder/sas.log") or die "Cannot retain SAS log: $!" if $log && -f $log;
my $expect=$s->{expected}//'COMPLETE';my $accepted=$r->{status} eq $expect && ($expect eq 'FAILED' || !$r->{is_error});
my $m={lesson=>$name,variant=>$variant,source=>$source,source_sha256=>sha256_hex($source_snapshot),submitted_sha256=>sha256_hex($submitted),client_sha256=>sha256_hex(bytes("$folder/client.pl")),job_pid=>$r->{pid},inputs=>\@inputs,dependencies=>\@dependencies,endpoint=>$url,target=>$target{saspy_cfgname},persistent=>JSON::PP::false,status=>$r->{status},exit_code=>$r->{exit_code},expected=>$expect,expected_outcome_observed=>$accepted?JSON::PP::true:JSON::PP::false,capture_complete=>JSON::PP::false,remote_output=>$remote,rpc=>\@rpc,(-f"$folder/sas.log"?(log_sha256=>sha256_hex(bytes("$folder/sas.log"))):()),scope=>'Actual plotting/teaching submission; no upstream study or native IDE transcript implied'};
save("$folder/manifest.json",$m);
if($r->{status} eq 'COMPLETE' && !$r->{is_error}) {
    make_path("$folder/artifacts"); my @dest=map{abs_path("$folder/artifacts")."/$_"}@files;
    make_path(dirname($_)) for @dest;
    my $d=job({download_file=>[map{"$remote/$_"}@files],download_local_path=>\@dest,%target});
    die "Download failed; existing files are not accepted as fresh results\n" unless $d->{status} eq 'COMPLETE' && !$d->{is_error};
    if($s->{collect_images}) {
        require Text::CSV;
        my $csv=Text::CSV->new({binary=>1});open my $images,'<',"$folder/artifacts/reader_images.csv" or die $!;
        $csv->getline($images);my @names;while(my $row=$csv->getline($images)){die 'Unsafe generated image filename' unless $row->[0]=~/^[A-Za-z0-9_.-]+\.(png|svg)$/i;push @names,$row->[0]}
        die 'Expected plot produced no owned images' unless @names;
        my @new=map{abs_path("$folder/artifacts")."/$_"}@names;
        my $images_job=job({download_file=>[map{"$remote/$_"}@names],download_local_path=>\@new,%target});
        die 'Image collection failed' unless $images_job->{status} eq 'COMPLETE' && !$images_job->{is_error};push @dest,@new;
    }
    my @artifacts;my $artifact_root=abs_path("$folder/artifacts");for my $p(@dest) { die "Missing/empty output $p\n" unless -s $p;push @artifacts,{name=>substr($p,length($artifact_root)+1),bytes=>-s$p,sha256=>sha256_hex(bytes($p))} }
    $m->{artifacts}=\@artifacts;$m->{rpc}=[@rpc];$m->{capture_complete}=JSON::PP::true;save("$folder/manifest.json",$m);
}
elsif($accepted && $expect eq 'FAILED') { $m->{capture_complete}=JSON::PP::true;save("$folder/manifest.json",$m) }
print "$name/$variant: $r->{status}; expected $expect; records in $folder\n";
exit($accepted?0:1);
