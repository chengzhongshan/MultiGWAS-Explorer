#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../DiffGWASDeps";
use SAS_Runtime_Selection;
use SAS_ODA_Runner;
use Test::More;
use File::Temp qw(tempdir);
use File::Spec;
use JSON::PP qw(encode_json decode_json);

my $dir=tempdir(CLEANUP=>1);
my $default=SAS_Runtime_Selection::resolve(env=>{});
is($default->{SASPY_CFGNAME},'oda','ODA stays the default');
is($default->{SAS_LOCAL_HOME_DIR},'','ODA has no local staging home');
my $local=SAS_Runtime_Selection::resolve(
    cli=>{cfgname=>'linuxlocal'}, spec=>{saspy_cfgname=>'winlocal'},
    env=>{SASPY_CFGNAME=>'oda'},output_dir=>$dir,artifact_stem=>'test');
is($local->{SASPY_CFGNAME},'linuxlocal','CLI profile overrides spec and environment');
like($local->{SAS_LOCAL_HOME_DIR},qr{SAS_Local_Home_test$},'local steps share a dedicated staging child directory');
is_deeply(SAS_Runtime_Selection::resolve(spec=>{saspy_cfgname=>'winlocal',saspy_cfgfile=>'local profile.py'},env=>{SASPY_CFGNAME=>'oda'})->{SASPY_CFGNAME},'winlocal','spec overrides environment');
is(SAS_Runtime_Selection::resolve(env=>{SASPY_CONFIG_NAME=>'linuxlocal'})->{SASPY_CFGNAME},'linuxlocal','legacy environment profile supported');
is(SAS_Runtime_Selection::resolve(cli=>{cfgname=>'oda'},env=>{SAS_LOCAL_HOME_DIR=>$dir})->{SAS_LOCAL_HOME_DIR},'','switching to ODA clears inherited local staging');
eval { SAS_Runtime_Selection::resolve(cli=>{cfgname=>'oda',home=>$dir}); };
like($@,qr/cannot be used/,'contradictory explicit ODA and local home rejected');
is(SAS_Runtime_Selection::resolve(cli=>{cfgname=>'my_custom_local',home=>$dir})->{SAS_LOCAL_HOME_DIR},$dir,'custom local profiles supported with explicit home');

# Exercise the real workflow config emitter, not just the resolver.
my $root=File::Spec->rel2abs("$Bin/..");
open my $in,'<',"$root/configs/spec_pgc_scz_sex_common_automation.json" or die $!;
my $spec=decode_json(do {local $/;<$in>});close $in;
$spec->{workdir}=$root;$spec->{configs_dir}=$dir;$spec->{output_dir}=$dir;
$spec->{input_dir}=$dir;$spec->{artifact_stem}='target_test';
$spec->{saspy_cfgname}='oda';
open my $out,'>',"$dir/spec.json" or die $!;print {$out} encode_json($spec);close $out;
{
    local $ENV{SASPY_CFGNAME}='oda';local $ENV{SASPY_CFGFILE}='';local $ENV{SAS_LOCAL_HOME_DIR}='';
    my @cmd=($^X,"$root/auto_prepare_and_run_diff_gwas.pl",'--spec',"$dir/spec.json",'--mode','configs',
        '--saspy-cfgname','winlocal','--saspy-cfgfile','C:/SAS config/local profile.py');
    open my $pipe,'-|',@cmd or die $!;my $text=do {local$/;<$pipe>};close $pipe;
    is($?>>8,0,'workflow accepts explicit local SAS in config-only mode without connecting');
    like($text,qr/Requested SASPy profile: winlocal/,'workflow reports the requested target');
}
open $in,'<',"$dir/auto_target_test_runner.json" or die $!;
my $cfg=decode_json(do {local$/;<$in>});close $in;
is($cfg->{SASPY_CFGNAME},'winlocal','saved plotting configuration pins profile');
is($cfg->{SASPY_CFGFILE},'C:/SAS config/local profile.py','configuration path with spaces retained');
like($cfg->{SAS_LOCAL_HOME_DIR},qr/SAS_Local_Home_target_test$/,'upload/submit/download prefixes share one home');
open my $pipe,'-|',$^X,"$root/DiffGWASDeps/emit_diff_gwas_runner_env.pl",'--config',"$dir/auto_target_test_runner.json" or die $!;
my $exports=do {local$/;<$pipe>};close $pipe;
like($exports,qr/export SASPY_CFGNAME='winlocal'/,'shell wrappers receive the selected target');
like($exports,qr/export SAS_LOCAL_HOME_DIR='/,'shell wrappers receive stable staging home');

# Exercise the actual framed Python handler using SASPy fixtures: no SAS account.
my $python=$ENV{PIPELINE_PYTHON_BIN} || 'python3';
my $script="$dir/target.py";
open $out,'>',$script or die $!;
print {$out} <<'PY';
import sys, os, json, types, socket, threading, runpy
from pathlib import Path
calls=[]
class Session:
    def __init__(self, **kw):
        calls.append(kw)
        if kw['cfgname']=='broken_local': raise OSError('local SAS unavailable')
        self.kw=kw
        self._io=types.SimpleNamespace(sascfg=types.SimpleNamespace(iomhost=None,ssh=None))
    def submit(self, code, **kw): return {'LOG':'NOTE: '+self.kw['cfgname']+' '+str(self.kw['cfgfile']), 'LST':'<html>fixture</html>'}
sys.modules['saspy']=types.SimpleNamespace(SASsession=Session)
os.environ['SASPY_CFGNAME']='oda'
os.environ['SASPY_CFGFILE']='server_first_client.py'
os.environ['SAS_ODA_SESSION_DEBUG_LOG']=str(Path(sys.argv[2])/'handler.log')
ns=runpy.run_path(sys.argv[1])
def request(name, cfgfile='test profile.py', home=''):
    a,b=socket.socketpair();a.settimeout(10)
    body=json.dumps({'cmd':'submit','session_id':'same_id','sas_config':{'cfgname':name,'cfgfile':cfgfile},'local_home':home,'load_macros':False,'code':'data _null_; run;'}).encode()
    t=threading.Thread(target=ns['handle_client'],args=(b,('127.0.0.1',1)))
    t.start();a.sendall(len(body).to_bytes(8,'big')+body)
    def receive(n):
        out=b''
        while len(out)<n:
            part=a.recv(n-len(out));assert part,'incomplete response';out+=part
        return out
    result=json.loads(receive(int.from_bytes(receive(8),'big')))
    a.close();t.join(10);assert not t.is_alive()
    return result
assert request('oda')['status']=='ok'
assert request('winlocal',home=sys.argv[2])['status']=='ok'
assert request('linuxlocal')['status']=='ok'
assert len(calls)==3 and [c['cfgname'] for c in calls]==['oda','winlocal','linuxlocal']
assert request('winlocal',home=sys.argv[2])['status']=='ok' and len(calls)==3
assert request('winlocal','second config.py')['status']=='ok' and len(calls)==4
error=request('broken_local');assert error['status']=='error' and 'local SAS unavailable' in error['error']
assert len(calls)==5 and calls[-1]['cfgname']=='broken_local'  # no ODA fallback
assert all(c['cfgfile']!='server_first_client.py' for c in calls)
assert len(ns['sessions'])==4
local=[s for s in ns['sessions'].values() if s.kw['cfgname']=='winlocal' and s.kw['cfgfile']=='test profile.py'][0]
assert local._pipeline_local_home==sys.argv[2].replace('\\','/')
print('Framed handler: ODA/Windows/Linux routing, config isolation, reuse and no fallback PASS')
PY
close $out;
my @args=map {SAS_ODA_Runner::_normalize_local_path_for_python($_,$python)}
    ($script,"$root/DiffGWASDeps/sas_oda_session_server.py",$dir);
is(system({$python} $python,@args),0,'actual session handler isolates explicit targets and preserves failures');

for my $path ('/home/test/figure.png','C:/SAS output/figure.png','/work/project/figure.svg') {
    open $out,'>',"$dir/plot.log" or die $!;print {$out} "NOTE: before\nThe final figure is put here:\n$path\nNOTE: after\n";close $out;
    open my $parse,'-|',$^X,"$root/DiffGWASDeps/extract_sas_figure_path.pl","$dir/plot.log" or die $!;
    my $actual=do {local$/;<$parse>};close $parse;
    is($actual,$path,'plot download recovers full absolute path including spaces');
}
done_testing;
