#!/usr/bin/env python3
"""Opt-in real SAS test of macro queries, copied HTML source and WORK reuse.

Use a configured SASPy interpreter and --live. Each invocation owns a unique
session/helper; it does not attach to or close another user's WORK session.
"""
import argparse
import json
import os
from pathlib import Path
import secrets
import shutil
import socket
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tests'))
from test_macroparas_html import validate_fixture


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--live',action='store_true',required=True)
    parser.add_argument('--runner',choices=('perl','python'),required=True)
    parser.add_argument('--perl',default=shutil.which('perl'))
    parser.add_argument('--saspy-cfgname',default=os.environ.get('SASPY_CFGNAME','oda'))
    parser.add_argument('--saspy-cfgfile',default=os.environ.get('SASPY_CFGFILE'))
    parser.add_argument('--macro-library',help='Complete existing SAS-side macro library for the optional bulk-import check.')
    args=parser.parse_args()
    token=secrets.token_hex(6);session='explorer_macroquery_'+token
    output=ROOT/('SAS_Output_RandID_macroquery_live_'+token);output.mkdir()
    env=dict(os.environ,OPEN_RESULT='0',SAS_ODA_AUTOLOAD_MACROS='0',
             SASPY_CFGNAME=args.saspy_cfgname,PIPELINE_PYTHON_BIN=sys.executable,
             MCP4SAS_PYTHON_STATE_DIR=str(output/'backend-state'))
    if args.saspy_cfgfile:env['SASPY_CFGFILE']=str(Path(args.saspy_cfgfile).resolve())
    with socket.socket() as reservation:
        reservation.bind(('127.0.0.1',0));port=reservation.getsockname()[1]
    env['SAS_ODA_SESSION_PORT']=str(port)
    cygwin=False
    if args.runner=='perl':
        if not args.perl:parser.error('Perl not found; pass --perl /path/to/perl')
        cygwin=subprocess.check_output([args.perl,'-e','print $^O'],env=env,text=True).strip()=='cygwin'
    def client_path(path):
        if cygwin:
            converter=str(Path(args.perl).with_name('cygpath.exe'))
            return subprocess.check_output([converter,'-u',str(path)],env=env,text=True).strip()
        return str(path)
    interpreter=args.perl if args.runner=='perl' else sys.executable
    script=ROOT/('run_sas_codes_or_files_in_ODA.'+('pl' if args.runner=='perl' else 'py'))
    common=[interpreter,client_path(script),'--persistent','--session-id',session,'--no-open-html']
    records=[]
    def submit(name,code):
        file=output/(name+'.sas');file.write_text(code,encoding='utf-8',newline='\n')
        directory=output/name
        start=time.monotonic()
        result=subprocess.run(common+['--file',client_path(file),'--output-prefix',client_path(directory)],
                              cwd=ROOT,env=env,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=240)
        (output/(name+'.console.log')).write_text(result.stdout+'\n'+result.stderr,encoding='utf-8')
        record=dict(name=name,exit_code=result.returncode,elapsed_seconds=round(time.monotonic()-start,3))
        records.append(record);print(json.dumps(record),flush=True)
        if result.returncode:raise AssertionError(name+' failed; inspect '+str(output))
        return directory
    try:
        directory=submit('fixture',(ROOT/'tests/macroparas_html_regression.sas').read_text())
        html=next(directory.glob('*.html')).read_text(encoding='utf-8')
        copied=validate_fixture(html)
        submit('seed-legacy',"""data work.mcp_macroparas_sentinel;sentinel=42;run;
%let MCP4SAS_MACRO_DIR=&_mp_html_dir;
%macro macroparas(macrorgx=);%put ERROR: OLD_MACROPARAS_MUST_BE_REPLACED;%mend;
""")
        directory=submit('query-fallback',"""%macroparas(macrorgx=CopyPreview,numlines2print=1000,print_outdsd=1);proc print;run;
data _null_;set work.mcp_macroparas_sentinel;if sentinel=42 then put 'PASS: WORK_PRESERVED';else put 'ERROR: WORK_CHANGED';run;
proc sql noprint;select count(*) into :_query_count trimmed from work.Macro_paras;quit;
%macro verify; %if &_query_count=2 %then %put PASS: TWO_MACROS_FOUND; %else %put ERROR: MACRO_COUNT_CHANGED; %mend; %verify;
""")
        log=next(directory.glob('*.info.txt')).read_text(encoding='utf-8')
        assert '\nPASS: WORK_PRESERVED' in log and '\nPASS: TWO_MACROS_FOUND' in log
        validate_fixture(next(directory.glob('*.html')).read_text(encoding='utf-8'))
        if args.macro_library:
            library=args.macro_library
            if any(c in library for c in '\r\n;'):raise ValueError('Invalid SAS macro-library path')
            submit('bundled-importer',f"""%importallmacros_ue(MacroDir={library},fileRgx=macroparas);
%macroparas(macrorgx=QueryGTEx4GeneID,numlines2print=0,print_outdsd=1);
%put MCP4SAS_BUNDLE_POST_IMPORT_OK;
""")
        directory=submit('run-copied-fixtures',copied+"""
data _null_;set work.mp_copied_a;if value=7 and label='A < B & C > D' then put 'PASS: COPIED_A_VALUE';else put 'ERROR: COPIED_A_CHANGED';run;
data _null_;set work.mp_copied_b;if flag=1 and note='It''s a <test> & "quote"' then put 'PASS: COPIED_B_VALUE';else put 'ERROR: COPIED_B_CHANGED';run;
""")
        log=next(directory.glob('*.info.txt')).read_text(encoding='utf-8')
        assert '\nPASS: COPIED_A_VALUE' in log and '\nPASS: COPIED_B_VALUE' in log
        print('PASS: real SAS HTML, copied programs, legacy replacement and persistent WORK.')
    finally:
        (output/'summary.json').write_text(json.dumps(dict(runner=args.runner,profile=args.saspy_cfgname,records=records),indent=2)+'\n')
        if args.runner=='python':
            subprocess.run(common+['--kill-saspy-sessions','--output-prefix',str(output/'close')],cwd=ROOT,env=env,
                           stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,timeout=40,check=False)
        else:
            def rpc(payload):
                raw=json.dumps(payload).encode()
                with socket.create_connection(('127.0.0.1',port),timeout=20) as connection:
                    connection.sendall(len(raw).to_bytes(8,'big')+raw)
                    connection.recv(65536)
            try:
                cfgfile=env.get('SASPY_CFGFILE','').replace('\\','/')
                rpc(dict(cmd='close_session',session_id=session,sas_config=dict(cfgname=args.saspy_cfgname,cfgfile=cfgfile)))
                rpc(dict(cmd='shutdown'))
            except OSError:pass


if __name__=='__main__':main()
