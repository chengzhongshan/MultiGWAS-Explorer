"""Exercise exact source staging and metadata queries with an inert SASPy stub."""
import importlib.util
from pathlib import Path
import re
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

root=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('explorer_helper',root/'DiffGWASDeps/sas_oda_session_server.py')
helper=importlib.util.module_from_spec(spec)
with patch.dict(sys.modules,{'saspy':types.SimpleNamespace(SASsession=lambda **kw:None)}):
    spec.loader.exec_module(helper)


class Session:
    def __init__(self,work,remote=False,upload_success=True):
        self.workpath=str(work)
        self._io=types.SimpleNamespace(sascfg=types.SimpleNamespace(iomhost=['oda.invalid'] if remote else []))
        self.calls=[];self.uploads=[];self.upload_success=upload_success;self.source=None

    def upload(self,source,target,overwrite):
        self.uploads.append((target,Path(source).read_bytes()))
        return {'Success':self.upload_success}

    def submit(self,code,results='html'):
        self.calls.append((code,results))
        if code.startswith('filename '):
            path=re.match(r"filename \w+ '((?:''|[^'])*)'",code).group(1).replace("''", "'")
            self.source=self.uploads[0][1] if self.uploads else Path(path).read_bytes()
        return {'LOG':'NOTE: user program completed','LST':'<html>kept</html>'}


class EfficientSubmitTests(unittest.TestCase):
    def test_small_program_remains_one_direct_submit(self):
        session=Session('unused')
        helper.submit_user_code(session,'proc print;run;')
        self.assertEqual(session.calls,[('proc print;run;','html')])
        self.assertEqual(session.uploads,[])

    def check_staging(self,remote):
        code='/* '+'variant label 日本語\n'*3000+' */\nproc print;run;\n'
        with tempfile.TemporaryDirectory(prefix='SAS source spaces ') as directory:
            session=Session(directory,remote=remote)
            result=helper.submit_user_code(session,code)
        self.assertEqual(session.source,code.encode('utf-8'))
        self.assertEqual(len(session.calls),1,'user code must not be replayed')
        self.assertEqual(len(session.uploads),int(remote))
        self.assertIn('%include ',session.calls[0][0])
        self.assertIn('staged source bytes=',result['LOG'])
        self.assertEqual(result['LST'],'<html>kept</html>')

    def test_large_local_source_preserves_exact_utf8(self):self.check_staging(False)
    def test_large_remote_source_uses_one_upload(self):self.check_staging(True)

    def test_failed_staging_does_not_submit(self):
        session=Session('/remote/work',remote=True,upload_success=False)
        with self.assertRaises(IOError):helper.submit_user_code(session,'/*'+'x'*40000+'*/')
        self.assertEqual(session.calls,[])

    def test_metadata_uses_one_text_round_trip(self):
        session=Session('/remote/work',remote=True)
        def submit(code,results='html'):
            session.calls.append((code,results))
            return {'LOG':'MCP4SAS_META_EXISTS=1\nMCP4SAS_META_SIZE=19\nMCP4SAS_META_CREATED=\nMCP4SAS_META_MODIFIED=\n'}
        session.submit=submit
        result=helper.run_fileinfo(session,'/remote/input.txt')
        self.assertTrue(result['exists']);self.assertEqual(result['size'],19)
        self.assertEqual(len(session.calls),1)
        self.assertEqual(session.calls[0][1],'text')


if __name__=='__main__':unittest.main(verbosity=2)
