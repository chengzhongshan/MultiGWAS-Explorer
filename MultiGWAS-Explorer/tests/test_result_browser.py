"""Cross-platform browser-launch regressions; no SAS login or browser required."""
import base64
import contextlib
import io
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'DiffGWASDeps'))
import result_browser as browser


class BrowserTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='SAS browser & labels ')
        self.addCleanup(self.temp.cleanup)
        self.html = Path(self.temp.name) / "result 'one' & two.html"
        self.html.write_text('<html><body>SAS output</body></html>', encoding='utf-8')
        env = patch.dict(os.environ)
        env.start()
        self.addCleanup(env.stop)
        os.environ.pop('OPEN_RESULT_BROWSER', None)

    def test_mac_uses_checked_native_open_without_firefox(self):
        with patch.object(browser.sys, 'platform', 'darwin'), patch.object(browser, '_find', return_value='/usr/bin/open'), patch.object(browser, '_checked') as run, patch.object(browser.webbrowser, 'open') as fallback:
            self.assertTrue(browser.open_result_html(self.html))
        run.assert_called_once_with(['/usr/bin/open', str(self.html.resolve())], None)
        fallback.assert_not_called()

    def test_failed_mac_launcher_is_not_reported_as_success(self):
        with patch.object(browser.sys, 'platform', 'darwin'), patch.object(browser, '_find', return_value='/usr/bin/open'), patch.object(browser, '_checked', side_effect=OSError('open exited 1')), patch.object(browser.webbrowser, 'open', return_value=False), contextlib.redirect_stderr(io.StringIO()) as errors, contextlib.redirect_stdout(io.StringIO()) as output:
            self.assertFalse(browser.open_result_html(self.html))
        self.assertIn('open exited 1', errors.getvalue())
        self.assertIn('Firefox is not required', errors.getvalue())
        self.assertNotIn('Requested', output.getvalue())

    def test_native_windows_uses_file_association(self):
        with patch.object(browser.sys, 'platform', 'win32'), patch.object(browser.os, 'startfile', create=True) as start, patch.object(browser.webbrowser, 'open') as fallback:
            self.assertTrue(browser.open_result_html(self.html))
        start.assert_called_once_with(str(self.html.resolve()))
        fallback.assert_not_called()

    def test_cygwin_without_cygstart_uses_converted_path_not_shell_code(self):
        native = "C:\\Users\\Reader & Co\\result 'one'.html"
        def find(name, system_path=None):
            return {'powershell.exe': 'powershell.exe', 'cygpath': '/usr/bin/cygpath'}.get(name)
        with patch.object(browser.sys, 'platform', 'cygwin'), patch.object(browser, '_find', side_effect=find), patch.object(browser, '_checked', side_effect=[native, '']) as run, patch.object(browser.webbrowser, 'open') as fallback:
            self.assertTrue(browser.open_result_html(self.html))
        command, env = run.call_args_list[1].args
        script = base64.b64decode(command[-1]).decode('utf-16le')
        self.assertNotIn(native, script)
        self.assertIn('$env:MCP4SAS_RESULT_FILE', script)
        self.assertEqual(env['MCP4SAS_RESULT_FILE'], native)
        fallback.assert_not_called()

    def test_cygstart_failure_falls_back_instead_of_stopping_after_fork(self):
        def find(name, system_path=None):
            return name
        with patch.object(browser.sys, 'platform', 'cygwin'), patch.object(browser, '_find', side_effect=find), patch.object(browser, '_checked', side_effect=[OSError('cygstart failed'), 'C:\\result.html', '']) as run:
            self.assertTrue(browser.open_result_html(self.html))
        self.assertEqual(run.call_count, 3)

    def test_headless_cygwin_does_not_send_unix_uri_to_windows(self):
        with patch.object(browser.sys, 'platform', 'cygwin'), patch.object(browser, '_find', return_value=None), patch.object(browser.webbrowser, 'open') as fallback, contextlib.redirect_stderr(io.StringIO()) as errors:
            self.assertFalse(browser.open_result_html(self.html))
        fallback.assert_not_called()
        self.assertIn('Open this saved HTML file manually', errors.getvalue())

    def test_linux_default_launcher_precedes_named_browser_fallback(self):
        with patch.dict(os.environ, DISPLAY=':0'), patch.object(browser.sys, 'platform', 'linux'), patch.object(browser, '_find', return_value='/usr/bin/xdg-open'), patch.object(browser, '_checked') as run:
            self.assertTrue(browser.open_result_html(self.html))
        self.assertEqual(run.call_args.args[0][0], '/usr/bin/xdg-open')

    def test_headless_linux_never_blocks_on_a_text_browser(self):
        with patch.dict(os.environ, {}, clear=True), patch.object(browser.sys, 'platform', 'linux'), patch.object(browser, '_checked') as run, patch.object(browser.webbrowser, 'open') as fallback, contextlib.redirect_stderr(io.StringIO()):
            self.assertFalse(browser.open_result_html(self.html))
        run.assert_not_called()
        fallback.assert_not_called()

    def test_checked_launcher_reports_exit_and_limits_wait(self):
        result = subprocess.CompletedProcess(['open'], 1, '', 'No application registered')
        with patch.object(browser.subprocess, 'run', return_value=result) as run:
            with self.assertRaisesRegex(OSError, 'exited 1'):
                browser._checked(['open', str(self.html)])
        self.assertEqual(run.call_args.kwargs['timeout'], 15)
        self.assertNotIn('shell', run.call_args.kwargs)

    def test_missing_file_is_advisory_and_diagnostic_opens_nothing(self):
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertFalse(browser.open_result_html(self.html.with_name('missing.html')))
        with patch.object(browser, 'open_result_html') as opener, contextlib.redirect_stdout(io.StringIO()) as output:
            self.assertEqual(browser.main(['--check']), 0)
        opener.assert_not_called()
        self.assertFalse(json.loads(output.getvalue())['firefox_required'])

    def test_false_environment_values_are_consistent(self):
        for value in ('0', 'false', 'FALSE', 'no', 'off'):
            self.assertFalse(browser.auto_open_enabled(value))
        self.assertTrue(browser.auto_open_enabled('1'))


if __name__ == "__main__":
    unittest.main()
