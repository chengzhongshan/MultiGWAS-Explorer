@echo off
setlocal
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0install_windows_portable_cygwin.ps1" %*
endlocal
