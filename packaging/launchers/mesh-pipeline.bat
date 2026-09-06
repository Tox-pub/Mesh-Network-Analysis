@echo off
REM  The pipeline, without the window. Mirrors the mesh-pipeline script in the
REM  Linux and macOS bundles, so the same commands work on all three.
REM
REM  Every argument is passed through unchanged:
REM      mesh-pipeline.bat --step viz
REM      mesh-pipeline.bat --step all --interactive
REM
REM  If this file will not run, invoke the module directly - it is the same
REM  command with nothing in front of it:
REM      "python\python.exe" -m mesh_aop.cli --step viz
setlocal
cd /d "%~dp0"

if not exist "python\python.exe" (
  echo [X] python\python.exe is missing.
  echo     Extract the whole zip, keeping the folder structure intact.
  exit /b 1
)

REM  python.exe, not pythonw.exe: the pipeline writes progress to stdout, and
REM  under pythonw stdout is None, so the first print ends the run.
"python\python.exe" -m mesh_aop.cli %*
exit /b %ERRORLEVEL%
