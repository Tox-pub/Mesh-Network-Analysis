@echo off
REM  MeSH Workbench - portable launcher
REM
REM  Runs the bundled Python, which is the official embeddable build signed by
REM  the Python Software Foundation. Nothing here needs installing and nothing
REM  is written outside this folder.
title MeSH Workbench
cd /d "%~dp0"
if not exist "python\python.exe" (
  echo [X] python\python.exe is missing.
  echo     Extract the whole zip, keeping the folder structure intact.
  echo.
  pause
  exit /b 1
)
start "" "python\pythonw.exe" "app\launch.py" %*
