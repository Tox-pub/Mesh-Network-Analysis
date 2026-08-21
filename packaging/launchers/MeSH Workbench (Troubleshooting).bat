@echo off
REM  Same launcher with a console attached, so a startup error stays readable.
title MeSH Workbench (console)
cd /d "%~dp0"
"python\python.exe" "app\launch.py" %*
echo.
echo Exited with code %ERRORLEVEL%
pause
