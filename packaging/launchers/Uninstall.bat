@echo off
REM  Removes the program and everything it downloaded or built - including the
REM  ETL working folder in the system temp directory, which deleting this folder
REM  would otherwise leave behind.
title MeSH Workbench - Uninstall
cd /d "%~dp0"
if not exist "python\python.exe" (
  echo [X] python\python.exe is missing - nothing to run.
  pause
  exit /b 1
)
"python\python.exe" -m mesh_aop.uninstall_cli --project "%~dp0."
echo.
set "FIN="
set /p FIN=Delete the remaining program folder too? [y/N]
if /i not "%FIN%"=="y" goto :done

REM  A folder cannot delete itself while a program inside it is running, so the
REM  removal is handed to a script in TEMP that waits for this one to exit.
> "%TEMP%\mw_cleanup.bat" echo @echo off
>>"%TEMP%\mw_cleanup.bat" echo timeout /t 3 /nobreak ^>nul
>>"%TEMP%\mw_cleanup.bat" echo rmdir /s /q "%~dp0."
>>"%TEMP%\mw_cleanup.bat" echo del "%%~f0"
cd /d "%TEMP%"
start "" /min cmd /c "%TEMP%\mw_cleanup.bat"
echo.
echo The folder will disappear in a few seconds.
timeout /t 2 /nobreak >nul
exit /b 0

:done
echo.
pause
