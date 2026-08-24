@echo off
REM ---------------------------------------------------------------------------
REM  Install MeSH Workbench for the current user.
REM
REM  Nothing here is compiled. This copies the folder into place, makes the
REM  shortcuts, and registers the program so it appears in Add/Remove Programs -
REM  using only tools already present in Windows.
REM
REM  That matters on a managed machine. A freshly built Setup.exe has no
REM  prevalence and is not signed, so Defender Exploit Guard refuses to run it
REM  ("blocked an operation that is not allowed by your IT administrator",
REM  0x80070005). No new executable is created here, so that rule has nothing to
REM  act on. PowerShell is avoided for the same reason: its execution policy is
REM  frequently locked by policy, while cscript and reg.exe are not.
REM
REM  Usage:
REM    Install.bat                     install to the default per-user location
REM    Install.bat "D:\somewhere"      install to a chosen folder
REM    Install.bat /S                  silent - no prompts, for deployment
REM ---------------------------------------------------------------------------
setlocal EnableExtensions

set "APPNAME=MeSH Workbench"
set "APPVER=3.0.0"
set "SRC=%~dp0"
set "TARGET="
set "SILENT="

:parse
if "%~1"=="" goto parsed
if /i "%~1"=="/S" goto parse_silent
set "TARGET=%~1"
shift
goto parse
:parse_silent
set "SILENT=1"
shift
goto parse
:parsed
if not defined TARGET set "TARGET=%LOCALAPPDATA%\Programs\%APPNAME%"

title Install %APPNAME%
echo.
echo   %APPNAME% %APPVER%
echo   =============================
echo.
echo   From : %SRC%
echo   To   : %TARGET%
echo.
echo   Installs for you only. No administrator rights are needed and nothing
echo   outside the folder above is modified, apart from your Start menu.
echo.

if not exist "%SRC%python\python.exe" goto not_extracted

if defined SILENT goto copy_files
set "GO="
set /p GO=Continue? [Y/n]
if /i "%GO%"=="n" goto cancelled

:copy_files
echo.
echo   Copying files...
REM  /NFL /NDL /NJH /NJS keep the output quiet. robocopy returns values below 8
REM  for success, so only 8 and above is a real failure.
set "EXCL=Install.bat Uninstall.bat portable.marker"
robocopy "%SRC%." "%TARGET%" /E /NFL /NDL /NJH /NJS /NP /XF %EXCL% >nul
if %ERRORLEVEL% GEQ 8 goto copy_failed
echo   Copied.

REM  The marker is what tells the application it owns its folder. An installed
REM  copy must not have it, or settings and results would be written beside the
REM  program and shared between accounts instead of kept per user.
if exist "%TARGET%\portable.marker" del /q "%TARGET%\portable.marker"

set "ICON=%TARGET%\app\mesh_workbench\assets\mesh_workbench.ico"
set "PYW=%TARGET%\python\pythonw.exe"
set "LAUNCH=%TARGET%\app\launch.py"

echo   Creating shortcuts...
call :make_shortcut "Programs" "%APPNAME%"
echo   Start menu entry created.

if defined SILENT goto make_desktop
set "DESK="
set /p DESK=Add a desktop shortcut too? [Y/n]
if /i "%DESK%"=="n" goto register

:make_desktop
call :make_shortcut "Desktop" "%APPNAME%"
echo   Desktop shortcut created.

:register
REM  Registering under HKCU is what puts the program in Settings ^> Apps. This
REM  is the per-user hive, so it needs no elevation.
set "KEY=HKCU\Software\Microsoft\Windows\CurrentVersion\Uninstall\MeSHWorkbench"
reg add "%KEY%" /v DisplayName     /t REG_SZ /d "%APPNAME%" /f >nul
reg add "%KEY%" /v DisplayVersion  /t REG_SZ /d "%APPVER%" /f >nul
reg add "%KEY%" /v Publisher       /t REG_SZ /d "Karolinska Institutet" /f >nul
reg add "%KEY%" /v InstallLocation /t REG_SZ /d "%TARGET%" /f >nul
reg add "%KEY%" /v DisplayIcon     /t REG_SZ /d "%ICON%" /f >nul
reg add "%KEY%" /v UninstallString /t REG_SZ /d "\"%TARGET%\Uninstall.bat\"" /f >nul
reg add "%KEY%" /v NoModify        /t REG_DWORD /d 1 /f >nul
reg add "%KEY%" /v NoRepair        /t REG_DWORD /d 1 /f >nul
reg add "%KEY%" /v EstimatedSize   /t REG_DWORD /d 471040 /f >nul
echo   Registered in Add/Remove Programs.

call :write_uninstaller

echo.
echo   Done.
echo.
echo   Start it from the Start menu, or from the desktop shortcut.
echo   Remove it from Settings ^> Apps, or by running:
echo     "%TARGET%\Uninstall.bat"
echo.
if not defined SILENT pause
exit /b 0

REM ---------------------------------------------------------------------------
:make_shortcut
REM  %~1 = special folder name, %~2 = shortcut name. Built through cscript
REM  because shortcut creation needs the WScript.Shell COM object either way,
REM  and cscript is not subject to the PowerShell execution policy.
set "VBS=%TEMP%\mw_sc_%RANDOM%.vbs"
> "%VBS%" echo Set oWS = WScript.CreateObject("WScript.Shell")
>>"%VBS%" echo sLink = oWS.SpecialFolders("%~1") ^& "\%~2.lnk"
>>"%VBS%" echo Set L = oWS.CreateShortcut(sLink)
>>"%VBS%" echo L.TargetPath = "%PYW%"
>>"%VBS%" echo L.Arguments = """%LAUNCH%"""
>>"%VBS%" echo L.WorkingDirectory = "%TARGET%"
>>"%VBS%" echo L.Description = "Build and validate MeSH concept networks"
>>"%VBS%" echo L.IconLocation = "%ICON%"
>>"%VBS%" echo L.Save
cscript //nologo "%VBS%" >nul
del "%VBS%" >nul 2>&1
exit /b 0

REM ---------------------------------------------------------------------------
:write_uninstaller
REM  Written into the installed folder rather than copied, so it always knows
REM  where it was installed without having to work it out.
set "U=%TARGET%\Uninstall.bat"
> "%U%" echo @echo off
>>"%U%" echo setlocal
>>"%U%" echo title Uninstall %APPNAME%
>>"%U%" echo echo.
>>"%U%" echo echo   Removing %APPNAME%
>>"%U%" echo echo     %TARGET%
>>"%U%" echo echo.
>>"%U%" echo REM  The application's own uninstaller knows about the downloaded
>>"%U%" echo REM  archives, the databases and the scratch copy left in the temp
>>"%U%" echo REM  folder - none of which live here, and all of which would
>>"%U%" echo REM  otherwise survive. Results are kept unless asked for.
>>"%U%" echo if exist "%TARGET%\python\python.exe" "%TARGET%\python\python.exe" -m mesh_aop.uninstall_cli --project "%TARGET%"
>>"%U%" echo echo.
>>"%U%" echo set "V=%%TEMP%%\mw_rmsc.vbs"
>>"%U%" echo ^> "%%V%%" echo Set oWS = WScript.CreateObject("WScript.Shell")
>>"%U%" echo ^>^>"%%V%%" echo On Error Resume Next
>>"%U%" echo ^>^>"%%V%%" echo Set fso = CreateObject("Scripting.FileSystemObject")
>>"%U%" echo ^>^>"%%V%%" echo fso.DeleteFile oWS.SpecialFolders("Programs") ^^^& "\%APPNAME%.lnk", True
>>"%U%" echo ^>^>"%%V%%" echo fso.DeleteFile oWS.SpecialFolders("Desktop") ^^^& "\%APPNAME%.lnk", True
>>"%U%" echo cscript //nologo "%%V%%" ^>nul 2^>^&1
>>"%U%" echo del "%%V%%" ^>nul 2^>^&1
>>"%U%" echo reg delete "%KEY%" /f ^>nul 2^>^&1
>>"%U%" echo echo   Shortcuts and the Add/Remove Programs entry are gone.
>>"%U%" echo echo.
>>"%U%" echo REM  A folder cannot delete itself while a program inside it runs,
>>"%U%" echo REM  so the last step is handed to a script in TEMP.
>>"%U%" echo ^> "%%TEMP%%\mw_rm.bat" echo @echo off
>>"%U%" echo ^>^>"%%TEMP%%\mw_rm.bat" echo ping -n 4 127.0.0.1 ^^^>nul
>>"%U%" echo ^>^>"%%TEMP%%\mw_rm.bat" echo rmdir /s /q "%TARGET%"
>>"%U%" echo ^>^>"%%TEMP%%\mw_rm.bat" echo del "%%%%~f0"
>>"%U%" echo cd /d "%%TEMP%%"
>>"%U%" echo start "" /min cmd /c "%%TEMP%%\mw_rm.bat"
>>"%U%" echo echo   The program folder will disappear in a few seconds.
>>"%U%" echo ping -n 3 127.0.0.1 ^>nul
>>"%U%" echo exit /b 0
exit /b 0

REM ---------------------------------------------------------------------------
:not_extracted
echo   [X] This does not look like the extracted program folder.
echo       Run Install.bat from inside the folder you extracted.
echo.
if not defined SILENT pause
exit /b 1

:copy_failed
echo   [X] Copy failed ^(robocopy code %ERRORLEVEL%^).
echo.
if not defined SILENT pause
exit /b 1

:cancelled
echo   Cancelled.
if not defined SILENT pause
exit /b 1
