@echo off
REM  Puts a MeSH Workbench shortcut on the Desktop, pointing back at this folder.
REM
REM  The program stays where it is - a shortcut is a pointer, not a copy - so
REM  moving or deleting this folder later breaks the shortcut, which is the
REM  intended behaviour for something that installs nothing.
title MeSH Workbench - Create desktop shortcut
setlocal

set "TARGET=%~dp0MeSH Workbench.bat"
set "ICON=%~dp0app\mesh_workbench\assets\mesh_workbench.ico"

if not exist "%TARGET%" (
  echo [X] "MeSH Workbench.bat" is not next to this file.
  echo     Run this from inside the extracted folder.
  echo.
  pause
  exit /b 1
)

REM  Built through WScript rather than PowerShell: shortcut creation needs the
REM  WScript.Shell COM object either way, and cscript is not affected by the
REM  PowerShell execution policy, which is commonly restricted on managed
REM  machines and would block a .ps1 outright.
set "VBS=%TEMP%\mw_shortcut.vbs"
> "%VBS%" echo Set oWS = WScript.CreateObject("WScript.Shell")
>>"%VBS%" echo sLink = oWS.SpecialFolders("Desktop") ^& "\MeSH Workbench.lnk"
>>"%VBS%" echo Set oLink = oWS.CreateShortcut(sLink)
>>"%VBS%" echo oLink.TargetPath = "%TARGET%"
>>"%VBS%" echo oLink.WorkingDirectory = "%~dp0"
>>"%VBS%" echo oLink.Description = "MeSH Workbench - concept network analysis"
>>"%VBS%" echo oLink.WindowStyle = 7
if exist "%ICON%" >>"%VBS%" echo oLink.IconLocation = "%ICON%"
>>"%VBS%" echo oLink.Save

cscript //nologo "%VBS%"
set "RC=%ERRORLEVEL%"
del "%VBS%" >nul 2>&1

if not "%RC%"=="0" (
  echo [X] The shortcut could not be created ^(code %RC%^).
  echo     Some managed machines block scripted shortcut creation. You can
  echo     right-click "MeSH Workbench.bat" and choose Send to ^> Desktop instead.
) else (
  echo Done - "MeSH Workbench" is on your Desktop.
)
echo.
pause
