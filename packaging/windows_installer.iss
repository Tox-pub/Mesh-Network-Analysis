; MeSH AOP Workbench - Windows installer (Inno Setup 6)
;
; Standard Program Files install with an uninstaller, Start Menu entry and an
; optional desktop icon - the arrangement Windows users already know. The data
; folder is deliberately NOT placed under Program Files: the pipeline writes
; tens of gigabytes and Program Files is both UAC-protected and usually on the
; smallest drive. The user picks that folder on first launch instead.
;
; Build:  "C:\Program Files (x86)\Inno Setup 6\ISCC.exe" installer.iss

#define AppName    "MeSH Workbench"
#define AppVersion "3.0.0"
#define AppExe     "MeshWorkbench.exe"

[Setup]
AppId={{7C3F1A62-9B4E-4E77-9E2D-5A1D2B8C4F10}
AppName={#AppName}
AppVersion={#AppVersion}
AppPublisher=Karolinska Institutet
DefaultDirName={autopf}\{#AppName}
DefaultGroupName={#AppName}
DisableProgramGroupPage=yes
OutputDir=dist
OutputBaseFilename=MeshWorkbench-{#AppVersion}-win64-setup
Compression=lzma2/max
SolidCompression=yes
WizardStyle=modern
ArchitecturesInstallIn64BitMode=x64compatible
ArchitecturesAllowed=x64compatible
PrivilegesRequired=admin
UninstallDisplayIcon={app}\{#AppExe}
; The bundle is ~313 MB; refuse to start an install that cannot finish.
ExtraDiskSpaceRequired=52428800

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "Create a &desktop shortcut"; \
  GroupDescription: "Additional shortcuts:"; Flags: unchecked

[Files]
Source: "dist\MeshWorkbench\{#AppExe}"; DestDir: "{app}"; Flags: ignoreversion
Source: "dist\MeshWorkbench\_internal\*"; DestDir: "{app}\_internal"; \
  Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{group}\{#AppName}"; Filename: "{app}\{#AppExe}"
Name: "{group}\Uninstall {#AppName}"; Filename: "{uninstallexe}"
Name: "{autodesktop}\{#AppName}"; Filename: "{app}\{#AppExe}"; Tasks: desktopicon

[Run]
Filename: "{app}\{#AppExe}"; Description: "Launch {#AppName}"; \
  Flags: nowait postinstall skipifsilent

[Messages]
WelcomeLabel2=This will install [name/ver] on your computer.%n%nThe application itself is about 313 MB. The PubMed data it works from is much larger - roughly 44 GB to download once, and about 11 GB kept afterwards. The application downloads that separately and lets you resume if interrupted.%n%nChoose a data folder on a drive with room when you first run it.

[Code]
{ Warn once, before any files are written, if the target drive is tight. The
  bundle needs ~350 MB but a user who installs here will point the data folder
  somewhere too, and finding that out after a 44 GB download is far worse. }
function NextButtonClick(CurPageID: Integer): Boolean;
var
  FreeMB: Cardinal;
  TotalMB: Cardinal;
begin
  Result := True;
  if CurPageID = wpSelectDir then
  begin
    if GetSpaceOnDisk(ExtractFileDrive(WizardDirValue), True, FreeMB, TotalMB) then
// GetSpaceOnDisk returns MB when the second argument is True
      if FreeMB < 1024 then
        Result := MsgBox('That drive has only ' + IntToStr(FreeMB) +
          ' MB free.' + #13#10#13#10 +
          'The application needs about 350 MB, and you will also need somewhere ' +
          'with roughly 44 GB for the PubMed data - though that can be a ' +
          'different drive, chosen when the application first runs.' + #13#10#13#10 +
          'Continue anyway?', mbConfirmation, MB_YESNO) = IDYES;
  end;
end;
