; MeSH Workbench - Windows installer (Inno Setup 6 or 7)
;
; Installs the SAME tree the portable zip contains. Nothing here is compiled:
; the only program that ever executes is python.exe from the official embeddable
; distribution, signed by the Python Software Foundation. An earlier version of
; this script packaged a PyInstaller build, and that unsigned binary was refused
; outright on a managed Windows machine with no error shown - which is the whole
; reason the portable route exists. This gives the conventional install
; experience without reintroducing that problem.
;
; Setup.exe itself is unsigned, so SmartScreen warns once on first run. Only a
; code-signing certificate removes that, and it applies to the installer alone -
; never to the installed application.
;
; Build:  python packaging/build_installer_windows.py
; (which compiles this with /DPortableDir=... pointing at a finished build)

#define AppName    "MeSH Workbench"
#define AppVersion "3.2.0"
#define AppIcon    "app\mesh_workbench\assets\mesh_workbench.ico"

; Where the finished portable tree lives. Passed in by the build script; the
; fallback only makes a hand-run compile fail with a clear message instead of
; silently producing an empty installer.
#ifndef PortableDir
  #define PortableDir "..\..\mesh_workbench_build\MeshWorkbench"
#endif

[Setup]
AppId={{7C3F1A62-9B4E-4E77-9E2D-5A1D2B8C4F10}
AppName={#AppName}
AppVersion={#AppVersion}
AppPublisher=Karolinska Institutet
AppSupportURL=https://github.com/Tox-pub/Mesh-Network-Analysis
DefaultDirName={autopf}\{#AppName}
DefaultGroupName={#AppName}
DisableProgramGroupPage=yes
OutputDir=dist
OutputBaseFilename=MeSH-Workbench-{#AppVersion}-win64-setup
Compression=lzma2/max
SolidCompression=yes
WizardStyle=modern
ArchitecturesInstallIn64BitMode=x64compatible
ArchitecturesAllowed=x64compatible
; Per-user by default, so no administrator prompt and no blocked install on a
; managed machine. The wizard still offers "for all users" to anyone who has the
; rights and wants it.
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog
UninstallDisplayIcon={app}\{#AppIcon}
UninstallDisplayName={#AppName} {#AppVersion}
SetupIconFile={#PortableDir}\{#AppIcon}
; The installed tree is ~460 MB. Refuse to start an install that cannot finish.
ExtraDiskSpaceRequired=104857600

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "Create a &desktop shortcut"; \
  GroupDescription: "Additional shortcuts:"

[Files]
; The .bat launchers and their README belong to the portable zip, where there
; are no shortcuts and no Add/Remove Programs entry. An installed copy starts
; from the Start menu and uninstalls through Windows, so shipping them here
; would only be more things to click by mistake - and the README would be
; telling the user to click files this install does not contain.
;
; portable.marker is excluded for a different reason: its presence is what tells
; the application to keep settings and results in its own folder. An installed
; copy must not have it, or every user of this machine would share one
; configuration and one results folder - and under Program Files it could not
; write either of them at all.
Source: "{#PortableDir}\*"; DestDir: "{app}"; \
  Excludes: "*.bat,README - Install and First Run.txt,portable.marker"; \
  Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
; Shortcuts run pythonw.exe directly rather than going through a .bat: no
; console window flashes, and the shortcut carries the application icon like any
; other installed program.
Name: "{group}\{#AppName}"; Filename: "{app}\python\pythonw.exe"; \
  Parameters: """{app}\app\launch.py"""; WorkingDir: "{app}"; \
  IconFilename: "{app}\{#AppIcon}"; Comment: "Build and validate MeSH concept networks"
Name: "{group}\Uninstall {#AppName}"; Filename: "{uninstallexe}"
Name: "{autodesktop}\{#AppName}"; Filename: "{app}\python\pythonw.exe"; \
  Parameters: """{app}\app\launch.py"""; WorkingDir: "{app}"; \
  IconFilename: "{app}\{#AppIcon}"; Tasks: desktopicon

[Run]
Filename: "{app}\python\pythonw.exe"; Parameters: """{app}\app\launch.py"""; \
  WorkingDir: "{app}"; Description: "Launch {#AppName}"; \
  Flags: nowait postinstall skipifsilent

[Messages]
WelcomeLabel2=This will install [name/ver] on your computer.%n%nThe application is about 460 MB. The PubMed data it works from is much larger - roughly 44 GB downloaded once, and about 8 GB kept afterwards. The application fetches that separately, and lets you choose which drive it goes on.

[Code]
{ Warn before any files are written if the target drive is tight. The
  application needs ~500 MB, but a user installing here will probably point the
  data folder somewhere too, and discovering that after a 44 GB download is far
  worse than being asked now. }
function NextButtonClick(CurPageID: Integer): Boolean;
var
  FreeMB, TotalMB: Cardinal;
begin
  Result := True;
  if CurPageID = wpSelectDir then
  begin
    { GetSpaceOnDisk returns MB when the second argument is True }
    if GetSpaceOnDisk(ExtractFileDrive(WizardDirValue), True, FreeMB, TotalMB) then
      if FreeMB < 1024 then
        Result := MsgBox('That drive has only ' + IntToStr(FreeMB) +
          ' MB free.' + #13#10#13#10 +
          'The application needs about 500 MB, and you will also need somewhere ' +
          'with roughly 52 GB for the PubMed data - though that can be a ' +
          'different drive, chosen inside the application.' + #13#10#13#10 +
          'Continue anyway?', mbConfirmation, MB_YESNO) = IDYES;
  end;
end;

{ Windows removes only what Setup installed. Everything the application later
  downloads or builds - the archive, the databases, and the scratch copy it
  leaves in the system temp folder - would survive an ordinary uninstall and is
  measured in tens of gigabytes. Ask, and hand the job to the application's own
  uninstaller so both routes agree on what "clean" means. Results are kept. }
procedure CurUninstallStepChanged(CurUninstallStep: TUninstallStep);
var
  ResultCode: Integer;
  PythonExe, Args: String;
begin
  if CurUninstallStep = usUninstall then
  begin
    PythonExe := ExpandConstant('{app}\python\python.exe');
    if not FileExists(PythonExe) then
      Exit;

    { Data first, and default to removing it: it is machine-generated, it is the
      bulk of what is on disk, and leaving tens of gigabytes behind after an
      uninstall is the behaviour users complain about. }
    if MsgBox('Also remove the downloaded PubMed data and the databases built' +
              ' from it?' + #13#10#13#10 +
              'This can be tens of gigabytes, and includes a working copy the' +
              ' build leaves in the Windows temp folder that removing the' +
              ' program would not otherwise clear.' + #13#10#13#10 +
              'This cannot be undone.',
              mbConfirmation, MB_YESNO) = IDNO then
      Exit;

    Args := ExpandConstant('-m mesh_aop.uninstall_cli --project "{app}" --yes');

    { Results are asked about separately, and separately declined by default.
      They are the user's own work, they live in their documents rather than in
      the install, and nothing else can reproduce them. }
    if MsgBox('Remove your results as well?' + #13#10#13#10 +
              'These are the figures, workbooks and reports this program' +
              ' produced for you, kept in your own documents folder.' + #13#10#13#10 +
              'Choose No to keep them.',
              mbConfirmation, MB_YESNO) = IDYES then
      Args := Args + ' --results';

    Exec(PythonExe, Args, ExpandConstant('{app}'), SW_SHOW,
         ewWaitUntilTerminated, ResultCode);
  end;
end;
