#!/usr/bin/env bash
#
# install.sh - install MeSH Workbench on macOS or Linux.
#
# There is no compiled binary and nothing to sign. The application is Python,
# so this creates a private virtual environment, installs the package into it,
# and puts a launcher where the desktop expects to find one:
#
#   Linux   ~/.local/share/applications/mesh-workbench.desktop   (appears in the
#           application menu, with the icon)
#   macOS   /Applications/MeSH Workbench.app                     (a real bundle,
#           launchable from Finder and Spotlight)
#
# Everything goes under the user's own home directory, so no administrator
# rights are needed and nothing is shared between accounts - the same layout the
# Windows build uses.
#
# Usage:
#   ./install.sh                     install for the current user
#   ./install.sh --prefix DIR        put the environment somewhere else
#   ./install.sh --uninstall         remove it again
#
set -euo pipefail

APP_NAME="MeSH Workbench"
SLUG="mesh-workbench"
MIN_MINOR=11          # requires-python in pyproject.toml is >=3.11,<3.14
MAX_MINOR=13

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PREFIX="${XDG_DATA_HOME:-$HOME/.local/share}/$SLUG"
UNINSTALL=0

while [ $# -gt 0 ]; do
    case "$1" in
        --prefix)    PREFIX="$2"; shift 2 ;;
        --uninstall) UNINSTALL=1; shift ;;
        -h|--help)   sed -n '3,25p' "$0"; exit 0 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done

case "$(uname -s)" in
    Darwin) PLATFORM=macos ;;
    Linux)  PLATFORM=linux ;;
    *)      echo "This installer supports macOS and Linux. On Windows use the" >&2
            echo "zip or the installer from the releases page." >&2; exit 1 ;;
esac

DESKTOP_FILE="$HOME/.local/share/applications/$SLUG.desktop"
APP_BUNDLE="/Applications/$APP_NAME.app"
[ -w /Applications ] || APP_BUNDLE="$HOME/Applications/$APP_NAME.app"

# ---------------------------------------------------------------- uninstall
if [ "$UNINSTALL" = "1" ]; then
    echo "Removing $APP_NAME"
    # The application's own uninstaller knows about downloaded data and caches,
    # which live outside this directory; run it before the environment goes.
    if [ -x "$PREFIX/venv/bin/mesh-uninstall" ]; then
        "$PREFIX/venv/bin/mesh-uninstall" --list || true
        echo
        echo "  To remove downloaded data and results as well, run:"
        echo "    $PREFIX/venv/bin/mesh-uninstall"
        echo "  (do that first if you want it - the command disappears with the next step)"
        echo
        printf "  Continue removing the application? [y/N] "
        read -r reply
        case "$reply" in [Yy]*) ;; *) echo "  cancelled"; exit 1 ;; esac
    fi
    rm -rf "$PREFIX"                       && echo "  removed $PREFIX"
    [ -e "$DESKTOP_FILE" ] && rm -f "$DESKTOP_FILE" && echo "  removed $DESKTOP_FILE"
    [ -e "$APP_BUNDLE" ]   && rm -rf "$APP_BUNDLE"  && echo "  removed $APP_BUNDLE"
    echo "Done."
    exit 0
fi

# ------------------------------------------------------------ find a python
# Newest first within the supported range: 3.13 is fine, 3.14 is not tested and
# the dependency pins in pyproject.toml would reject it anyway.
PYTHON=""
for minor in $(seq $MAX_MINOR -1 $MIN_MINOR); do
    if command -v "python3.$minor" >/dev/null 2>&1; then
        PYTHON="python3.$minor"; break
    fi
done
if [ -z "$PYTHON" ] && command -v python3 >/dev/null 2>&1; then
    v=$(python3 -c 'import sys; print(sys.version_info[1])')
    if [ "$v" -ge "$MIN_MINOR" ] && [ "$v" -le "$MAX_MINOR" ]; then PYTHON=python3; fi
fi
if [ -z "$PYTHON" ]; then
    echo "No suitable Python found." >&2
    echo "Need Python 3.$MIN_MINOR to 3.$MAX_MINOR. Install one, then run this again:" >&2
    if [ "$PLATFORM" = macos ]; then
        echo "    brew install python@3.12          # or python.org/downloads" >&2
    else
        echo "    sudo apt install python3.12 python3.12-venv python3-tk" >&2
        echo "    (or the equivalent for your distribution)" >&2
    fi
    exit 1
fi
echo "Using $PYTHON ($("$PYTHON" -c 'import sys; print(sys.version.split()[0])'))"

# tkinter ships separately on many Linux distributions and its absence is only
# discovered when the window fails to open - check now, while it can be fixed.
if ! "$PYTHON" -c 'import tkinter' >/dev/null 2>&1; then
    echo "This Python has no tkinter, which the application's window needs." >&2
    if [ "$PLATFORM" = linux ]; then
        echo "    sudo apt install python3-tk        # Debian/Ubuntu" >&2
        echo "    sudo dnf install python3-tkinter   # Fedora" >&2
    else
        echo "    brew install python-tk" >&2
    fi
    exit 1
fi

# ------------------------------------------------------------------ install
echo "Installing into $PREFIX"
mkdir -p "$PREFIX"
"$PYTHON" -m venv "$PREFIX/venv"
"$PREFIX/venv/bin/python" -m pip install --quiet --upgrade pip
echo "  resolving dependencies (this takes a few minutes the first time)"
"$PREFIX/venv/bin/python" -m pip install --quiet "$REPO"

BIN="$PREFIX/venv/bin/mesh-workbench"
[ -x "$BIN" ] || { echo "install finished but $BIN is missing" >&2; exit 1; }

ICON_SRC="$REPO/src/mesh_workbench/assets/mesh_workbench.png"
ICON="$PREFIX/mesh_workbench.png"
[ -f "$ICON_SRC" ] && cp "$ICON_SRC" "$ICON"

if [ "$PLATFORM" = linux ]; then
    mkdir -p "$(dirname "$DESKTOP_FILE")"
    cat > "$DESKTOP_FILE" <<EOF
[Desktop Entry]
Type=Application
Name=$APP_NAME
Comment=Build and validate MeSH concept networks
Exec=$BIN
Icon=$ICON
Terminal=false
Categories=Science;Education;
EOF
    chmod +x "$DESKTOP_FILE"
    command -v update-desktop-database >/dev/null 2>&1 && \
        update-desktop-database "$(dirname "$DESKTOP_FILE")" >/dev/null 2>&1 || true
    echo "  menu entry: $DESKTOP_FILE"
else
    # A minimal .app: Finder and Spotlight only need the bundle layout and a
    # plist, and the "executable" may be a shell script.
    mkdir -p "$APP_BUNDLE/Contents/MacOS" "$APP_BUNDLE/Contents/Resources"
    cat > "$APP_BUNDLE/Contents/Info.plist" <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0"><dict>
  <key>CFBundleName</key><string>$APP_NAME</string>
  <key>CFBundleDisplayName</key><string>$APP_NAME</string>
  <key>CFBundleIdentifier</key><string>se.ki.meshworkbench</string>
  <key>CFBundleVersion</key><string>3.1.0</string>
  <key>CFBundleShortVersionString</key><string>3.1.0</string>
  <key>CFBundleExecutable</key><string>run</string>
  <key>CFBundlePackageType</key><string>APPL</string>
  <key>NSHighResolutionCapable</key><true/>
</dict></plist>
EOF
    printf '#!/bin/sh\nexec "%s" "$@"\n' "$BIN" > "$APP_BUNDLE/Contents/MacOS/run"
    chmod +x "$APP_BUNDLE/Contents/MacOS/run"
    [ -f "$ICON" ] && cp "$ICON" "$APP_BUNDLE/Contents/Resources/icon.png"
    echo "  application: $APP_BUNDLE"
fi

cat <<EOF

Installed.

  Launch it        from your application menu, or run:
                     $BIN
  Command line     $PREFIX/venv/bin/mesh-pipeline --help
  Remove it        $0 --uninstall

Settings and logs live under your user profile; results and data are chosen
the first time the application runs.
EOF
