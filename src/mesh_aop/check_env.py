# -*- coding: utf-8 -*-
"""
check_env.py - environment and dependency verification (the `mesh-check-env` command).

Validates that the active interpreter and installed packages can run the
pipeline before a long job is started.

It checks the Python version against the supported range, verifies that every
dependency declared in pyproject.toml is present and satisfies its version
constraint, installs the package in editable mode when requested, and resolves
the common `community` / `python-louvain` name collision. It also scans for
cross-package conflicts and, on Linux, provisions the OS-level libraries needed
for headless figure rendering.
"""

import sys
import subprocess
import os
import re
from importlib.metadata import version, PackageNotFoundError

# tomllib is native to Python 3.11+
try:
    import tomllib
except ImportError:
    print("[!] CRITICAL ERROR: tomllib not found. Python 3.11+ is strictly required.")
    sys.exit(1)

# packaging ships with pip/setuptools; used to verify that installed versions
# actually satisfy the pyproject specifiers. Degrade gracefully if it's absent.
try:
    from packaging.requirements import Requirement
    _HAS_PACKAGING = True
except ImportError:
    _HAS_PACKAGING = False

def check_environment(toml_path="pyproject.toml", auto_install=False):
    """
    Checks Python version and utilizes importlib.metadata to verify dependencies safely
    against the pyproject.toml specifications.
    """
    print("\n" + "<"*30 + ">"*30)
    print("   SYSTEM ENVIRONMENT & DEPENDENCY CHECK")
    print("<"*30 + ">"*30 + "\n")

    # 1. CHECK PYTHON VERSION
    current_py = sys.version_info
    if current_py >= (3, 13):
        print(f"[!] WARNING: Python {sys.version.split()[0]} is NOT supported by this pipeline.")
        print("    A core dependency (node2vec) requires numpy<2.0, and numpy<2.0 has NO")
        print("    prebuilt wheel for Python 3.13+. Installation will try to COMPILE numpy")
        print("    from source and fail unless a C/C++ compiler (MSVC) is installed.")
        print("    >>> Fix: install Python 3.11 or 3.12 and rebuild the virtual environment.")
        if not auto_install:
            try:
                input("    Press Enter to attempt anyway (likely to fail), or Ctrl+C to exit...")
            except KeyboardInterrupt:
                print("\n\n[!] Exiting: Environment check cancelled by user.")
                return False
    elif current_py < (3, 11):
        print(f"[!] WARNING: Python {sys.version.split()[0]} is too old; Python 3.11 or 3.12 is required.")
        if not auto_install:
            try:
                input("    Press Enter to attempt to continue, or Ctrl+C to exit...")
            except KeyboardInterrupt:
                print("\n\n[!] Exiting: Environment check cancelled by user.")
                return False

    # 2. CHECK PYPROJECT FILE
    if not os.path.exists(toml_path):
        print(f"\n[!] ERROR: '{toml_path}' not found.")
        print("    Ensure you are running this script from the project root.")
        return False

    try:
        with open(toml_path, "rb") as f:
            project_data = tomllib.load(f)
        required_lines = project_data.get("project", {}).get("dependencies", [])
    except Exception as e:
        print(f"\n[!] ERROR: Failed to parse '{toml_path}'.")
        print(f"    Details: {e}")
        return False

    if not required_lines:
        print(f"[!] Notice: No dependencies found in '{toml_path}'. Skipping dependency check.")
        return True

    # 3. VERIFY INDIVIDUAL DEPENDENCIES (presence AND version satisfaction)
    to_install = []
    version_conflicts = []
    print(f"Analyzing {len(required_lines)} dependencies...\n")

    for req in required_lines:
        # Resolve the distribution name and version specifier from the requirement.
        specifier = None
        if _HAS_PACKAGING:
            try:
                req_obj = Requirement(req)
                pkg_name = req_obj.name
                specifier = req_obj.specifier
            except Exception:
                pkg_name = re.split(r'[=><~!\[;]', req, maxsplit=1)[0].strip()
        else:
            pkg_name = re.split(r'[=><~!\[;]', req, maxsplit=1)[0].strip()

        try:
            installed_ver = version(pkg_name)
        except PackageNotFoundError:
            print(f"    [!] Missing completely: {req}")
            to_install.append(req)
            continue
        except Exception as e:
            print(f"    [!] Error checking {pkg_name}: {e}")
            to_install.append(req)
            continue

        # Installed, but verify the version actually satisfies our constraint.
        if specifier is not None and len(specifier) and installed_ver not in specifier:
            print(f"    [!] VERSION CONFLICT: {pkg_name} {installed_ver} is installed, "
                  f"but this project requires {pkg_name}{specifier}.")
            version_conflicts.append(f"{pkg_name} {installed_ver} (needs {specifier})")
            to_install.append(req)

    if version_conflicts:
        print(f"\n[!] {len(version_conflicts)} installed package(s) do not satisfy the project's version constraints:")
        for c in version_conflicts:
            print(f"        - {c}")
        print("    These will be reinstalled to a compatible version below.")

    # 3.5 TRAP CHECK: COMMUNITY VS PYTHON-LOUVAIN
    try:
        import community
        if not hasattr(community, 'best_partition'):
            print("\n    [!] NAME COLLISION DETECTED: The generic 'community' package is shadowing 'python-louvain'.")
            print("        Uninstalling the conflicting package now...")
            subprocess.run([sys.executable, "-m", "pip", "uninstall", "-y", "community"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if "python-louvain" not in to_install:
                to_install.append("python-louvain")
    except ImportError:
        pass # Handled by the standard requirements check

    # 4. TARGETED PACKAGE INSTALLATION
    if to_install:
        print(f"\n[!] Found {len(to_install)} missing or unverified dependencies.")

        if not auto_install:
            try:
                choice = input("    Run package installation (pip install -e .) now? (y/n): ").strip().lower()
                if choice != 'y':
                    print("\n[!] Exiting: User declined to install required dependencies.")
                    return False
            except KeyboardInterrupt:
                print("\n\n[!] Exiting: Installation prompt cancelled by user.")
                return False
        else:
            print("    Auto-install flag detected. Proceeding with installation.")

        print("\nStarting package installation...")
        try:
            # Install the local directory as an editable package, resolving all dependencies
            result = subprocess.run(
                [sys.executable, "-m", "pip", "install", "-e", "."],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
            )
            if result.returncode != 0:
                print(f"\n[!] CRITICAL ERROR: Package installation failed.")
                print("--- PIP ERROR LOG ---")
                print(result.stderr.strip())
                print("---------------------")
                return False
        except Exception as e:
            print(f"\n[!] CRITICAL ERROR: Unexpected system failure installing package.")
            print(f"    Details: {e}")
            return False
        print("\n[+] Package and dependencies installed successfully.")

    # 5. CROSS-PACKAGE CONFLICT SCAN (catches conflicts from ANY installed library)
    _report_pip_conflicts()

    print("\n[+] Environment and dependencies are fully verified.")
    return True

def _report_pip_conflicts():
    """Runs `pip check` to surface dependency conflicts among ALL installed packages.

    This catches the general case the user cares about: any library whose
    installed version is incompatible with another package's requirements (the
    same class of problem as node2vec pinning numpy<2). pip check inspects the
    whole environment, not just this project's direct dependencies.
    """
    print("\n<<< Scanning for cross-package dependency conflicts (pip check) >>>")
    try:
        result = subprocess.run(
            [sys.executable, "-m", "pip", "check"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        output = (result.stdout + result.stderr).strip()
        if result.returncode == 0:
            print("    [+] No dependency conflicts detected.")
        else:
            print("    [!] WARNING: dependency conflicts detected among installed libraries:")
            for line in output.splitlines():
                if line.strip():
                    print(f"        - {line.strip()}")
            print("    These can cause import errors or incorrect behaviour at runtime.")
            print("    Resolve by installing compatible versions (e.g. on Python 3.13 several")
            print("    scientific packages have no wheels — use Python 3.11/3.12 and reinstall).")
    except Exception as e:
        print(f"    [!] Could not run conflict scan: {e}")

def check_windows_activation_policy():
    """On Windows, warn if the PowerShell execution policy blocks venv activation.

    A 'Restricted'/'AllSigned' policy makes `Activate.ps1` fail with
    "running scripts is disabled on this system", which stops users from
    activating the virtual environment for the shorthand `python ...` commands.
    """
    if sys.platform != "win32":
        return
    try:
        result = subprocess.run(
            ["powershell", "-NoProfile", "-Command", "Get-ExecutionPolicy"],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=15
        )
        policy = (result.stdout or "").strip()
    except Exception:
        return  # PowerShell unavailable; nothing to warn about

    if policy.lower() in ("restricted", "allsigned"):
        print("\n" + "<"*30 + ">"*30)
        print(f"[!] WINDOWS NOTE: PowerShell execution policy is '{policy}'.")
        print("<"*30 + ">"*30)
        print("    This BLOCKS virtual-environment activation. If you see")
        print("    'Activate.ps1 cannot be loaded ... running scripts is disabled',")
        print("    enable local scripts once (no admin required):")
        print("        Set-ExecutionPolicy -Scope CurrentUser -ExecutionPolicy RemoteSigned")
        print("    Then: & \"$env:USERPROFILE\\mesh_env\\Scripts\\Activate.ps1\"")
        print("    Or skip activation entirely and call the venv Python by full path:")
        print("        & \"$env:USERPROFILE\\mesh_env\\Scripts\\python.exe\" -m mesh_aop.cli --step all --interactive")

def provision_kaleido_dependencies():
    """Installs required OS-level shared libraries and the headless Chrome binary."""
    if not sys.platform.startswith('linux'):
        print("\n[!] Non-Linux OS detected. Skipping OS-level Kaleido dependencies.")
        return

    print("\nProvisioning system dependencies for headless Chrome rendering...")

    # 1. Install Linux shared libraries
    try:
        env = os.environ.copy()
        env["DEBIAN_FRONTEND"] = "noninteractive"

        subprocess.run(["apt-get", "update"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env)
        libs = [
            "libnss3", "libatk-bridge2.0-0", "libxcomposite1",
            "libxdamage1", "libxrandr2", "libgbm1", "libasound2"
        ]
        subprocess.run(["apt-get", "install", "-y"] + libs, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, env=env)
        print("    [+] System libraries provisioned.")
    except subprocess.CalledProcessError as e:
        print(f"    [!] Failed to install system libraries: {e}")
        print("        Note: If not running as root/Colab, this requires sudo.")

    # 2. Install Headless Chrome via Plotly
    try:
        print("    Downloading headless Chromium binary (auto-confirming prompt)...")
        
        subprocess.run(
            ["plotly_get_chrome"],
            check=True,
            input="y\n",
            text=True
        )
        print("    [+] Headless Chromium binary installed.")
    except subprocess.CalledProcessError:
        print(f"    [!] Warning: Failed to install headless Chromium.")
        print("        Static HTML files will be generated as fallbacks.")
    except FileNotFoundError:
        print("    [!] Error: 'plotly_get_chrome' command not found.")
        print("        Ensure 'plotly >= 6.1.1' is installed via pyproject.toml first.")

def main():
    """Entry point for `mesh-check-env`: verify the environment, then provision rendering libraries."""
    auto_mode = "--auto" in sys.argv
    try:
        success = check_environment(auto_install=auto_mode)
        if not success:
            sys.exit(1)
    except KeyboardInterrupt:
        print("\n\n[!] Exiting: Script manually interrupted by user.")
        sys.exit(130)

    check_windows_activation_policy()
    provision_kaleido_dependencies()

if __name__ == "__main__":
    main()