# -*- coding: utf-8 -*-
"""
Launch the Workbench.

Resolves the two things the window needs and cannot infer reliably: the project
directory holding mesh_config.json, and the interpreter that should run pipeline
steps. Both can be overridden by environment variable, which is how the portable
build points them at its own bundled Python.
"""

import os
import sys


def _project_dir():
    """Where mesh_config.json lives.

    Prefer an explicit override, then the working directory if it looks like a
    project, then the repository this file was installed from - so running the
    command from anywhere still finds the right configuration.
    """
    env = os.environ.get('MESH_REPO')
    if env and os.path.isdir(env):
        return env
    if os.path.exists(os.path.join(os.getcwd(), 'mesh_config.json')):
        return os.getcwd()
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(os.path.dirname(here))          # src/mesh_workbench -> repo
    return root if os.path.exists(os.path.join(root, 'mesh_config.json')) else os.getcwd()


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv and argv[0] in ('-h', '--help'):
        print(__doc__.strip())
        print('\n  --project DIR   directory containing mesh_config.json')
        return 0
    project = _project_dir()
    if '--project' in argv:
        project = argv[argv.index('--project') + 1]
    python_exe = os.environ.get('MESH_PYTHON') or sys.executable
    from .app import main as run_app
    run_app(project, python_exe)
    return 0


if __name__ == '__main__':
    sys.exit(main())
