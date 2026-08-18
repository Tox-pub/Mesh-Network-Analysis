# -*- coding: utf-8 -*-
"""
mesh_workbench - desktop front-end for the mesh_aop pipeline.

A window over the same configuration and the same pipeline steps the
`mesh-pipeline` command exposes. Settings are read from and written back to
mesh_config.json, and each step runs in a child process so a long analysis
cannot block the interface.

Entry points
    mesh-workbench                 console script (see pyproject.toml)
    python -m mesh_workbench       equivalent
"""

__version__ = '3.0.0'
__all__ = ['app', 'runner', 'settings_schema']
