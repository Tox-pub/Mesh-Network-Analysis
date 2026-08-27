# -*- coding: utf-8 -*-
"""
app.py - the MeSH AOP Workbench window.

Four screens in one fixed-size window: Data setup, Settings, Running, Results.
The settings form is generated from schema.TABS, so the UI never hard-codes a
configuration key; the description pane underneath reports what the focused
control does, its default, and the caveat that matters.

The pipeline runs in a child process (see runner.py). The UI polls a queue on a
timer rather than blocking, so the window stays live for the whole run - which
for the benchmark step is over an hour.

Styling is deliberately plain Win32: grey face, ridged bevels, no themed widgets.
That is the look the tool was asked for, and it happens to be what tkinter gives
natively on Windows.
"""

import json
import os
import shutil
import subprocess
import sys
import tkinter as tk
import tkinter.font as tkfont
from tkinter import filedialog, messagebox

from . import settings_schema as schema
from . import __version__
from .runner import (AT_CHECKPOINT, DONE, LOG, PAUSED, PHASE, PROGRESS,
                     TRANSIENT, PipelineRunner)

FACE = '#c0c0c0'
FIELD = '#ffffff'
FOCUS = '#ffffe1'
INK = '#000000'
DIM = '#404040'
OK = '#006400'
ERR = '#800000'
WARN = '#7a5c00'
NAVY = '#000080'
CONSOLE_BG = '#000000'
# The description pane at the foot of Settings has to hold the longest note in
# the schema without clipping it - a caveat that is cut off mid-sentence is
# worse than absent, because the reader cannot tell there was more. HELP_H is
# measured against every field by the test suite rather than guessed.
W, H = 880, 720
HELP_H = 150


class Workbench(tk.Tk):
    # Folder fields worth showing a resolved path in rather than a blank. These
    # are the two the first-run dialog asks about, and the two that decide where
    # tens of gigabytes land - so leaving them empty hides the one fact the user
    # most needs. input_dir, output_dir and the ETL workspace are deliberately
    # NOT here: they are overrides, and filling one in pins it, so a later change
    # to the data folder would silently fail to move anything.
    _PREFILL_DIRS = {}

    def __init__(self, repo_dir, python_exe):
        super().__init__()
        self.repo_dir = repo_dir
        self.python_exe = python_exe
        # A portable copy keeps its settings in its own folder; an installed one
        # keeps them under the user's profile, because the install directory is
        # shared and often read-only. paths.py decides, and the CLI asks it too.
        from mesh_aop import paths as _paths
        self.paths = _paths
        self._PREFILL_DIRS = {
            'directories.results_dir': _paths.default_user_results_dir,
            'directories.data_dir': _paths.default_user_data_dir,
        }
        self.cfg_path = str(_paths.config_path(repo_dir))
        self.cfg = self._load_cfg()
        self.vars = {}
        self.runner = PipelineRunner(repo_dir, python_exe)
        self._elapsed = 0
        self._tick_job = None
        # Set when a run stops for AOP annotation, so _finish can tell an
        # intentional pause apart from an ordinary success - both exit zero.
        self._paused_annotation = None
        self._scan_rows = []

        self.title('MeSH Workbench')
        self.configure(bg=FACE)
        self.geometry(f'{W}x{H}')
        self.minsize(W, H)
        self.protocol('WM_DELETE_WINDOW', self._on_close)
        self._set_icon()

        self.f_ui = tkfont.Font(family='MS Sans Serif', size=8)
        self.f_bold = tkfont.Font(family='MS Sans Serif', size=8, weight='bold')
        self.f_mono = tkfont.Font(family='Consolas', size=9)
        self.option_add('*Font', self.f_ui)

        self._menu()
        self.body = tk.Frame(self, bg=FACE)
        self.body.pack(fill='both', expand=True)
        self.screens = {}
        for name, build in (('setup', self._build_setup), ('settings', self._build_settings),
                            ('running', self._build_running), ('results', self._build_results),
                            ('integrity', self._build_integrity),
                            ('uninstall', self._build_uninstall)):
            fr = tk.Frame(self.body, bg=FACE)
            self.screens[name] = fr
            build(fr)
        self.show('setup')
        # After the window exists, so the dialog has a parent to sit over.
        self.after(150, self.first_run_locations)

    # ---------------------------------------------------------------- config
    def _load_cfg(self):
        """Read the settings file, creating it from the defaults on first run.

        A missing file is the normal state of a fresh install, not a failure -
        reporting it as one greeted every new user with "[Errno 2] No such file
        or directory". A file that exists but cannot be parsed is a real problem
        and still says so.
        """
        # Remembered, because this is what tells first_run_locations that this
        # is a first run. Testing for the file later cannot work: seeding the
        # defaults here creates it, so by the time the dialog asked, the answer
        # was always "it exists" and the dialog never appeared at all - nobody
        # was ever asked where to put the ~52 GB it is there to ask about.
        self.first_run = not os.path.exists(self.cfg_path)
        if self.first_run:
            return self._write_default_cfg()
        try:
            with open(self.cfg_path, encoding='utf-8') as fh:
                return json.load(fh)
        except (OSError, ValueError) as exc:
            messagebox.showerror(
                'Configuration',
                f'{self.cfg_path} could not be read, so the defaults are in '
                f'use. Your settings are not lost - fix or delete the file and '
                f'restart.\n\n{exc}')
            return {}

    def _write_default_cfg(self):
        """Seed a settings file from the pipeline's own factory defaults.

        Derived rather than shipped as a second copy: a checked-in defaults file
        would drift from MeshConfig the first time a default changed, and the
        two front-ends would disagree again. Transient run switches (the
        underscore-prefixed ones) are never written.
        """
        try:
            from mesh_aop.config_parser import MeshConfig
            defaults = MeshConfig(config_path='\0').factory_defaults
        except Exception:
            return {}
        cfg = {k: v for k, v in defaults.items() if not k.startswith('_')}
        try:
            os.makedirs(os.path.dirname(self.cfg_path) or '.', exist_ok=True)
            with open(self.cfg_path, 'w', encoding='utf-8') as fh:
                json.dump(cfg, fh, indent=4)
        except OSError:
            pass          # read-only location: run from the defaults in memory
        return cfg

    def save_cfg(self):
        for key, var in self.vars.items():
            fld = self._field(key)
            raw = var.get()
            try:
                if fld.kind == 'int':
                    val = int(float(raw))
                elif fld.kind == 'float':
                    val = float(raw)
                elif fld.kind == 'bool':
                    val = bool(raw)
                else:
                    val = str(raw)
            except (TypeError, ValueError):
                messagebox.showerror('Invalid value',
                                     f'"{fld.label}" needs a {fld.kind}; got "{raw}".')
                return False
            schema.put(self.cfg, key, val)
        bad = self._weight_sum()
        if abs(bad - 1.0) > 1e-6:
            if not messagebox.askyesno(
                    'Node weights do not total 1.00',
                    f'The four node weight factors total {bad:.2f}, not 1.00.\n\n'
                    'Node weights will not be normalised and scores will not be '
                    'comparable with other runs.\n\nSave anyway?'):
                return False
        try:
            # The per-user settings folder does not exist until something writes
            # to it, and on a first run that something is this.
            os.makedirs(os.path.dirname(self.cfg_path) or '.', exist_ok=True)
            with open(self.cfg_path, 'w', encoding='utf-8') as fh:
                json.dump(self.cfg, fh, indent=4)
        except Exception as exc:
            messagebox.showerror('Configuration', f'Could not write:\n\n{exc}')
            return False
        self.status.config(text='Saved to mesh_config.json')
        self._refresh_effective_paths()
        return True

    def _field(self, key):
        for _, fields in schema.TABS:
            for f in fields:
                if f.key == key:
                    return f
        raise KeyError(key)

    def _weight_sum(self):
        tot = 0.0
        for k in schema.WEIGHT_KEYS:
            try:
                tot += float(self.vars[k].get())
            except (KeyError, TypeError, ValueError):
                pass
        return tot

    # ------------------------------------------------------------------ chrome
    def _menu(self):
        bar = tk.Menu(self)
        m = tk.Menu(bar, tearoff=0)
        m.add_command(label='Save configuration', command=self.save_cfg)
        m.add_command(label='Restore defaults', command=self.restore_defaults)
        m.add_separator()
        m.add_command(label='Exit', command=self._on_close)
        bar.add_cascade(label='File', menu=m)
        r = tk.Menu(bar, tearoff=0)
        for step, label, _ in schema.STEPS:
            r.add_command(label=label, command=lambda s=step: self.start_run(s))
        bar.add_cascade(label='Run', menu=r)
        d = tk.Menu(bar, tearoff=0)
        d.add_command(label='Data setup…', command=lambda: self.show('setup'))
        bar.add_cascade(label='Database', menu=d)
        # The results screen was built and filled after every successful run but
        # nothing ever displayed it - there was no menu entry and no button.
        bar.add_command(label='Results', command=self.open_results)
        t = tk.Menu(bar, tearoff=0)
        # Held on the instance so their state can follow the run. An entry that
        # would fail if clicked is greyed rather than left enabled to raise an
        # error box - which is the difference between a program that tells you
        # where you are and one that lets you find out the hard way.
        t.add_command(label='Annotate AOP levels...', command=self.open_annotation)
        t.add_command(label='Show the annotation guide', command=self.open_annotation_guide)
        t.add_separator()
        t.add_command(label='Check and repair files...', command=self.open_integrity)
        t.add_separator()
        t.add_command(label='Uninstall...', command=self.open_uninstall)
        self.menu_tools = t
        self._tools_idx = {'annotate': 0, 'guide': 1, 'check': 3}
        bar.add_cascade(label='Tools', menu=t)

        h = tk.Menu(bar, tearoff=0)
        # The manual documents every part of what ships here - the pipeline
        # steps, the settings, the annotation workflow, the run controls, the
        # repair tool, the ledger and the PRISMA report - so it is first and
        # named for what it is rather than for its filename.
        h.add_command(label='The manual - everything this package does',
                      command=self.open_manual)
        h.add_separator()
        h.add_command(label='Installing and updating (INSTALL.md)',
                      command=lambda: self.open_doc('INSTALL.md'))
        h.add_command(label='Read me first (README.md)',
                      command=lambda: self.open_doc('README.md'))
        h.add_separator()
        h.add_command(label='Reference figures - what they are',
                      command=lambda: self.open_doc(
                          os.path.join('results', 'reference_figures',
                                       'PROVENANCE.md')))
        h.add_command(label='How to annotate AOP levels',
                      command=self.open_annotation_guide)
        h.add_separator()
        h.add_command(label='About', command=lambda: messagebox.showinfo(
            'About', f'MeSH Workbench {__version__}\n\nDesktop front-end for the '
                     'mesh_aop pipeline.'))
        bar.add_cascade(label='Help', menu=h)
        self.config(menu=bar)
        self._sync_tools_state()

    # ------------------------------------------------------------ tools state
    def _annotation_file(self):
        """The run-annotations file for the current prefix, or None if absent."""
        try:
            results = self.vars.get('directories.results_dir')
            results = results.get() if results else ''
            if not results:
                results = str(self.paths.default_user_results_dir(self.repo_dir))
            prefix = self._current_prefix()
            path = os.path.join(results, f'{prefix}_run_annotations.csv')
            return path if os.path.exists(path) else None
        except Exception:
            return None

    def _current_prefix(self):
        try:
            if self.vars['control_flags.use_reference_data'].get():
                return 'DAC_Mesh'
            return self.vars['control_flags.custom_file_prefix'].get() or 'DAC_Mesh'
        except Exception:
            return 'DAC_Mesh'

    def _sync_tools_state(self):
        """Grey out what cannot be done right now.

        The annotation entries only mean anything once a run has produced a file
        to annotate, which is exactly when the user starts looking for them.
        """
        menu = getattr(self, 'menu_tools', None)
        if menu is None:
            return
        anno = self._annotation_file()
        for name, present in (('annotate', bool(anno)), ('guide', bool(anno))):
            try:
                menu.entryconfig(self._tools_idx[name],
                                 state='normal' if present else 'disabled')
            except Exception:
                pass

    def open_annotation(self):
        anno = self._annotation_file()
        if not anno:
            messagebox.showinfo(
                'Nothing to annotate yet',
                'The annotation file is written once the network has been built '
                '(Step 3). Run that first, then come back here.')
            return
        self._open_path(anno)

    def open_annotation_guide(self):
        anno = self._annotation_file()
        if not anno:
            # Reachable from the Help menu before any run has produced a file,
            # so show the instructions rather than doing nothing.
            from mesh_aop import guides
            tmp = os.path.join(self._results_dir(),
                               f'{self._current_prefix()}_run_annotations.csv')
            path = guides.write_annotation_guide(tmp)
            if path:
                self._open_path(path)
            else:
                messagebox.showinfo(
                    'How to annotate',
                    'The annotation file is written once the network has been '
                    'built (Step 3). It lists every MeSH term in the network, '
                    'and you fill in the "aop_level" column with one of:\n\n'
                    + '\n'.join('    ' + lvl for lvl in guides.AOP_LEVELS)
                    + '\n\nKeep the file semicolon-delimited - MeSH headings '
                      'contain commas.')
            return
        from mesh_aop import guides
        guide = os.path.join(os.path.dirname(anno), guides.ANNOTATION_GUIDE_NAME)
        if not os.path.exists(guide):
            # Written when a run pauses; produce it on demand for a run that did
            # not pause but whose annotations the user wants to edit anyway.
            guides.write_annotation_guide(anno)
        self._open_path(guide) if os.path.exists(guide) else self._reveal(anno)

    def open_manual(self):
        self.open_doc('HELP.md')

    def open_doc(self, name):
        """Open one of the shipped documents, wherever this copy keeps them."""
        candidates = [
            os.path.join(self.repo_dir, name),
            os.path.join(os.path.dirname(self.repo_dir), name),
            os.path.join(self.repo_dir, 'docs', name),
        ]
        for path in candidates:
            if os.path.exists(path):
                self._open_path(path)
                return
        messagebox.showinfo(
            'Not found', f'{name} is not in this copy of the program.\n\n'
                         f'Looked in:\n' + '\n'.join('  ' + c for c in candidates))

    def first_run_locations(self):
        """Ask where results and data should live, the first time only.

        Shown when no settings file exists yet. Only the two locations the user
        genuinely owns are asked about - settings, logs and scratch are
        machinery and go to the profile without discussion.

        A portable copy is self-contained by definition, so it is never asked.
        """
        if self.paths.is_portable(self.repo_dir) or not getattr(self, 'first_run', False):
            return

        win = tk.Toplevel(self)
        # Only leaves via a button. It used to swallow the close box in silence,
        # which looks like a hung window; _make_modal answers instead.
        self._make_modal(win, 'MeSH Workbench - first run')

        tk.Label(win, text='Where should your files go?', bg=FACE, font=self.f_bold,
                 anchor='w').pack(fill='x', padx=12, pady=(12, 2))
        tk.Label(win, bg=FACE, anchor='w', justify='left', wraplength=520,
                 text='These can be changed later on the Folders tab. Settings and '
                      'logs are kept privately under your user profile and are not '
                      'asked about.').pack(fill='x', padx=12, pady=(0, 10))

        rows = (
            ('Results', 'Figures, workbooks and reports.',
             'directories.results_dir', self.paths.default_user_results_dir(self.repo_dir)),
            ('Data', 'Downloaded archives and the databases built from them. '
                     'This is the large one.',
             'directories.data_dir', self.paths.default_user_data_dir(self.repo_dir)),
        )
        chosen = {}
        for label, blurb, key, default in rows:
            box = self._sunken(win)
            box.pack(fill='x', padx=12, pady=4)
            tk.Label(box, text=label, bg=FACE, font=self.f_bold, anchor='w'
                     ).pack(fill='x', padx=8, pady=(6, 0))
            tk.Label(box, text=blurb, bg=FACE, fg=DIM, anchor='w', justify='left',
                     wraplength=500).pack(fill='x', padx=8)
            line = tk.Frame(box, bg=FACE)
            line.pack(fill='x', padx=8, pady=(4, 8))
            var = tk.StringVar(value=str(default))
            chosen[key] = var
            tk.Entry(line, textvariable=var, bg=FIELD, relief='sunken', bd=2
                     ).pack(side='left', fill='x', expand=True)
            tk.Button(line, text='…', width=3,
                      command=lambda v=var: self._browse_dir(v)).pack(side='left', padx=(4, 0))

        act = tk.Frame(win, bg=FACE)
        act.pack(fill='x', padx=12, pady=(4, 12))

        def apply(values):
            """Record the choice everywhere it has to appear, and make it real.

            Writing only to self.cfg was a silent data loss: the Folders tab is
            bound to self.vars, and save_cfg writes every var back into the
            config - so the empty field overwrote the chosen path with "" the
            first time the user pressed Run. The answer has to land in the
            form as well, or it does not survive.
            """
            for key, value in values.items():
                value = str(value).strip()
                if not value:
                    continue
                value = os.path.abspath(os.path.expanduser(value))
                if not self._ensure_writable(value):
                    return False
                schema.put(self.cfg, key, value)
                if key in self.vars:
                    self.vars[key].set(value)
            # Writing now means the dialog does not reappear, and the rest of
            # the application can assume a settings file exists.
            self.save_cfg()
            self._refresh_effective_paths()
            return True

        def accept():
            if apply({k: v.get() for k, v in chosen.items()}):
                win.destroy()

        def use_defaults():
            # The resolved paths, not blanks: the point of showing them is that
            # the user can see afterwards where their files went.
            if apply({key: default for _, _, key, default in rows}):
                win.destroy()

        tk.Button(act, text='Continue', width=14, font=self.f_bold, command=accept
                  ).pack(side='right')
        tk.Button(act, text='Use defaults', width=14,
                  command=use_defaults).pack(side='right', padx=6)

        win.update_idletasks()
        x = self.winfo_rootx() + (self.winfo_width() - win.winfo_width()) // 2
        y = self.winfo_rooty() + 80
        win.geometry(f'+{max(x, 0)}+{max(y, 0)}')
        win.grab_set()
        self.wait_window(win)
        self.scan_data()

    def _ensure_writable(self, path):
        """Create a folder the user named, and prove it can be written to.

        The pipeline creates its folders anyway - mkdir(parents=True) - but not
        until a run starts, so a path on a drive that is missing, full or
        read-only was accepted here in silence and failed an hour later inside
        a step that had nothing to do with it. Creating it now turns that into
        a sentence, while the user is still looking at the field they typed.
        """
        try:
            os.makedirs(path, exist_ok=True)
            probe = os.path.join(path, '.mesh_write_test')
            with open(probe, 'w', encoding='utf-8') as fh:
                fh.write('')
            os.remove(probe)
            return True
        except OSError as exc:
            messagebox.showerror(
                'That folder cannot be used',
                f'{path}\n\n{exc}\n\n'
                'The folder was created if it did not exist, but it could not '
                'be written to. Check the drive is present and that you have '
                'permission, or choose somewhere else.')
            return False

    def _effective_paths(self):
        """Where things will actually be written, with every default resolved.

        Read from MeshConfig rather than reconstructed, so this cannot drift
        from what the pipeline does.
        """
        try:
            from mesh_aop.config_parser import MeshConfig
            cfg = MeshConfig(config_path=self.cfg_path)
            return [
                ('Downloaded archives and databases', str(cfg.active_raw_dir)),
                ('Working files (networks, scores)', str(cfg.active_source_dir)),
                ('Results', str(cfg.results_dir)),
                ('Figures', str(cfg.figures_dir)),
                ('Logs', str(cfg.log_dir)),
                ('Settings file', self.cfg_path),
            ]
        except Exception as exc:
            return [('(could not be resolved)', str(exc))]

    def _refresh_effective_paths(self):
        """Update the read-only summary at the foot of the Folders tab."""
        box = getattr(self, 'eff_text', None)
        if box is None:
            return
        try:
            box.config(state='normal')
            box.delete('1.0', 'end')
            for label, path in self._effective_paths():
                box.insert('end', f'{label:<36}{path}\n')
            box.config(state='disabled')
        except Exception:
            pass

    def _browse_dir(self, var):
        d = filedialog.askdirectory(initialdir=var.get() or os.path.expanduser('~'))
        if d:
            var.set(os.path.normpath(d))

    def _set_icon(self):
        """Give the window and taskbar entry the application icon.

        iconbitmap takes the .ico on Windows; elsewhere Tk wants a PhotoImage,
        so the PNG is used instead. The reference is kept on the instance because
        Tk does not own the image and it would otherwise be collected, leaving a
        blank icon. Never fatal: a missing icon is not a reason to fail to start.
        """
        assets = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'assets')
        try:
            if os.name == 'nt':
                ico = os.path.join(assets, 'mesh_workbench.ico')
                if os.path.exists(ico):
                    self.iconbitmap(default=ico)
                    return
            png = os.path.join(assets, 'mesh_workbench.png')
            if os.path.exists(png):
                self._icon_img = tk.PhotoImage(file=png)
                self.iconphoto(True, self._icon_img)
        except tk.TclError:
            pass

    def show(self, name):
        for fr in self.screens.values():
            fr.pack_forget()
        self.screens[name].pack(fill='both', expand=True)

    def _sunken(self, parent, **kw):
        return tk.Frame(parent, bg=FACE, relief='sunken', bd=2, **kw)

    @staticmethod
    def _setup_cols(frame):
        """One column geometry for the data table, shared by header and rows.

        Column 3 is reserved on every row even though only the archive row puts a
        button there - without it that row's Status and Size shift left and the
        table looks ragged.
        """
        frame.columnconfigure(0, minsize=330, weight=1)
        frame.columnconfigure(1, minsize=80)
        frame.columnconfigure(2, minsize=80)
        frame.columnconfigure(3, minsize=130)

    # ------------------------------------------------------------- screen: setup
    # ------------------------------------------------------------ screen: setup
    #
    # Every path here comes from MeshConfig, never from a repo-relative guess.
    # The previous version hard-coded `data/raw/...` beneath the program folder
    # and the DAC_Mesh prefix; on an installed copy the pipeline keeps its data
    # under the user profile, so the screen reported "Missing" for files that
    # were plainly there - and the one thing this screen exists to answer is
    # whether a file is present.

    def _setup_config(self):
        """The live configuration, or None if it cannot be built."""
        try:
            from mesh_aop.config_parser import MeshConfig
            return MeshConfig(config_path=self.cfg_path)
        except Exception:
            return None

    def _build_setup(self, root):
        tk.Label(root, text='Local data sources', bg=FACE, font=self.f_bold,
                 anchor='w').pack(fill='x', padx=10, pady=(10, 2))
        tk.Label(root, bg=FACE, anchor='w', justify='left',
                 text='The pipeline runs offline once these are in place. '
                      'Nothing is downloaded or removed unless you ask.'
                 ).pack(fill='x', padx=10, pady=(0, 4))
        self.setup_where = tk.Label(root, bg=FACE, anchor='w', fg=DIM,
                                    font=self.f_mono)
        self.setup_where.pack(fill='x', padx=10, pady=(0, 6))

        box = self._sunken(root)
        box.pack(fill='both', expand=True, padx=10)
        hdr = tk.Frame(box, bg=FACE)
        hdr.pack(fill='x', padx=4, pady=(3, 0))
        self._setup_cols(hdr)
        for c, (text, anchor) in enumerate((('Source', 'w'), ('Status', 'w'),
                                            ('Size', 'e'), ('', 'e'))):
            tk.Label(hdr, text=text, bg=FACE, font=self.f_bold, anchor=anchor
                     ).grid(row=0, column=c, sticky='we',
                            padx=(0, 8) if c == 2 else 0)
        tk.Frame(box, bg='#808080', height=1).pack(fill='x', padx=4, pady=(2, 0))
        canvas = tk.Canvas(box, bg=FACE, highlightthickness=0)
        sb = tk.Scrollbar(box, command=canvas.yview)
        canvas.config(yscrollcommand=sb.set)
        sb.pack(side='right', fill='y')
        canvas.pack(side='left', fill='both', expand=True)
        self.setup_rows = tk.Frame(canvas, bg=FACE)
        canvas.create_window((0, 0), window=self.setup_rows, anchor='nw')
        self.setup_rows.bind('<Configure>',
                             lambda e: canvas.config(scrollregion=canvas.bbox('all')))

        foot = tk.Frame(root, bg=FACE)
        foot.pack(fill='x', padx=10, pady=9)

        # The daily-update choice, in the open. It used to be the third dialog
        # in a chain that only started once you had already agreed to a build,
        # so nobody could see whether it was on without committing to one.
        opts = tk.Frame(foot, bg=FACE)
        opts.pack(fill='x', pady=(0, 6))
        self.var_with_updates = tk.BooleanVar(value=False)
        tk.Checkbutton(opts, variable=self.var_with_updates, bg=FACE,
                       text='Also fetch the daily update files when building',
                       anchor='w').pack(side='left')
        tk.Label(opts, bg=FACE, fg=DIM, anchor='w',
                 text='  (the baseline is a yearly snapshot; the updates carry '
                      'everything published since)').pack(side='left')

        self.setup_total = tk.Label(foot, bg=FACE, relief='sunken', bd=2, anchor='w')
        self.setup_total.pack(fill='x', ipady=2, pady=(0, 6))
        act = tk.Frame(foot, bg=FACE)
        act.pack(fill='x')
        tk.Button(act, text='Continue', width=11, font=self.f_bold,
                  command=lambda: self.show('settings')).pack(side='right', padx=2)
        tk.Button(act, text='Refresh', width=11, command=self.scan_data
                  ).pack(side='right', padx=2)
        self.after(120, self.scan_data)

    def scan_data(self):
        """Measure what is on disk, at the paths the pipeline actually uses."""
        for w in self.setup_rows.winfo_children():
            w.destroy()

        cfg = self._setup_config()
        if cfg is None:
            self.setup_where.config(text='  (the configuration could not be read)')
            self.setup_total.config(text='   Nothing could be measured.')
            return

        raw, proc = cfg.active_raw_dir, cfg.active_source_dir
        prefix = cfg.prefix
        self.setup_where.config(
            text=f'  data  {raw}\n  work  {proc}')

        f = cfg.files
        items = [
            dict(key='master', label='Master annotation database',
                 path=str(f['master_db']),
                 note='Built from the archive. Required for every run. Hours to rebuild.',
                 build='Build', rebuild='Rebuild', delete=True, cost='hours'),
            dict(key='baseline', label='PubMed baseline archive',
                 path=str(raw / 'pubmed_baseline'),
                 note='The yearly snapshot, ~44 GB. Only needed to build the database.',
                 build='Download', rebuild='Re-download', delete=True, cost='hours'),
            dict(key='updates', label='PubMed daily updates',
                 path=str(raw / 'pubmed_updates'),
                 note='Records published since the baseline snapshot.',
                 build='Download', rebuild='Re-download', delete=True, cost='minutes'),
            dict(key='desc', label='MeSH descriptor file',
                 path=str(raw / 'desc2025.xml'),
                 note='Defines the stop-word vocabulary. Downloaded automatically.',
                 delete=True, cost='minutes'),
            dict(key='pmids', label=f'Retrieved PMIDs  ({prefix})',
                 path=str(f['pmids_db']),
                 note='Everything the search and its citation hops returned.',
                 delete=True, cost='hours'),
            dict(key='cleaned', label=f'Citation database  ({prefix})',
                 path=str(f['cleaned_db']),
                 note='Query results and citation generations, with MeSH attached.',
                 delete=True, cost='minutes'),
            dict(key='relevance', label=f'Relevance database  ({prefix})',
                 path=str(f['relevance_db']),
                 note='Per-article scores for this project.',
                 delete=True, cost='minutes'),
        ]

        total = reclaim = 0.0
        for it in items:
            gb = _size_gb(it['path'])
            present = gb is not None
            total += gb or 0.0
            if it['key'] in ('baseline', 'updates'):
                reclaim += gb or 0.0

            row = tk.Frame(self.setup_rows, bg=FACE)
            row.pack(fill='x', pady=2, padx=4)
            self._setup_cols(row)
            tk.Label(row, text=it['label'], bg=FACE, font=self.f_bold, anchor='w'
                     ).grid(row=0, column=0, sticky='w')
            tk.Label(row, text=it['note'], bg=FACE, fg=DIM, anchor='w'
                     ).grid(row=1, column=0, sticky='w')
            tk.Label(row, text=('Present' if present else 'Missing'), bg=FACE,
                     fg=(OK if present else ERR), font=self.f_bold, anchor='w'
                     ).grid(row=0, column=1, sticky='w')
            tk.Label(row, text=(f'{gb:,.2f} GB' if present else '-'), bg=FACE,
                     anchor='e').grid(row=0, column=2, sticky='e', padx=(0, 8))

            buttons = tk.Frame(row, bg=FACE)
            buttons.grid(row=0, column=3, rowspan=2, sticky='e', padx=4)
            label = it.get('rebuild' if present else 'build')
            if label:
                tk.Button(buttons, text=label, width=12,
                          command=lambda i=it, p=present: self._setup_build(i, p)
                          ).pack(side='left', padx=2)
            if it.get('delete') and present:
                tk.Button(buttons, text='Delete', width=8,
                          command=lambda i=it, g=gb: self._setup_delete(i, g)
                          ).pack(side='left', padx=2)
            tk.Frame(row, bg='#b0b0b0', height=1).grid(
                row=2, column=0, columnspan=4, sticky='we', pady=(3, 0))

        self.setup_total.config(
            text=f'   On disk now {total:,.2f} GB' +
                 (f'   -   {reclaim:,.2f} GB of that is downloaded archive, '
                  f'removable once the database is built'
                  if reclaim else ''))

    # -- typed confirmation ------------------------------------------------
    #
    # Every action on this screen either destroys something that took hours to
    # produce or starts something that will take hours to finish. A yes/no box
    # is one mis-aimed click; typing the word is a deliberate act, and it forces
    # a reading of what the word is attached to.

    def _make_modal(self, win, title, on_close=None):
        """Centre a dialog on the window and decide what closing it means.

        `on_close=None` refuses the close box outright: the dialog has to be
        answered with one of the buttons on it. That is deliberate for anything
        destructive or expensive - the whole point of the dialog is that the
        reader takes in what is at stake, and a title-bar X is the one way out
        that requires reading nothing. Cancel is always present, so nobody is
        trapped; they are only asked to say which they meant.

        Passing a callable instead routes the close box (and Escape) to it, for
        dialogs where dismissing genuinely is one of the choices.
        """
        win.title(title)
        win.configure(bg=FACE)
        win.transient(self)
        win.resizable(False, False)
        if on_close is None:
            def refuse():
                win.bell()
                self._flash(win)
            win.protocol('WM_DELETE_WINDOW', refuse)
            win.bind('<Escape>', lambda _e: refuse())
        else:
            win.protocol('WM_DELETE_WINDOW', on_close)
            win.bind('<Escape>', lambda _e: on_close())
        return win

    def _flash(self, win):
        """Draw the eye back to the dialog when its close box is refused."""
        try:
            original = win.cget('bg')
            win.configure(bg=WARN)
            win.after(90, lambda: win.configure(bg=original))
        except Exception:
            pass

    def _place_modal(self, win, offset=60):
        win.update_idletasks()
        win.geometry('+%d+%d' % (self.winfo_rootx() + offset,
                                 self.winfo_rooty() + offset))
        win.grab_set()

    def confirm_typed(self, title, word, message, detail=''):
        """Return True only if the user types `word` exactly.

        There is no way out of this dialog except its own two buttons - see
        _make_modal. The action it guards either destroys something that took
        hours to produce or starts something that will take hours to finish.
        """
        win = tk.Toplevel(self)
        self._make_modal(win, title)
        body = tk.Frame(win, bg=FACE)
        body.pack(fill='both', expand=True, padx=16, pady=14)
        tk.Label(body, text=title, bg=FACE, font=self.f_bold, anchor='w'
                 ).pack(fill='x')
        tk.Label(body, text=message, bg=FACE, anchor='w', justify='left',
                 wraplength=520).pack(fill='x', pady=(6, 6))
        if detail:
            tk.Label(body, text=detail, bg=FACE, anchor='w', justify='left',
                     fg=ERR, wraplength=520, font=self.f_bold).pack(fill='x',
                                                                    pady=(0, 6))
        tk.Label(body, bg=FACE, anchor='w', fg=DIM, justify='left', wraplength=520,
                 text='Answer with one of the buttons below - the close box is '
                      'disabled here on purpose.').pack(fill='x', pady=(0, 8))
        tk.Label(body, text=f'Type {word} to confirm:', bg=FACE, anchor='w'
                 ).pack(fill='x', pady=(4, 2))
        typed = tk.StringVar()
        entry = tk.Entry(body, bg=FIELD, relief='sunken', bd=2, textvariable=typed)
        entry.pack(fill='x')
        entry.focus_set()

        state = {'ok': False}
        row = tk.Frame(body, bg=FACE)
        row.pack(fill='x', pady=(12, 0))
        go = tk.Button(row, text=word.title(), width=12, state='disabled',
                       command=lambda: (state.update(ok=True), win.destroy()))
        go.pack(side='left')
        tk.Button(row, text='Cancel', width=10,
                  command=win.destroy).pack(side='right')

        # Traced on the variable rather than bound to <KeyRelease>: a word put
        # in by a right-click paste produces no key event, and the button would
        # have stayed disabled with the correct word sitting in the box.
        def check(*_):
            go.config(state='normal' if typed.get().strip() == word else 'disabled')
        typed.trace_add('write', check)
        entry.bind('<Return>',
                   lambda e: go.invoke() if str(go['state']) == 'normal' else None)

        self._place_modal(win)
        self.wait_window(win)
        return state['ok']

    def _setup_delete(self, item, gb):
        cost = {'hours': 'hours', 'minutes': 'minutes'}.get(item['cost'], 'time')
        detail = ('Rebuilding this takes HOURS.' if item['cost'] == 'hours' else '')
        if not self.confirm_typed(
                f"Delete: {item['label']}", 'DELETE',
                f"{item['path']}\n\n"
                f"This frees {gb:,.2f} GB. Nothing else is touched, and the "
                f"pipeline rebuilds it when it next needs it - which costs "
                f"{cost}.", detail):
            return
        try:
            if os.path.isdir(item['path']):
                shutil.rmtree(item['path'])
            else:
                os.remove(item['path'])
                # A database's write-ahead log and health record describe a file
                # that no longer exists; left behind, SQLite will try to replay
                # the log into whatever is built next.
                for suffix in ('-wal', '-shm', '-journal', '.health.json'):
                    side = item['path'] + suffix
                    if os.path.exists(side):
                        os.remove(side)
        except OSError as exc:
            messagebox.showerror('Could not delete', f"{item['path']}\n\n{exc}")
            return
        self.scan_data()

    def _setup_build(self, item, present):
        """Build or download one item. Step 0 covers the first three."""
        if item['key'] == 'desc':
            messagebox.showinfo(
                'Downloaded automatically',
                'The MeSH descriptor file is fetched from the NLM whenever the '
                'support files are rebuilt. Run Step 1 to refresh it.')
            return
        if item['key'] in ('pmids', 'cleaned', 'relevance'):
            messagebox.showinfo(
                'Produced by a run',
                f"{item['label']} is written by the pipeline itself.\n\n"
                'Delete it here if you want it rebuilt, then run the step that '
                'produces it.')
            return

        cfg = self._setup_config()
        extra = ['--build-database']
        wants_updates = bool(self.var_with_updates.get())

        if item['key'] == 'updates':
            # Updates alone: compile from what is on disk, plus the update files.
            extra += ['--skip-baseline-download', '--with-updates']
            if not self.confirm_typed(
                    'Fetch the daily updates', 'REBUILD',
                    'The daily update files will be downloaded and compiled into '
                    'the master database.\n\nThe baseline archive already on '
                    'disk is not downloaded again.'):
                return
            self.start_run('baseline', extra, title='fetching the daily updates')
            return

        archive = str(cfg.active_raw_dir / 'pubmed_baseline') if cfg else ''
        have_archive = os.path.isdir(archive) and any(os.scandir(archive))

        if item['key'] == 'baseline':
            if present and not self.confirm_typed(
                    'Download the baseline again', 'REBUILD',
                    'The PubMed baseline archive is already here. Downloading it '
                    'again fetches roughly 44 GB and takes hours.',
                    'You almost certainly do not need this. To rebuild the '
                    'database from the archive you already have, use Rebuild on '
                    'the master database row instead.'):
                return
            if wants_updates:
                extra.append('--with-updates')
            self.start_run('baseline', extra, title='downloading the PubMed baseline')
            return

        # The master database.
        if present:
            if not self.confirm_typed(
                    'Rebuild the master annotation database', 'REBUILD',
                    'The existing database is deleted and built again from the '
                    'archive.\n\nNothing else in the pipeline can run while it '
                    'is missing.',
                    'This takes HOURS.'):
                return
            extra.append('--rebuild-corrupt')
        elif not have_archive and not self.confirm_typed(
                'Download and build', 'REBUILD',
                'No archive is present, so it is downloaded first: roughly 44 GB, '
                'and several hours on a fast connection.\n\nThe download resumes '
                'if interrupted, and the archive can be deleted afterwards.'):
            return

        if have_archive:
            extra.append('--skip-baseline-download')
        if wants_updates:
            extra.append('--with-updates')
        self.start_run('baseline', extra, title='building the master database')

    # ---------------------------------------------------------- screen: settings
    def _build_settings(self, root):
        top = self._sunken(root)
        top.pack(fill='x', padx=10, pady=(10, 8))
        tk.Label(top, text='Pipeline step:', bg=FACE).pack(side='left', padx=(8, 6), pady=6)
        self.step_var = tk.StringVar(value='all')
        om = tk.OptionMenu(top, self.step_var, *[s for s, _, _ in schema.STEPS],
                           command=lambda *_: self._step_changed())
        om.config(bg=FACE, activebackground=FACE, highlightthickness=1, width=12)
        om.pack(side='left', pady=4)
        self.step_desc = tk.Label(top, bg=FACE, fg=DIM, anchor='w')
        self.step_desc.pack(side='left', padx=8)
        self.step_eta = tk.Label(top, bg=FACE, fg=OK, font=self.f_bold)
        self.step_eta.pack(side='right', padx=10)

        strip = tk.Frame(root, bg=FACE)
        strip.pack(fill='x', padx=10)
        self.sheet = self._sunken(root)
        self.sheet.pack(fill='both', expand=True, padx=10)
        self.panes = {}
        self.tab_btns = {}
        for name, fields in schema.TABS:
            b = tk.Button(strip, text=name, relief='raised', bd=2, bg=FACE,
                          padx=6, pady=2, command=lambda n=name: self._tab(n))
            b.pack(side='left', padx=1)
            self.tab_btns[name] = b
            pane = tk.Frame(self.sheet, bg=FACE)
            self.panes[name] = pane
            self._build_fields(pane, fields, tab_name=name)

        self.help = tk.Frame(root, bg=FACE, relief='sunken', bd=2, height=HELP_H)
        self.help.pack(fill='x', padx=10, pady=(8, 0))
        self.help.pack_propagate(False)
        self.help_title = tk.Label(self.help, bg=FACE, font=self.f_bold, anchor='w')
        self.help_title.pack(fill='x', padx=8, pady=(5, 0))
        self.help_what = tk.Label(self.help, bg=FACE, anchor='w', justify='left',
                                  wraplength=W - 60)
        self.help_what.pack(fill='x', padx=8)
        self.help_def = tk.Label(self.help, bg=FACE, fg=DIM, anchor='w',
                                 justify='left', wraplength=W - 60)
        self.help_def.pack(fill='x', padx=8)
        self.help_note = tk.Label(self.help, bg=FACE, fg=ERR, anchor='w',
                                  justify='left', wraplength=W - 60)
        self.help_note.pack(fill='x', padx=8)
        self._help(None, 'Settings',
                   'Click any field to see what it controls, its default, and '
                   'anything worth knowing before changing it.', '', '')

        act = tk.Frame(root, bg=FACE)
        act.pack(fill='x', padx=10, pady=8)
        tk.Button(act, text='Run', width=12, font=self.f_bold,
                  command=lambda: self.start_run(self.step_var.get())).pack(side='left')
        tk.Button(act, text='Save', width=12, command=self.save_cfg).pack(side='left', padx=4)
        tk.Button(act, text='Defaults', width=12,
                  command=self.restore_defaults).pack(side='left', padx=4)
        tk.Button(act, text='Data setup…', width=12,
                  command=lambda: self.show('setup')).pack(side='left', padx=4)
        self.status = tk.Label(act, bg=FACE, relief='sunken', bd=2, anchor='w',
                               text='Ready.')
        self.status.pack(side='left', fill='x', expand=True, padx=(8, 0), ipady=2)
        self._tab(schema.TABS[0][0])
        self._step_changed()
        self._refresh_effective_paths()

    def _build_fields(self, pane, fields, tab_name=None):
        grid = tk.Frame(pane, bg=FACE)
        grid.pack(fill='x', padx=10, pady=10)
        grid.columnconfigure(1, weight=1)
        for r, f in enumerate(fields):
            cur = schema.get(self.cfg, f.key, f.default)
            # An empty folder field means "work it out from the defaults", which
            # is correct behaviour and useless to look at: the user cannot see
            # where their files are going. Show the resolved path instead, for
            # the two locations the first-run dialog asks about.
            #
            # Not for a portable copy: there the defaults are relative to the
            # program folder, and writing an absolute path in would break the
            # copy the moment it was moved to another machine or a USB stick.
            if (f.key in self._PREFILL_DIRS and not str(cur or '').strip()
                    and not self.paths.is_portable(self.repo_dir)):
                cur = str(self._PREFILL_DIRS[f.key](self.repo_dir))
            if f.kind == 'bool':
                var = tk.BooleanVar(value=bool(cur))
                w = tk.Checkbutton(grid, text=f.label, variable=var, bg=FACE,
                                   activebackground=FACE, anchor='w')
                w.grid(row=r, column=0, columnspan=2, sticky='we', pady=3)
            else:
                var = tk.StringVar(value='' if cur is None else str(cur))
                tk.Label(grid, text=f.label, bg=FACE, anchor='w'
                         ).grid(row=r, column=0, sticky='w', padx=(0, 10), pady=3)
                w = tk.Entry(grid, textvariable=var, bg=FIELD, relief='sunken', bd=2)
                w.grid(row=r, column=1, sticky='we', pady=3)
                if f.kind == 'path':
                    tk.Button(grid, text='…', width=3,
                              command=lambda v=var: self._browse(v)
                              ).grid(row=r, column=2, padx=(4, 0))
            self.vars[f.key] = var
            for ev in ('<FocusIn>', '<Button-1>', '<Enter>'):
                w.bind(ev, lambda _e, fl=f: self._help(fl))
            if f.kind in ('int', 'float') and not f.kind == 'bool':
                w.bind('<FocusIn>', lambda e: e.widget.config(bg=FOCUS), add='+')
                w.bind('<FocusOut>', lambda e: e.widget.config(bg=FIELD), add='+')

        # A standing note for the tab, if it has one. The description pane below
        # only ever shows the field you are on, so anything you need to know
        # before choosing a value has nowhere else to be said.
        if tab_name == 'Folders':
            # Where things actually land, with every default resolved. The
            # fields above are inputs and several of them mean "decide for me";
            # this is the answer, and it is what the user came to find out.
            box = tk.Frame(pane, bg=FACE, relief='ridge', bd=2)
            box.pack(fill='both', expand=True, padx=10, pady=(0, 10))
            tk.Label(box, text='Where your files will actually go', bg=FACE,
                     font=self.f_bold, anchor='w').pack(fill='x', padx=8, pady=(6, 0))
            tk.Label(box, bg=FACE, fg=DIM, anchor='w', justify='left', wraplength=760,
                     text='Resolved from the settings above. A blank field means '
                          'the default, which is shown here in full. Press Save '
                          'to update this after a change.'
                     ).pack(fill='x', padx=8, pady=(0, 4))
            self.eff_text = tk.Text(box, bg=FACE, relief='flat', height=7,
                                    font=self.f_mono, wrap='none', cursor='arrow',
                                    highlightthickness=0, padx=8)
            self.eff_text.pack(fill='both', expand=True, pady=(0, 6))
            self.eff_text.config(state='disabled')

        note = schema.TAB_NOTES.get(tab_name)
        if note:
            title, text = note
            box = tk.Frame(pane, bg=FACE, relief='ridge', bd=2)
            box.pack(fill='both', expand=True, padx=10, pady=(0, 10))
            tk.Label(box, text=title, bg=FACE, font=self.f_bold, anchor='w'
                     ).pack(fill='x', padx=8, pady=(6, 2))
            body = tk.Text(box, bg=FACE, relief='flat', wrap='word',
                           font=self.f_ui, height=11, cursor='arrow',
                           highlightthickness=0, padx=8, pady=0)
            body.insert('1.0', text)
            # Read-only, but still selectable: these are instructions someone
            # may well want to copy.
            body.config(state='disabled')
            body.pack(fill='both', expand=True, padx=0, pady=(0, 6))

    def _browse(self, var):
        p = filedialog.askopenfilename(
            initialdir=self.repo_dir,
            filetypes=[('Ground truth', '*.xlsx *.xls *.csv'), ('All files', '*.*')])
        if p:
            try:
                p = os.path.relpath(p, self.repo_dir).replace('\\', '/')
            except ValueError:
                pass
            var.set(p)

    def _tab(self, name):
        for n, p in self.panes.items():
            p.pack_forget()
            self.tab_btns[n].config(font=self.f_ui, relief='raised')
        self.panes[name].pack(fill='both', expand=True)
        self.tab_btns[name].config(font=self.f_bold, relief='sunken')

    def _help(self, fld, title=None, what=None, deflt=None, note=None):
        self.help_title.config(text=title or (fld.label if fld else ''))
        self.help_what.config(text=what if what is not None else (fld.what if fld else ''))
        self.help_def.config(text=deflt if deflt is not None else (fld.deflt if fld else ''))
        n = note if note is not None else ((fld.note or '') if fld else '')
        self.help_note.config(text=('Note: ' + n) if n else '')

    def _step_changed(self):
        s = self.step_var.get()
        for step, label, eta in schema.STEPS:
            if step == s:
                self.step_desc.config(text=label.split(' - ', 1)[-1])
                self.step_eta.config(text='est. ' + eta)

    def restore_defaults(self):
        if not messagebox.askyesno(
                'Restore default settings?',
                'Every value on every tab returns to its shipped default, '
                'including your NCBI credentials and folder choices.\n\n'
                'Settings already written to mesh_config.json are not changed '
                'until you save again.\n\nRestore defaults?', icon='warning'):
            return
        for _, fields in schema.TABS:
            for f in fields:
                self.vars[f.key].set(f.default if f.kind != 'bool' else bool(f.default))
        self.status.config(text='Defaults restored - not yet saved.')

    # ----------------------------------------------------------- screen: running
    def _build_running(self, root):
        pad = dict(padx=10)
        self.run_title = tk.Label(root, bg=FACE, font=self.f_bold, anchor='w')
        self.run_title.pack(fill='x', pady=(10, 6), **pad)

        self.phase_lab = tk.Label(root, bg=FACE, anchor='w')
        self.phase_lab.pack(fill='x', **pad)
        self.bar_overall = _Bar(root)
        self.bar_overall.pack(fill='x', pady=(2, 8), **pad)

        self.sub_lab = tk.Label(root, bg=FACE, anchor='w', fg=DIM)
        self.sub_lab.pack(fill='x', **pad)
        self.bar_sub = _Bar(root)
        self.bar_sub.pack(fill='x', pady=(2, 8), **pad)

        strip = tk.Frame(root, bg=FACE)
        strip.pack(fill='x', pady=(0, 6), **pad)
        self.lab_spin = tk.Label(strip, bg=FACE, relief='sunken', bd=2, width=14, anchor='w')
        self.lab_spin.pack(side='left', ipady=2)
        self.lab_elapsed = tk.Label(strip, bg=FACE, relief='sunken', bd=2, width=18, anchor='w')
        self.lab_elapsed.pack(side='left', padx=4, ipady=2)
        self.lab_transient = tk.Label(strip, bg=FACE, relief='sunken', bd=2, anchor='w')
        self.lab_transient.pack(side='left', fill='x', expand=True, ipady=2)

        self.console = tk.Text(root, bg=CONSOLE_BG, fg='#c0c0c0', font=self.f_mono,
                               relief='sunken', bd=2, wrap='word', height=14,
                               insertbackground='#c0c0c0')
        self.console.pack(fill='both', expand=True, **pad)
        for tag, col in (('ok', '#2ecc40'), ('warn', '#ffdc00'),
                         ('err', '#ff4136'), ('dim', '#7a7a7a')):
            self.console.tag_config(tag, foreground=col)
        self.console.config(state='disabled')

        act = tk.Frame(root, bg=FACE)
        act.pack(fill='x', pady=9, **pad)
        # Pause and Abort are deliberately separate buttons, not one control
        # with a confirmation: they do entirely different things to three hours
        # of work, and a user reaching for "stop" under time pressure should not
        # have to read a dialog to discover which one they got.
        self.btn_pause = tk.Button(act, text='Pause', width=10, command=self.pause_run)
        self.btn_pause.pack(side='left')
        self.btn_cancel = tk.Button(act, text='Abort run', width=11, command=self.cancel_run)
        self.btn_cancel.pack(side='left', padx=4)
        self.btn_back = tk.Button(act, text='Back to settings', width=16,
                                  command=lambda: self.show('settings'))
        self.btn_back.pack(side='left', padx=(0, 4))
        self.btn_results = tk.Button(act, text='View results', width=14,
                                     state='disabled', command=self.open_results)
        self.btn_results.pack(side='left', padx=(0, 4))
        self.run_status = tk.Label(act, bg=FACE, relief='sunken', bd=2, anchor='w')
        self.run_status.pack(side='left', fill='x', expand=True, padx=(8, 0), ipady=2)

    def start_run(self, step, extra=None, title=None):
        if self.runner.is_running():
            messagebox.showinfo('Already running',
                                'A pipeline step is already in progress.')
            return
        if not self.save_cfg():
            return
        step = step or 'all'
        title = title or step
        if not self._confirm_overwrite(step):
            return
        self.show('running')
        self.run_title.config(text=f'Running: {title}')
        self._log_clear()
        self._log(f'--- starting {title} ---', 'dim')
        self.bar_overall.set(0)
        self.bar_sub.set(0)
        self.phase_lab.config(text='Starting…')
        self.sub_lab.config(text='')
        self._elapsed = 0
        self._spin_i = 0
        self.btn_cancel.config(state='normal')
        self.btn_pause.config(state='normal', text='Pause')
        self._paused_annotation = None
        self.runner.start(step, extra)
        self._pump()

    def _confirm_overwrite(self, step):
        """Warn before replacing results this prefix already has.

        Silent when there is nothing to replace, which is the whole point: a
        prefix that has not been used for this step raises no question, so a
        fresh project never meets a dialog it has no reason to read.
        """
        # Step 0 is driven from the Data Setup screen, which asks for a typed
        # confirmation of its own; asking twice would train people to click past.
        if step == 'baseline':
            return True
        try:
            from mesh_aop import integrity
            from mesh_aop.config_parser import MeshConfig
            config = MeshConfig(config_path=self.cfg_path)
            at_risk = integrity.would_overwrite(config, step)
        except Exception:
            return True                 # never block a run over a warning
        if not at_risk:
            return True

        total = sum(size for _, _, size in at_risk)
        lines = []
        for label, path, size in at_risk[:12]:
            lines.append(f'    {label}   ({size/1e6:,.1f} MB)')
        if len(at_risk) > 12:
            lines.append(f'    ... and {len(at_risk) - 12} more')

        return messagebox.askyesno(
            'Existing results will be replaced',
            f'The prefix "{config.prefix}" has already been used for this part '
            f'of the pipeline. Running it again replaces:\n\n'
            + '\n'.join(lines) +
            f'\n\n{len(at_risk)} item(s), {total/1e6:,.1f} MB in total.\n\n'
            'To keep them, cancel and change the project prefix on the Search '
            'tab - two prefixes never touch each other\'s files.\n\n'
            'Steps that are already complete are skipped rather than redone, so '
            'this is also how you resume an interrupted run.\n\n'
            'Continue?', icon='warning', default='no')

    def pause_run(self):
        """Ask the run to stop at its next safe point, or let a stopped one go on."""
        if not self.runner.is_running():
            return
        if self.runner.is_paused() or self.runner.pause_pending():
            self.runner.resume()
            self.btn_pause.config(text='Pause')
            self.run_status.config(text='Running.', fg=INK)
            return
        if not self.runner.can_pause():
            messagebox.showinfo(
                'Pause unavailable',
                'This run was not started with pause control, so it can only be '
                'left to finish or aborted.')
            return
        if self.runner.pause():
            self.btn_pause.config(text='Resume')
            # Say what is actually true: the request is in, and the run stops
            # when it reaches a point where stopping is safe. Claiming it is
            # already paused while the progress bar keeps moving would be a
            # small lie that costs the user their trust in every other message.
            self.run_status.config(
                text='Pausing - it will stop at the next safe point (this can '
                     'take a few minutes on a long step).', fg=WARN)

    def cancel_run(self):
        if not self.runner.is_running():
            self.show('settings')
            return
        # Say what is actually lost. "Partial outputs may be left on disk" told
        # the user nothing they could act on; the elapsed time and the offer to
        # pause instead are the two facts that change the decision.
        mins = self.runner.elapsed() / 60.0
        spent = (f'{mins:.0f} minutes' if mins >= 1 else 'less than a minute')
        paused = self.runner.is_paused()
        if messagebox.askyesno(
                'Abort this run?',
                f'The step is stopped for good and its work is lost - '
                f'{spent} so far.\n\n'
                'Finished steps are kept, so a later run resumes from the last '
                'one that completed rather than starting over. Any half-written '
                'file is removed by Tools > Check and repair files.\n\n'
                + ('The run is paused at a safe point; it will be released and '
                   'then stopped.\n\n' if paused
                   else 'To stop it only temporarily, use Pause instead.\n\n')
                + 'Abort it?',
                icon='warning', default='no'):
            self.runner.cancel()
            self.btn_pause.config(state='disabled')
            self.run_status.config(text='Aborting…', fg=WARN)

    def _annotation_pause(self, anno_path):
        """The run stopped on purpose and is waiting for the user.

        A console message is not enough here: the pause can last days, this is
        the one point where the pipeline needs a person, and the file it needs
        editing is somewhere the user has never had reason to look.
        """
        self._paused_annotation = anno_path
        self._sync_tools_state()
        self.bar_overall.set(100)
        self.lab_spin.config(text='  Waiting for you')
        self.run_status.config(text='Paused for AOP annotation.', fg=WARN)

        from mesh_aop.guides import AOP_LEVELS
        levels = ', '.join(AOP_LEVELS)
        win = tk.Toplevel(self)
        # Here dismissing IS one of the choices - "Later" is on the dialog - so
        # the close box is allowed, and routed to exactly what that button does
        # rather than to some undefined fourth outcome.
        self._make_modal(win, 'The run is paused - it needs your annotation',
                         on_close=lambda: win.destroy())

        body = tk.Frame(win, bg=FACE)
        body.pack(fill='both', expand=True, padx=16, pady=14)
        tk.Label(body, text='The network is built. It now needs you.',
                 bg=FACE, font=self.f_bold, anchor='w').pack(fill='x')
        tk.Label(
            body, bg=FACE, justify='left', anchor='w', wraplength=560,
            text=('Nothing has gone wrong. Every MeSH term in the network has to '
                  'be placed on the adverse outcome pathway, and that is a '
                  'judgement the pipeline cannot make for you.\n\n'
                  'Open the file below and fill in the "aop_level" column. It is '
                  'semicolon-delimited - keep it that way, because MeSH headings '
                  'contain commas.')
        ).pack(fill='x', pady=(6, 8))

        tk.Label(body, text='Allowed levels, in pathway order:', bg=FACE,
                 anchor='w', font=self.f_bold).pack(fill='x')
        tk.Label(body, text='    ' + levels, bg=FACE, anchor='w',
                 wraplength=560, justify='left').pack(fill='x', pady=(0, 8))
        tk.Label(body, bg=FACE, anchor='w', justify='left', wraplength=560, fg=DIM,
                 text=('Leave a term as "Uncategorized" if it genuinely is not on '
                       'the pathway. It stays in the network and in the '
                       'topological figures, and is left out of the AOP ones.')
                 ).pack(fill='x', pady=(0, 8))

        tk.Label(body, text='The file to edit:', bg=FACE, anchor='w',
                 font=self.f_bold).pack(fill='x')
        box = tk.Entry(body, font=self.f_mono, relief='sunken', bd=2)
        box.insert(0, anno_path)
        box.config(state='readonly')
        box.pack(fill='x', pady=(2, 10))

        tk.Label(body, bg=FACE, anchor='w', justify='left', wraplength=560,
                 text=('When you have saved it, come back and run Step 4 - '
                       'Figures. Nothing earlier is recomputed, so it is quick, '
                       'and you can edit and re-run as often as you like.')
                 ).pack(fill='x', pady=(0, 12))

        row = tk.Frame(body, bg=FACE)
        row.pack(fill='x')
        tk.Button(row, text='Open the file', width=15,
                  command=lambda: self._open_path(anno_path)).pack(side='left')
        tk.Button(row, text='Show me the folder', width=18,
                  command=lambda: self._reveal(anno_path)).pack(side='left', padx=4)
        tk.Button(row, text='Run Step 4 now', width=16,
                  command=lambda: (win.destroy(), self.start_run('viz', title='Step 4 - Figures'))
                  ).pack(side='left', padx=(0, 4))
        tk.Button(row, text='Later', width=9, command=win.destroy).pack(side='right')

        self._place_modal(win, offset=40)

    def _open_path(self, path):
        """Open a file with whatever the system uses for it."""
        try:
            if sys.platform == 'win32':
                os.startfile(path)                      # noqa: S606
            elif sys.platform == 'darwin':
                subprocess.call(['open', path])
            else:
                subprocess.call(['xdg-open', path])
        except Exception as exc:
            messagebox.showerror('Could not open it', f'{path}\n\n{exc}')

    def _reveal(self, path):
        """Open the containing folder with the file selected, where that is possible."""
        folder = os.path.dirname(os.path.abspath(path))
        try:
            if sys.platform == 'win32':
                subprocess.Popen(['explorer', '/select,', os.path.abspath(path)])
            elif sys.platform == 'darwin':
                subprocess.call(['open', '-R', path])
            else:
                subprocess.call(['xdg-open', folder])
        except Exception:
            self._open_path(folder)

    def _pump(self):
        for ev in self.runner.drain():
            kind = ev[0]
            if kind == LOG:
                self._log(ev[1], ev[2])
            elif kind == TRANSIENT:
                self.lab_transient.config(text=' ' + ev[1][:110])
            elif kind == PHASE:
                self.phase_lab.config(text=ev[1])
                self.bar_sub.set(0)
            elif kind == PROGRESS:
                done, total = ev[1], ev[2]
                if total:
                    self.bar_sub.set(100.0 * done / total)
                    self.sub_lab.config(text=f'{done:,} / {total:,}')
            elif kind == PAUSED:
                self._annotation_pause(ev[1])
            elif kind == AT_CHECKPOINT:
                self.runner._mark_paused(ev[1])
                if ev[1]:
                    self.run_status.config(
                        text='Paused at a safe point - nothing is being computed.',
                        fg=WARN)
                    self.lab_spin.config(text='  Paused')
                else:
                    self.run_status.config(text='Running.', fg=INK)
                    self.btn_pause.config(text='Pause')
            elif kind == DONE:
                self._finish(ev[1], ev[2])
                return
        if self.runner.is_running():
            # From the runner, which excludes time spent paused - an elapsed
            # clock that counts a lunch break makes every duration a lie.
            self._elapsed = self.runner.elapsed()
            if self.runner.is_paused():
                self.lab_spin.config(text='  Paused')
            elif int(self._elapsed * 5) % 3 == 0:
                self._spin_i = (self._spin_i + 1) % 4
                self.lab_spin.config(text=' ' + '|/-\\'[self._spin_i] + '  Working…')
            self.lab_elapsed.config(text=' Elapsed ' + _hms(int(self._elapsed)))
            self._tick_job = self.after(200, self._pump)
        else:
            self._tick_job = self.after(200, self._pump)

    def _finish(self, rc, elapsed):
        self.lab_spin.config(text='  Finished')
        self.bar_overall.set(100 if rc == 0 else self.bar_overall.value)
        self.btn_cancel.config(state='disabled')
        self.btn_pause.config(state='disabled', text='Pause')
        self._sync_tools_state()
        if self._paused_annotation:
            # The pipeline exits 0 when it pauses on purpose. Reporting that as
            # "completed successfully" would be true of the process and quite
            # wrong about the analysis, which is unfinished and waiting.
            self._log('--- paused for AOP annotation ---', 'warn')
            self.run_status.config(text='Paused for annotation - see the dialog.', fg=WARN)
            return
        if rc == 0:
            self._log(f'--- completed in {elapsed/60:.1f} min ---', 'ok')
            self.run_status.config(text='Completed successfully.', fg=OK)
            self.refresh_results()
            self.btn_results.config(state='normal')
        elif rc == -1:
            self._log('--- aborted ---', 'warn')
            self.run_status.config(text='Aborted.', fg=WARN)
        elif rc == 130:
            # The pipeline saw the stop request at a checkpoint and exited
            # tidily. Worth distinguishing from a crash: nothing is half-written.
            self._log('--- stopped at a safe point ---', 'warn')
            self.run_status.config(
                text='Stopped at a safe point - no file was left half-written.',
                fg=WARN)
        else:
            self._log(f'--- failed, exit code {rc} ---', 'err')
            self.run_status.config(text=f'Failed (exit {rc}).', fg=ERR)

    def _log_clear(self):
        self.console.config(state='normal')
        self.console.delete('1.0', 'end')
        self.console.config(state='disabled')

    def _log(self, text, tag=''):
        self.console.config(state='normal')
        self.console.insert('end', text + '\n', tag or ())
        self.console.see('end')
        self.console.config(state='disabled')

    # ----------------------------------------------------------- screen: results
    def _build_results(self, root):
        self.res_title = tk.Label(root, bg=FACE, font=self.f_bold, anchor='w')
        self.res_title.pack(fill='x', padx=10, pady=(10, 4))

        # The overview first. After a run this is what someone actually wants:
        # a single account of what went in, what was excluded at each stage and
        # why, and what came out. The file listing answers a different and much
        # less common question, so it goes underneath.
        tk.Label(root, text='Run overview', bg=FACE, font=self.f_bold,
                 anchor='w').pack(fill='x', padx=10, pady=(4, 2))
        obox = self._sunken(root)
        obox.pack(fill='both', expand=True, padx=10)
        osb = tk.Scrollbar(obox)
        osb.pack(side='right', fill='y')
        self.res_overview = tk.Text(obox, bg=FIELD, font=self.f_mono,
                                    relief='flat', wrap='none', height=18,
                                    yscrollcommand=osb.set)
        osb.config(command=self.res_overview.yview)
        oxsb = tk.Scrollbar(obox, orient='horizontal',
                            command=self.res_overview.xview)
        oxsb.pack(side='bottom', fill='x')
        self.res_overview.config(xscrollcommand=oxsb.set)
        self.res_overview.pack(fill='both', expand=True)

        tk.Label(root, text='Files', bg=FACE, font=self.f_bold,
                 anchor='w').pack(fill='x', padx=10, pady=(8, 2))
        box = self._sunken(root)
        box.pack(fill='both', expand=True, padx=10)
        sb = tk.Scrollbar(box)
        sb.pack(side='right', fill='y')
        self.res_text = tk.Text(box, bg=FIELD, font=self.f_mono, relief='flat',
                                wrap='none', height=8, yscrollcommand=sb.set)
        sb.config(command=self.res_text.yview)
        self.res_text.pack(fill='both', expand=True)

        act = tk.Frame(root, bg=FACE)
        act.pack(fill='x', padx=10, pady=9)
        tk.Button(act, text='Open results folder', width=18, font=self.f_bold,
                  command=self._open_results).pack(side='left')
        self.btn_res_report = tk.Button(
            act, text='Open the full report', width=18, state='disabled',
            command=lambda: self._open_path(self._prisma_report_path()))
        self.btn_res_report.pack(side='left', padx=4)
        self.btn_res_figure = tk.Button(
            act, text='Open the diagram', width=15, state='disabled',
            command=self._open_prisma_figure)
        self.btn_res_figure.pack(side='left', padx=(0, 4))
        tk.Button(act, text='Refresh', width=10,
                  command=self.refresh_results).pack(side='left', padx=4)
        tk.Button(act, text='Back to settings', width=16,
                  command=lambda: self.show('settings')).pack(side='right')

    def _results_dir(self):
        """The configured results folder, not a repo-relative guess."""
        try:
            from mesh_aop.config_parser import MeshConfig
            return str(MeshConfig(config_path=self.cfg_path).results_dir)
        except Exception:
            var = self.vars.get('directories.results_dir')
            chosen = var.get().strip() if var else ''
            return chosen or str(self.paths.default_user_results_dir(self.repo_dir))

    def _prisma_report_path(self):
        return os.path.join(self._results_dir(),
                            f'{self._current_prefix()}_prisma_flow_report.txt')

    def _open_prisma_figure(self):
        """Prefer the SVG: its text can be selected and copied out."""
        base = os.path.join(self._results_dir(), 'figures',
                            f'{self._current_prefix()}_prisma_flow')
        for ext in ('svg', 'png', 'jpeg', 'tif', 'pdf'):
            if os.path.exists(base + '.' + ext):
                self._open_path(base + '.' + ext)
                return
        messagebox.showinfo('Not there yet',
                            'The flow diagram is written when a run finishes.')

    def refresh_results(self):
        d = self._results_dir()
        prefix = self._current_prefix()

        # -- the overview -------------------------------------------------
        report = self._prisma_report_path()
        self.res_overview.config(state='normal')
        self.res_overview.delete('1.0', 'end')
        if os.path.exists(report):
            try:
                with open(report, encoding='utf-8') as fh:
                    self.res_overview.insert('end', fh.read())
                self.btn_res_report.config(state='normal')
            except OSError as exc:
                self.res_overview.insert('end', f'The report could not be read:\n{exc}')
        else:
            self.btn_res_report.config(state='disabled')
            self.res_overview.insert('end',
                '  No run overview yet for the prefix "%s".\n\n'
                '  One is written at the end of every run, into\n'
                '      %s\n\n'
                '  It accounts for the whole search: records retrieved per\n'
                '  citation generation, what the MeSH screen excluded and why,\n'
                '  what each optimiser kept, what the consensus and the largest\n'
                '  connected component discarded, and what was finally included.\n'
                % (prefix, report))
        self.res_overview.config(state='disabled')

        figures = os.path.join(d, 'figures')
        has_fig = any(os.path.exists(os.path.join(figures, f'{prefix}_prisma_flow.{e}'))
                      for e in ('svg', 'png', 'jpeg', 'tif', 'pdf'))
        self.btn_res_figure.config(state='normal' if has_fig else 'disabled')

        # -- the file listing ---------------------------------------------
        rows = []
        for base, _, files in os.walk(d):
            for f in sorted(files):
                path = os.path.join(base, f)
                try:
                    rows.append((os.path.getmtime(path),
                                 os.path.relpath(path, d).replace('\\', '/'),
                                 os.path.getsize(path)))
                except OSError:
                    continue
        rows.sort(reverse=True)
        self.res_title.config(text=f'Results  -  {len(rows)} files in {d}')
        self.res_text.config(state='normal')
        self.res_text.delete('1.0', 'end')
        if not rows:
            self.res_text.insert('end', '  Nothing here yet.\n')
        else:
            import datetime
            self.res_text.insert('end', f'{"modified":<18}{"size":>12}  file\n')
            self.res_text.insert('end', '-' * 92 + '\n')
            for mt, rel, sz in rows[:400]:
                ts = datetime.datetime.fromtimestamp(mt).strftime('%Y-%m-%d %H:%M')
                self.res_text.insert('end', f'{ts:<18}{sz/1024:>10,.0f} KB  {rel}\n')
            if len(rows) > 400:
                self.res_text.insert('end', f'... and {len(rows)-400} more\n')
        self.res_text.config(state='disabled')

    def open_results(self):
        self.refresh_results()
        self.show('results')

    # -------------------------------------------------------- screen: integrity
    def _build_integrity(self, root):
        tk.Label(root, text='Check and repair files', bg=FACE, font=self.f_bold,
                 anchor='w').pack(fill='x', padx=10, pady=(10, 2))
        tk.Label(root, bg=FACE, anchor='w', justify='left', wraplength=760,
                 text='A run that was interrupted - a machine that slept, a disk '
                      'that filled, a sync client mid-write - can leave a file '
                      'that looks complete and is not. Everything the pipeline '
                      'depends on is opened and read here.\n'
                      'Nothing listed is your own work: all of it is rebuilt from '
                      'what is left. Your results are never touched.'
                 ).pack(fill='x', padx=10, pady=(0, 8))

        box = self._sunken(root)
        box.pack(fill='both', expand=True, padx=10)
        canvas = tk.Canvas(box, bg=FACE, highlightthickness=0)
        bar = tk.Scrollbar(box, command=canvas.yview)
        canvas.config(yscrollcommand=bar.set)
        bar.pack(side='right', fill='y')
        canvas.pack(side='left', fill='both', expand=True)
        self.ck_rows = tk.Frame(canvas, bg=FACE)
        canvas.create_window((0, 0), window=self.ck_rows, anchor='nw')
        self.ck_rows.bind('<Configure>',
                          lambda e: canvas.config(scrollregion=canvas.bbox('all')))

        foot = tk.Frame(root, bg=FACE)
        foot.pack(fill='x', padx=10, pady=9)
        self.ck_summary = tk.Label(foot, bg=FACE, relief='sunken', bd=2, anchor='w',
                                   justify='left')
        self.ck_summary.pack(fill='x', ipady=2, pady=(0, 6))
        act = tk.Frame(foot, bg=FACE)
        act.pack(fill='x')
        self.btn_ck_scan = tk.Button(act, text='Scan again', width=13,
                                     command=lambda: self.scan_integrity(deep=False))
        self.btn_ck_scan.pack(side='left')
        self.btn_ck_deep = tk.Button(act, text='Thorough check', width=15,
                                     command=lambda: self.scan_integrity(deep=True))
        self.btn_ck_deep.pack(side='left', padx=4)
        self.btn_ck_fix = tk.Button(act, text='Remove ticked', width=15,
                                    font=self.f_bold, state='disabled',
                                    command=self.repair_integrity)
        self.btn_ck_fix.pack(side='left', padx=4)
        tk.Button(act, text='Back to settings', width=16,
                  command=lambda: self.show('settings')).pack(side='right')

    def open_integrity(self):
        self.show('integrity')
        self.after(60, lambda: self.scan_integrity(deep=False))

    def scan_integrity(self, deep=False):
        """Examine every artefact and list what needs a decision."""
        from mesh_aop import integrity
        from mesh_aop.config_parser import MeshConfig

        for w in self.ck_rows.winfo_children():
            w.destroy()
        self._scan_rows = []
        self.ck_summary.config(text='Checking...' + (
            ' A thorough check reads the whole database and can take minutes.'
            if deep else ''), fg=INK)
        self.update_idletasks()

        if not self.save_cfg():
            return
        try:
            config = MeshConfig(config_path=self.cfg_path)
            artifacts = integrity.scan(config, deep=deep)
        except Exception as exc:
            self.ck_summary.config(text=f'The check could not run: {exc}', fg=ERR)
            return

        colour = {integrity.OK: OK, integrity.MISSING: DIM,
                  integrity.SUSPECT: WARN}
        for a in artifacts:
            row = tk.Frame(self.ck_rows, bg=FACE)
            row.pack(fill='x', pady=1)
            needs_action = a.status in (integrity.EMPTY, integrity.CORRUPT,
                                        integrity.ORPHAN, integrity.SUSPECT)
            if needs_action:
                var = tk.BooleanVar(value=a.broken)   # suspect starts unticked
                tk.Checkbutton(row, variable=var, bg=FACE,
                               command=self._ck_recount).pack(side='left')
                self._scan_rows.append((var, a))
            else:
                tk.Label(row, text='   ', bg=FACE).pack(side='left')
            mark = {integrity.OK: '[ok]', integrity.MISSING: '[--]'}.get(a.status, '[!!]')
            tk.Label(row, text=mark, bg=FACE, font=self.f_mono,
                     fg=colour.get(a.status, ERR), width=5, anchor='w').pack(side='left')
            tk.Label(row, text=a.label, bg=FACE, anchor='w', width=34,
                     font=self.f_bold if needs_action else self.f_ui).pack(side='left')
            detail = 'not built yet' if a.status == integrity.MISSING else a.detail
            if a.status != integrity.MISSING and a.cost == integrity.COST_HOURS \
                    and a.status != integrity.OK:
                detail += '   [rebuilding this takes hours]'
            tk.Label(row, text=detail, bg=FACE, anchor='w', wraplength=420,
                     justify='left', fg=colour.get(a.status, ERR)).pack(side='left')

        text = integrity.summarise(artifacts)
        bad = integrity.problems(artifacts)
        self.ck_summary.config(text=text.split('\n\n')[0] if not bad else text,
                               fg=OK if not bad else WARN)
        self._ck_recount()
        self._scan_artifacts = artifacts

    def _ck_recount(self):
        chosen = [a for var, a in self._scan_rows if var.get()]
        self.btn_ck_fix.config(state='normal' if chosen else 'disabled',
                               text=f'Remove {len(chosen)} file(s)' if chosen
                               else 'Remove ticked')

    def repair_integrity(self):
        from mesh_aop import integrity
        chosen = [a for var, a in self._scan_rows if var.get()]
        if not chosen:
            return
        costly = [a for a in chosen if a.cost == integrity.COST_HOURS]
        warning = ''
        if costly:
            warning = ('\n\nThese take HOURS to rebuild:\n'
                       + '\n'.join(f'    {a.label}  ({a.size_mb:,.0f} MB)'
                                   for a in costly))
        total = sum(a.size for a in chosen) / 1e6
        if not messagebox.askyesno(
                'Remove these files?',
                f'{len(chosen)} file(s) will be deleted, freeing {total:,.0f} MB.'
                f'{warning}\n\nThe pipeline rebuilds whatever is missing on the '
                'next run. Remove them?', icon='warning', default='no'):
            return
        removed, freed, failures = integrity.remove(chosen)
        step = integrity.resume_step(getattr(self, '_scan_artifacts', []))
        msg = f'Removed {len(removed)} file(s), {freed/1e6:,.0f} MB freed.'
        if failures:
            msg += ('\n\nCould not remove:\n'
                    + '\n'.join(f'    {a.label}: {err}' for a, err in failures)
                    + '\n\nA file in use cannot be deleted - close anything that '
                      'has it open, or restart, and try again.')
        if step:
            msg += (f'\n\nResume from:\n    {integrity.STEP_LABEL.get(step, step)}')
        messagebox.showinfo('Repair', msg)
        self.scan_integrity(deep=False)

    # -------------------------------------------------------- screen: uninstall
    def _build_uninstall(self, root):
        tk.Label(root, text='Uninstall MeSH Workbench', bg=FACE, font=self.f_bold,
                 anchor='w').pack(fill='x', padx=10, pady=(10, 2))
        tk.Label(root, bg=FACE, anchor='w', justify='left',
                 text='Everything this program downloaded, built or installed is '
                      'listed below - including files it keeps outside this '
                      'folder. Choose what to remove.'
                 ).pack(fill='x', padx=10, pady=(0, 8))

        box = self._sunken(root)
        box.pack(fill='both', expand=True, padx=10)
        canvas = tk.Canvas(box, bg=FACE, highlightthickness=0)
        bar = tk.Scrollbar(box, command=canvas.yview)
        canvas.config(yscrollcommand=bar.set)
        bar.pack(side='right', fill='y')
        canvas.pack(side='left', fill='both', expand=True)
        self.un_rows = tk.Frame(canvas, bg=FACE)
        canvas.create_window((0, 0), window=self.un_rows, anchor='nw')
        self.un_rows.bind('<Configure>',
                          lambda e: canvas.config(scrollregion=canvas.bbox('all')))
        self.un_canvas = canvas

        foot = tk.Frame(root, bg=FACE)
        foot.pack(fill='x', padx=10, pady=9)
        self.un_total = tk.Label(foot, bg=FACE, relief='sunken', bd=2, anchor='w')
        self.un_total.pack(fill='x', ipady=2, pady=(0, 6))

        act = tk.Frame(foot, bg=FACE)
        act.pack(fill='x')
        tk.Label(act, text='Type REMOVE to confirm:', bg=FACE).pack(side='left')
        self.un_confirm = tk.StringVar()
        ent = tk.Entry(act, textvariable=self.un_confirm, width=12, bg=FIELD,
                       relief='sunken', bd=2)
        ent.pack(side='left', padx=6)
        self.un_go = tk.Button(act, text='Remove', width=12, font=self.f_bold,
                               state='disabled', command=self._do_uninstall)
        self.un_go.pack(side='left')
        tk.Button(act, text='Cancel', width=11,
                  command=lambda: self.show('settings')).pack(side='right')
        tk.Button(act, text='Rescan', width=11,
                  command=self.scan_uninstall).pack(side='right', padx=4)
        # The button stays dead until the word is typed exactly - this deletes
        # tens of gigabytes and there is no undo.
        self.un_confirm.trace_add('write', lambda *_: self._un_gate())

    def _un_gate(self):
        ready = (self.un_confirm.get().strip() == 'REMOVE'
                 and any(v.get() for v in getattr(self, 'un_vars', {}).values()))
        self.un_go.config(state='normal' if ready else 'disabled')

    def open_uninstall(self):
        self.show('uninstall')
        self.after(60, self.scan_uninstall)

    def scan_uninstall(self):
        from mesh_aop import uninstall as U
        for w in self.un_rows.winfo_children():
            w.destroy()
        self.un_total.config(text='   Scanning…')
        self.update_idletasks()

        # No MeshConfig here: constructing one calls refresh_paths(), which
        # creates every directory it resolves - so the uninstaller would
        # recreate data/ and results/ on the rescan that follows a removal.
        # Sorted once and kept in that order. Enumerating a sorted copy while
        # indexing the unsorted list maps each checkbox to the wrong item, which
        # is how a tick on "results" ends up deleting the archive.
        self.un_items = sorted(U.inventory(self.repo_dir),
                               key=lambda x: (x.category, -x.bytes))
        self.un_vars = {}

        current = None
        for i, it in enumerate(self.un_items):
            if it.category != current:
                current = it.category
                tk.Label(self.un_rows, text=U.CATEGORY_LABEL[current], bg=FACE,
                         font=self.f_bold, anchor='w').pack(fill='x', pady=(8, 1), padx=6)
            row = tk.Frame(self.un_rows, bg=FACE)
            row.pack(fill='x', padx=6)
            if it.removable:
                var = tk.BooleanVar(value=it.category in U.DEFAULT_SELECTED)
                self.un_vars[i] = var
                tk.Checkbutton(row, text=it.label, variable=var, bg=FACE,
                               activebackground=FACE, anchor='w',
                               command=self._un_recount).pack(side='left')
            else:
                tk.Label(row, text=f'{it.label}  (kept)', bg=FACE, fg=DIM,
                         anchor='w').pack(side='left')
            tk.Label(row, text=f'{it.gb:,.2f} GB', bg=FACE, anchor='e'
                     ).pack(side='right')
            if it.note:
                tk.Label(self.un_rows, text=f'      {it.note}', bg=FACE, fg=DIM,
                         anchor='w', justify='left').pack(fill='x', padx=6)
        self._un_recount()

    def _un_recount(self):
        sel = sum(self.un_items[i].bytes for i, v in self.un_vars.items() if v.get())
        tot = sum(i.bytes for i in self.un_items)
        self.un_total.config(
            text=f'   Selected {sel/1e9:,.2f} GB of {tot/1e9:,.2f} GB found')
        self._un_gate()

    def _do_uninstall(self):
        from mesh_aop import uninstall as U
        chosen = [self.un_items[i] for i, v in self.un_vars.items() if v.get()]
        if not chosen:
            return
        gb = sum(c.bytes for c in chosen) / 1e9
        if not messagebox.askyesno(
                'Remove these files?',
                f'{len(chosen)} item(s) will be deleted, freeing {gb:,.2f} GB.\n\n'
                'This cannot be undone.\n\nContinue?', icon='warning'):
            return

        portable = U.is_portable(self.repo_dir)
        removed, freed, failures = U.remove(chosen)
        msg = f'Removed {removed} of {len(chosen)} item(s), freeing {freed/1e9:,.2f} GB.'
        if failures:
            msg += (f'\n\n{len(failures)} could not be removed, usually because '
                    'another program still has them open:\n\n' +
                    '\n'.join(f'{p}: {e}' for p, e in failures[:5]))
        if not portable:
            msg += ('\n\nThe package itself is still installed. To remove it, run:\n\n'
                    f'    {U.pip_hint()}')
        else:
            msg += ('\n\nTo finish, close this window and delete the program '
                    'folder:\n\n' + self.repo_dir)
        messagebox.showinfo('Uninstall', msg)
        self.scan_uninstall()

    def _open_results(self):
        self._open_path(self._results_dir())

    def _on_close(self):
        if self.runner.is_running():
            if not messagebox.askyesno('Quit?', 'A pipeline step is still '
                                                'running. Stop it and quit?',
                                       icon='warning'):
                return
            self.runner.cancel()
        if self._tick_job:
            try:
                self.after_cancel(self._tick_job)
            except Exception:
                pass
        self.destroy()


class _Bar(tk.Canvas):
    """Blocky Win95 progress bar - discrete cells, not a smooth fill."""

    def __init__(self, parent, height=18):
        super().__init__(parent, height=height, bg=FIELD, relief='sunken', bd=2,
                         highlightthickness=0)
        self.value = 0.0
        self.bind('<Configure>', lambda _e: self._draw())

    def set(self, pct):
        self.value = max(0.0, min(100.0, float(pct)))
        self._draw()

    def _draw(self):
        self.delete('all')
        w = self.winfo_width() or 1
        h = self.winfo_height() or 1
        cell, gap = 9, 2
        n = int((w - 4) / (cell + gap))
        lit = int(round(n * self.value / 100.0))
        for i in range(lit):
            x = 2 + i * (cell + gap)
            self.create_rectangle(x, 2, x + cell, h - 2, fill=NAVY, outline='')


def _size_gb(path):
    if not os.path.exists(path):
        return None
    if os.path.isfile(path):
        return os.path.getsize(path) / 1e9
    tot = 0
    for base, _, files in os.walk(path):
        for f in files:
            try:
                tot += os.path.getsize(os.path.join(base, f))
            except OSError:
                pass
    return tot / 1e9


def _hms(s):
    return '%02d:%02d:%02d' % (s // 3600, (s // 60) % 60, s % 60)


def main(repo_dir=None, python_exe=None):
    repo_dir = repo_dir or os.environ.get('MESH_REPO') or os.getcwd()
    Workbench(repo_dir, python_exe or sys.executable).mainloop()
