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
from .runner import DONE, LOG, PHASE, PROGRESS, TRANSIENT, PipelineRunner

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
W, H = 880, 610


class Workbench(tk.Tk):
    def __init__(self, repo_dir, python_exe):
        super().__init__()
        self.repo_dir = repo_dir
        self.python_exe = python_exe
        # A portable copy keeps its settings in its own folder; an installed one
        # keeps them under the user's profile, because the install directory is
        # shared and often read-only. paths.py decides, and the CLI asks it too.
        from mesh_aop import paths as _paths
        self.paths = _paths
        self.cfg_path = str(_paths.config_path(repo_dir))
        self.cfg = self._load_cfg()
        self.vars = {}
        self.runner = PipelineRunner(repo_dir, python_exe)
        self._elapsed = 0
        self._tick_job = None

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
        if not os.path.exists(self.cfg_path):
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
        t.add_command(label='Uninstall…', command=self.open_uninstall)
        bar.add_cascade(label='Tools', menu=t)
        h = tk.Menu(bar, tearoff=0)
        h.add_command(label='About', command=lambda: messagebox.showinfo(
            'About', f'MeSH Workbench {__version__}\n\nDesktop front-end for the '
                     'mesh_aop pipeline.'))
        bar.add_cascade(label='Help', menu=h)
        self.config(menu=bar)

    def first_run_locations(self):
        """Ask where results and data should live, the first time only.

        Shown when no settings file exists yet. Only the two locations the user
        genuinely owns are asked about - settings, logs and scratch are
        machinery and go to the profile without discussion.

        A portable copy is self-contained by definition, so it is never asked.
        """
        if self.paths.is_portable(self.repo_dir) or os.path.exists(self.cfg_path):
            return

        win = tk.Toplevel(self)
        win.title('MeSH Workbench - first run')
        win.configure(bg=FACE)
        win.transient(self)
        win.resizable(False, False)
        win.protocol('WM_DELETE_WINDOW', lambda: None)   # only leaves via a button

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

        def accept():
            for key, var in chosen.items():
                value = var.get().strip()
                if value:
                    schema.put(self.cfg, key, value)
            # Writing now means the dialog does not reappear, and the rest of
            # the application can assume a settings file exists.
            self.save_cfg()
            win.destroy()

        tk.Button(act, text='Continue', width=14, font=self.f_bold, command=accept
                  ).pack(side='right')
        tk.Button(act, text='Use defaults', width=14,
                  command=lambda: (self.save_cfg(), win.destroy())).pack(side='right', padx=6)

        win.update_idletasks()
        x = self.winfo_rootx() + (self.winfo_width() - win.winfo_width()) // 2
        y = self.winfo_rooty() + 80
        win.geometry(f'+{max(x, 0)}+{max(y, 0)}')
        win.grab_set()
        self.wait_window(win)
        self.scan_data()

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
    def _build_setup(self, root):
        tk.Label(root, text='Local data sources', bg=FACE, font=self.f_bold,
                 anchor='w').pack(fill='x', padx=10, pady=(10, 2))
        tk.Label(root, bg=FACE, anchor='w', justify='left',
                 text='The pipeline runs offline once these are in place. '
                      'Nothing is re-downloaded unless you ask.'
                 ).pack(fill='x', padx=10, pady=(0, 8))

        box = self._sunken(root)
        box.pack(fill='both', expand=True, padx=10)
        # header uses the same grid geometry as the rows, or the columns drift
        hdr = tk.Frame(box, bg=FACE)
        hdr.pack(fill='x', padx=4, pady=(3, 0))
        self._setup_cols(hdr)
        for c, (text, anchor) in enumerate((('Source', 'w'), ('Status', 'w'),
                                            ('Size', 'e'), ('', 'e'))):
            tk.Label(hdr, text=text, bg=FACE, font=self.f_bold, anchor=anchor
                     ).grid(row=0, column=c, sticky='we',
                            padx=(0, 8) if c == 2 else 0)
        tk.Frame(box, bg='#808080', height=1).pack(fill='x', padx=4, pady=(2, 0))
        self.setup_rows = tk.Frame(box, bg=FACE)
        self.setup_rows.pack(fill='both', expand=True)

        foot = tk.Frame(root, bg=FACE)
        foot.pack(fill='x', padx=10, pady=9)
        self.setup_total = tk.Label(foot, bg=FACE, relief='sunken', bd=2, anchor='w')
        self.setup_total.pack(side='left', fill='x', expand=True, ipady=2, padx=(0, 8))
        tk.Button(foot, text='Refresh', width=11, command=self.scan_data
                  ).pack(side='right', padx=2)
        tk.Button(foot, text='Continue', width=11, font=self.f_bold,
                  command=lambda: self.show('settings')).pack(side='right', padx=2)
        self.after(120, self.scan_data)

    def scan_data(self):
        """Measure what is actually on disk - no assumed sizes."""
        for w in self.setup_rows.winfo_children():
            w.destroy()
        d = self.repo_dir
        items = [
            ('Master annotation database', os.path.join(d, 'data/raw/master_mesh_database.db'),
             'Built from the archive. Required for every run.'),
            ('PubMed baseline archive', os.path.join(d, 'data/raw/pubmed_baseline'),
             'Only needed to build or rebuild the database.'),
            ('MeSH descriptor file', os.path.join(d, 'data/raw/desc2025.xml'),
             'Defines the stop-word vocabulary.'),
            ('Relevance database', os.path.join(d, 'data/processed/DAC_Mesh_mean_relevancy.db'),
             'Per-article scores for the current project.'),
            ('Citation database', os.path.join(d, 'data/processed/DAC_Mesh_cleaned_pmids.db'),
             'Query results and citation generations.'),
        ]
        total = 0.0
        reclaim = 0.0
        for label, path, note in items:
            gb = _size_gb(path)
            total += gb or 0.0
            if 'baseline' in label.lower():
                reclaim = gb or 0.0
            # grid, not nested packs: a fixed-width frame with pack_propagate
            # off collapses to zero height and clips its labels to slivers.
            row = tk.Frame(self.setup_rows, bg=FACE)
            row.pack(fill='x', pady=2, padx=4)
            self._setup_cols(row)
            present = gb is not None
            tk.Label(row, text=label, bg=FACE, font=self.f_bold, anchor='w'
                     ).grid(row=0, column=0, sticky='w')
            tk.Label(row, text=note, bg=FACE, fg=DIM, anchor='w'
                     ).grid(row=1, column=0, sticky='w')
            tk.Label(row, text=('Present' if present else 'Missing'), bg=FACE,
                     fg=(OK if present else ERR), font=self.f_bold, anchor='w'
                     ).grid(row=0, column=1, sticky='w')
            tk.Label(row, text=(f'{gb:,.1f} GB' if present else '—'), bg=FACE,
                     anchor='e').grid(row=0, column=2, sticky='e', padx=(0, 8))
            if 'baseline' in label.lower() and present:
                tk.Button(row, text='Delete & reclaim',
                          command=lambda p=path, g=gb: self.reclaim(p, g)
                          ).grid(row=0, column=3, rowspan=2, sticky='e', padx=4)
            # Step 0 is the only way to produce the master database, and nothing
            # else in the pipeline runs without it.
            if label.startswith('Master annotation'):
                tk.Button(row, text=('Rebuild' if present else 'Build'),
                          command=lambda p=present: self.build_database(p)
                          ).grid(row=0, column=3, rowspan=2, sticky='e', padx=4)
            tk.Frame(row, bg='#b0b0b0', height=1).grid(
                row=2, column=0, columnspan=4, sticky='we', pady=(3, 0))
        self.setup_total.config(
            text=f'   On disk now {total:,.1f} GB' +
                 (f'   ·   after reclaiming the archive {total - reclaim:,.1f} GB'
                  if reclaim else ''))

    def reclaim(self, path, gb):
        if not messagebox.askyesno(
                'Delete the PubMed baseline archive?',
                f'This removes the archive and frees {gb:,.1f} GB.\n\n'
                'The master annotation database is unaffected and the pipeline '
                'keeps working.\n\nYou would need to download the archive again '
                '- several hours - only if you rebuild the database from '
                'scratch.\n\nDelete it?', icon='warning'):
            return
        try:
            shutil.rmtree(path)
        except OSError as exc:
            messagebox.showerror('Could not delete',
                                 f'The archive was not removed:\n\n{exc}')
            return
        self.scan_data()

    def build_database(self, present):
        """Run Step 0: fetch the PubMed baseline and compile the master database.

        The switches travel on the command line, not through mesh_config.json -
        a rebuild recorded in the config would repeat on every later run.
        """
        archive = os.path.join(self.repo_dir, 'data/raw/pubmed_baseline')
        have_archive = os.path.isdir(archive) and any(os.scandir(archive))
        extra = ['--build-database']

        if present and not messagebox.askyesno(
                'Rebuild the master database?',
                'The existing master annotation database will be deleted and '
                'built again from the archive.\n\nThis takes hours. Nothing '
                'else in the pipeline can run while it is missing.\n\nRebuild?',
                icon='warning'):
            return
        if present:
            extra.append('--rebuild-corrupt')

        if have_archive:
            if messagebox.askyesno(
                    'Use the archive already on disk?',
                    'A PubMed baseline archive is already here.\n\n'
                    'Yes  - compile from it, no download.\n'
                    'No   - download the archive again first (~44 GB).',):
                extra.append('--skip-baseline-download')
        elif not messagebox.askyesno(
                'Download the PubMed baseline?',
                'No archive is present, so it has to be downloaded first: '
                'roughly 44 GB, and several hours on a fast connection.\n\n'
                'The download resumes if it is interrupted, and the archive can '
                'be deleted afterwards.\n\nStart?', icon='warning'):
            return

        if messagebox.askyesno(
                'Include the daily updates?',
                'The baseline is a yearly snapshot. The daily update files carry '
                'records published since.\n\nFetch them too?'):
            extra.append('--with-updates')

        self.start_run('baseline', extra, title='building the master database')

    # ---------------------------------------------------------- screen: settings
    def _build_settings(self, root):
        top = self._sunken(root)
        top.pack(fill='x', padx=10, pady=(10, 8))
        tk.Label(top, text='Pipeline step:', bg=FACE).pack(side='left', padx=(8, 6), pady=6)
        self.step_var = tk.StringVar(value='benchmark')
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
            self._build_fields(pane, fields)

        self.help = tk.Frame(root, bg=FACE, relief='sunken', bd=2, height=88)
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

    def _build_fields(self, pane, fields):
        grid = tk.Frame(pane, bg=FACE)
        grid.pack(fill='both', expand=True, padx=10, pady=10)
        grid.columnconfigure(1, weight=1)
        for r, f in enumerate(fields):
            cur = schema.get(self.cfg, f.key, f.default)
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
        self.btn_cancel = tk.Button(act, text='Cancel', width=12, command=self.cancel_run)
        self.btn_cancel.pack(side='left')
        self.btn_back = tk.Button(act, text='Back to settings', width=16,
                                  command=lambda: self.show('settings'))
        self.btn_back.pack(side='left', padx=4)
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
        self.runner.start(step, extra)
        self._pump()

    def cancel_run(self):
        if not self.runner.is_running():
            self.show('settings')
            return
        if messagebox.askyesno('Cancel the run?',
                               'The pipeline stops where it is. Partial outputs '
                               'may be left on disk.\n\nCancel it?', icon='warning'):
            self.runner.cancel()
            self.run_status.config(text='Cancelling…')

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
            elif kind == DONE:
                self._finish(ev[1], ev[2])
                return
        if self.runner.is_running():
            self._elapsed += 0.2
            if int(self._elapsed * 5) % 3 == 0:
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
        if rc == 0:
            self._log(f'--- completed in {elapsed/60:.1f} min ---', 'ok')
            self.run_status.config(text='Completed successfully.', fg=OK)
            self.refresh_results()
            self.btn_results.config(state='normal')
        elif rc == -1:
            self._log('--- cancelled ---', 'warn')
            self.run_status.config(text='Cancelled.', fg=WARN)
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
        self.res_title.pack(fill='x', padx=10, pady=(10, 6))
        box = self._sunken(root)
        box.pack(fill='both', expand=True, padx=10)
        self.res_text = tk.Text(box, bg=FIELD, font=self.f_mono, relief='flat',
                                wrap='none')
        sb = tk.Scrollbar(box, command=self.res_text.yview)
        self.res_text.config(yscrollcommand=sb.set)
        sb.pack(side='right', fill='y')
        self.res_text.pack(fill='both', expand=True)
        act = tk.Frame(root, bg=FACE)
        act.pack(fill='x', padx=10, pady=9)
        tk.Button(act, text='Open results folder', width=18, font=self.f_bold,
                  command=self._open_results).pack(side='left')
        tk.Button(act, text='Refresh', width=12,
                  command=self.refresh_results).pack(side='left', padx=4)
        tk.Button(act, text='Back to settings', width=16,
                  command=lambda: self.show('settings')).pack(side='left', padx=4)

    def open_results(self):
        self.refresh_results()
        self.show('results')

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

    def refresh_results(self):
        d = os.path.join(self.repo_dir, 'results')
        rows = []
        for base, _, files in os.walk(d):
            for f in sorted(files):
                p = os.path.join(base, f)
                try:
                    sz = os.path.getsize(p)
                    mt = os.path.getmtime(p)
                except OSError:
                    continue
                rel = os.path.relpath(p, d).replace('\\', '/')
                rows.append((mt, rel, sz))
        rows.sort(reverse=True)
        self.res_title.config(text=f'Results — {len(rows)} files in {d}')
        self.res_text.config(state='normal')
        self.res_text.delete('1.0', 'end')
        self.res_text.insert('end', f'{"modified":<18}{"size":>12}  file\n')
        self.res_text.insert('end', '-' * 92 + '\n')
        import datetime
        for mt, rel, sz in rows[:400]:
            ts = datetime.datetime.fromtimestamp(mt).strftime('%Y-%m-%d %H:%M')
            self.res_text.insert('end', f'{ts:<18}{sz/1024:>10,.0f} KB  {rel}\n')
        self.res_text.config(state='disabled')

    def _open_results(self):
        d = os.path.join(self.repo_dir, 'results')
        try:
            if os.name == 'nt':
                os.startfile(d)                       # noqa: S606
            else:
                subprocess.Popen(['xdg-open', d])
        except Exception as exc:
            messagebox.showerror('Open folder', str(exc))

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
