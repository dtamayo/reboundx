from codecs import open
import os
import inspect
import sys 
import sysconfig

def get_reb_paths(sitepackagesdir):
    try:
        import rebound
        rebdir = os.path.dirname(inspect.getfile(rebound))
        version = rebound.__version__ 
    except:
        raise AttributeError("REBOUND was not installed.")

    try: # try to get local rebound directory if using editable pip installs
        with open(sitepackagesdir+'rebound-'+version+".dist-info/direct_url.json") as f:
            lines = f.readlines()
            for l in lines:
                blocks = l.split('"')
                if 'url' in blocks:
                    for block in blocks:
                        if block.startswith('file://'):
                            path = block.strip('file:')
        return rebdir, path+'/'
    except:
        return rebdir, ""

try:
    from setuptools import setup, Extension
    from setuptools.command.build_ext import build_ext as _build_ext
except ImportError:
    print("Installing REBOUNDx requires setuptools.  Do 'pip install setuptools'.")
    sys.exit(1)

suffix = sysconfig.get_config_var('EXT_SUFFIX')
if suffix is None:
    suffix = ".so"

# Try to get git hash
try:
    import subprocess
    ghash = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode("ascii")
    ghash_arg = "-DREBXGITHASH="+ghash.strip()
except:
    ghash_arg = "-DREBXGITHASH=49fd028521361fecd07728faa0b0c829efb2f908" #GITHASHAUTOUPDATE

def _find_dumpbin(compiler=None):
    """Locate dumpbin.exe. It ships alongside cl.exe in every MSVC install.

    Order of attempts:
      1. Alongside the cl.exe used by `compiler` (if given and initialized).
      2. Anywhere on PATH (after Developer Command Prompt activation).
      3. Glob the standard Visual Studio install trees.
    """
    import glob
    import shutil

    if compiler is not None and getattr(compiler, 'cc', None):
        sibling = os.path.join(os.path.dirname(compiler.cc), 'dumpbin.exe')
        if os.path.exists(sibling):
            return sibling

    on_path = shutil.which('dumpbin.exe')
    if on_path:
        return on_path

    program_files_roots = [
        os.environ.get('ProgramFiles', r'C:\Program Files'),
        os.environ.get('ProgramFiles(x86)', r'C:\Program Files (x86)'),
    ]
    candidates = []
    for root in program_files_roots:
        pattern = os.path.join(
            root, 'Microsoft Visual Studio', '*', '*',
            'VC', 'Tools', 'MSVC', '*', 'bin', 'Host*', 'x64', 'dumpbin.exe',
        )
        candidates.extend(glob.glob(pattern))
    if candidates:
        # glob.glob() does not promise a stable order across filesystems;
        # sort lexicographically so the picked dumpbin is reproducible
        # across machines and over time. The last entry is, in practice,
        # the highest-version MSVC install since path components are
        # version-numbered (e.g., 14.39.33519).
        candidates.sort()
        return candidates[-1]
    raise RuntimeError(
        "dumpbin.exe not found. Install Visual Studio Build Tools or run "
        "this inside a Developer Command Prompt."
    )


def _write_def_from_objs(obj_files, def_path, module_name, compiler=None):
    """Generate a .def file listing every public, defined symbol found in
    the compiled .obj files. Same approach CMake uses for its
    WINDOWS_EXPORT_ALL_SYMBOLS target property.

    MSVC doesn't auto-export DLL symbols (unlike gcc); without a .def
    or __declspec(dllexport) annotations everywhere, reboundx's ctypes
    layer can't find rebx_version_str or any of the rebx_* functions.
    """
    import subprocess
    dumpbin = _find_dumpbin(compiler)
    symbols = set()
    for obj in obj_files:
        out = subprocess.run(
            [dumpbin, "-SYMBOLS", obj],
            capture_output=True, text=True, check=True,
        ).stdout
        for line in out.splitlines():
            # dumpbin "External | name" lines; skip UNDEF (imports) and keep
            # only symbols that look like C identifiers. Internal MSVC-mangled
            # names start with $ or ?, which we filter out.
            if "External" not in line or "UNDEF" in line:
                continue
            parts = line.split("|")
            if len(parts) < 2:
                continue
            name = parts[-1].strip().split()[0]
            if not name or not (name[0].isalpha() or name[0] == "_"):
                continue
            if name.startswith("$") or name.startswith("?") or name.startswith("__"):
                continue
            # Don't export the PyInit stub; Python already sees it via the
            # module's export table, and listing it in a .def confuses things.
            if name.startswith("PyInit_"):
                continue
            symbols.add(name)

    with open(def_path, "w") as f:
        f.write(f"LIBRARY {module_name}\n")
        f.write("EXPORTS\n")
        for sym in sorted(symbols):
            f.write(f"    {sym}\n")


def _ensure_rebound_import_lib(rebdir, compiler=None):
    """Ensure an MSVC-compatible import library exists for rebound's DLL.

    rebound's Windows wheel ships librebound.<suffix>.pyd (the DLL) but no
    corresponding .lib (import library), which MSVC's linker needs. We
    generate one on the fly using dumpbin + lib.exe, placing it next to
    the DLL so setup.py's default library search path finds it.

    No-op if the .lib already exists (user or earlier build produced it).
    Unix is untouched — .so linking doesn't need an import library.
    """
    import glob
    import subprocess

    # Find librebound.<suffix>.pyd — the DLL. Prefer the .pyd whose ABI
    # tag matches the active interpreter's EXT_SUFFIX so we don't import
    # the wrong CPython build's binary in mixed environments.
    import sysconfig as _sc
    ext_suffix = _sc.get_config_var('EXT_SUFFIX') or '.pyd'
    pattern = os.path.join(rebdir, '..', 'librebound*.pyd')
    dll_candidates = glob.glob(pattern)
    if not dll_candidates:
        # Maybe installed as a namespace package with the .pyd inside the folder.
        pattern = os.path.join(rebdir, 'librebound*.pyd')
        dll_candidates = glob.glob(pattern)
    if not dll_candidates:
        raise RuntimeError(
            f"Could not locate librebound*.pyd near {rebdir}. "
            "Is rebound installed?"
        )
    # Prefer ABI-matching candidates; fall back to a deterministic sort.
    matching = [p for p in dll_candidates if p.endswith(ext_suffix)]
    dll_candidates = sorted(matching) if matching else sorted(dll_candidates)
    dll_path = dll_candidates[0]
    dll_dir = os.path.dirname(dll_path)
    dll_name = os.path.basename(dll_path)
    # librebound.cp312-win_amd64.pyd -> librebound.cp312-win_amd64.lib
    lib_name = dll_name.rsplit('.', 1)[0] + '.lib'
    lib_path = os.path.join(dll_dir, lib_name)
    if os.path.exists(lib_path):
        return lib_path  # already there

    dumpbin = _find_dumpbin(compiler)
    # Find lib.exe — it's next to dumpbin.exe.
    lib_exe = os.path.join(os.path.dirname(dumpbin), 'lib.exe')
    if not os.path.exists(lib_exe):
        raise RuntimeError(f"lib.exe not found next to dumpbin at {dumpbin}")

    # Step 1: extract the DLL's exports into a .def file.
    out = subprocess.run(
        [dumpbin, '-exports', dll_path],
        capture_output=True, text=True, check=True,
    ).stdout
    def_path = lib_path.rsplit('.', 1)[0] + '.def'
    with open(def_path, 'w') as f:
        f.write(f"LIBRARY {dll_name}\n")
        f.write("EXPORTS\n")
        for line in out.splitlines():
            # dumpbin -exports lines look like:
            #   1    0 000465D0 PyInit_librebound
            # Grab the 4th whitespace-separated field when the first three
            # are ordinal / hint / RVA (all hex digits).
            parts = line.split()
            if len(parts) < 4:
                continue
            ord_, hint, rva, *name_parts = parts
            if not (ord_.isdigit() and all(c in '0123456789ABCDEFabcdef' for c in hint)
                    and all(c in '0123456789ABCDEFabcdef' for c in rva)):
                continue
            name = name_parts[0]
            if not (name and (name[0].isalpha() or name[0] == '_')):
                continue
            f.write(f"    {name}\n")

    # Step 2: produce the import .lib from the .def. Pick the /machine
    # value from the running interpreter so cross-arch builds (32-bit
    # CPython, ARM64 native) work without a hand-edit. x64 covers the
    # 99% case; the others are wired in for when someone needs them.
    import platform
    machine = platform.machine().upper()
    machine_flag = {
        'AMD64': '/machine:x64',
        'X86_64': '/machine:x64',
        'X86': '/machine:x86',
        'I386': '/machine:x86',
        'ARM64': '/machine:arm64',
    }.get(machine, '/machine:x64')
    subprocess.run(
        [lib_exe, f'/def:{def_path}', f'/out:{lib_path}', machine_flag],
        capture_output=True, text=True, check=True,
    )
    return lib_path


class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)
        if "PYODIDE" in os.environ:
            return None
       
        # get site-packages dir to add to paths in case reb & rebx installed simul in tmp dir
        sitepackagesdir = sysconfig.get_path('platlib')+'/'
        rebdir, editable_rebdir = get_reb_paths(sitepackagesdir)

        print("***", rebdir, "***", sitepackagesdir, "***", editable_rebdir, "***")
        self.include_dirs.append(rebdir)
        #self.include_dirs.append(editable_rebdir)
        sources = [ 'src/central_force.c', 'src/core.c', 'src/exponential_migration.c', 'src/gas_damping_timescale.c', 'src/gas_dynamical_friction.c', 'src/gr.c', 'src/gravitational_harmonics.c', 'src/gr_full.c', 'src/gr_potential.c', 'src/inner_disk_edge.c', 'src/input.c', 'src/integrate_force.c', 'src/integrator_euler.c', 'src/integrator_implicit_midpoint.c', 'src/integrator_rk2.c', 'src/integrator_rk4.c', 'src/interpolation.c', 'src/lense_thirring.c', 'src/linkedlist.c', 'src/modify_mass.c', 'src/modify_orbits_direct.c', 'src/modify_orbits_forces.c', 'src/output.c', 'src/radiation_forces.c', 'src/rebxtools.c', 'src/steppers.c', 'src/stochastic_forces.c', 'src/tides_constant_time_lag.c', 'src/tides_dynamical.c', 'src/tides_spin.c', 'src/track_min_distance.c', 'src/type_I_migration.c', 'src/yarkovsky_effect.c'],
        
        self.library_dirs.append(rebdir+'/../')
        self.library_dirs.append(sitepackagesdir)
        for ext in self.extensions:
            # MSVC's linker doesn't understand rpath or runtime_library_dirs
            # (Windows DLLs use PATH / same-dir resolution instead). Skip them
            # on Windows; on Unix, keep the original behavior.
            if sys.platform != 'win32':
                ext.runtime_library_dirs.append(rebdir+'/../')
                ext.extra_link_args.append('-Wl,-rpath,'+rebdir+'/../')
                ext.runtime_library_dirs.append(sitepackagesdir)
                ext.extra_link_args.append('-Wl,-rpath,'+sitepackagesdir)
        if editable_rebdir:
            self.library_dirs.append(editable_rebdir)
            for ext in self.extensions:
                if sys.platform != 'win32':
                    ext.runtime_library_dirs.append(editable_rebdir)
                    ext.extra_link_args.append('-Wl,-rpath,'+editable_rebdir)

    def build_extension(self, ext):
        # On Windows, pre-compile the objects so we can scan them for
        # exported symbols, write a .def file, and pass /DEF:file.def to
        # the linker. Then let super link normally (it re-checks .obj
        # timestamps and does the link step). Without this, the resulting
        # DLL exports nothing and reboundx/__init__.py fails on the very
        # first ctypes lookup of rebx_version_str.
        if sys.platform == 'win32':
            # Make sure rebound's import library exists (rebound's wheel
            # ships the .pyd but no .lib; MSVC's linker needs both). The
            # helper is a no-op when the .lib is already present.
            import sysconfig as _sc
            sitepackagesdir = _sc.get_path('platlib') + '/'
            rebdir, _ = get_reb_paths(sitepackagesdir)
            _ensure_rebound_import_lib(rebdir, self.compiler)
            sources = list(ext.sources or [])
            extra_postargs = list(ext.extra_compile_args or [])
            macros = list(ext.define_macros or [])
            objects = self.compiler.compile(
                sources, output_dir=self.build_temp,
                macros=macros, include_dirs=ext.include_dirs or [],
                debug=self.debug, extra_postargs=extra_postargs,
                depends=ext.depends or [],
            )
            def_path = os.path.join(self.build_temp, ext.name + '.def')
            _write_def_from_objs(objects, def_path, ext.name + suffix, self.compiler)
            ext.extra_link_args = list(ext.extra_link_args or []) + ['/DEF:' + def_path]
        _build_ext.build_extension(self, ext)


extra_link_args=[]
if sys.platform == 'darwin':
    config_vars = sysconfig.get_config_vars()
    config_vars['LDSHARED'] = config_vars['LDSHARED'].replace('-bundle', '-shared')
    extra_link_args.append('-Wl,-install_name,@rpath/libreboundx'+suffix)
if sys.platform == 'win32':
    # /GL- disables whole-program optimization, which produces "anonymous"
    # LTCG object files whose symbols dumpbin can't read. We need real
    # symbol tables so setup.py can auto-generate the .def export list.
    extra_compile_args=['/GL-', ghash_arg, '-DLIBREBOUNDX', '-D_GNU_SOURCE']
else:
    # Default compile args
    extra_compile_args=['-fstrict-aliasing', '-O3','-std=c99','-Wno-unknown-pragmas', ghash_arg, '-DLIBREBOUNDX', '-D_GNU_SOURCE', '-fPIC']

# Option to disable FMA in CLANG. 
FFP_CONTRACT_OFF = os.environ.get("FFP_CONTRACT_OFF", None)
if FFP_CONTRACT_OFF:
    extra_compile_args.append('-ffp-contract=off')

libreboundxmodule = Extension('libreboundx',
        sources = [ 'src/central_force.c', 'src/core.c', 'src/exponential_migration.c', 'src/gas_damping_timescale.c', 'src/gas_dynamical_friction.c', 'src/gr.c', 'src/gravitational_harmonics.c', 'src/gr_full.c', 'src/gr_potential.c', 'src/inner_disk_edge.c', 'src/input.c', 'src/integrate_force.c', 'src/integrator_euler.c', 'src/integrator_implicit_midpoint.c', 'src/integrator_rk2.c', 'src/integrator_rk4.c', 'src/interpolation.c', 'src/lense_thirring.c', 'src/linkedlist.c', 'src/modify_mass.c', 'src/modify_orbits_direct.c', 'src/modify_orbits_forces.c', 'src/output.c', 'src/radiation_forces.c', 'src/rebxtools.c', 'src/steppers.c', 'src/stochastic_forces.c', 'src/tides_constant_time_lag.c', 'src/tides_dynamical.c', 'src/tides_spin.c', 'src/track_min_distance.c', 'src/type_I_migration.c', 'src/yarkovsky_effect.c'],
                    include_dirs = ['src'],
                    library_dirs = [],
                    runtime_library_dirs = [] if sys.platform == 'win32' else ["."],
                    # On Windows, MSVC's linker doesn't auto-drop the "lib"
                    # prefix the way GCC/ld does, so we must name the full
                    # file on disk. REBOUND ships as librebound<ext-suffix>.pyd;
                    # we also depend on an import library with the same basename.
                    libraries=['librebound'+suffix[:suffix.rfind('.')]]
                              if sys.platform == 'win32'
                              else ['rebound'+suffix[:suffix.rfind('.')]],
                    define_macros=[ ('LIBREBOUNDX', None) ],
                    extra_compile_args=extra_compile_args,
                    extra_link_args=extra_link_args,
                    )

here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='reboundx',
    version='4.6.2',
    description='A library for including additional forces in REBOUND',
    long_description=long_description,
    url='https://github.com/dtamayo/reboundx',
    author='Daniel Tamayo',
    author_email='dtamayo@hmc.edu',
    license='GPL',
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 5 - Production/Stable',

        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',
        'Topic :: Scientific/Engineering :: Astronomy',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
    ],
    keywords='astronomy astrophysics nbody integrator',
    packages=['reboundx'],
    cmdclass={'build_ext':build_ext},
    setup_requires=['rebound>=4.6.0,<5.0.0'],
    install_requires=['rebound>=4.6.0,<5.0.0'],
    tests_require=['rebound>=4.6.0,<5.0.0','numpy'],
    test_suite="reboundx.test",
    ext_modules = [libreboundxmodule],
    zip_safe=False)
