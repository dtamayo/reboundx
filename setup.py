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
    ghash_arg = "-DREBXGITHASH=d52096a49427f25aa6d5f6c1a3f83fe98389c559" #GITHASHAUTOUPDATE

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
        sources = [ 'src/central_force.c', 'src/core.c', 'src/exponential_migration.c', 'src/gas_damping_timescale.c', 'src/gas_dynamical_friction.c', 'src/gr.c', 'src/gr_full.c', 'src/gr_potential.c', 'src/gravitational_harmonics.c', 'src/inner_disk_edge.c', 'src/input.c', 'src/integrate_force.c', 'src/integrator_euler.c', 'src/integrator_implicit_midpoint.c', 'src/integrator_rk2.c', 'src/integrator_rk4.c', 'src/interpolation.c', 'src/lense_thirring.c', 'src/linkedlist.c', 'src/modify_mass.c', 'src/modify_orbits_direct.c', 'src/modify_orbits_forces.c', 'src/output.c', 'src/radiation_forces.c', 'src/rebxtools.c', 'src/steppers.c', 'src/stochastic_forces.c', 'src/tides_constant_time_lag.c', 'src/tides_spin.c', 'src/track_min_distance.c', 'src/type_I_migration.c', 'src/yarkovsky_effect.c'],
        
        self.library_dirs.append(rebdir+'/../')
        self.library_dirs.append(sitepackagesdir)
        for ext in self.extensions:
            ext.runtime_library_dirs.append(rebdir+'/../')
            ext.extra_link_args.append('-Wl,-rpath,'+rebdir+'/../')
            ext.runtime_library_dirs.append(sitepackagesdir)
            ext.extra_link_args.append('-Wl,-rpath,'+sitepackagesdir)
        if editable_rebdir:
            self.library_dirs.append(editable_rebdir)
            for ext in self.extensions:
                ext.runtime_library_dirs.append(editable_rebdir)
                ext.extra_link_args.append('-Wl,-rpath,'+editable_rebdir)


extra_link_args=[]
if sys.platform == 'darwin':
    config_vars = sysconfig.get_config_vars()
    config_vars['LDSHARED'] = config_vars['LDSHARED'].replace('-bundle', '-shared')
    extra_link_args.append('-Wl,-install_name,@rpath/libreboundx'+suffix)
if sys.platform == 'win32':
    extra_compile_args=[ghash_arg, '-DLIBREBOUNDX', '-D_GNU_SOURCE']
else:
    # Default compile args
    extra_compile_args=['-fstrict-aliasing', '-O3','-std=c99','-Wno-unknown-pragmas', ghash_arg, '-DLIBREBOUNDX', '-D_GNU_SOURCE', '-fPIC']

# Option to disable FMA in CLANG. 
FFP_CONTRACT_OFF = os.environ.get("FFP_CONTRACT_OFF", None)
if FFP_CONTRACT_OFF:
    extra_compile_args.append('-ffp-contract=off')

libreboundxmodule = Extension('libreboundx',
        sources = [ 'src/central_force.c', 'src/core.c', 'src/exponential_migration.c', 'src/gas_damping_timescale.c', 'src/gas_dynamical_friction.c', 'src/gr.c', 'src/gr_full.c', 'src/gr_potential.c', 'src/gravitational_harmonics.c', 'src/inner_disk_edge.c', 'src/input.c', 'src/integrate_force.c', 'src/integrator_euler.c', 'src/integrator_implicit_midpoint.c', 'src/integrator_rk2.c', 'src/integrator_rk4.c', 'src/interpolation.c', 'src/lense_thirring.c', 'src/linkedlist.c', 'src/modify_mass.c', 'src/modify_orbits_direct.c', 'src/modify_orbits_forces.c', 'src/output.c', 'src/radiation_forces.c', 'src/rebxtools.c', 'src/steppers.c', 'src/stochastic_forces.c', 'src/tides_constant_time_lag.c', 'src/tides_spin.c', 'src/track_min_distance.c', 'src/type_I_migration.c', 'src/yarkovsky_effect.c'],
                    include_dirs = ['src'],
                    library_dirs = [],
                    runtime_library_dirs = ["."],
                    libraries=['rebound'+suffix[:suffix.rfind('.')]],
                    define_macros=[ ('LIBREBOUNDX', None) ],
                    extra_compile_args=extra_compile_args,
                    extra_link_args=extra_link_args,
                    )

here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='reboundx',
    version='4.4.2',
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
    setup_requires=['rebound>=4.0.0'],
    install_requires=['rebound>=4.0.0'],
    tests_require=['rebound>=4.0.0','numpy'],
    test_suite="reboundx.test",
    ext_modules = [libreboundxmodule],
    zip_safe=False)
