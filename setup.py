from codecs import open
import os
import inspect
import sys 
from distutils import sysconfig
from distutils.sysconfig import get_python_lib 

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
    ghash_arg = "-DREBXGITHASH="+ghash
except:
    ghash_arg = "-DREBXGITHASH=86015f5b20b1fa1a6622d3ca6aeef10b56f1e032" #GITHASHAUTOUPDATE

class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)
        if "PYODIDE" in os.environ:
            return None

        try:
            import rebound
        except ImportError:
            print("REBOUNDx did not automatically install REBOUND.  Please let me know if this happens (tamayo.daniel@gmail.com), and try first installing REBOUND (https://rebound.readthedocs.org/en/latest/python_quickstart.html")
            sys.exit(1)
        try:
            version = rebound.__version__ # Added in 2.12.1
        except AttributeError:
            print("REBOUNDx did not automatically install a recent enough version of REBOUND.  Please let me know if this happens (tamayo.daniel@gmail.com), and try upgrading REBOUND.  See 5.3 in https://rebound.readthedocs.org/en/latest/python_quickstart.html")
            sys.exit(1)

        rebdir = os.path.dirname(inspect.getfile(rebound))
        # get site-packages dir to add to paths in case reb & rebx installed simul in tmp dir
        rebdirsp = get_python_lib()+'/'#[p for p in sys.path if p.endswith('site-packages')][0]+'/'
        self.include_dirs.append(rebdir)
        sources = [ 'src/modify_mass.c', 'src/integrator_euler.c', 'src/modify_orbits_forces.c', 'src/integrator_rk2.c', 'src/track_min_distance.c', 'src/tides_spin.c', 'src/rebxtools.c', 'src/inner_disk_edge.c', 'src/gravitational_harmonics.c', 'src/gr_potential.c', 'src/core.c', 'src/integrator_rk4.c', 'src/input.c', 'src/central_force.c', 'src/stochastic_forces.c', 'src/gr.c', 'src/modify_orbits_direct.c', 'src/tides_constant_time_lag.c', 'src/yarkovsky_effect.c', 'src/gr_full.c', 'src/steppers.c', 'src/integrate_force.c', 'src/interpolation.c', 'src/type_I_migration.c', 'src/output.c', 'src/radiation_forces.c', 'src/integrator_implicit_midpoint.c', 'src/exponential_migration.c', 'src/linkedlist.c'],
        
        self.library_dirs.append(rebdir+'/../')
        self.library_dirs.append(rebdirsp)
        for ext in self.extensions:
            ext.runtime_library_dirs.append(rebdir+'/../')
            ext.extra_link_args.append('-Wl,-rpath,'+rebdir+'/../')
            ext.runtime_library_dirs.append(rebdirsp)
            ext.extra_link_args.append('-Wl,-rpath,'+rebdirsp)

from distutils.version import LooseVersion

extra_link_args=[]
if sys.platform == 'darwin':
    from distutils import sysconfig
    vars = sysconfig.get_config_vars()
    vars['LDSHARED'] = vars['LDSHARED'].replace('-bundle', '-shared')
    extra_link_args.append('-Wl,-install_name,@rpath/libreboundx'+suffix)

libreboundxmodule = Extension('libreboundx',
        sources = [ 'src/modify_mass.c', 'src/integrator_euler.c', 'src/modify_orbits_forces.c', 'src/integrator_rk2.c', 'src/track_min_distance.c', 'src/tides_spin.c', 'src/rebxtools.c', 'src/inner_disk_edge.c', 'src/gravitational_harmonics.c', 'src/gr_potential.c', 'src/core.c', 'src/integrator_rk4.c', 'src/input.c', 'src/central_force.c', 'src/stochastic_forces.c', 'src/gr.c', 'src/modify_orbits_direct.c', 'src/tides_constant_time_lag.c', 'src/yarkovsky_effect.c', 'src/gr_full.c', 'src/steppers.c', 'src/integrate_force.c', 'src/interpolation.c', 'src/type_I_migration.c', 'src/output.c', 'src/radiation_forces.c', 'src/integrator_implicit_midpoint.c', 'src/exponential_migration.c', 'src/linkedlist.c'],
                    include_dirs = ['src'],
                    library_dirs = [],
                    runtime_library_dirs = ["."],
                    libraries=['rebound'+suffix[:suffix.rfind('.')]],
                    define_macros=[ ('LIBREBOUNDX', None) ],
                    extra_compile_args=['-fstrict-aliasing', '-D_GNU_SOURCE', '-O3','-std=c99', '-fPIC', '-Wpointer-arith', ghash_arg],
                    extra_link_args=extra_link_args,
                    )

here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='reboundx',
    version='3.9.1',
    description='A library for including additional forces in REBOUND',
    long_description=long_description,
    url='https://github.com/dtamayo/reboundx',
    author='Daniel Tamayo',
    author_email='tamayo.daniel@gmail.com',
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
    setup_requires=['rebound>=3.23.0', 'numpy'],
    install_requires=['rebound>=3.23.0', 'numpy'],
    tests_require=["numpy","matplotlib"],
    test_suite="reboundx.test",
    ext_modules = [libreboundxmodule],
    zip_safe=False)
