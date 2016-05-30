from codecs import open
import os
import inspect
import sys 
import sysconfig

try:
    from setuptools import setup, Extension
    from setuptools.command.build_ext import build_ext as _build_ext
except ImportError:
    print("Installing REBOUNDx requires setuptools.  Do 'pip install setuptools'.")
    sys.exit(1)

suffix = sysconfig.get_config_var('EXT_SUFFIX')
if suffix is None:
    suffix = ".so"

class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)

        try:
            import rebound
        except ImportError:
            print("REBOUNDx did not automatically install REBOUND.  Please let me know if this happens (tamayo.daniel@gmail.com), and try first installing REBOUND (http://rebound.readthedocs.org/en/latest/python_quickstart.html")
            sys.exit(1)
        try:
            version = rebound.__version__ # Added in 2.12.1
        except AttributeError:
            print("REBOUNDx did not automatically install a recent enough version of REBOUND.  Please let me know if this happens (tamayo.daniel@gmail.com), and try upgrading REBOUND.  See 5.3 in http://rebound.readthedocs.org/en/latest/python_quickstart.html")
            sys.exit(1)

        rebdir = os.path.dirname(inspect.getfile(rebound))
        self.include_dirs.append(rebdir)
        self.library_dirs.append(rebdir+'/../')
        for ext in self.extensions:
            ext.runtime_library_dirs.append(rebdir+'/../')
            ext.extra_link_args.append('-Wl,-rpath,'+rebdir+'/../')

from distutils.version import LooseVersion

extra_link_args=[]
if sys.platform == 'darwin':
    from distutils import sysconfig
    vars = sysconfig.get_config_vars()
    vars['LDSHARED'] = vars['LDSHARED'].replace('-bundle', '-shared')
    extra_link_args.append('-Wl,-install_name,@rpath/libreboundx'+suffix)

libreboundxmodule = Extension('libreboundx',
                    sources = [ 'src/core.c', 'src/gr.c', 'src/gr_full.c', 'src/gr_potential.c', 'src/modify_mass.c', 'src/modify_orbits_direct.c', 'src/modify_orbits_forces.c', 'src/radiation_forces.c', 'src/rebxtools.c'],
                    include_dirs = ['src'],
                    library_dirs = [],
                    runtime_library_dirs = ["."],
                    libraries=['rebound'+suffix[:-3]], #take off .so from the suffix
                    define_macros=[ ('LIBREBOUNDX', None) ],
                    extra_compile_args=['-fstrict-aliasing', '-O3','-std=c99','-march=native', '-fPIC', '-Wpointer-arith'],
                    extra_link_args=extra_link_args,
                    )

here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='reboundx',
    version='2.8.4',
    description='A library for including additional forces in REBOUND',
    long_description=long_description,
    url='http://github.com/dtamayo/reboundx',
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
    setup_requires=['rebound>=2.13.6'],
    install_requires=['rebound>=2.13.6'],
    tests_require=["numpy","matplotlib"],
    test_suite="reboundx.test",
    ext_modules = [libreboundxmodule],
    zip_safe=False)
