from codecs import open
import os
import inspect
import sys 
from glob import glob
import sysconfig
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext

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

suffix = sysconfig.get_config_var('EXT_SUFFIX')
if suffix is None:
    suffix = ".so"

# Try to get git hash
try:
    import subprocess
    ghash = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode("ascii")
    ghash_arg = "-DREBXGITHASH="+ghash.strip()
except:
    ghash_arg = "-DREBXGITHASH=3a00b45c45de6271a322d22593e61204bd84a56c" #GITHASHAUTOUPDATE

class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)
        if "PYODIDE" in os.environ:
            return None
       
        # get site-packages dir to add to paths in case reb & rebx installed simul in tmp dir
        sitepackagesdir = sysconfig.get_path('platlib')+'/'
        rebdir, editable_rebdir = get_reb_paths(sitepackagesdir)
        rebdir_parent = os.path.dirname(rebdir)
        print("***", rebdir, "***", sitepackagesdir, "***", editable_rebdir, "***")
        ## include both as one seems required for editable installs, the other for wheels.
        self.include_dirs.append(rebdir_parent+"/src")
        self.include_dirs.append(rebdir+"/include")
        #self.include_dirs.append(editable_rebdir)
        
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
    extra_compile_args=[ghash_arg, '-D_GNU_SOURCE']
else:
    # Default compile args
    extra_compile_args=['-fstrict-aliasing', '-O3','-std=c99','-Wno-unknown-pragmas', ghash_arg, '-D_GNU_SOURCE', '-fPIC']

# Option to disable FMA in CLANG. 
FFP_CONTRACT_OFF = os.environ.get("FFP_CONTRACT_OFF", None)
if FFP_CONTRACT_OFF:
    extra_compile_args.append('-ffp-contract=off')

libreboundxmodule = Extension(
    'libreboundx',
    sources = sorted(glob("src/*.c")),
    include_dirs = ['src'],
    libraries=['rebound'+suffix[:suffix.rfind('.')]],
    extra_compile_args=extra_compile_args,
    extra_link_args=extra_link_args,
    )

setup(ext_modules=[libreboundxmodule],
    cmdclass={'build_ext':build_ext},
    )
