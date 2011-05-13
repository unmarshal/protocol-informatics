from distutils.core import setup
from distutils.extension import Extension

import sys, os.path

# Use Pyrex
#import distutils.sysconfig
#if distutils.sysconfig.get_config_var('CC').startswith("gcc"):
#    pyrex_compile_options = []
#else:
pyrex_compile_options = []

if sys.platform == "win32" and len(sys.argv) < 2:
    sys.argv[1:] = ["bdist_wininst"]

# Compiling Pyrex modules to .c and .so
try:
    import Pyrex.Distutils
except ImportError:
    distutils_extras = {}
    pyrex_suffix = ".c"
else:
    class pyrex_build_ext(Pyrex.Distutils.build_ext):
        def pyrex_compile(self, source):
            from Pyrex.Compiler.Main import CompilationOptions, default_options
            options = CompilationOptions(default_options)
            result = Pyrex.Compiler.Main.compile(source, options)
            if result.num_errors <> 0:
                sys.exit(1)
    distutils_extras = {
        "cmdclass": {
            'build_ext': pyrex_build_ext}}
    pyrex_suffix = ".pyx"

def PIExtension(module_name):
    path = module_name.replace('.', '/')
    return Extension(module_name, [path + pyrex_suffix],
            extra_compile_args = pyrex_compile_options)

setup(
    name = "PI",
    version = "0.01",
    url = "http://www.baselineresearch.net/PI",
    author_email = "mbeddoe@baselineresearch.net",
    description = "Protocol analysis toolkit using bioinformatics algorithms",
    packages = ["PI"],
    ext_modules=[
        PIExtension("PI.align")],**distutils_extras
)
