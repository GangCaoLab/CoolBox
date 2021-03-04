from setuptools import setup, find_packages
from setuptools import Extension
from setuptools.command.build_ext import build_ext
import setuptools
import sys
import re


classifiers = [
    "Development Status :: 3 - Alpha",
    "Operating System :: POSIX",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    "Intended Audience :: Science/Research",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "Topic :: Scientific/Engineering :: Visualization",
]


keywords = [
    'genomics',
    'bioinformatics',
    'visualization',
    'Jupyter',
]


def get_version():
    with open("coolbox/__init__.py") as f:
        for line in f.readlines():
            m = re.match("__version__ = '([^']+)'", line)
            if m:
                ver = m.group(1)
                return ver
        raise IOError("Version information can not found.")


def get_long_description():
    #with open("README.md") as f:
    #    desc = f.read()
    #return desc
    return "https://github.com/GangCaoLab/CoolBox"


def get_install_requires():
    requirements = []
    with open('requirements.txt') as f:
        for line in f:
            requirements.append(line.strip())
    return requirements

### For compile strawC

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)


ext_modules = [
    Extension(
        'coolbox.utilities.hic.strawC',
        ['coolbox/utilities/hic/strawC_src/straw.cpp'],
        include_dirs=[
            # Path to pybind11 headers
            get_pybind_include(),
            get_pybind_include(user=True)
        ],
        language='c++'
    ),
]


# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def cpp_flag(compiler):
    """Return the -std=c++[11/14/17] compiler flag.
    The newer version is prefered over c++11 (when it is available).
    """
    flags = ['-std=c++14', '-std=c++11']

    for flag in flags:
        if has_flag(compiler, flag): return flag

    raise RuntimeError('Unsupported compiler -- at least C++11 support '
                       'is needed!')


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }
    l_opts = {
        'msvc': [],
        'unix': ['-lcurl', '-lz'],
    }

    if sys.platform == 'darwin':
        darwin_opts = ['-stdlib=libc++', '-mmacosx-version-min=10.7']
        c_opts['unix'] += darwin_opts
        l_opts['unix'] += darwin_opts

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        link_opts = self.l_opts.get(ct, [])
        if ct == 'unix':
            opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args = opts
            ext.extra_link_args = link_opts
        build_ext.build_extensions(self)


###

setup(
    name='coolbox',
    author='Weize Xu',
    author_email='vet.xwz@gmail.com',
    version=get_version(),
    license='GPLv3',
    description='Jupyter notebook based genomic data visulization toolkit.',
    long_description=get_long_description(),
    keywords=keywords,
    url='https://github.com/GangCaoLab/CoolBox',
    packages=find_packages(),
    scripts=['scripts/coolbox'],
    include_package_data=True,
    zip_safe=False,
    classifiers=classifiers,
    install_requires=get_install_requires(),
    python_requires='>=3.6, <4',
    setup_requires=['pybind11>=2.4'],
    ext_modules=ext_modules,
    cmdclass={'build_ext': BuildExt},
)
