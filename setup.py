from setuptools import setup, find_packages
import re


classifiers = [
    "Development Status :: 3 - Alpha",
    "Operating System :: POSIX",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
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
    return "See https://github.com/GangCaoLab/CoolBox"


def get_install_requires():
    requirements = []
    with open('requirements.txt') as f:
        for line in f:
            requirements.append(line.strip())
    return requirements


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
    python_requires='>=3.7, <4',
)
