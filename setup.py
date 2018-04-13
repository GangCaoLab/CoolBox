from setuptools import setup, find_packages

import  re


classifiers = [
    "Development Status :: 3 - Alpha",
    "Operating System :: OS Independent",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.4",
    "Programming Language :: Python :: 3.5",
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
    with open("README.md") as f:
        desc = f.read()
    return desc


def get_install_requires():
    requirements = [
        "numpy >= 1.12.*",
        "matplotlib >= 2.0.0",
        "pandas >= 0.22.0",
        "intervaltree >= 2.1.0",
        "pybigwig >= 0.3.7",
        "cooler >= 0.7.6",
        "ipywidgets >= 7.1.2",
    ]
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
    url='https://github.com/Nanguage/CoolBox',
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    classifiers=classifiers,
    install_requires=get_install_requires(),
    python_requires='>=3.4, <4',
)