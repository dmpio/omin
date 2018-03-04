"""A setuptools based setup module.
"""
from setuptools import setup, find_packages

import os
here = os.getcwd()

with open(os.path.join(here, 'omin', '__version__')) as f:
    __version__ = f.read().strip()

setup(
    name='omin',
    version=__version__,
    description='Tools for omics analysis',
    url='https://github.com/draperjames/skunkworks',
    author='James Draper',
    author_email='james.draper@duke.edu',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Topic :: Software Development :: Build Tools',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    keywords='pandas',
    packages=find_packages(),
	package_data = {'omin': ['__version__',
                             'databases/mitocarta/*.ipynb',# FIXME : Remove all but pickle or bz.
                             'databases/mitocarta/*.p',
                             'databases/mitocarta/*.pickle',
                             'databases/mitocarta/*.xls',
                             'databases/mitocarta/*.xlsx',
                             'databases/mitocarta/*.gz'],},
    install_requires=['pandas',
                      'xlrd',# Needed for pandas export of DataFrames to xlsx
                      'numpy',
                      'guipyter',
                      'dominate',
                      'pandomics',
                      'matplotlib_venn'],
    entry_points = {
        'console_scripts': [
            'omin=omin.cli.cli:main',
            ],
        },
)
