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
	package_data = {'omin': ['__version__', 'databases/mitocarta/*.p']},
    install_requires=['pandas',
                      'xlrd',# Needed for pandas export of DataFrames to xlsx
                      'numpy',
                      'scipy',
                      'dominate',
                      'statsmodels',
                      'intermine',
                      'lxml',
                      'beautifulsoup4',
                      'matplotlib_venn',
                      'guipyter',
                      'pandomics'],
    entry_points = {
        'console_scripts': [
            'omin=omin.cli.cli:main',
            ],
        },
)
