"""A setuptools based setup module.
"""
from setuptools import setup, find_packages

from os import path

from . import omin

# here = path.abspath(path.dirname(__file__))
__version__ = "0.0.19"

setup(
    name='omin',
    version=omin.__version__,
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
	package_data = {'omin': ['databases/*.txt',
	                         'databases/*.p',
							 'databases/*.pickle',
	                         'databases/*.xls',
							 'databases/*.xlsx',
							 'databases/mitocarta/*.ipynb',
							 'databases/mitocarta/*.p',
							 'databases/mitocarta/*.pickle',
							 'databases/mitocarta/*.xls',
							 'databases/mitocarta/*.xlsx',
							 'databases/mitocarta/*.gz'],},
    install_requires=['pandas',
                      'xlrd',
                      'numpy',
                      'guipyter',
                      'dominate',
                      'matplotlib_venn'],
    entry_points = {
        'console_scripts': [
            'omin=omin.cli.cli:main',
            ],
        },
)
