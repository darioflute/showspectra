#!/usr/bin/env python
"""Setup script fo installing showspectra."""
from distutils.core import setup
import json

with open('sospex/version.json') as fp:
    _info = json.load(fp)

with open("README.md", "r") as fh:
    long_description = fh.read()

config = {
    'name': 'showspectra',
    'version': _info['version'],
    'description': 'Show Spectra',
    'long_description': 'The package shows spectra and fits lines',
    'author': 'Dario Fadda',
    'author_email': 'darioflute@gmail.com',
    'url': 'https://github.com/darioflute/showspectra.git',
    'download_url': 'https://github.com/darioflute/showspectra',
    'license': 'GPLv3+',
    'packages': ['showspectra'],
    'scripts': ['bin/showspectra'],
    'include_package_data': True,
    'package_data': {'showspectra': ['icons/*.png', 'icons/*.gif',
                                     'stylesheet.css', 'Templates/*.fit.gz', 'version.json']},
    'classifiers':[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: GPLv3+ License",
            "Operating System :: OS Independent",
            "Intended Audience :: Science/Research", 
            "Development Status :: 4 - Beta",
            "Topic :: Scientific/Engineering :: Visualization",
            ],
}

setup(**config)
