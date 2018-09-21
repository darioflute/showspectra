#!/usr/bin/env python
"""Setup script fo installing showspectra."""

from setuptools import setup

config = {
    'name': 'showspectra',
    'version': '0.07',
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
                                     'stylesheet.css', 'Templates/*.fit.gz']}
}

setup(**config)
