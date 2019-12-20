# showspectra
Display and fit optical spectra from multi-object facilities

The code allows one to display spectra from multi-object instruments (fibers, MOS).
For the moment, it is limited to WIYN and VIMOS data which have to be formatted in three files containing flux, errors, and sky. 

The code is written in Python 3 and it is available on anaconda.

## Installation

The current release requires the installation of [anaconda with Python 3.7](https://www.anaconda.com/distribution):

Once anaconda (Python 3.7) is installed, check if you have 'limfit' installed:

conda install -c astropy limfit

and then install showspectra:

conda install -c darioflute showspectra

## Usage

Start the code by writing:

showspectra

You will be asked to select the instrument, and then the files containing fluxes, errors, and sky.

## Reading the analysis

The analysis is saved in JSON files, which are plain ASCII files.
An utility to read them in available in the showspectra package.
You can read some examples in the notebook: [Read JSON output](notebooks/ReadLines.ipynb)
