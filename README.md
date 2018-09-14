# showspectra
Display and fit optical spectra from multi-object facilities

The code allows one to display spectra from multi-object instruments (fibers, MOS).
For the moment, it is limited to WIYN and VIMOS data which have to be formatted in three files containing flux, errors, and sky. 

The code is written in Python 3 and it is available on anaconda.

## Installation

Once anaconda (Python 3.x) is installed, do:

conda install -c darioflute showspectra

## Usage

Start the code by writing:

showspectra

You will be asked to select the instrument, and then the files containing fluxes, errors, and sky.
