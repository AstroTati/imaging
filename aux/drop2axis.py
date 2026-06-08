"""
Created on Thu Jan 27 08:29:01 2022
Tatiana M. Rodriguez

============================================================================
drop2axis.py
============================================================================
Drops 2 axis in a VLA continuum fits file. This is needed to make it 
compatible with astropy plotting packages. 

Input:
    path to file(s).
Output:
    fits file(s) with only 2 axis in their header. Overwrites, so be careful!
    
Note:
    I used copilot to optimize this script I wrote a few years back.
    I tested it on my data, but let me know if you encounter any issues.

"""

import astropy
from astropy.coordinates import Angle
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.wcs import WCS
import os

def drop2axis(filename, outname, zeroes=False):
    fmt = "PC{:02d}_{:02d}" if zeroes else "PC{}_{}"
    # In older header versions there are zeroes. 
    # Check your header to decided whether zeroes should be true or false.

    with astropy.io.fits.open(filename) as hdul:
        hdu = hdul[0]
        header = hdu.header.copy()

        keywords = (
            [f"{kw}{n}" for kw in ('CTYPE','CRVAL','CRPIX','CDELT','CUNIT','NAXIS') for n in (3, 4)]
            + [fmt.format(i, j) for i in (3, 4) for j in (1, 2, 3, 4)]
            + [fmt.format(i, j) for i in (1, 2) for j in (3, 4)]
        )

        for kw in keywords:
            header.remove(kw, ignore_missing=True)

        astropy.io.fits.writeto(outname, hdu.data[0, 0], header, overwrite=True)

def load(filename):
    file = get_pkg_data_filename(filename)
    hdu = fits.open(file)[0]
    # print(hdu.header) # for checks
    data = hdu.data
    wcs = WCS(hdu.header)
    
    return(data, wcs, hdu)

# ============================================================================
# Inputs and run
# ============================================================================

files_path = '/path/to/fits_files'

for filename in os.listdir(files_path):
    if filename.endswith(".fits"):
        file = "{}/{}".format(files_path,filename)
        drop2axis(file,file,zeroes = False)   # Modify outpout file path if you don't want it to overwrite


### another option: a single file
# file = '/path/to/file.fits'
# drop2axis(file,file, zeroes = False)










