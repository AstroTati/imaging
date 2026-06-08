#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ONE MAY NEED TO HAVE ONLY 2 AXIS IN THEIR FITS
FILE TO USE AT LEAST SOME ASTROPY PACKAGES.
THIS SCRIPT GETS RID OF 2 AXIS IN VLA DATA HEADER.

YOU CAN EITHER GIVE IT ONLY 1 FILE OR A PATH AND 
RUN FOR ALL FILES INSIDE A DIRECTORY.

INPUT:
    - original fits file or path to files.
OUTPUT
    - fits file(s) w/ only 2 axis.


Created on Thu Jan 27 08:29:01 2022
@author: trodriguez
"""

import astropy
from astropy.coordinates import Angle
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.wcs import WCS
import os

def drop2axis(filename, outname, zeroes=False):
    fmt = "PC{:02d}_{:02d}" if zeroes else "PC{}_{}"
    # In older header versions there are zeroes. You can check your header and check if
    # zeroes should be true or false.

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

####### INPUTS AND RUN
files_path = '/path/to/fits_files'

for filename in os.listdir(files_path):
    if filename.endswith(".fits"):
        file = "{}/{}".format(files_path,filename)
        drop2axis(file,file,zeroes = False)   # This is because I want it to overwrite. 
                                            # One could also define a file_output name ofc.

### another option: a single file
# file = '/path/to/file.fits'
# drop2axis(file,file, zeroes = False)










