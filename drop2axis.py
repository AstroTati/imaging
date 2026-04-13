#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 08:29:01 2022

@author: trodriguez
"""
import astropy
from astropy.coordinates import Angle
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.wcs import WCS


def drop2axis(filename, outname):
    hdu = astropy.io.fits.open(filename)[0]
    try:    
        for kw in 'CTYPE', 'CRVAL','CRPIX','CDELT', 'CUNIT','NAXIS':
            for n in 3,4:
                hdu.header.remove(f"{kw}{n}")
        for kw in 'PC3_1','PC3_2','PC3_3','PC3_4','PC4_1','PC4_2','PC4_3','PC4_4',:
            hdu.header.remove(f"{kw}")
        for kw in 'PC1_3','PC1_4','PC2_3','PC2_4':
            hdu.header.remove(f"{kw}")
        astropy.io.fits.writeto(outname, hdu.data[0,0], hdu.header,overwrite=True)
    except KeyError:
        pass

def drop2axis_wzeroes(filename, outname):
    hdu = astropy.io.fits.open(filename)[0]
    try:
        for kw in 'CTYPE', 'CRVAL','CRPIX','CDELT', 'CUNIT','NAXIS':
            for n in 3,4:
                hdu.header.remove(f"{kw}{n}")
        for kw in 'PC03_01','PC03_02','PC03_03','PC03_04','PC04_01','PC04_02','PC04_03','PC04_04',:
            hdu.header.remove(f"{kw}")
        for kw in 'PC01_03','PC01_04','PC02_03','PC02_04':
            hdu.header.remove(f"{kw}")
        astropy.io.fits.writeto(outname, hdu.data[0,0], hdu.header,overwrite=True)
    except KeyError:
        pass


# def drop2axis_wzeroes(filename, outname):
    # hdu = astropy.io.fits.open(filename)[0]
    # for kw in 'CTYPE', 'CRVAL','CRPIX','CDELT', 'CUNIT','NAXIS':
    #     for n in 3,4:
    #         hdu.header.remove(f"{kw}{n}")
    # for kw in 'PC003001','PC003003','PC003004','PC004001','PC004002','PC004003','PC004004':
    #     hdu.header.remove(f"{kw}")
    # for kw in 'PC001003','PC001004','PC002003','PC002004':
    #     hdu.header.remove(f"{kw}")
    # for kw in 'PC003002','PC001003':
    #     hdu.header.remove(f"{kw}")
    # astropy.io.fits.writeto(outname, hdu.data[0,0], hdu.header,overwrite=True)
    # except KeyError:
    #     print('error!',KeyError)
    #     pass
   
    # 'PC00302',
    

def load(filename):
    file = get_pkg_data_filename(filename)
    hdu = fits.open(file)[0]
    # print(hdu.header)
    data = hdu.data
    wcs = WCS(hdu.header)
    
    return(data, wcs, hdu)


import os
files_path = '/home/tatush/Desktop/Projects/VLA-G023'

for filename in os.listdir(files_path):
    if filename.endswith("wmask.fits"):
        file = "{}/{}".format(files_path,filename)
        drop2axis_wzeroes(file,file)
        drop2axis(file,file)
        print(file)
    
    
# drop2axis_wzeroes('/home/tatush/Desktop/VLA22A092/Final_cont_fits/18151_rosC.fits',
            # '/home/tatush/Desktop/VLA22A092/Final_cont_fits/18151_rosC.fits')

# load('/home/trodriguez/Desktop/pyscripts+txt/VLA22A-092/fits_files/19035Cband-Aarray-full.fits')



# 18470_rosK.fits
# 18151_rosC.fits
# 18182_rosC.fits
# 18566_rosC.fits
# 19012_rosC.fits


# import numpy
# import pyfits

# filename = '/home/tatush/Desktop/VLA22A092/Final_cont_fits/18151_rosC.fits'
# outfile  = '/home/tatush/Desktop/VLA22A092/Final_cont_fits/18151_rosC_test.fits'

# # Get the shape of the file
# fitsfile=pyfits.open(filename)
# image = fitsfile[0].data
# header =fitsfile[0].header
# z = image.shape[0]  # No. channels                  
# y = image.shape[1]  # No. x pixels
# # x = image.shape[3]  # No. y pixels

# newimage = numpy.reshape(image,[z,y])

# pyfits.core.writeto(outfile,newimage,header, clobber=True)













