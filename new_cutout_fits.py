#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 11:10:49 2023

@author: tatush
"""

from astropy.io import fits
import numpy as np

def cut_fits(input_file, output_file, center_x, center_y, size=350):
    # Open the input FITS file
    with fits.open(input_file) as hdul:
        # Get the data and header from the original file
        data = hdul[0].data
        header = hdul[0].header
        
        print(data.shape)
        # Determine the starting and ending positions for the cutout
        start_x = max(center_x - size // 2, 0)
        end_x = min(start_x + size, header['NAXIS1'])
        start_y = max(center_y - size // 2, 0)
        end_y = min(start_y + size, header['NAXIS2'])
        
        # Create the cutout by slicing the data
        cutout = data[start_y:end_y, start_x:end_x]
        
        # Update header information for the cutout
        header['NAXIS1'] = end_x - start_x
        header['NAXIS2'] = end_y - start_y
        header['CRPIX1'] -= start_x
        header['CRPIX2'] -= start_y
        
        # Write the cutout data and updated header to a new FITS file
        fits.writeto(output_file, cutout, header, overwrite=True)
        
        
path = '/home/tatush/Desktop/Projects/VLA-22A092/Final_cont_fits/'
input_file_name = path + '18345_spw0_15.pbcor.fits'
output_file_name = path + '18345_spw0-15_rosK.fits'


center_x_position = 1218  # Replace with desired X center position
center_y_position = 1284  # Replace with desired Y center position

cut_fits(input_file_name, output_file_name, center_x_position, center_y_position)


# # 18470_rosK.fits
# 885, 1297