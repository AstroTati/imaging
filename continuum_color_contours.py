#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 13:26:56 2023

@author: tatush

=== 
This script is to plot a radio continuum image using VLA data.
It will run only if you have 2 axis (check drop2axis.py)


I will plot one data set in contours (data_1) and one data set in both contours (data_2). 
The data have different angular resolution, the latter having a smaller beam size. 
===

"""

from astropy.wcs import WCS
import astropy
from astropy.io import fits
from astropy.coordinates import Angle
from astropy.utils.data import get_pkg_data_filename
import numpy as np
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.nddata import Cutout2D
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator

## I don't like warnings :3
import warnings
warnings.filterwarnings("ignore")

## ====================================================================
## ===================FUNCTION DEFINITION==============================

## load function
def load(filename):
    file = get_pkg_data_filename(filename)
    hdu = fits.open(file)[0]
    data = hdu.data
    wcs = WCS(hdu.header)
    
    return(data, wcs, hdu)


def find_center_and_scale(hdr):
    obsra, obsdec = hdr["OBSRA"], hdr["OBSDEC"]
    scale_ra, scale_dec = hdr["CDELT1"], hdr["CDELT2"] # deg/pixel
    npxl_ra, npxl_dec = hdr["NAXIS1"], hdr["NAXIS2"]     
    
    return obsra, obsdec, scale_ra, scale_dec, npxl_ra, npxl_dec  # in deg

def cut(data, wcs, hdu, center, box):
    ### I want to plot only a small area around my source, not the whole image.
    size = u.Quantity(box, u.pix) # units could be also pixels
    cutout = Cutout2D(data,position=center,size=size,wcs=wcs)
    hdu.header.update(cutout.wcs.to_header())
    hdu.data = cutout.data
    wcs_cut = WCS(hdu.header)
    
    return(hdu.data,wcs_cut,hdu)

def find_beam(hdr):
    ### Finds the beam info in the header
    major = hdr["BMAJ"] # in degrees from VLA
    minor = hdr["BMIN"]
    pa = hdr["BPA"] # position angle in degrees from VLA
    
    return (major_data_2 minor_data_2, pa)

## ====================================================================
## ======================DATA DEFINITION===============================

file_1 = '/the-path/continuum_contours.fits' ## I will use only the contours of these data
data_1 , wcs_1 , hdu_1 = load(file_1)
data_1 = np.ma.masked_invalid(data_1) # Mask blanks/NaN
data_1_deg = SkyCoord(286.5067083333333,6.776722222222222, frame='fk5', unit='deg')
data_1_pix = wcs_1.world_to_pixel(data_1_deg)
# print(data_1_pix) # check
data_1 , wcs_1 , hdu_1 = cut(data_1 , wcs_1 , hdu_1, (1048,1224), box=(55,55))  ## This is the center of the area I want to plot
                                                                                ## It is not necessary the center of the image,  
                                                                                ## hence I give it the precise pixel. Box size in pix.


data_2 , wcs_2 , hdu_2 = load('/the-path/continuum_image_color.fits')
data_2 = np.ma.masked_invalid(data_2) # Mask blanks/NaN
rms=data_2.std()
center_deg = SkyCoord(286.5067083333333,6.776722222222222, frame='fk5', unit='deg')
center_pix = wcs_2.world_to_pixel(center_deg)
# print(center_pix)
data_2 , wcs_2 , hdu_2 = cut(data_2 , wcs_2 , hdu_2, center=(935,1020),box=(300,300))


## ====================================================================
## ======================PLOTTING======================================

fig = plt.figure(figsize=(15,15))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=wcs_1) # The lower res data determines the wcs
    

## axis label and ticks
ra = ax.coords[0]
dec = ax.coords[1]
ra.set_axislabel("RA (J2000)", minpad=0.8, fontsize=25)
dec.set_axislabel("Dec (J2000)", minpad=-1.0, fontsize=25)
ra.set_major_formatter('hh:mm:ss.ss')
ax.set_title('source name',fontsize=25)
# I'll plot only 40 pixels
ax.set_xlim(7,47)
ax.set_ylim(5,45)
    

ra.display_minor_ticks(True)
dec.display_minor_ticks(True)    
ax.tick_params(which='both',direction='in',color='white',length=10,width=2, labelsize=20)
ax.tick_params(which='minor', length=5)

## Color plot
im=ax.imshow(data_1, vmin=-1e-7, vmax=data_1.max(), origin='lower', cmap='plasma', transform=ax.get_transform(wcs_1))

## Color bar
cbar = fig.colorbar(im,shrink=0.8)
for t in cbar.ax.get_yticklabels():
    t.set_fontsize(18)
cbar.ax.set_ylabel(r'$\mu$Jy/beam', fontsize=20)
cbar.ax.set_yticklabels(['0','50','100','150','200','250','300','350','400'])

## Contour plot
plt.rcParams["lines.linewidth"] = 2
levels = [rms*-3,rms*3,rms*5,rms*10,rms*15]      
ax.contour(data_2, levels=levels, colors='k', transform=ax.get_transform(wcs_2))



## Scale bar
d = 2.3 # kpc
sb_length = (500) /(d * 1000)     # I want a 500 au scale bar
pix_scale = 0.05                  # pixel size
bsize = sb_length / pix_scale     # size of bar in pix

xlim1, xlim2 = ax.get_xlim()
ylim1, ylim2 = ax.get_ylim()

aux = (xlim2-2-bsize+xlim2-2)/2.
ax.plot([xlim2-2-bsize,xlim2-2], [ylim2-1,ylim2-1], linewidth=3, c='w', alpha=1.0)

scale_bar_fontsize=16
scale_bar_text = '500 au'

ax.text(aux, ylim2-2.5, scale_bar_text, fontsize=22, horizontalalignment='center', color='w')

## Plot beam
obsra, obsdec, scale_ra, scale_dec, npxl_ra, npxl_dec = find_center_and_scale(hdu_1.header)
major_data_2, minor_data_2, pa_data_2 = find_beam(hdu_2.header)

plt.rcParams["hatch.linewidth"] = 3
obsra, obsdec, scale_ra, scale_dec, npxl_ra, npxl_dec = find_center_and_scale(hdu_1.header)
major_data_2, minor_data_2, pa = find_beam(hdu_2.header)
# beam_ra  = obsra  - 1 *scale_ra*npxl_ra + major_data_2
# beam_dec = obsdec - 1 *scale_dec*npxl_dec + major_data_2
# I was having issues with the method above so I forced the beam position manually.
beam_ra = Angle('19h06m01.663s').deg
beam_dec = Angle('6d46m35.35s').deg

major_data_1, minor_data_1, PA_data_1 = find_beam(hdu_1.header)
beam_data_1 = Ellipse((beam_ra, beam_dec), major_data_1, minor_data_1, 90.0-PA_data_1, # BPA relative pos
                    transform=ax.get_transform('fk5'),
                    edgecolor='w', facecolor='none',lw=3)
ax.add_patch(beam_data_1)

beam_data_2 = Ellipse((beam_ra, beam_dec), major_data_2, minor_data_2, 90.0-pa_data_2, # BPA relative pos
                    transform=ax.get_transform('fk5'),
                    edgecolor='w', facecolor='w')
ax.add_patch(beam_data_2)


## optional for saving plot
# plot_name = 'plot_name.pdf' # or .png
# plt.savefig(plot_name, bbox_inches='tight') 
