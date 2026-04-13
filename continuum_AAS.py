#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 13:26:56 2023

@author: tatush
"""

from astropy.wcs import WCS
import astropy
from astropy.io import fits
from astropy.coordinates import Angle
from astropy.utils.data import get_pkg_data_filename
# from regions import CirclePixelRegion, PixCoord
import numpy as np
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.nddata import Cutout2D
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator

import warnings
warnings.filterwarnings("ignore")

#%% load function
def load(filename):
    file = get_pkg_data_filename(filename)
    hdu = fits.open(file)[0]

    data = hdu.data
    wcs = WCS(hdu.header)
    
    return(data, wcs, hdu)


def find_center_and_scale(hdr):
    obsra, obsdec = hdr["OBSRA"], hdr["OBSDEC"]
    scale_ra, scale_dec = hdr["CDELT1"], hdr["CDELT2"] # deg/pixel
    # center_pixel_ra, center_pixel_dec = hdr["CRPIX1"], hdr["CRPIX2"]
    npxl_ra, npxl_dec = hdr["NAXIS1"], hdr["NAXIS2"]     
    return obsra, obsdec, scale_ra, scale_dec, npxl_ra, npxl_dec  # in deg

#%% cut function
def cut(data, wcs, hdu,center,box):
    size = u.Quantity(box, u.pix)
    cutout = Cutout2D(data,position=center,size=size,wcs=wcs)
    hdu.header.update(cutout.wcs.to_header())
    hdu.data = cutout.data

    wcs_cut = WCS(hdu.header)
    
    return(hdu.data,wcs_cut,hdu)


def cut2(data, wcs, hdu,center,box):
    size = u.Quantity(box, u.deg)
    cutout = Cutout2D(data,position=center,size=size,wcs=wcs)
    hdu.header.update(cutout.wcs.to_header())
    hdu.data = cutout.data

    wcs_cut = WCS(hdu.header)
    
    return(hdu.data,wcs_cut,hdu)

    
#%% find beam function
def find_beam(hdr):
    major = hdr["BMAJ"] # in degree
    minor = hdr["BMIN"]
    pa = hdr["BPA"] # position angle    
    return (major, minor, pa)


#%%


file2 = '/home/tatush/Desktop/19035/fits_files/19035-Kband-Barray-full.fits' # C band
RoseroK_data , RoseroK_wcs , RoseroK_hdu = load(file2)
centerK_deg = SkyCoord(286.5067083333333,6.776722222222222, frame='fk5', unit='deg')
centerK_pix = RoseroK_wcs.world_to_pixel(centerK_deg)
# print(centerK_pix)
RoseroK_data , RoseroK_wcs , RoseroK_hdu = cut(RoseroK_data , RoseroK_wcs , 
                                               RoseroK_hdu, (1048,1224),box=(55,55))


mydata , mywcs , myhdu = load('/home/tatush/Desktop/19035/fits_files/19035-cont-0.013-im2000-pI-r0.5.pbcor.fits')
mydata = np.ma.masked_invalid(mydata)
rms=mydata.std()
center_deg = SkyCoord(286.5067083333333,6.776722222222222, frame='fk5', unit='deg')
center_pix = mywcs.world_to_pixel(center_deg)
# print(center_pix)
mydata , mywcs , myhdu = cut(mydata , mywcs , myhdu, center=(935,1020),box=(300,300))



## plot
fig = plt.figure(figsize=(15,15))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=RoseroK_wcs)
    

    # -- axis label and ticks
ra = ax.coords[0]
dec = ax.coords[1]
ra.set_axislabel("RA (J2000)", minpad=0.8, fontsize=25)
dec.set_axislabel("Dec (J2000)", minpad=-1.0, fontsize=25)
ra.set_major_formatter('hh:mm:ss.ss')
ax.set_title('IRAS 19035+0641 A',fontsize=25)
ax.set_xlim(7,47)
ax.set_ylim(5,45)
    

ra.display_minor_ticks(True)
dec.display_minor_ticks(True)    
ax.tick_params(which='both',direction='in',color='white',length=10,width=2,
                   labelsize=20)
ax.tick_params(which='minor', length=5)

    # -- color plot
im=ax.imshow(RoseroK_data, vmin=-1e-7, vmax=RoseroK_data.max(), origin='lower', 
               cmap='plasma', transform=ax.get_transform(RoseroK_wcs))


plt.rcParams["lines.linewidth"] = 2
# rms = 13.5e-6
levels = [rms*-3,rms*3,rms*5,rms*10,rms*15]      
ax.contour(mydata, levels=levels, colors='k', transform=ax.get_transform(mywcs))

# plt.plot(pixel_pos.item(0),pixel_pos.item(1),'r+', markersize=100)

    # -- color bar
cbar = fig.colorbar(im,shrink=0.8)
for t in cbar.ax.get_yticklabels():
    t.set_fontsize(18)
cbar.ax.set_ylabel(r'$\mu$Jy/beam', fontsize=20)
cbar.ax.set_yticklabels(['0','50','100','150','200','250','300','350','400'])

## -- Scale bar
d = 2.3 # kpc
sb_length = (500) /(d * 1000) # I want a 1000 au scale bar
print(sb_length)
pix_scale = 0.05 # cell size
bsize = sb_length / pix_scale # size of bar in pix

xlim1, xlim2 = ax.get_xlim()
ylim1, ylim2 = ax.get_ylim()

aux = (xlim2-2-bsize+xlim2-2)/2.
ax.plot([xlim2-2-bsize,xlim2-2], 
                  [ylim2-1,ylim2-1], 
                  linewidth=3, c='w', alpha=1.0)

scale_bar_fontsize=16
scale_bar_text = '500 au'

ax.text(aux, ylim2-2.5, scale_bar_text, 
                  fontsize=22,
                  horizontalalignment='center', 
                  color='w')



obsra, obsdec, scale_ra, scale_dec, npxl_ra, npxl_dec = find_center_and_scale(RoseroK_hdu.header)
major, minor, pa = find_beam(myhdu.header)


plt.rcParams["hatch.linewidth"] = 3
# Find image center and scale
obsra, obsdec, scale_ra, scale_dec, npxl_ra, npxl_dec = find_center_and_scale(RoseroK_hdu.header)
major, minor, pa = find_beam(myhdu.header)
# beam_ra  = obsra  - 1 *scale_ra*npxl_ra + major 
# beam_dec = obsdec - 1 *scale_dec*npxl_dec + major
beam_ra = Angle('19h06m01.663s').deg
beam_dec = Angle('6d46m35.35s').deg

majorC, minorC, paC = find_beam(RoseroK_hdu.header)
beam3 = Ellipse((beam_ra, beam_dec), majorC, minorC, 90.0-paC, # BPA relative pos
                    transform=ax.get_transform('fk5'),
                    edgecolor='w', facecolor='none',lw=3)
ax.add_patch(beam3)

beam1 = Ellipse((beam_ra, beam_dec), major, minor, 90.0-pa, # BPA relative pos
                    transform=ax.get_transform('fk5'),
                    edgecolor='w', facecolor='w')
ax.add_patch(beam1)
