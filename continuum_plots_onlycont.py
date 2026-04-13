#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  8 15:34:53 2022

@author: trodriguez
"""

## CONTINUUM IMAGING
# ------------------
# ------------------
# This script will generate a plot where the contours of the Rosero+16 data
# will be overlaid to our VLA/22A-092 continuum data (in colors).
# ------------------


## --- CALLING PACKAGES ---
# --------------------------
from astropy.wcs import WCS
from astropy.visualization import wcsaxes
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
import pandas as pd
import ast
from astropy import units as u
# from astropy.visualization.wcsaxes import add_scalebar

import warnings
warnings.filterwarnings("ignore")


## --- FUNCTION DEFINING ---
# --------------------------

# LOADING DATA
def load(filename):
    file = get_pkg_data_filename(filename)
    hdu = fits.open(file)[0]

    data = hdu.data

    # Ensure the WCS is limited to 2 dimensions
    wcs = WCS(hdu.header, naxis=2)

    return data, wcs, hdu
# --------------------------

# FIND BEAM
def find_beam(hdr):
    major = hdr["BMAJ"] # in degree
    minor = hdr["BMIN"]
    pa = hdr["BPA"] # position angle    
    return (major, minor, pa)
# --------------------------

# FIND IMAGE CENTER AND SCALE
def find_center_and_scale(hdr):
    # obsra, obsdec = hdr["OBSRA"], hdr["OBSDEC"]
    scale_ra, scale_dec = hdr["CDELT1"], hdr["CDELT2"] # deg/pixel
    # center_pixel_ra, center_pixel_dec = hdr["CRPIX1"], hdr["CRPIX2"]
    npxl_ra, npxl_dec = hdr["NAXIS1"], hdr["NAXIS2"]     
    return scale_ra, scale_dec, npxl_ra, npxl_dec  # in deg
# --------------------------
# READ AND CONVERT RA/DEC COLUMNS
def read_ra_dec(file):
    df = pd.read_csv(file, sep='\t')
    ra = df['RA'].apply(ast.literal_eval).explode().astype(float).to_numpy()
    dec = df['DEC'].apply(ast.literal_eval).explode().astype(float).to_numpy()
    return ra, dec


all_masers_ra, all_masers_dec= read_ra_dec('CIIMMs/all.csv')
# CIIMMs_RA_Beu02, CIIMMs_Dec_Beu02= read_ra_dec('methanol_masers_degrees_Beuther02.csv')
# 


# ----------------------------------------------------

df = pd.read_csv('continuum_info.csv',sep='\t')
# src = ['20293', '19413']
src = ['18345','18151','G34','G53.11B', '20293', 'UYSO1','IRDC18223',\
    '18553','G53.25mm2','18470','18566','18182','19012', 'G11','G35',\
        'G23','G53.25mm4','G53.11mm2','20343','18521','18517','18440',\
            '18264', '19035','19413']
       
# d_far = [2,6.8]


#%%

for i in range(len(src)):
    
## LOAD FILES
    path = '/home/tatush/Desktop/Projects/VLA-22A092/Final_cont_fits/'
    
    our_fits = path + src[i] + '_ours.fits'
    rosK_fits = path + src[i] + '_rosK.fits'
    
    
    mydata , mywcs , myhdu = load(our_fits)
    mydata = np.ma.masked_invalid(mydata)
    rms = mydata.std()
    
    
    rosK_data, rosK_wcs, rosK_hdu = load(rosK_fits)
    rosK_data = np.ma.masked_invalid(rosK_data)
    if src[i] == '18345':
        rmsK=3e-5
    else:
        rmsK = rosK_data.std()
    
    
## IMAGE SETTINGS  
    center_str = df.loc[i, 'center']  # Get the string representation of the list
    center = eval(center_str) 
    

    center_deg = SkyCoord(center[0],center[1], frame='fk5', unit='deg')
    center_pix = mywcs.world_to_pixel(center_deg)
    
    box = eval(df.loc[i,'box'])

    size = u.Quantity((box[0],box[0]), u.pix)
    cutout = Cutout2D(mydata,position=center_pix,size=size,wcs=mywcs)
    myhdu.header.update(cutout.wcs.to_header())
    mydata = cutout.data
    mywcs = WCS(myhdu.header)
    
## SET UP FIGURE
    
    # fig = plt.figure(figsize=(15,15)) # ONE PLOT
    fig = plt.figure(figsize=(36, 15)) # 2 PLOTS SIDE BY SIDE
    # plt.figure(figsize=(15,15))
    # ax1 = plt.axes(projection=mywcs)

    ax1 = fig.add_subplot(1, 2, 1, projection=mywcs)    
    plt.subplots_adjust(hspace=0.2)
    
## TITLE
    title = df.loc[i, 'Source']

    ax1.set_title(title,fontsize=24)
    
## COORDINATES & AXIS
    ra = ax1.coords[0]
    dec = ax1.coords[1]
    ra.set_axislabel("RA (J2000)", minpad=.8, fontsize=25)
    dec.set_axislabel("Dec (J2000)", minpad=-1.0, fontsize=25)
    ra.set_major_formatter('hh:mm:ss.ss')
    ra.display_minor_ticks(True)
    dec.display_minor_ticks(True)    
    ax1.tick_params(which='both',direction='in',color='k',length=15,width=2,
                        labelsize=24)
    ax1.tick_params(which='minor', length=5)

    
    ## CONTOURS    
    our_cont_str = df.loc[i, 'ours_cont']  # Get the string representation of the list
    our_cont_list = eval(our_cont_str)  # Convert the string to a list using eval
    ax1.contour(mydata, levels=our_cont_list, colors='k',
                    transform=ax1.get_transform(mywcs))
    
    # IMAGE LIMITS
    ax1.set_xlim(-.5,box[0]-.5) 
    ax1.set_ylim(-.5,box[0]-.5)
    
    
    rosK_str = df.loc[i, 'rosK_cont']  # Get the string representation of the list
    rosK_list = eval(rosK_str)  # Convert the string to a list using eval
    ax1.contour(rosK_data, levels=rosK_list, colors='g',
                    transform=ax1.get_transform(rosK_wcs))

## BEAM    
    beam_pos  = mywcs.pixel_to_world(box[0]/12,box[0]/12) # 1/15 of the cut box
    
    ## Rosero K-band beam
    scale_raK, scale_decK, npxl_raK, npxl_decK = find_center_and_scale(rosK_hdu.header)
    majorK, minorK, paK = find_beam(rosK_hdu.header)
    
    
    beamK = Ellipse(([beam_pos.ra.deg, beam_pos.dec.deg]), majorK, minorK, 90.0-paK, # BPA relative pos
                    transform=ax1.get_transform('icrs'),
                    edgecolor='g', facecolor='none',lw=3)#, hatch="//",lw=3)
    ax1.add_patch(beamK)
        
   
    ## Our beam
    scale_ra, scale_dec, npxl_ra, npxl_dec = find_center_and_scale(myhdu.header)
    major, minor, pa = find_beam(myhdu.header)
    
    beam = Ellipse(([beam_pos.ra.deg, beam_pos.dec.deg]), major, minor, 90.0-pa, # BPA relative pos
                    transform=ax1.get_transform('icrs'),
                    edgecolor='k', facecolor='k',lw=3)
    ax1.add_patch(beam)
    

   
## SCALE BAR

    d = df.loc[i, 'd'] # input distance in kpc
    d = d*1000 # make it pc

    d = df.loc[i, 'd']
    pix_scale = 0.013 # PIX SIZE
    
    # Define the scale bar length you want in AU
    physical_scale_length = 1000  # in AU
    
    # Get the distance to the source in parsecs (my input is in kpc in this case)
    distance_to_source_pc = d * 1000  
    # dist_far = (d_far[i]*1000*physical_scale_length)/ distance_to_source_pc

    # Calculate the scale bar length in pixels
    scale_bar_length_pixels = physical_scale_length / distance_to_source_pc /  pix_scale
    # print(src[i],scale_bar_length_pixels)

    # Convert the scale bar length to the data coordinates
    scale_bar_length_data = scale_bar_length_pixels * mywcs.wcs.cdelt[0]
    
    # Calculate the scale bar text based on the scale bar length (e.g., in AU)
    scale_bar_text = f"{physical_scale_length} au"
    # scale_bar_text = f"{physical_scale_length}/{int(dist_far)} au"
    
    image_width = box[0]
    image_height = box[0]
    
    aux = box[0]/10
    aux_y = box[0]/25
    
    pos_x = [image_width - scale_bar_length_pixels - aux, image_width - aux]
    pos_y = [image_height - aux_y, image_height - aux_y]
    
    # Draw the scale bar
    ax1.plot(pos_x,pos_y , color='k', lw=8)
    
    ax1.text((pos_x[0]+pos_x[1])/2, pos_y[0]-aux_y, scale_bar_text, color='k',
            ha='center', va='center', fontsize=24)

    
## CONTINUUM LABELS    
    labels_og = df.loc[i, 'labels'] 
    labels = eval(labels_og)
    
    lab_RA_str = df.loc[i, 'lab_RA']  # Get the string representation of the list
    lab_RA = eval(lab_RA_str)  # Convert the string to a list using eval
    
    lab_Dec_str = df.loc[i, 'lab_Dec']  # Get the string representation of the list
    lab_Dec = eval(lab_Dec_str)  # Convert the string to a list using eval
    
    for j in range(len(lab_Dec)):
        plt.text(lab_RA[j],lab_Dec[j],labels[j],size=35,color='k', 
                  transform=ax1.get_transform('fk5'))
    
## WATER MASERS
    # Some sources don't have water masers detected, so the water maser plotting
    # is activated with a boolean
    wm = df.loc[i,'WM_if']
    if wm:
        wm_ra_str = df.loc[i, 'WM_RA'] 
        wm_dec_str = df.loc[i, 'WM_Dec'] 
        
        wm_ra = eval(wm_ra_str)
        wm_dec = eval(wm_dec_str)
    
        for k in range(len(wm_ra)):
            kk = k + 2
            plt.plot(wm_ra[k],wm_dec[k],'r+', fillstyle='full', markersize=30, mew=3, alpha=1, 
                      markerfacecolor='r',transform=ax1.get_transform('fk5'))
            # plt.text(wm_ra[k]+0.00003,wm_dec[k]-0.000035,str(kk),size=35,color='r', 
            #          transform=ax.get_transform('fk5'))
            # THE LABELS WEREN'T LOOKING GOOD, I DECIDED NOT TO PLOT THEM
    
## CIIMMs 
    
    aux_ra = ax1.get_xlim()
    aux_dec = ax1.get_ylim()
    aux1 = mywcs.pixel_to_world(aux_ra[0]+(box[0]/2),aux_dec[0]+(box[0]/2))
    
        
    for k in range(len(all_masers_ra)):
        pos =  SkyCoord(all_masers_ra[k],all_masers_dec[k], frame='fk5', unit='deg')
        delta_ra = (pos.ra - aux1.ra).to(u.arcsec)
        if abs(delta_ra.value) <= 20:              
            plt.plot(all_masers_ra[k],all_masers_dec[k],'ko', fillstyle='full', markersize=10, 
                      mew=1, alpha=.8,markerfacecolor='r',  transform=ax1.get_transform('fk5'))
            
    
    # plt.show()

# DRAW BOX TO SHOW THE AREA ZOOMED-IN ON THE SECOND FIGURE
    from matplotlib.patches import Rectangle
    center_im = box[0] / 2

    x_lower_left = center_im - (box[1] / 2)
    y_lower_left = center_im - (box[1] / 2)
    
    r = Rectangle((x_lower_left, y_lower_left), box[1], box[1], edgecolor='k', facecolor='none', lw=1.5)
    
    ax1.add_patch(r)
    
## ADD TEXT OF CONTOUR LEVELS
    x = ax1.get_xlim()
    text_y_pos1 = -(box[0]/6)
    
    
    K_levels = [x / rmsK for x in rosK_list]
    formatted_Klevels_aux = ['{:.1f}'.format(a) for a in K_levels]
    formatted_Klevels= ', '.join(formatted_Klevels_aux)    
    
    
    rmsK_aux = rmsK*1e6
    formatted_rmsK = '{:.2f}'.format(rmsK_aux)  # Change precision as needed
    
       
    plt.text(min(x), text_y_pos1, f'[{formatted_Klevels}] x {formatted_rmsK} $\mu$Jy beam$^{{-1}}$', ha='left',fontsize=22,color='g')
    

    # text_y_pos = -(box[1]/6)
    # levels = [x / rms for x in our_cont_list]
    # formatted_levels_aux = ['{:.1f}'.format(a) for a in levels]
    # formatted_levels= ', '.join(formatted_levels_aux)
    # rms_aux = rms*1e6
    # formatted_rms = '{:.2f}'.format(rms_aux)  # Change precision as needed
    # plt.text(min(x),text_y_pos, f'[{formatted_levels}] x {formatted_rms} $\mu$Jy$\,$beam$^{{-1}}$', ha='left',fontsize=22)
#     # plt.text(min(x), text_y_pos2, f' ', ha='left',fontsize=22)
    
#%% SECOND PLOT - SUBFIGURE 2
# LOAD FILES
            
    mydata , mywcs , myhdu = load(our_fits)
    
## IMAGE SETTINGS  
    center_str = df.loc[i, 'center']  # Get the string representation of the list
    center = eval(center_str) 
    
    center_deg = SkyCoord(center[0],center[1], frame='fk5', unit='deg')
    center_pix = mywcs.world_to_pixel(center_deg)
    
    size = u.Quantity((box[1],box[1]), u.pix)
    cutout = Cutout2D(mydata,position=center_pix,size=size,wcs=mywcs)
    myhdu.header.update(cutout.wcs.to_header())
    mydata = cutout.data
    mywcs = WCS(myhdu.header)
    
## SET UP FIGURE
    ax2 = fig.add_subplot(1, 2, 2, projection=mywcs)

    # FOR 1 PLOT:
    # fig = plt.figure(figsize=(15,15))
    # ax = plt.subplot(1,2,2)
    # ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], projection=mywcs) # add projection
    
## TITLE
    title = df.loc[i, 'Source']

    ax2.set_title(title,fontsize=24)
    
## COORDINATES & AXIS
    ra = ax2.coords[0]
    dec = ax2.coords[1]
    ra.set_axislabel("RA (J2000)", minpad=0.8, fontsize=25)
    dec.set_axislabel("Dec (J2000)", minpad=-1.0, fontsize=25)
    ra.set_major_formatter('hh:mm:ss.ss')
    ra.display_minor_ticks(True)
    dec.display_minor_ticks(True)    
    ax2.tick_params(which='both',direction='in',color='k',length=15,width=2,
                        labelsize=24)
    ax2.tick_params(which='minor', length=5)

    
## CONTOURS    
    our_cont_str = df.loc[i, 'ours_cont']  # Get the string representation of the list
    our_cont_list = eval(our_cont_str)  # Convert the string to a list using eval
    ax2.contour(mydata, levels=our_cont_list, colors='k',
                    transform=ax2.get_transform(mywcs))

        # IMAGE LIMITS
    ax2.set_xlim(-.5,box[1]-.5) 
    ax2.set_ylim(-.5,box[1]-.5)
    
## SCALE BAR
    d = df.loc[i, 'd']
    pix_scale = 0.013
    
    # Define the physical scale bar length in AU
    physical_scale_length = 500  # in AU
    
    # Get the distance to the source in parsecs (replace this value with the actual distance)
    distance_to_source_pc = d * 1000  # Replace with the actual distance in parsecs
    # dist_far = (d_far[i]*1000*physical_scale_length)/ distance_to_source_pc

    # Calculate the scale bar length in pixels
    scale_bar_length_pixels = physical_scale_length / distance_to_source_pc / pix_scale

    # Convert the scale bar length to the data coordinates
    scale_bar_length_data = scale_bar_length_pixels * mywcs.wcs.cdelt[0]
    
    # Calculate the scale bar text based on the scale bar length (e.g., in AU)
    # scale_bar_text = f"{physical_scale_length}/{int(dist_far)} au"
    scale_bar_text = f"{physical_scale_length} au"
    
    image_width = box[1]
    image_height = box[1]
    
    aux = box[1]/10
    aux_y = box[1]/25
    
    pos_x = [image_width - scale_bar_length_pixels - aux, image_width - aux]
    pos_y = [image_height - aux_y, image_height - aux_y]
    
    # Draw the scale bar
    ax2.plot(pos_x,pos_y , color='k', lw=8)
    
    ax2.text((pos_x[0]+pos_x[1])/2, pos_y[0]-aux_y, scale_bar_text, color='k',
            ha='center', va='center', fontsize=24)
    
    
    
## CONTINUUM LABELS    
    labels_og = df.loc[i, 'labels'] 
    labels = eval(labels_og)    
    
    lab_RA_str = df.loc[i, 'lab_RA']  # Get the string representation of the list
    lab_RA = eval(lab_RA_str)  # Convert the string to a list using eval
    
    lab_Dec_str = df.loc[i, 'lab_Dec']  # Get the string representation of the list
    lab_Dec = eval(lab_Dec_str)  # Convert the string to a list using eval
    
    xlim1_pix, xlim2_pix = ax2.get_xlim()
    ylim1_pix, ylim2_pix = ax2.get_ylim()
    
    blc = mywcs.pixel_to_world(xlim1_pix, ylim1_pix)
    trc = mywcs.pixel_to_world(xlim2_pix, ylim2_pix)
    
    # for j in range(len(lab_Dec)):
    for j, (xi, yi) in enumerate(zip(lab_RA, lab_Dec)):
        if  trc.ra.value <= xi <= blc.ra.value and blc.dec.value <= yi <= trc.dec.value:
            plt.text(lab_RA[j],lab_Dec[j],labels[j],size=35,color='k', 
                      transform=ax2.get_transform('fk5'))
    
## WATER MASERS
    wm = df.loc[i,'WM_if']
    if wm:
        wm_ra_str = df.loc[i, 'WM_RA'] 
        wm_dec_str = df.loc[i, 'WM_Dec'] 
        
        wm_ra = eval(wm_ra_str)
        wm_dec = eval(wm_dec_str)
    
        for k in range(len(wm_ra)):
            kk = k + 2
            
            plt.plot(wm_ra[k],wm_dec[k],'r+', fillstyle='full', markersize=30, mew=4, alpha=1, 
                      markerfacecolor='r',transform=ax2.get_transform('fk5'))
    
    
    aux_ra = ax2.get_xlim()
    aux_dec = ax2.get_ylim()
    aux2 = mywcs.pixel_to_world(aux_ra[1]+(box[1]/2),aux_dec[1]+(box[1]/2))
    
        
    for k in range(len(all_masers_ra)):
        pos =  SkyCoord(all_masers_ra[k],all_masers_dec[k], frame='fk5', unit='deg')
        delta_ra = (pos.ra - aux2.ra).to(u.arcsec)
        if abs(delta_ra.value) <= 20:              
            plt.plot(all_masers_ra[k],all_masers_dec[k],'ko', fillstyle='full', markersize=12, 
                      mew=1, alpha=.8,markerfacecolor='r',  transform=ax2.get_transform('fk5'))
            
            
    
    
## BEAM    
    beam_pos  = mywcs.pixel_to_world(box[1]/14,box[1]/15) 
    
    scale_ra, scale_dec, npxl_ra, npxl_dec = find_center_and_scale(myhdu.header)
    major, minor, pa = find_beam(myhdu.header)
    
    beam = Ellipse(([beam_pos.ra.deg, beam_pos.dec.deg]), major, minor, 90.0-pa, # BPA relative pos
                    transform=ax2.get_transform('icrs'),
                    edgecolor='k', facecolor='k',lw=3)
    ax2.add_patch(beam)
    
    
    text_y_pos = -(box[1]/6)
    levels = [x / rms for x in our_cont_list]
    formatted_levels_aux = ['{:.1f}'.format(a) for a in levels]
    formatted_levels= ', '.join(formatted_levels_aux)
    rms_aux = rms*1e6
    formatted_rms = '{:.2f}'.format(rms_aux)  # Change precision as needed
    plt.text(min(x),text_y_pos, f'[{formatted_levels}] x {formatted_rms} $\mu$Jy$\,$beam$^{{-1}}$', ha='left',fontsize=22)


    plt.text(min(x), text_y_pos, f' ', ha='left',fontsize=22)
        
    # path_to_plot = '/home/tatush/Desktop/Projects/VLA-22A092/Continuum_pdf_images/'
    # plot_name = path_to_plot + src[i] + '_1plot.pdf'
    # plt.savefig(plot_name, bbox_inches='tight')
    # plt.close()

    # print('Done with',src[i])
    
    
    
    
#     plt.show()
    
    path_to_plot = '/home/tatush/Desktop/Projects/VLA-22A092/Continuum_pdf_images/'
    plot_name = path_to_plot + src[i] + '_1plot_cont.pdf'
    plt.savefig(plot_name, bbox_inches='tight')
    plt.close()

    print('Done with',src[i])
    
    
    