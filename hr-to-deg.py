#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 13:05:50 2022

@author: trodriguez
"""

from astropy.coordinates import Angle
import numpy as np
import pandas as pd

# df = pd.read_csv('/home/tatush/Desktop/VLA22A092/aux.csv',sep='\t')
# a = df['RA'] # center small
# b = df['Dec']
a = ['18h21m08.9907038240s']
b = ['-14d31m46.6854447519s']

ra=[]
dec=[]
for i in range(len(a)):
    ra = np.append(ra,Angle(a[i]).degree)
    dec=np.append(dec,Angle(b[i]).degree)
ra = ', '.join(str(element) for element in ra)
dec = ', '.join(str(element) for element in dec)
print(ra,dec)



