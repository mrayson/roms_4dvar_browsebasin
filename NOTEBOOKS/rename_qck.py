#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 15:01:49 2020

@author: ivica
"""

import os, glob
from netCDF4 import Dataset, num2date, date2num
import numpy as np
    
def read_time(filein):
    nc = Dataset(filein,'r')
    time = num2date(nc.variables['ocean_time'][0], nc.variables['ocean_time'].units)
    time_day = date2num(time, 'days since 2000-1-1')
    nc.close()
    return time.strftime('%Y%m%d_%H'), time_day

type = 'qck'
files = sorted(glob.glob('%s*.nc' %type))
for i in range(0, len(files)):
   fin = files[i]
   time, time_day = read_time(fin)
#   fout = 'roms_nws_%s__%s.nc' %(time, type)
   fout = 'roms_nws_%d_%s.nc' %(np.int(time_day), type)
   print(fout)
   os.system('mv %s %s' %(fin, fout))
   #os.system('ncks -x -vwetdry_mask_psi,wetdry_mask_rho,wetdry_mask_u,wetdry_mask_v --ppc default=.3 %s -O -o %s' %(fout, fout))
