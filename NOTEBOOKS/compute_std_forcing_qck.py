# -*- coding: utf-8 -*-
"""
Ugly python script for creating standard deviation and average files for 4DVAR ROMS 
and ATMOSPHERIC boundary layer corrections (this is first file we need)

I am using just subset from daily history (and quick) files for testing 

Template for std_frc netCDF file is at magnus, with prefilled cordinates etc:
/scratch/pawsey0219/ijanekovic/NWS_ASSIM/std_frc_template.nc
and locally in /mnt/qnap/AUSTRALIA/SHELL/4DVAR/

Note, it is the same for avg_frc, just file name and variable values differ (std -> avg) 

IJ

"""

from netCDF4 import Dataset, num2date, date2num
from datetime import datetime
import numpy as np
from os import system as sys
from bunch import Bunch

# files are sorted with start day in name (days starting from 2000-1-1)
 
def data_read(file, var):
    try:
        nc = Dataset(file, 'r')
    except:
        print('Error opening file %s' %file)
        sys.exit()
    return nc.variables[var][:]
    nc.close()
    

vars_frc = ['ocean_time', 'shflux', 'ssflux', 'sustr', 'svstr']
t1 = datetime(2017,2,1)     # this is our start of 4DVar 
t2 = datetime(2017,5,31)
t1_day = np.int(date2num(t1, 'days since 2000-1-1'))
t2_day = np.int(date2num(t2, 'days since 2000-1-1'))
day_recs = 24
# using overlap of 2 days, 1 day before and 1 day after the assim cycle so using 6 files/cycle with 24 recs
cycles = np.arange(t1_day, t2_day, 4)
for c in cycles:
    date = num2date(c, 'days since 2000-1-1')
    for i in np.arange(c-1, np.min(c+5, t2_day)):
        fn = "roms_nws_%d_qck.nc" %i
        if i==c-1:
            data = Bunch()
            for v in vars_frc:
                data[v] = data_read(fn, v)
        else:
            for v in vars_frc:
                data[v] = np.append(data[v],data_read(fn, v), axis=0)
    print('done with reading data')
    # prepare files
    outfile = 'std_f_%d.nc' %c
    sys('cp std_f_template.nc %s' %outfile)
    nc_std = Dataset(outfile, 'a')
    outfile = 'avg_f_%d.nc' %c
    sys('cp avg_f_template.nc %s' %outfile)
    nc_avg = Dataset(outfile, 'a')  
    for v in vars_frc:
        if v == 'ocean_time':
            std = 0
            avg = 0
            #units = nc_std[v].units
            #std = date2num(date, '%s since 2000-1-1' %units)   # actually std vars are references since init time so no need
            #avg = date2num(date, '%s since 2000-1-1' %units)
        elif v == 'shflux':
            ## have to make for each hour separate std and avg and then take mean of them as diurnal cycle
            tmps=[]; tmpa=[]            # surely there is better way, more pythonic
            for j in np.arange(day_recs):
                tmps.append(np.std(data[v][j:-1:day_recs,:,:], axis = 0))
                tmpa.append(np.mean(data[v][j:-1:day_recs,:,:], axis = 0))
            std = np.ma.masked_greater(np.std(tmps, axis = 0, keepdims=True), 1e9)
            avg = np.ma.masked_greater(np.mean(tmpa, axis = 0, keepdims=True), 1e9)
            del(tmps); del(tmpa)
        else:
            # compute std and avg
            std = np.std(data[v], axis = 0, keepdims=True)
            avg = np.mean(data[v], axis = 0, keepdims=True)
        nc_std.variables[v][:] = np.ma.masked_greater(std, 1e9)
        nc_avg.variables[v][:] = np.ma.masked_greater(avg, 1e9)
        del(std); del(avg)

        print('done for %s' %v)   
    nc_std.close()
    nc_avg.close()
