#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Merge observation files into one file for given date range and step.
    For file merging function merg_files from seapy.roms.obs is used

    usage:
        obs_merge_files.py drange --step --obs

    drange: yyyymmdd yyyymmdd
    step: number of days within one file, this is assimilation cycle window lenght    
    obs: list of obs Type of data, e.g. ADCP, HF, if not specified list is a dict. 
    key from nws_obsgen name_dict

    Examples
    --------
    obs_merge_files.py 20170201 20170301 --step 4
    obs_merge_files.py 20140101 20140112 --step 4 --obs SST_NOAA ADCP TEMP

    Written by B.Grcic
"""

from seapy.roms.obs import merge_files
from seapy.model import asgrid
from netCDF4 import date2num
from datetime import datetime
import os, glob

import warnings 
warnings.filterwarnings('ignore')

#            folder          file
name_dict = {'HF_RADAR':    'HF*', 
             'SST_HIMAWARI':'SST*',
             'GLIDER':      'GLIDER*',
             'CTD':         'CTD*',
             'TEMP':        'TEMP*',
             'ADCP':        'ADCP*',
             'ARGO':        'ARGO*',
             'AVISO':       'AVISO*'
            }

path = './'
dt = 600/86400.
    
if __name__=='__main__':
    
    from argparse import ArgumentParser
    
    parser = ArgumentParser()
    parser.add_argument('drange', nargs = '+', type = str, help = 'date range example: 20140101 20140401')
    parser.add_argument('--step', help = 'number of days')
    parser.add_argument('--obs', default = None, nargs = '+', type = str, help = 'list of obs files')
    args = parser.parse_args()
    

    start, end = args.drange
    step = int(args.step)
    units = 'days since 2000-01-01'
    
    sdate = int(date2num(datetime.strptime(start,'%Y%m%d'), units=units))
    edate =int(date2num(datetime.strptime(end,'%Y%m%d'), units=units))
    drange = [(i, i+step) for i in range(sdate,edate,step)]
    
    # remove obs on the edge of the model grid within 5 cells
    grid_file = os.path.join(path, 'grid.nc')
    grid = asgrid(grid_file)
    pad = 5    
    limits = {'north':grid.eta_rho-pad, 'south':pad, 'east':grid.xi_rho-pad, 'west':pad}

    obs_list = args.obs    
    if args.obs is None:
        obs_list = name_dict.keys()
        
    obs_files = []
    for name in obs_list:
        data_dir = os.path.join(path, 'GRIDDED', name)
        fname = name_dict[name] + '.nc'
        flist = glob.glob(os.path.join(data_dir, fname))
        if flist == []:
            continue
        for f in flist:            
            obs_files.append(f)

    print('Run: ', args.drange, 'step: ', args.step)
    print('List of files to merged:')    
    for x in obs_files: print(x)        
    new_file = str(os.path.join(path, 'GRIDDED/Merged_obs', 'obs_#.nc' ))
    merge_files(obs_files, new_file, drange, dt=dt, limits=limits)         
    
