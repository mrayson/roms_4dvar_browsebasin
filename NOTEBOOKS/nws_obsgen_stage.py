# -*- coding: utf-8 -*-
"""
    Script that collects files in specified folders and calls nws_obsgen gridding for each file.
    In case of ADCP/HF RADAR data, scripts will rotate east and north components and then call gridding.

    usage:

    nws_obsgen_stage.py type --drange 20140101 20140104 --type sst_himawari --daily True

    type : sst_himawari, adcp, argo, glider, hf_radar etc...
    drange: yyyymmdd yyyymmdd
    daily: true if gridded obs need to be saved in daily files

    Examples
    --------

    nws_obsgen_stage.py --type sst_himawari
    nws_obsgen_stage.py --type adcp --drange 20140101 20140401 --daily True
    
    
    It is hardcoded to look for grid in grid.nc file and in RAW folder for raw data with subfolder named as obs_type.  
    Outputs is in GRIDDED folder with the same filename prefix. 
    If you need you can change that.
    
    Written by B.Grcic
"""

import os, glob, sys
from nws_obsgen import call_gridder
import warnings, warnings.filterwarnings('ignore')

wdir = os.path.dirname(os.path.realpath('__file__'))

def _makedirs(dirs):
    try:
        os.makedirs(dirs)
    except OSError:
        pass

def _exist(file):
    if glob.glob(file) == []:
        return False
    return True    

if __name__=='__main__':
    
    from argparse import ArgumentParser
    
    parser = ArgumentParser()
    parser.add_argument('--type', default = 'sst_himawari', 
    help = 'data type to grid; adcp, argo, sst_hmawari, sst_bom, aviso, glider, hf_radar, ctd ..., default is --type sst_himawari')
    parser.add_argument('--drange', default = 'None', nargs = '+', type = str, 
    help = 'process data for specific date range, --drange 20140101 20140401 in yyyymmdd format, no default')
    parser.add_argument('--daily', default = False, help = 'create daily obs files --daily True/False, default is False')
    args = parser.parse_args()

    # look for the folder where to store gridded obs files 
    obs_type = args.type.upper()
    new_dir = os.path.join(wdir, 'GRIDDED/%s' %obs_type)
    _makedirs(new_dir)
    
    # look for the folder where are raw data
    folder = os.path.join(wdir, 'RAW', obs_type)
    flist = glob.glob(os.path.join(folder, '*.nc'))
    if flist == []:
        print('%s Folder is empty' %os.path.join(folder, '*.nc'))        
        sys.exit()
    
    # grid file that will be used for mapping 
    grid_file = os.path.join(wdir, 'grid.nc')
    
    for f in flist:
        path, name = os.path.split(f)    
        if (obs_type == 'ADCP' or obs_type == 'HF_RADAR'):        
            print('Will rotate components into the grid cooradinate') 
        
        new_file = os.path.join(new_dir, name)
        # check if already exist
        if _exist(new_file):
            print('File %s exist, remove it to continue' %new_file)            
            continue
        
        print('Gridding file: ', f)
        call_gridder(grid_file, f, new_file, obs_type, args.drange, args.daily)
            
    
