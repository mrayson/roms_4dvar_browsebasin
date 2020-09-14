#!/usr/bin/env python3
"""
  
  Module to process observations:
    obsgen : class to convert from raw data to ROMS observations using
             specific subclasses
    
  epoch = datetime(2000,1,1)

  Examples
  --------
  nws_obsgen.py SST_HIMAWARI_2017.nc --obs sst_himawari


I.J.

"""
import os, sys
from scipy import signal
import numpy as np
from netCDF4 import Dataset, date2num, num2date
import seapy
from seapy.roms.obsgen import obsgen
from datetime import datetime
import warnings, warnings.filterwarnings('ignore')

# define provenance for each observation platform so we can track them in assimilation
seapy.roms.obs.obs_provenance.update({ 190:"GLIDER_SG_150",
                                       200:"CTD_KISSME_1",
                                       240:"ARGO_CTS_D59",
                                       300:"SST_HIMAWARI",
                                       316:"SST_NOAA16",
                                       318:"SST_NOAA18",
                                       319:"SST_NOAA19",
                                       340:"SST_BOM",
                                       400:"AVISO",
                                       600:"RADAR",
                                       700:"ADCP_KISSME_1",
                                       800:"TZ_KISSME_1",
                                       })

def _makedirs(dirs):
    try:
        os.makedirs(dirs)
    except OSError:
        pass
    
def open_nc(fname):
    try:
        nc_data = Dataset(fname)
    except:
        print('Error opening file %s' %fname)
        sys.exit()
    return nc_data

def slice_time(time, drange):
    d1 = date2num(datetime.strptime(drange[0],'%Y%m%d'),time.units)
    d2 = date2num(datetime.strptime(drange[1],'%Y%m%d'),time.units)
    idx = np.argwhere(np.logical_and(time[:] >= d1, time[:] <= d2 ))
    return idx.flatten()
    

def new_dim(time, lat, lon, mask):  # Build the time, lon and lat arrays of appropriate size
    nt, nlat, nlon = len(time), lat.size, lon.size        
    time = np.resize(time, (nlat, nlon, nt))
    time = np.squeeze(np.transpose(time, (2,0,1)))[mask]
    lon = np.resize(lon, (nlat, nt, nlon))
    lon = np.squeeze(np.transpose(lon, (1,0,2)))[mask]
    lat = np.resize(lat, (nlon, nt, nlat))
    lat = np.squeeze(np.transpose(lat, (1,2,0)))[mask]
    return time, lat, lon           

def make_daily_obs(obs, fname):
    start = int(obs.time[0])
    end = int(obs.time[-1])
    head, tail = os.path.split(fname)
    tail = tail.split('.')[0]
    for t in range(start, end+1):
        ind = np.where(np.logical_and(obs.time >= t, obs.time < t+1))
        new_file = os.path.join(head, '%s_%s.nc' %(tail, t))
        obs[ind].to_netcdf(new_file)
        print('Grided file is %s' %new_file)
    return

class sst_himawari(obsgen):
   """
   class to process SST map netcdf files into ROMS observation
   files. This handles temperature fields.
   """
   def __init__(self, grid, dt, reftime=seapy.default_epoch, temp_error=0.4,
                temp_limits=None):
       self.temp_error = temp_error
       if temp_limits is None:
           self.temp_limits = (15,36)
       else:
           self.temp_limits = temp_limits
       super().__init__(grid, dt, reftime)

   def convert_file(self, file, drange=None, title="Himawari SST Obs"):
       """
       Load an SST file and convert into an obs structure for given time range (idx)

       """
       # Load SST Data
       nc = open_nc(file)
       lon = nc.variables["lon"][:]
       lat = nc.variables["lat"][:]
       time = nc.variables["time"]
       if drange is None:
            idx = np.arange(0,len(time))
       else:    
            idx = slice_time(time, drange)
            print('Using slice in time from %d till %d record' %(idx[0],idx[-1]))
       dtime = num2date(time[idx], time.units) - self.epoch
       time = list(map(lambda x: np.round(x.total_seconds()*seapy.secs2day*24*6)/24/6, dtime))
#       if temp.units in ('kelvin','Kelvin'):
#           dat += - 273.15       # hack for Matt's data
       dat = nc.variables["sst_night_mu"][idx] -273.15
       dat = np.ma.masked_outside(np.squeeze(dat),self.temp_limits[0], self.temp_limits[1])
       err = np.squeeze(nc.variables["sst_night_var"][idx])    
       err = np.maximum(np.sqrt(err), self.temp_error)	# we have variance instead of std and instrument error limit
       bad = np.where(err >= 2)	# those are really bad and do not consider them (makes smaller obs files)
       dat[bad] = np.ma.masked
       mask = np.ma.getmaskarray(dat)
       err[bad] = np.ma.masked
       nc.close()
       # Build the time, lon and lat arrays of appropriate size
       time, lat, lon = new_dim(time, lat, lon, ~mask)
       data = [seapy.roms.obs.raw_data("TEMP", "SST_HIMAWARI", dat.compressed(),
                                       err.compressed(), self.temp_error)]
       # call obs.gridder
       return seapy.roms.obs.gridder(self.grid, time, lon, lat, None,
                                    data, self.dt, title)

def call_gridder(grid_file, filename, new_file, obs_type, drange=None, daily=False):
    """
    Convenience function to call gridder.

    grid_file:  path to grid file
    filename:   netCDF file to grid
    obs:        data type
    drange:     date range ex 20140101 20140401
    """
    
    mygrid = seapy.model.asgrid(grid_file)
    # time-step of the model in days * 2 at least
    # 10 minutes time step for time window obs average
    dt = 600/86400.                
    
    gdict = {'SST_HIMAWARI': sst_himawari, 
             }    
    gen = gdict[obs_type](mygrid, dt)
    obs = gen.convert_file(filename, drange=drange)
    
    # return obs
    if daily and len(obs.time) >1:
        make_daily_obs(obs, new_file)
    else:         
        obs.to_netcdf(new_file)
        print('File %s has been gridded. New file is %s' %(filename, new_file))    
        
if __name__=='__main__':
    
    from argparse import ArgumentParser
    
    parser = ArgumentParser()
    parser.add_argument('filename', help='full path to raw file')
    parser.add_argument('--obs', default= 'sst_himawari', help='observation')
    parser.add_argument('--drange', default= None, nargs='+', type=str, help='date range ex 20140101 20140401')
    args = parser.parse_args()

    obs_type = args.obs.upper()
    
    wdir = os.path.dirname(os.path.realpath('__file__'))
    grid_file = os.path.join(wdir, 'grid.nc')
    fpath, file = os.path.split(args.filename)               
    # define new file filename   
    new_file = os.path.join(wdir, 'GRIDDED', file)        
    call_gridder(grid_file, args.filename, new_file, obs_type, args.drange, args.daily)
         
