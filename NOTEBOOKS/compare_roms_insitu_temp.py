#!/usr/bin/env python
# coding: utf-8

# # Plot validation between SUNTANS, ROMS and mooring data
# 

# In[1]:


# import sys
# sys.path.append('/home/mrayson/group/mrayson/code')
# !export PYTHONPATH="$PYTHONPATH:home/mrayson/code/mycurrents"
get_ipython().system('echo $PYTHONPATH')


# In[2]:


import numpy as np 
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
from datetime import datetime
import os
from scipy import signal
from scipy.interpolate import interp1d

import soda.utils.mysignal as sp
from soda.utils.timeseries import timeseries
from soda.dataio.suntans.sunprofile import Profile
from soda.utils.maptools import ll2utm
from soda.utils.modvsobs import ModVsObs
from soda.dataio import netcdfio
import soda.dataio as io
from soda.utils.myproj import MyProj
from soda.dataio.roms.romsio import roms_timeseries, get_depth

#from octant.tools import isoslice
from mycurrents import oceanmooring as om

from matplotlib import rcParams

#rcParams['font.family'] = 'sans-serif'
#rcParams['font.sans-serif'] = ['Bitstream Vera Sans']
#rcParams['font.serif'] = ['Bitstream Vera Sans']
rcParams["font.size"] = "12"
rcParams['axes.labelsize']='medium'


# In[3]:


# %matplotlib notebook


# In[4]:


#!ls /scratch/pawsey0219/ijanekovic/NWS_ASSIM/CYCLE
get_ipython().system('ls /scratch/pawsey0219/ijanekovic/ROMS_CYCLE/outputs/')
# romsfile = '/scratch/pawsey0219/ijanekovic/NWS_ASSIM/CYCLE/kp150_0.nc'
# ds = xr.open_dataset(romsfile)
# ds


# In[11]:


romsfile = '/scratch/pawsey0219/ijanekovic/ROMS_CYCLE/outputs/sta_6341.nc'
ds = xr.open_dataset(romsfile)
ds


# In[6]:


def get_comparison(lonstation, latstation, P, scenario, uv_ds ):
    # Load the suntans profile
    #sunfile = '../SCENARIOS/TEST_{}/ROMS_SUNTANS_{}_Profile.nc'.format(scenario, scenario)
    sunfile = '../SCENARIOS/TEST_{}/ROMS_SUNTANS_{}_Profile.nc'.format(scenario, scenario)

    sunts = Profile(sunfile)

    # Set the project
    #P = MyProj(None, utmzone=51, isnorth=False)
    # lonstation = uv_ds.X
    # latstation = uv_ds.Y
    xstation, ystation = P.to_xy(lonstation,latstation)


    #print 'Grabbing suntans points...'
    mo_u, Umod = get_suntans_mo(xstation, ystation, uv_ds.u.values,                     uv_ds.time.values, np.array([uv_ds.Z]), sunts, 'uc', 'm/s', )

    mo_v, Vmod = get_suntans_mo(xstation, ystation, uv_ds.v.values,                     uv_ds.time.values, np.array([uv_ds.Z]), sunts, 'vc', 'm/s')
    
    return mo_u, mo_v, Umod, Vmod


# In[7]:


def compare_imos_roms_temp(sunfile, Tfile, stationT, P, icdatastr, t1,t2, mo_T=None, plot=True):
    print(72*'#')
    print('#\t Station: %s'%stationT)
    #print 'Loading the suntans profile...'
    if mo_T is None:
        sunTS = Profile(sunfile)
        sunTS
        ##########
        # Load the 400 and 200 m velocity and temperature data
        ##########
        #print 'Loading mooring data...'

        Tobs = om.from_netcdf(Tfile, 'temperature', group=stationT).clip(t1,t2)
        #Tobs.X, Tobs.Y

        # Set the project
        #     P = MyProj('merc')
        xstation, ystation = P.to_xy(Tobs.X,Tobs.Y)
        #xstation, ystation

        #print 'Grabbing suntans points...'
        mo_T, Umod = get_roms_mo_temp(xstation, ystation, Tobs, sunTS, 'temp', 'degC')
    
    if plot:
        mo_T.printStatsZ()

        # Filter and  Convert to oceanmooring objects
        #dsobs = om.OceanMooring(mo_T.TSobs.t, mo_T.TSobs.filt_uneven(34*3600.,order=3), -mo_T.TSobs.Z, positive='down')
        #dsmod = om.OceanMooring(mo_T.TSmod.t, mo_T.TSmod.filt_uneven(34*3600.,order=3), -mo_T.TSmod.Z, positive='down')
        dsobs = om.OceanMooring(mo_T.TSobs.t, mo_T.TSobs.y, -mo_T.TSobs.Z, positive='up')
        dsmod = om.OceanMooring(mo_T.TSmod.t, mo_T.TSmod.y, -mo_T.TSmod.Z, positive='up')
        ####
        # Plot
        plt.figure(figsize=(8,9))
        ax1=plt.subplot2grid((3,3),(0,0),colspan=3)
        dsobs.contourf(np.arange(7,32,1))
        dsobs.contourf(np.arange(7,32,1), filled=False, cbar=False)

        plt.ylabel('Depth')
        ax1.set_xticklabels([])
        
        plt.title(stationT)


        ax2=plt.subplot2grid((3,3),(1,0),colspan=3)
        dsmod.contourf(np.arange(7,32,1))
        dsmod.contourf(np.arange(7,32,1), filled=False, cbar=False)

        plt.ylabel('Depth [m]')
        
        plt.title('%s'%icdatastr)

        # Mean temp
        ax3=plt.subplot2grid((3,3),(2,0),colspan=1)
        plt.plot(mo_T.meanMod, -mo_T.Z,'ro')
        plt.plot(mo_T.meanObs, -mo_T.Z,'bo')
        plt.xlabel('Temperature [$^{\circ}$C]')
        plt.ylabel('Depth [m]')
        plt.grid(b=True)
        plt.xlim([7,30])
        plt.legend(('Mod.','Obs.'), loc='upper left')

        # Bias
        ax4=plt.subplot2grid((3,3),(2,1),colspan=1)
        plt.plot(mo_T.bias, -mo_T.Z,'ko')
        plt.xlim(-3.25,3.25)
        ax4.set_xticks([-2,-1,0,1,2])
        ax4.set_yticklabels([])
        plt.xlabel('Bias [$^{\circ}$C]')
        plt.grid(b=True)


        # RMSE
        ax5=plt.subplot2grid((3,3),(2,2),colspan=1)
        plt.plot(mo_T.rmse, -mo_T.Z,'ko')
        plt.xlim(-3,3)
        ax5.set_yticklabels([])
        plt.xlabel('RMSE [$^{\circ}$C]')
        plt.grid(b=True)

        plt.tight_layout()
        
        #plt.savefig('../FIGURES/IMOS_v_SUNTANS_%s_%s.png'%(icdatastr,stationT), dpi=150)
        plt.show()
    
    return mo_T


# In[15]:


# Load a ROMS timeseries object

#tsroms = roms_timeseries(romfile, xyin, varname='ueastward')
def get_roms_station_old(romsfile, romsvar, xyin, zin):
    dsroms = xr.open_dataset(romsfile)
    #print(dsroms)
    zroms = get_depth(dsroms.s_rho.values, dsroms.Cs_r.values, dsroms.hc.values, dsroms.h.values)# , \
    #         Vtransform=dsroms.Vtransform.values)
    xroms = dsroms.lon_rho.values
    yroms = dsroms.lat_rho.values

    # Get the point
    #dist = np.abs( (xyin[0]-xroms) + 1j*(xyin[1]-yroms))
    #idx = np.argwhere(dist==dist.min())[0,0]
    #xyin, xroms[idx], yroms[idx], zroms[:,idx]

    data = dsroms[romsvar][:]
    #print(zroms[None,...].shape, data.shape)
    #print(zroms.ravel())

    Fi = interp1d(zroms[None,...].squeeze(), data.squeeze(), axis=1)
    
    return om.OceanMooring(dsroms.ocean_time.values, Fi(zin), zin)

def get_roms_station(romsfile, romsvar, xyin, zin):
    dsroms = xr.open_dataset(romsfile)
    #print(dsroms)
    zroms = get_depth(dsroms.s_rho.values, dsroms.Cs_r.values, dsroms.hc.values, dsroms.h.values)# , \
    #         Vtransform=dsroms.Vtransform.values)
    xroms = dsroms.lon_rho.values
    yroms = dsroms.lat_rho.values

    # Get the point
    dist = np.abs( (xyin[0]-xroms) + 1j*(xyin[1]-yroms))
    idx = np.argwhere(dist==dist.min())[0,0]
    print(xyin, xroms[idx], yroms[idx], zroms[:,idx])

    data = dsroms[romsvar][:,idx,:]
    print(zroms[...].shape, data.shape)
    #print(zroms.ravel())
    
    Fi = interp1d(zroms[...,idx].squeeze(), data.squeeze(), axis=1, fill_value='extrapolate')
    
    return om.OceanMooring(dsroms.ocean_time.values, Fi(zin), zin)

def get_roms_station_da(romsfile, romsvar, xyin, zin, cycle=2, ncycles=3):
    dsroms = xr.open_dataset(romsfile)
    #print(dsroms)
    zroms = get_depth(dsroms.s_rho.values, dsroms.Cs_r.values, dsroms.hc.values, dsroms.h.values)# , \
    #         Vtransform=dsroms.Vtransform.values)
    xroms = dsroms.lon_rho.values
    yroms = dsroms.lat_rho.values
    
    nt = dsroms.ocean_time.size//ncycles

    # Get the point
    dist = np.abs( (xyin[0]-xroms) + 1j*(xyin[1]-yroms))
    idx = np.argwhere(dist==dist.min())[0,0]
    xyin, xroms[idx], yroms[idx], zroms[:,idx]
    
    t0 = cycle*nt
    t1 = t0 + nt

    data = dsroms[romsvar][t0:t1,idx,:]
    print(zroms[...].shape, data.shape)
    #print(zroms.ravel())

    Fi = interp1d(zroms[...,idx].squeeze(), data.squeeze(), axis=1, fill_value='extrapolate')
    
    return om.OceanMooring(dsroms.ocean_time.values[t0:t1], Fi(zin), zin)


def get_roms_mo(X,Y, Uobs, tobs, zobs, romsfile, sunvar, units):
    """
    Get the suntans evuivalent points
    """
    Umod=get_roms_station(romsfile, sunvar, [X,Y], zobs)

    return ModVsObs(Umod.t, Umod.y, tobs, Uobs,             Z=zobs, units=units, varname=sunvar), Umod


# In[12]:




# Add the site location
P = MyProj(None, utmzone=51,isnorth=False)

lon_prelude, lat_prelude = P.to_ll(534322.7, 8475878.4)

sites = pd.DataFrame(np.array([[lon_prelude, 123.3535, 123.2803, 123.346383,127.56335],            [lat_prelude, -13.7621, -13.8170, -13.75895, -9.8122]]).T,
            index=['Prelude','NP250','DWR','KP150_phs2','ITFTIS'], columns=['lon','lat'])

# Load the u/v data as an xarray Dataset
# uv_ds = xr.open_dataset(uvfile)
# print(uv_ds, sites)


# In[9]:


days = range(6241,6321,4)
#romsfile = '/scratch/pawsey0219/ijanekovic/ROMS_CYCLE/outputs/sta_ver_6241.nc'
romsfiles = ['/scratch/pawsey0219/ijanekovic/ROMS_CYCLE/outputs/sta_{}.nc'.format(dd) for dd in days]

###
romsfile = romsfiles[5]
# IMOS station file
Tfile = '../DATA/KP150_Gridded_TP.nc'
stationT = 'KP150_phs2_T'
###

Tobsraw = om.from_netcdf(Tfile, 'watertemp', group=stationT)
zobs = Tobsraw.Z

xyin = np.array([sites['lon']['KP150_phs2'], sites['lat']['KP150_phs2']])
romsvar = 'temp'
Tmod = get_roms_station_da(romsfile, romsvar, xyin, zobs, cycle=0 )

t1, t2, dt = Tmod.t[0],Tmod.t[-1], (Tmod.t[1] - Tmod.t[0]).astype(float)*1e-9

Tobs = Tobsraw.clip(t1,t2).resample(dt,ndt=1)


mo_T = ModVsObs(Tmod.t, Tmod.y, Tobs.t, Tobs.y,             Z=-zobs, units='degC', varname='romsvar', interpmodel=False)

compare_imos_roms_temp(None, None, stationT, None, 'ROMS Phase', t1,t2, mo_T=mo_T, plot=True)

# plt.savefig('../FIGURES/temp_comparison_{}_ROMS.png'.format(stationT),dpi=150)


# In[10]:


romsfile = romsfiles[3]

stationT = 'ITFTIS'
xyin = np.array([sites['lon'][stationT], sites['lat'][stationT]])
romsvar = 'temp'

#Tfile = '../DATA/IMOS_Stack_Mooring_ITFTIS_201607_TEMP_60s.nc' # The first few cycles are in this folder
Tfile = '../DATA/IMOS_Stack_Mooring_ITFTIS_201702_TEMP_60s.nc'
dso = xr.open_dataset(Tfile)
Tobsraw = om.OceanMooring(dso.TIME.values, dso.TEMP.values, dso.NOMINAL_DEPTH.values)
zobs = -Tobsraw.Z

Tmod = get_roms_station_da(romsfile, romsvar, xyin, zobs, cycle=0)
t1, t2, dt = Tmod.t[0],Tmod.t[-1], (Tmod.t[1] - Tmod.t[0]).astype(float)*1e-9

#Tobsraw = om.from_netcdf(Tfile, 'TEMP', group=None).clip(t1,t2)
Tobs = Tobsraw.clip(t1,t2).resample(dt,ndt=1)


mo_T = ModVsObs(Tmod.t, Tmod.y, Tobs.t, Tobs.y,             Z=-zobs, units='degC', varname='romsvar', interpmodel=False)

compare_imos_roms_temp(None, None, stationT, None, 'ROMS Phase 0 ', t1,t2, mo_T=mo_T, plot=True)


# # Compare ROMS w/ and w/out 4DVAR

# In[13]:


def plot_da(mo_T, mo_T_da, titlestr):
    # Mean temp
    ax3=plt.subplot2grid((1,3),(0,0),colspan=1)
    plt.plot(mo_T.meanMod, -mo_T.Z,'ro')
    plt.plot(mo_T_da.meanMod, -mo_T.Z,'ko')
    plt.plot(mo_T.meanObs, -mo_T.Z,'bo')
    plt.xlabel('Temperature [$^{\circ}$C]')
    plt.ylabel('Depth [m]')
    plt.grid(b=True)
    plt.xlim([7,30])
    plt.legend(('ROMS','ROMS-4DVAR','Obs.'), loc='upper left')

    # Bias
    ax4=plt.subplot2grid((1,3),(0,1),colspan=1)
    plt.plot(mo_T.bias, -mo_T.Z,'ro')
    plt.plot(mo_T_da.bias, -mo_T_da.Z,'ko')
    plt.xlim(-3.25,3.25)
    ax4.set_xticks([-2,-1,0,1,2])
    ax4.set_yticklabels([])
    plt.xlabel('Bias [$^{\circ}$C]')
    plt.grid(b=True)
    plt.title(titlestr)


    # RMSE
    ax5=plt.subplot2grid((1,3),(0,2),colspan=1)
    plt.plot(mo_T.rmse, -mo_T.Z,'ro')
    plt.plot(mo_T_da.rmse, -mo_T_da.Z,'ko')
    plt.xlim(0,3)
    ax5.set_yticklabels([])
    plt.xlabel('RMSE [$^{\circ}$C]')
    plt.grid(b=True)
    plt.legend(('ROMS','ROMS-4DVAR'), loc='upper right')

    plt.tight_layout()

    #plt.savefig('../FIGURES/IMOS_v_SUNTANS_%s_%s.png'%(icdatastr,stationT), dpi=150)


# In[22]:


days = range(6241,6345,4)
romsfiles = ['/scratch/pawsey0219/ijanekovic/ROMS_CYCLE/outputs/sta_ver_{}.nc'.format(dd) for dd in days]
romsfiles_da = ['/scratch/pawsey0219/ijanekovic/ROMS_CYCLE/outputs/sta_{}.nc'.format(dd) for dd in days]

###
fileno = 25
###

romsfile = romsfiles[fileno]
romsfile_da = romsfiles_da[fileno]

stationT = 'ITFTIS'
xyin = np.array([sites['lon'][stationT], sites['lat'][stationT]])

### Data loading varies between sites
#Tfile = '../DATA/IMOS_Stack_Mooring_ITFTIS_201607_TEMP_60s.nc' # The first few cycles are in this folder
Tfile = '../DATA/IMOS_Stack_Mooring_ITFTIS_201702_TEMP_60s.nc'
dso = xr.open_dataset(Tfile)
Tobsraw = om.OceanMooring(dso.TIME.values, dso.TEMP.values, dso.NOMINAL_DEPTH.values)
zobs = -Tobsraw.Z

# Load the file w/out DA
romsvar = 'temp'
Tmod = get_roms_station(romsfile, romsvar, xyin, zobs)


t1, t2, dt = Tmod.t[0],Tmod.t[-1], (Tmod.t[1] - Tmod.t[0]).astype(float)*1e-9

#Tobsraw = om.from_netcdf(Tfile, 'TEMP', group=None).clip(t1,t2)
Tobs = Tobsraw.clip(t1,t2).resample(dt,ndt=1)

# 
Tmod_da = get_roms_station_da(romsfile_da, romsvar, xyin, zobs, cycle=0)



mo_T = ModVsObs(Tmod.t, Tmod.y, Tobs.t, Tobs.y,             Z=-zobs, units='degC', varname='romsvar', interpmodel=False)

mo_T_da = ModVsObs(Tmod_da.t, Tmod_da.y, Tobs.t, Tobs.y,             Z=-zobs, units='degC', varname='romsvar', interpmodel=False)

####
plt.figure(figsize=(12,4))
titlestr = '{} - Start Time {}'.format(stationT, days[fileno])
plot_da(mo_T, mo_T_da, titlestr)

outfile = '../FIGURES/InSitu_T_validation_4Dvar_{}_{}.png'.format(stationT, days[fileno])
plt.savefig(outfile, dpi=150)
print(t1,t2)


# In[29]:


fileno = 23
romsfile = romsfiles[fileno]
romsfile_da = romsfiles_da[fileno]

stationT = 'KP150_phs2'
xyin = np.array([sites['lon'][stationT], sites['lat'][stationT]])

### Data loading varies between sites
#Tfile = '../DATA/IMOS_Stack_Mooring_ITFTIS_201607_TEMP_60s.nc' # The first few cycles are in this folder
Tfile = '../DATA/IMOS_Stack_Mooring_ITFTIS_201702_TEMP_60s.nc'
dso = xr.open_dataset(Tfile)
Tobsraw = om.OceanMooring(dso.TIME.values, dso.TEMP.values, dso.NOMINAL_DEPTH.values)
zobs = -Tobsraw.Z

Tfile = '../DATA/KP150_Gridded_TP.nc'
filegroup = 'KP150_phs2_T'
Tobsraw = om.from_netcdf(Tfile, 'watertemp', group=filegroup)
zobs = Tobsraw.Z

# Load the file w/out DA
romsvar = 'temp'
Tmod = get_roms_station(romsfile, romsvar, xyin, zobs)


t1, t2, dt = Tmod.t[0],Tmod.t[-1], (Tmod.t[1] - Tmod.t[0]).astype(float)*1e-9

#Tobsraw = om.from_netcdf(Tfile, 'TEMP', group=None).clip(t1,t2)
Tobs = Tobsraw.clip(t1,t2).resample(dt,ndt=1)

# 
Tmod_da = get_roms_station_da(romsfile_da, romsvar, xyin, zobs, cycle=0)



mo_T = ModVsObs(Tmod.t, Tmod.y, Tobs.t, Tobs.y,             Z=-zobs, units='degC', varname='romsvar', interpmodel=False)

mo_T_da = ModVsObs(Tmod_da.t, Tmod_da.y, Tobs.t, Tobs.y,             Z=-zobs, units='degC', varname='romsvar', interpmodel=False)

####
plt.figure(figsize=(12,4))
titlestr = 'Site: {}, Start Time: {}'.format(stationT, days[fileno])
plot_da(mo_T, mo_T_da, titlestr)

outfile = '../FIGURES/InSitu_T_validation_4Dvar_{}_{}.png'.format(stationT, days[fileno])
plt.savefig(outfile, dpi=150)


# # Old Stuff

# In[36]:


romsfile = '/scratch/pawsey0219/ijanekovic/NWS_ASSIM/CYCLE/kp150_1.nc'

t1,t2 = datetime(2017,2,1), datetime(2017,2,5)
stationT = 'KP150_phs2_T'
Tobsraw = om.from_netcdf(Tfile, 'watertemp', group=stationT).clip(t1,t2)
Tobs = Tobsraw.resample(2*3600,ndt=1)

zobs = Tobs.Z
#zobs = Tobs.zvar.mean(axis=1)

xyin = np.array([sites['lon']['KP150_phs2'], sites['lat']['KP150_phs2']])
romsvar = 'temp'
Tmod = get_roms_station(romsfile, romsvar, xyin, zobs )


mo_T = ModVsObs(Tmod.t, Tmod.y, Tobs.t, Tobs.y,             Z=-zobs, units='degC', varname='romsvar', interpmodel=False)

compare_imos_suntans_temp(None, None, stationT, None, 'ROMS Phase 1 ', t1,t2, mo_T=mo_T, plot=True)

# plt.savefig('../FIGURES/temp_comparison_{}_ROMS.png'.format(stationT),dpi=150)


# In[37]:


romsfile = '/scratch/pawsey0219/ijanekovic/NWS_ASSIM/CYCLE/kp150_2.nc'

t1,t2 = datetime(2017,2,1), datetime(2017,2,5)
stationT = 'KP150_phs2_T'
Tobsraw = om.from_netcdf(Tfile, 'watertemp', group=stationT).clip(t1,t2)
Tobs = Tobsraw.resample(2*3600,ndt=1)

zobs = Tobs.Z
#zobs = Tobs.zvar.mean(axis=1)

xyin = np.array([sites['lon']['KP150_phs2'], sites['lat']['KP150_phs2']])
romsvar = 'temp'
Tmod = get_roms_station(romsfile, romsvar, xyin, zobs )


mo_T = ModVsObs(Tmod.t, Tmod.y, Tobs.t, Tobs.y,             Z=-zobs, units='degC', varname='romsvar', interpmodel=False)

compare_imos_suntans_temp(None, None, stationT, None, 'ROMS Phase 2 ', t1,t2, mo_T=mo_T, plot=True)


# In[21]:


get_ipython().system('ncdump -h ../DATA/IMOS_Stack_Mooring_ITFTIS_201607_TEMP_60s.nc')


# In[39]:


# Load the TIMOR sea data
Tfile = '../DATA/IMOS_Stack_Mooring_ITFTIS_201607_TEMP_60s.nc'
t1,t2 = datetime(2017,2,1), datetime(2017,2,5)

dso = xr.open_dataset(Tfile)
Tobsraw = om.OceanMooring(dso.TIME.values, dso.TEMP.values, dso.NOMINAL_DEPTH.values).clip(t1,t2)
#Tobsraw = om.from_netcdf(Tfile, 'TEMP', group=None).clip(t1,t2)
Tobs = Tobsraw.resample(2*3600,ndt=1)


# In[40]:


romsfile = '/scratch/pawsey0219/ijanekovic/NWS_ASSIM/CYCLE/itftis_0.nc'

# t1,t2 = datetime(2017,2,1), datetime(2017,2,5)
# stationT = 'KP150_phs2_T'
# Tobsraw = om.from_netcdf(Tfile, 'watertemp', group=stationT).clip(t1,t2)
# Tobs = Tobsraw.resample(2*3600,ndt=1)

stationT = 'ITFTIS'

zobs = -Tobs.Z
#zobs = Tobs.zvar.mean(axis=1)

xyin = np.array([sites['lon']['ITFTIS'], sites['lat']['ITFTIS']])
romsvar = 'temp'
Tmod = get_roms_station(romsfile, romsvar, xyin, zobs )


mo_T = ModVsObs(Tmod.t, Tmod.y, Tobs.t, Tobs.y,             Z=-zobs, units='degC', varname='romsvar', interpmodel=False)

compare_imos_suntans_temp(None, None, stationT, None, 'ROMS Phase 0 ', t1,t2, mo_T=mo_T, plot=True)


# In[41]:


romsfile = '/scratch/pawsey0219/ijanekovic/NWS_ASSIM/CYCLE/itftis_1.nc'

# t1,t2 = datetime(2017,2,1), datetime(2017,2,5)
# stationT = 'KP150_phs2_T'
# Tobsraw = om.from_netcdf(Tfile, 'watertemp', group=stationT).clip(t1,t2)
# Tobs = Tobsraw.resample(2*3600,ndt=1)

stationT = 'ITFTIS'

zobs = -Tobs.Z
#zobs = Tobs.zvar.mean(axis=1)

xyin = np.array([sites['lon']['ITFTIS'], sites['lat']['ITFTIS']])
romsvar = 'temp'
Tmod = get_roms_station(romsfile, romsvar, xyin, zobs )


mo_T = ModVsObs(Tmod.t, Tmod.y, Tobs.t, Tobs.y,             Z=-zobs, units='degC', varname='romsvar', interpmodel=False)

compare_imos_suntans_temp(None, None, stationT, None, 'ROMS Phase 1', t1,t2, mo_T=mo_T, plot=True)


# In[43]:


romsfile = '/scratch/pawsey0219/ijanekovic/NWS_ASSIM/CYCLE/itftis_2.nc'

# t1,t2 = datetime(2017,2,1), datetime(2017,2,5)
# stationT = 'KP150_phs2_T'
# Tobsraw = om.from_netcdf(Tfile, 'watertemp', group=stationT).clip(t1,t2)
# Tobs = Tobsraw.resample(2*3600,ndt=1)

stationT = 'ITFTIS'

zobs = -Tobs.Z
#zobs = Tobs.zvar.mean(axis=1)

xyin = np.array([sites['lon']['ITFTIS'], sites['lat']['ITFTIS']])
romsvar = 'temp'
Tmod = get_roms_station(romsfile, romsvar, xyin, zobs )


mo_T = ModVsObs(Tmod.t, Tmod.y, Tobs.t, Tobs.y,             Z=-zobs, units='degC', varname='romsvar', interpmodel=False)

compare_imos_suntans_temp(None, None, stationT, None, 'ROMS Phase 2', t1,t2, mo_T=mo_T, plot=True)


# In[27]:


plt.figure()
#plt.subplot(211)
#Tmod.contourf(np.arange(12,31))

plt.subplot(212)
Tobs.contourf(np.arange(6,31))


# In[ ]:




