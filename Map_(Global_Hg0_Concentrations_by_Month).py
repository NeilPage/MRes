#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 14:39:31 2017

THIS SCRIPT PLOTS THE GLOBAL HG0 CONCENTRATIONS FOR EACH MONTH OF THE YEAR

@author: ncp532
"""
# File system packages
from netCDF4 import Dataset				# function used to open single netcdf file
from netCDF4 import MFDataset				# function used to open multiple netcdf files

# Drawing packages
from mpl_toolkits.basemap import Basemap	# function to plot maps
import matplotlib.pyplot as plt			# import package as shorter nickname
import matplotlib as mpl				    # import package as shorter nickname
import matplotlib.colors as colors			# import package as shorter nickname
from matplotlib import cm

# Date and Time handling package
from datetime import datetime,timedelta		# functions to handle date and time

# Data handing packages
import math						         # package to do mathematical calculations (such as log, exponential)
import numpy as np

# Define the dataset
dataset1 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_InvOcean_Files/trac_avg.201*.nc')
#dataset2 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Default_Files/trac_avg.201*.nc')


# Get the GEOS-Chem data for the stations
airden1=dataset1.variables['BXHGHT_S__AIRNUMDE'][:] 
#airden2=dataset2.variables['BXHGHT_S__AIRNUMDE'][:] 
D1 =(dataset1.variables['IJ_AVG_S__Hg0'][:])*(1e-3)*200.59/(6.0221*1e23)*airden1*1e6 
#D2 =(dataset2.variables['IJ_AVG_S__Hg0'][:])*(1e-3)*200.59/(6.0221*1e23)*airden2*1e6 

time=dataset1.variables['time'][:]                           # in minutes
lat=dataset1.variables['lat'][:]                             # in degrees
lon=dataset1.variables['lon'][:]                             # in degrees

# Convert GEOS-Chem time into standard time (GEOS-Chem time is hours since 1985/01/01/0000)
d0=datetime(1985,1,1,0)
date = []
for t in time:
    hrs=timedelta(hours=int(t))
    date.append(d0+hrs)

def base(x,tit):
    #x=x[-1,0,:,:]
    #Create the map and plot
    lons, lats = np.meshgrid(lon,lat)
    m=Basemap(llcrnrlat=-85,  urcrnrlat=85,
              llcrnrlon=-180, urcrnrlon=179,
              resolution='c',projection='merc')
    #np.transpose because it mixes the lat lon order
    cs=m.pcolor(lons, lats, x, latlon=True,vmin=0.5, vmax=4.5, cmap=cm.plasma)#,
                #vmin=0,vmax=1300)
    m.drawcoastlines()
    m.drawparallels([0],labels=[0,0,0,0]) # draw equator, no label
    # Plot location markers
    lons = [144.6832, 131.0447, 151.1018, 18.4897, 77.3423]
    lats = [-40.6832, -12.2491, -32.4777, -34.3535, -37.4748]
    x,y = m(lons, lats)
    m.plot(x, y, 'bo', markersize=3)
    #add title, colorbar
    q=m.colorbar(cs,"right",size="5%", pad="2%")
    #plt.clim(-1200, 1200) #use just for the total
    q.set_label('$ng/m^3$', fontsize=10)
    plt.title(tit, fontsize=10)
    plt.xlabel('Longitude [$^\circ$]', fontsize=10)
    plt.ylabel('Latitude [$^\circ$]', fontsize=10)
    plt.suptitle('Mean Hg$^0$ Concentrations 2013', fontsize=15)
    return m
# Plot the axis for each graph

fig = plt.figure()
plt.subplots_adjust(hspace=0.3)
ax1=plt.subplot(431)
ax1=base(D1[0,0,:,:], tit='Jan')

ax1=plt.subplot(432)
ax1=base(D1[1,0,:,:], tit='Feb')

ax1=plt.subplot(433)
ax1=base(D1[2,0,:,:], tit='Mar')

ax1=plt.subplot(434)
ax1=base(D1[3,0,:,:], tit='Apr')

ax1=plt.subplot(435)
ax1=base(D1[4,0,:,:], tit='May')

ax1=plt.subplot(436)
ax1=base(D1[5,0,:,:], tit='Jun')

ax1=plt.subplot(437)
ax1=base(D1[6,0,:,:], tit='Jul')

ax1=plt.subplot(438)
ax1=base(D1[7,0,:,:], tit='Aug')

ax1=plt.subplot(439)
ax1=base(D1[8,0,:,:], tit='Sep')

ax1=plt.subplot(4,3,10)
ax1=base(D1[9,0,:,:], tit='Oct')

ax1=plt.subplot(4,3,11)
ax1=base(D1[10,0,:,:], tit='Nov')

ax1=plt.subplot(4,3,12)
ax1=base(D1[11,0,:,:], tit='Dec')