#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 16:19:13 2017

THIS SCRIPT PLOTS A MAP OF THE REGIONAL GLOBAL HG0 CONCENTRATION

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
import numpy as np					    # import package as shorter nickname - Numpy as great at handling multidimensional data arrays.

# Define the dataset
#dataset1 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Default_Files/trac_avg.201*.nc')	#define the dataset
dataset1 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_InvOcean_Files/trac_avg.201*.nc')	#define the dataset

#MAP 1
# Get the GEOS-Chem data for the stations
#airden1=dataset1.variables['BXHGHT_S__AIRNUMDE'][:] # for Hg0
airden1=dataset1.variables['BXHGHT_S__AIRNUMDE'][:,0,:,:] # for anthro, biomass, soil emissions
#wet1=(dataset1.variables['IJ_AVG_S__Hg0'][:])*(1e-3)*200.59/(6.0221*1e23)*airden1*1e6 # Hg0
#wet1=(dataset1.variables['HG_SRCE__Hg0_an'][:])*(1e-3)*200.59/(6.0221*1e23)*airden1*1e6    # anthro emissions
#wet1=(dataset1.variables['HG_SRCE__Hg_bb'][:])*(1e-3)*200.59/(6.0221*1e23)*airden1*1e6   # biomass emissions
wet1=(dataset1.variables['HG_SRCE__Hg_so'][:])*(1e-3)*200.59/(6.0221*1e23)*airden1*1e6   # soil emissions
time1=dataset1.variables['time'][:]                           # in minutes
lat1=dataset1.variables['lat'][:]                             # in degrees
lon1=dataset1.variables['lon'][:]                             # in degrees

# Average the time for wet1
avg=wet1
wet1_avg=np.mean(avg,axis=0)

# Convert GEOS-Chem time into standard time (GEOS-Chem time is hours since 1985/01/01/0000)
d0=datetime(1985,1,1,0)
date1 = []
for t in time1:
    hrs=timedelta(hours=int(t))
    date1.append(d0+hrs)


def plotkg(x, blah):
    lons, lats = np.meshgrid(lon1,lat1)
    m=Basemap(llcrnrlat=-47,  urcrnrlat=-2,        # set the latitude limits
              llcrnrlon=10, urcrnrlon=160,         # set the longitude limits
              resolution='c',projection='merc')

    # np.transpose because it mixes the lat lon order
    cs=m.pcolor(lons, lats, x, latlon=True, vmin=-100, vmax=450, cmap=cm.plasma)   # define the colorbar attributes (colors used, add "_r" to flip the color scale)
            #vmin=1e13,vmax=1e18,norm = LogNorm())
    # plt.clim(-1.5e7, 1.5e7)
    m.drawcoastlines()
    
    if blah=='Y':
        cb=m.colorbar(cs,"bottom",size="10%", pad="20%")
        cb.set_label(r"$ng/m^3$", fontsize =10)
        cb.ax.tick_params(labelsize=10)
    parallels = np.arange(-45,-2, 5.)
    # Plot parallels and meridians
    m.drawparallels(parallels,labels=[True,False,True,False], fontsize=10)
    meridians = np.arange(15,160,15.)
    m.drawmeridians(meridians,labels=[True,False,False,True], fontsize=10)
    # Plot location markers
    lons = [144.6832, 131.0447, 151.1018, 18.4897, 77.3423]
    lats = [-40.6832, -12.2491, -32.4777, -34.3535, -37.4748]
    #obs = [0.909, 0.974, 0.808]
    x,y = m(lons, lats)
    m.scatter(x, y, c='green', vmin=-100, vmax=450, cmap=cm.plasma, edgecolors='black', s=100)
    
    
    # Plot location labels
    labels = ['Cape Grim', 'Gunn Point', 'Glenville', 'Cape Point', 'Amsterdam Island']
    x_offsets = [-400000, -50000, -300000, -30000, -30000]
    y_offsets = [100000, -150000, 100000, -250000, 100000]
    for label, xpt, ypt, x_offset, y_offset in zip(labels, x, y, x_offsets, y_offsets):
        plt.text(xpt+x_offset, ypt+y_offset, label)
    # labels = [left,right,top,bottom]
    plt.xlabel('Longitude [$^\circ$]', fontsize=10)
    plt.ylabel('Latitude [$^\circ$]', fontsize=10)
    plt.title("Annual Hg emissions from soil sources", fontsize=15)
    #ax1.xaxis.set_label_coords(0.5, -0.1)
    #ax1.yaxis.set_label_coords(-0.07, 0.5)
    return m

fig = plt.figure()

fig=plt.subplot(111) # options graph 1 (vertical no, horizontal no, graph no)
#ax1 = plotkg(wet1_avg[0,:,:])
ax1 = plotkg(wet1_avg[:,:],'Y')
