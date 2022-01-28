#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May  2 10:50:22 2017

@author: ncp532
"""
# File system packages
from netCDF4 import Dataset				# function used to open single netcdf file
from netCDF4 import MFDataset				# function used to open multiple netcdf files

# Drawing packages
import matplotlib.pyplot as plt             # import package as shorter nickname
from matplotlib.cm import get_cmap
import matplotlib.ticker as ticker
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from wrf import to_np, getvar, CoordPair, vertcross

# Data handing packages
import numpy as np                          # import package as shorter nickname - Numpy as great at handling multidimensional data arrays.
import pandas as pd
from scipy import signal, stats

# Date and Time handling package
from datetime import datetime,timedelta		# functions to handle date and time

#------------------------------------------------------------------------------
# Define the dataset
dataset1 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Default_Files/trac_avg.201*.nc')       # Default
dataset2 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_InvOcean_Files/trac_avg.201*.nc')      # InvOcean

# Default
# Get the GEOS-Chem data for the stations
airdena=dataset1.variables['BXHGHT_S__AIRNUMDE'][:]                           # dry air number density in (molecules air/m3)
wet1=(dataset1.variables['IJ_AVG_S__Hg0'][:])*(1e-3)*200.59/(6.0221*1e23)*airdena*1e6 # Hg0
# Convert the dataset from ppt) to ng/m3 = ppt*(10^-12)*(molecular mass Hg/molar volume Hg)*(molar volume Hg/Avogadro's number)*(10^9)*(air density in molecules/m3)
timea=dataset1.variables['time'][:]                           # in minutes
leva=dataset1.variables['lev'][:]                             # atmospheric level (0 to 46)
lata=dataset1.variables['lat'][:]                             # latitude 

# InvOcean
# Get the GEOS-Chem data for the stations
airdenb=dataset2.variables['BXHGHT_S__AIRNUMDE'][:]                           # dry air number density in (molecules air/m3)
wet2=(dataset2.variables['IJ_AVG_S__Hg0'][:])*(1e-3)*200.59/(6.0221*1e23)*airdenb*1e6
wetHg2=(dataset1.variables['IJ_AVG_S__Hg2'][:])*(1e-3)*200.59/(6.0221*1e23)*airdenb*1e6
wetHgP=(dataset1.variables['IJ_AVG_S__HgP'][:])*(1e-3)*200.59/(6.0221*1e23)*airdenb*1e6
timeb=dataset2.variables['time'][:]                           # in minutes
levb=dataset2.variables['lev'][:]                             # atmospheric level (0 to 46)
latb=dataset2.variables['lat'][:]                             # latitude 

#------------------------------------------------------------------------------
# Calculate the overall mean

OM_1 = wet2[:,:,:,:]   # Hg0 Average
OM_1Ag = np.mean(OM_1,axis=0)
OM_1Avg = np.mean(OM_1Ag,axis=2)
#------------------------------------------------------------------------------
# Plot the axis for each graph
fig = plt.figure()
plt.subplots_adjust(hspace=0.5)

# Graph 1
ax=plt.subplot(111) # options graph 2 (vertical no, horizontal no, graph no)

wspd_contours = ax.contourf(to_np(wspd_cross), cmap=get_cmap("jet"))

# Add the color bar
plt.colorbar(wspd_contours, ax=ax)












# SIMULATIONS
# Plot the monthly mean for the location
plt.plot(lata, OM_1Avg, ls='--', color='red', linewidth=1, label="Cape Grim")               # Cape Grim
plt.plot(lata, OM_2Avg, ls='--', color='blue', linewidth=1, label="Gunn Point")             # Gunn Point
plt.plot(lata, OM_3Avg, ls='--', color='green', linewidth=1, label="Glenville")             # Glenville
plt.plot(lata, OM_4Avg, ls='--', color='darkorange', linewidth=1, label="Amsterdam Island") # Amsterdam Island
plt.plot(lata, OM_5Avg, ls='--', color='yellow', linewidth=1, label="Cape Point")           # Cape Point
plt.plot(lata, OM_6Avg, ls='--', color='black', linewidth=1, label="Hg$^0$ Average")               # Average

#plt.plot(lata, OM_1Avg, ls='--', color='red', linewidth=1, label="Default")                  # Normal
#plt.plot(lata, OM_2Avg, ls='--', color='blue', linewidth=1, label="InvOcean")                # InvOcean
#plt.plot(latg, OM_3Avg, ls='--', color='green', linewidth=1, label="NewChem/NewOcean")       # NewChem/NewOcean
#plt.plot(latg, OM_4Avg, ls='--', color='darkorange', linewidth=1, label="NewChem/SlabOcean") # NewChem/SlabOcean

plt.margins(1, 0.1)
#plt.xticks(rotation=15)
ax.yaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.1)) #World/SH
#ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05)) #Study Area
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.01))
ax.set_xlim(-90, 90) # World
ax.set_ylim(0.75, 1.65) # World


# Plot the axis labels, legend and title
plt.ylabel('Hg$^0$ (ng/m$^3$)', fontsize=10)
plt.xlabel('Latitude (degrees)', fontsize=10)
plt.title("Latitudinal variation of Hg$^0$ concentrations in the surface boundary layer", fontsize=15)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

# Graph 2
#ax=plt.subplot(212) # options graph 2 (vertical no, horizontal no, graph no)



