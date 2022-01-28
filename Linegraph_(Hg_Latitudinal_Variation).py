#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May  2 10:50:22 2017

THIS SCRIPT PLOTS A GRAPH OF THE LATITUDINAL VARIATION IN ATMOSPHERIC HG CONCENTRATIONS

THERE ARE 3 GRAPH OPTIONS:
    1) GLOBAL SCALE
    2) SOUTHERN HEMISPHERE SCALE
    3) STUDY AREA SCALE

@author: ncp532
"""
# File system packages
from netCDF4 import Dataset				# function used to open single netcdf file
from netCDF4 import MFDataset				# function used to open multiple netcdf files

# Drawing packages
import matplotlib.pyplot as plt             # import package as shorter nickname
import matplotlib.dates as mdates           # 
import matplotlib.ticker as ticker

# Data handing packages
import numpy as np                          # import package as shorter nickname - Numpy as great at handling multidimensional data arrays.
import pandas as pd
from scipy import signal, stats

# Date and Time handling package
from datetime import datetime,timedelta		# functions to handle date and time

#------------------------------------------------------------------------------
# Define the dataset
dataset1 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Default_Files/trac_avg.201*.nc')    # Default
dataset2 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_InvOcean_Files/trac_avg.201*.nc')   # InvOcean
dataset7 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Horowitz/Offline_MITgcm_ocean*.nc') # NewChem/NewOcean                    
dataset8 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Horowitz/SlabOcean.201*.nc')        # NewChem/SlabOcean
dataset9 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Horowitz/AIRDENS.2010*.nc')         # AIRDENS

# OBSERVATION
# Retrieve the observation data
yCG = np.genfromtxt("/Users/ncp532/Documents/Hg_data/Hg0SiteDaily_CG_2013-2016_Daily_YearlyMean.csv", delimiter=",", skip_header=1, usecols=(1), invalid_raise=False) # Hg Concentration
yGP = np.genfromtxt("/Users/ncp532/Documents/Hg_data/Hg0SiteDaily_GP_2013-2016_Daily_YearlyMean.csv", delimiter=",", skip_header=1, usecols=(1), invalid_raise=False) # Hg Concentration
yGL = np.genfromtxt("/Users/ncp532/Documents/Hg_data/Hg0SiteDaily_GL_2013-2016_Daily_YearlyMean.csv", delimiter=",", skip_header=1, usecols=(1), invalid_raise=False) # Hg Concentration
yAMS = np.genfromtxt("/Users/ncp532/Documents/Hg_data/Hg0SiteDaily_AMS_2013-2016_Daily_YearlyMean.csv", delimiter=",", skip_header=1, usecols=(1), invalid_raise=False) # Hg Concentration
yCP = np.genfromtxt("/Users/ncp532/Documents/Hg_data/Hg0SiteDaily_CP_2013-2016_Daily_YearlyMean.csv", delimiter=",", skip_header=1, usecols=(1), invalid_raise=False) # Hg Concentration

# Default
# Get the GEOS-Chem data for the stations
airdena=dataset1.variables['BXHGHT_S__AIRNUMDE'][:]                           # dry air number density in (molecules air/m3)
wet1=(dataset1.variables['IJ_AVG_S__Hg0'][:])*(1e-3)*200.59/(6.0221*1e23)*airdena*1e6 # Hg0
# Convert the dataset from ppt) to ng/m3 = ppt*(10^-12)*(molecular mass Hg/molar volume Hg)*(molar volume Hg/Avogadro's number)*(10^9)*(air density in molecules/m3)
lata=dataset1.variables['lat'][:]                             # latitude 

# InvOcean
# Get the GEOS-Chem data for the stations
airdenb=dataset2.variables['BXHGHT_S__AIRNUMDE'][:]                           # dry air number density in (molecules air/m3)
wet2=(dataset2.variables['IJ_AVG_S__Hg0'][:])*(1e-3)*200.59/(6.0221*1e23)*airdenb*1e6
wetHg2=(dataset1.variables['IJ_AVG_S__Hg2'][:])*(1e-3)*200.59/(6.0221*1e23)*airdenb*1e6
wetHgP=(dataset1.variables['IJ_AVG_S__HgP'][:])*(1e-3)*200.59/(6.0221*1e23)*airdenb*1e6

# NewChem/NewOcean
# Get the GEOS-Chem data for the stations
airden3=dataset9.variables['BXHGHT_S__AIRNUMDE'][:]                           # dry air number density in (molecules air/m3)
wet7=(dataset7.variables['IJ_AVG_S__Hg0'][:])*(1e-3)*200.59/(6.0221*1e23)*airden3 #cm3 to m3 is a factor of 1e-6 #(2.5*1e25)
# Convert the dataset from ppt) to ng/m3 = ppt*(10^-12)*(molecular mass Hg/molar volume Hg)*(molar volume Hg/Avogadro's number)*(10^9)*(air density in molecules/m3)

# NewChem/SlabOcean
# Get the GEOS-Chem data for the stations
wet8=(dataset8.variables['IJ_AVG_S__Hg0'][:])*(1e-3)*200.59/(6.0221*1e23)*airden3 #cm3 to m3 is a factor of 1e-6 #(2.5*1e25)
# Convert the dataset from ppt) to ng/m3 = ppt*(10^-12)*(molecular mass Hg/molar volume Hg)*(molar volume Hg/Avogadro's number)*(10^9)*(air density in molecules/m3)

#------------------------------------------------------------------------------
# Calculate the overall mean
#MODEL
OM_1 = wet2[:,0,:,130]   # Cape Grim (25,130)
OM_1Avg = np.mean(OM_1,axis=0)
OM_2 = wet2[:,0,:,124]   # Gunn Point (39,124)
OM_2Avg = np.mean(OM_2,axis=0)
OM_3 = wet2[:,0,:,132]   # Glenville (29,132)
OM_3Avg = np.mean(OM_3,axis=0)
OM_4 = wet2[:,0,:,103]   # Amsterdam Island (26,103)
OM_4Avg = np.mean(OM_4,axis=0)
OM_5 = wet2[:,0,:,79]   # Cape Point (28,79)
OM_5Avg = np.mean(OM_5,axis=0)
OM_6 = wet2[:,0,:,:]   # Hg0 Average
OM_6Ag = np.mean(OM_6,axis=0)
OM_6Avg = np.mean(OM_6Ag,axis=1)
OM_7 = wetHg2[:,0,:,:]   # HgII Average
OM_7Ag = np.mean(OM_7,axis=0)
OM_7Avg = np.mean(OM_7Ag,axis=1)
OM_8 = wetHgP[:,0,:,:]   # HgP Average
OM_8Ag = np.mean(OM_8,axis=0)
OM_8Avg = np.mean(OM_8Ag,axis=1)
#OM_2 = wet2[:,0,:,67]   # InvOcean
#OM_2Avg = np.mean(OM_2,axis=0)
#OM_3 = wet7[:,0,:,67] # NewChem/NewOcean [:,0,12,67]
#OM_3Avg = np.mean(OM_3,axis=0)
#OM_4 = wet8[:,0,:,67] # NewChem/SlabOcean
#OM_4Avg = np.mean(OM_4,axis=0)

#OBSERVATION
OM_CG = np.nanmean(yCG)             # Cape Grim
OM_GP = np.nanmean(yGP)             # Gunn Point
OM_GL = np.nanmean(yGL)             # Glenville
OM_AMS = np.nanmean(yAMS)           # Amsterdam Island
OM_CP = np.nanmean(yCP)             # Cape Point

#SIMULATED
S_CG = np.mean(wet2[:,0,25,130]) # Cape Grim
S_GP = np.mean(wet2[:,0,39,124]) # Gunn Point
#------------------------------------------------------------------------------
# Plot the axis for each graph
fig = plt.figure()
plt.subplots_adjust(hspace=0.5)

# Graph 1
ax=plt.subplot(211) # options graph 2 (vertical no, horizontal no, graph no)

# SIMULATIONS
# Plot the monthly mean for the location
plt.plot(lata, OM_1Avg, ls='--', color='red', linewidth=1, label="Cape Grim (Longitude: 144.6899$^\circ$)")               # Cape Grim
plt.plot(lata, OM_2Avg, ls='--', color='blue', linewidth=1, label="Gunn Point (Longitude: 131.0447$^\circ$)")             # Gunn Point
plt.plot(lata, OM_3Avg, ls='--', color='green', linewidth=1, label="Glenville (Longitude: 151.1018$^\circ$)")             # Glenville
plt.plot(lata, OM_4Avg, ls='--', color='darkorange', linewidth=1, label="Amsterdam Island (Longitude: 77.3423$^\circ$)") # Amsterdam Island
plt.plot(lata, OM_5Avg, ls='--', color='brown', linewidth=1, label="Cape Point (Longitude: 18.4897$^\circ$)")             # Cape Point
plt.plot(lata, OM_6Avg, ls='--', color='black', linewidth=1, label="Hg$^0$ Average")               # Average

# OBSERVATIONS
# Plot the observed mean for the location
plt.plot(-40.6832,OM_CG, "o", color='red', label="Cape Grim Observations \n  (Mean: "+str("%7.3f"%(OM_CG))+" ng/m$^3$)") # Cape Grim
plt.plot(-12.2491,OM_GP, "o", color='blue', label="Gunn Point Observations \n  (Mean: "+str("%7.3f"%(OM_GP))+" ng/m$^3$)") # Gunn Point
plt.plot(-32.4777,OM_GL, "o", color='green', label="Glenville Observations \n  (Mean: "+str("%7.3f"%(OM_GL))+" ng/m$^3$)") # Glenville
plt.plot(-37.4748,OM_AMS, "o", color='darkorange', label="Amsterdam Island Observations \n  (Mean: "+str("%7.3f"%(OM_AMS))+" ng/m$^3$)") # Amsterdam Island
plt.plot(-34.3535,OM_CP, "o", color='brown', label="Cape Point Observations \n  (Mean: "+str("%7.3f"%(OM_CP))+" ng/m$^3$)") # Cape Point
# Plot the observed gradient between GP and CG
x1, y1 = [-40.6832,-12.2491],[OM_CG,OM_GP]
grad1 = (OM_GP-OM_CG)/((-12.2491)-(-40.6832))
#plt.plot(x1, y1, ls='--', color='darkslategray', linewidth=2, label="Observed Australian Gradient:\n  ("+str("%7.4f"%(grad1))+" ng/m$^3$/degree)")

# Plot the simulated gradient between GP and CG
x2, y2 = [-40.6832,-12.2491],[S_CG,S_GP]
grad2 = (S_GP-S_CG)/((-12.2491)-(-40.6832))
#plt.plot(x2, y2, ls='--', color='darkgreen', linewidth=2, label="Simulated Australian Gradient:\n  ("+str("%7.4f"%(grad2))+" ng/m$^3$/degree)")

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
#ax.yaxis.set_major_locator(ticker.MultipleLocator(0.10)) #Study Area
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.01))
ax.set_xlim(-90, 90) # World
#ax.set_xlim(-90, 0) # SH
#ax.set_xlim(-50, 0) # Study Area 
#ax.set_ylim(0.75, 1.65) # World
#ax.set_ylim(0.78, 1.46) # Southern Hemisphere
#ax.set_ylim(0.80, 1.30) # Study Area

# Plot the axis labels, legend and title
plt.ylabel('Hg$^0$ (ng/m$^3$)', fontsize=10)
plt.xlabel('Latitude (degrees)', fontsize=10)
plt.title("Latitudinal variation of Hg$^0$ concentrations in the surface boundary layer", fontsize=15)
plt.legend(bbox_to_anchor=(1.06, 1), loc=2, borderaxespad=0.)

# Graph 2
ax1=plt.subplot(212) # options graph 2 (vertical no, horizontal no, graph no)
ax2 = ax1.twinx()

# SIMULATIONS
# Plot the monthly mean for the location
a, = ax1.plot(lata, OM_6Avg, ls='--', color='black', linewidth=1, label="Hg$^0$ Average")      # Hg0
b, = ax2.plot(lata, OM_7Avg, ls='--', color='cyan', linewidth=1, label="Hg$^I$$^I$ Average") # HgII
c, = ax2.plot(lata, OM_8Avg, ls='--', color='magenta', linewidth=1, label="Hg$^P$ Average")    # HgP
p = [a, b, c]

plt.margins(1, 0.1)
#plt.xticks(rotation=15)
#ax.yaxis.set_ticks_position('both')

# SET AX1
ax1.xaxis.set_major_locator(ticker.MultipleLocator(10))
ax1.xaxis.set_minor_locator(ticker.MultipleLocator(5))
#ax1.yaxis.set_major_locator(ticker.MultipleLocator(0.1)) #World/SH
ax1.yaxis.set_major_locator(ticker.MultipleLocator(0.10)) #Study Area
ax1.yaxis.set_minor_locator(ticker.MultipleLocator(0.01))
ax1.set_xlim(-90, 90) # World
#ax1.set_xlim(-90, 0) # SH
#ax1.set_xlim(-50, 0) # Study Area
#ax1.set_ylim(0.75, 1.65) # World
#ax1.set_ylim(0.78, 1.22) # Southern Hemisphere
#ax1.set_ylim(0.95, 1.35) # Study Area

# SET AX2
ax2.yaxis.set_major_locator(ticker.MultipleLocator(0.01)) #World/SH
#ax2.yaxis.set_major_locator(ticker.MultipleLocator(0.001)) #Study Area
ax2.yaxis.set_minor_locator(ticker.MultipleLocator(0.001))
ax2.set_ylim(-0.001, 0.052) # World
#ax2.set_ylim(-0.001, 0.052) # Southern Hemisphere
#ax2.set_ylim(-0.001, 0.011) # Study Area

# Plot the axis labels, legend and title
ax1.set_ylabel('Hg$^0$ (ng/m$^3$)', fontsize=10)
ax1.set_xlabel('Latitude (degrees)', fontsize=10)
ax2.set_ylabel('Hg$^I$$^I$ and Hg$^P$ (ng/m$^3$)', fontsize=10)
plt.title("Latitudinal variation of Hg concentrations in the surface boundary layer", fontsize=15)
plt.legend(p, [p_.get_label() for p_ in p], bbox_to_anchor=(1.06, 1), loc=2, borderaxespad=0.)
#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

