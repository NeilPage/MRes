#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May  2 10:50:22 2017

THIS SCRIPT PLOTS A SCATTERPLOT OF GEM CONCENTRATIONS AGAINST RELATIVE HUMIDITY

GEOS-CHEM COORDINATES FOR CAPE GRIM ISLAND ARE [:,0,6,56]

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
dataset1 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Default_Files/coards.201*.nc')  # Default
dataset2 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_InvOcean_Files/coards.201*.nc') # InvOcean

# OBSERVATION
# Retrieve the observation data
x1 = np.genfromtxt("/Users/ncp532/Documents/Hg_data/Hg0SiteDaily_AMS_2013-2016_Daily_YearlyMean.csv", delimiter=",", skip_header=1, usecols=(0), invalid_raise=False) # day number
y1 = np.genfromtxt("/Users/ncp532/Documents/Hg_data/Hg0SiteDaily_AMS_2013-2016_Daily_YearlyMean.csv", delimiter=",", skip_header=1, usecols=(1), invalid_raise=False) # Hg Concentration

# Default
# Get the GEOS-Chem data for the stations
airdena=dataset1.variables['TIME_SER__AIRDEN'][:,0,6,56]                           # dry air number density in (molecules air/m3)
wet1=(dataset1.variables['IJ_AVG_S__Hg0'][:,0,6,56])*(1e-3)*200.59/(6.0221*1e23)*airdena*1e6 #cm3 to m3 is a factor of 1e-6 #(2.5*1e25)
# Convert the dataset from ppt) to ng/m3 = ppt*(10^-12)*(molecular mass Hg/molar volume Hg)*(molar volume Hg/Avogadro's number)*(10^9)*(air density in molecules/m3)
timea=dataset1.variables['time'][:]                           # in minutes

# InvOcean
# Get the GEOS-Chem data for the stations
airdenb=dataset2.variables['TIME_SER__AIRDEN'][:,0,6,56] # Air Density
rh=dataset2.variables['TIME_SER__RH'][:,0,6,56] # Relative Humidity 
wet2=(dataset2.variables['IJ_AVG_S__Hg0'][:,0,6,56])*(1e-3)*200.59/(6.0221*1e23)*airdenb*1e6
#------------------------------------------------------------------------------
# DATETIME

# Default
# Convert GEOS-Chem time into standard time (GEOS-Chem time is hours since 1985/01/01/0000)
d0=datetime(1985,1,1,0)
datea = []
for t in timea:
    hrs=timedelta(hours=int(t))
    datea.append(d0+hrs)
    
#------------------------------------------------------------------------------
# keep only the Good Values in the observations (remove outliers)
OM_OBS = np.nanmean(y1)  # mean of Observations
OSTD_OBS = np.nanstd(y1) # Std of Observations
STDDN = (OM_OBS - 2*OSTD_OBS)
STDUP = (OM_OBS + 2*OSTD_OBS)

# 1st filter 
goodvals =np.where([i > (OM_OBS - 2*OSTD_OBS) for i in y1])[0] #(x > mean - 2 * sd)
Y1_filtered=[y1[i] for i in goodvals]
X1clean =[wet1[i] for i in goodvals]
X2clean =[wet2[i] for i in goodvals]
rhclean =[rh[i] for i in goodvals]
dateclean =[datea[i] for i in goodvals]

# 2nd filter
goodvals =np.where([i < (OM_OBS + 2*OSTD_OBS) for i in Y1_filtered])[0] #(x < mean + 2 * sd)
Y1clean =np.array([Y1_filtered[i] for i in goodvals])
X1clean =np.array([X1clean[i] for i in goodvals])
X2clean =np.array([X2clean[i] for i in goodvals])
rhclean =np.array([rhclean[i] for i in goodvals])
dateclean =np.array([dateclean[i] for i in goodvals])

STD_DN=np.empty(len(y1))
STD_DN.fill(STDDN)
STD_UP=np.empty(len(y1))
STD_UP.fill(STDUP)
#------------------------------------------------------------------------------
# Calculate the overall mean
OM_1 = np.nanmean(wet1)   # Default
OM_2 = np.nanmean(wet2)   # InvOcean
OM_OBS = np.nanmean(y1)             # Observations
#------------------------------------------------------------------------------
# Calculate the Coefficient of Correlation (r)

# Default
r_row1, p_value1 = stats.pearsonr(rhclean,X1clean)
slope1, intercept1, r1, p1, std_err1= stats.linregress(rhclean,X1clean)

# InvOcean
r_row2, p_value2 = stats.pearsonr(rhclean,X2clean)
slope2, intercept2, r2, p2, std_err2= stats.linregress(rhclean,X2clean)

# Observations
r_row3, p_value3 = stats.pearsonr(rhclean,Y1clean)
slope3, intercept3, r3, p3, std_err3= stats.linregress(rhclean,Y1clean)
#------------------------------------------------------------------------------
# Plot the axis for each graph
fig = plt.figure()
plt.subplots_adjust(hspace=0.5)

# Graph 1
ax=plt.subplot(111) # options graph 2 (vertical no, horizontal no, graph no)

# OBSERVATIONS
line1, = plt.plot(rh,y1, "o",color='black', linewidth=1, label="Observations:")

# SIMULATIONS
line2, = plt.plot(rh, wet1, "o", color='red', label="Default:")
line3, = plt.plot(rh, wet2, "o", color='blue', label="InvOcean") 

# Plot the regression line
line4, = plt.plot(rh, intercept1 + slope1*rh, color='red', label="Default:\n (slope: "+str("%7.4f"%(p_value1))+"$\pm$ "+str("%7.4f"%(std_err1))+"%, r: "+str("%7.4f"%(r_row1))+")")
line5, = plt.plot(rh, intercept2 + slope2*rh, color='blue',label="InvOcean:\n (slope: "+str("%7.4f"%(p_value2))+"$\pm$ "+str("%7.4f"%(std_err2))+"%, r: "+str("%7.4f"%(r_row2))+")")
line6, = plt.plot(rh, intercept3 + slope3*rh, color='black', label="Observations:\n (slope: "+str("%7.4f"%(p_value3))+"$\pm$ "+str("%7.4f"%(std_err3))+"%, r: "+str("%7.4f"%(r_row3))+")")
line7 = plt.plot(rh, STD_DN, color='black', alpha=0.5)
line8 = plt.plot(rh, STD_UP, color='black', alpha=0.5)

ax.yaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.01))
ax.set_xlim(49.75, 100.25) # Daily
ax.set_ylim(0.542, 1.558)  # Daily

# Plot the axis labels, legend and title
plt.ylabel('Hg$^0$ (ng/m$^3$)', fontsize=10) # Daily and Monthly
plt.xlabel('Relative Humidity (%)', fontsize=10)
plt.title("Relationship between Hg$^0$ and Relative Humidity (Cape Grim, Tas)", fontsize=15)
#plt.legend(handles=[line6,line4,line5], loc=1, handlelength=0, handletextpad=0)
plt.annotate("Observations:\n (slope: "+str("%7.4f"%(slope3))+"$\pm$ "+str("%7.4f"%(std_err3))+" ng m$^-$$^3$ %$^-$$^1$, r: "+str("%7.4f"%(r3))+")", xy=(50.1,1.51), color='black', fontweight='bold')
plt.annotate("Default:\n (slope: "+str("%7.4f"%(slope1))+"$\pm$ "+str("%7.4f"%(std_err1))+" ng m$^-$$^3$ %$^-$$^1$, r: "+str("%7.4f"%(r1))+")", xy=(50.1,1.46), color='red', fontweight='bold')
plt.annotate("InvOcean:\n (slope: "+str("%7.4f"%(slope2))+"$\pm$ "+str("%7.4f"%(std_err2))+" ng m$^-$$^3$ %$^-$$^1$, r: "+str("%7.4f"%(r2))+")", xy=(50.1,1.41), color='blue', fontweight='bold')
plt.annotate("Upper limit: "+str("%7.4f"%(STDUP))+" ng m$^-$$^3$ %$^-$$^1$", xy=(50.1,1.39), color='black', fontweight='bold', alpha=0.5)
plt.annotate("Lower limit: "+str("%7.4f"%(STDDN))+" ng m$^-$$^3$ %$^-$$^1$", xy=(50.1,1.36), color='black', fontweight='bold', alpha=0.5)

