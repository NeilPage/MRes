#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May  2 10:50:22 2017

THIS SCRIPT PLOTS A SCATTERPLOT OF SIMULATED GEM CONCENTRATIONS AGAINST OBSERVED GEM CONCENTRATIONS

GEOS-CHEM COORDINATES FOR CAPE GRIM AT 2x2.5 RESOLUTION ARE [:,0,6,56]

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
# CONVERT DATETIME TO MONTHTIME

# Default
# Convert datea to monthsa
monthsa = [int(d.strftime('%Y%m')) for d in datea]
monthsa = set(monthsa)
monthsa = np.sort([str(m)+"01" for m in monthsa])
monthsa = [datetime.strptime(m,"%Y%m%d") for m in monthsa]

#------------------------------------------------------------------------------
# CONVERT DATETIME TO YEARTIME

# Default
# Convert datea to yearsa
yearsa = [int(d.strftime('%Y')) for d in datea]
yearsa = set(yearsa)
yearsa = np.sort([str(y)+"01" for y in yearsa])
yearsa = [datetime.strptime(y,"%Y%m") for y in yearsa]

#------------------------------------------------------------------------------    
# Calculate the monthly mean
def monthly(x, date):
    df = pd.DataFrame({'X':x}, index=date) 
    df = df.resample('MS').mean()
    #Reset the index
    df =df.reset_index()
    #extract the values
    x=df['X']
    date=df['index']  
    #convert the pandas series date to list
    date = date.tolist()
    return x,date 
# Cape Grim
wet1_Mavg, date_avg=monthly(wet1,datea) # Default
wet2_Mavg, date_avg=monthly(wet2,datea) # InvOcean
obs_Mavg, date_avg=monthly(y1,datea)              # Observations
#------------------------------------------------------------------------------
# Calculate the annual mean
def yearly(x, date):
    df = pd.DataFrame({'X':x}, index=date) 
    df = df.resample('AS').mean()
    #Reset the index
    df =df.reset_index()
    #extract the values
    x=df['X']
    date=df['index']  
    #convert the pandas series date to list
    date = date.tolist()
    return x,date 

# Cape Grim
wet1_Yavg, date_avg=yearly(wet1,datea) # Default
wet2_Yavg, date_avg=yearly(wet2,datea) # InvOcean
obs_Yavg, date_avg=yearly(y1,datea)              # Observations
#------------------------------------------------------------------------------
# Calculate the overall mean
OM_1 = np.nanmean(wet1)   # Default
OM_2 = np.nanmean(wet2)   # InvOcean
OM_OBS = np.nanmean(y1)             # Observations
#------------------------------------------------------------------------------
# Calculate the seasonal variation
# need to average all the values for januray, feburary, etc...
# also plot the spatial standard deviation.

# Default
monthavg1=np.zeros([12]) # Set the size
monthvals1= [ dt.month for dt in monthsa ] # Extract the month values from the date 
for monthind in np.arange(12):
    dataind=monthvals1 == monthind+1
    monthavg1[monthind] = np.mean(wet1_Mavg[dataind]) # find the mean of months with same month value

# InvOcean
monthavg2=np.zeros([12]) # Set the size
monthvals2= [ dt.month for dt in monthsa ] # Extract the month values from the date 
for monthind in np.arange(12):
    dataind=monthvals2 == monthind+1
    monthavg2[monthind] = np.mean(wet2_Mavg[dataind]) # find the mean of months with same month value

# OBSERVATIONS
monthavg3=np.zeros([12]) # Set the size
monthvals3= [ dt.month for dt in monthsa ] # Extract the month values from the date 
for monthind in np.arange(12):
    dataind=monthvals3 == monthind+1
    monthavg3[monthind] = np.nanmean(obs_Mavg[dataind]) # find the mean of months with same month value
#------------------------------------------------------------------------------
# Calculate the decimal logarithm of the daily average
dl1 = np.log10(wet1) # Default
dl2 = np.log10(wet2) # InvOcean
dlObs = np.log10(y1) # Observations
#------------------------------------------------------------------------------
# Calculate the Coefficient of Correlation (r)

# DAILY AVERAGES
finiteY1mask = np.isfinite(y1) # Scan for NaN values
Y1clean = y1[finiteY1mask]      # Remove NaN values from the obs

# Default
X1clean = wet1[finiteY1mask] # Remove values for wet1 corresponding to obs NaN values
r_rowD1, p_valueD1 = stats.pearsonr(X1clean,Y1clean)
slopeD1, interceptD1, rD1, pD1, std_errD1= stats.linregress(X1clean,Y1clean)

# InvOcean
X2clean = wet2[finiteY1mask] # Remove values for wet2 corresponding to obs NaN values
r_rowD2, p_valueD2 = stats.pearsonr(X1clean,Y1clean)
slopeD2, interceptD2, rD2, pD2, std_errD2= stats.linregress(X2clean,Y1clean)

# MONTHLY AVERAGES
finiteY2mask = np.isfinite(obs_Mavg) # Scan for NaN values
Y2clean = obs_Mavg[finiteY2mask]      # Remove NaN values from the obs

# Default
X3clean = wet1_Mavg[finiteY2mask]
r_rowM1, p_valueM1 = stats.pearsonr(X3clean,Y2clean)
slopeM1, interceptM1, rM1, pM1, std_errM1= stats.linregress(X3clean,Y2clean)

# InvOcean
X4clean = wet2_Mavg[finiteY2mask]
r_rowM2, p_valueM2 = stats.pearsonr(X4clean,Y2clean)
slopeM2, interceptM2, rM2, pM2, std_errM2= stats.linregress(X4clean,Y2clean)

# Decimal Logarithm
finiteY3mask = np.isfinite(dlObs)
Y3clean = dlObs[finiteY3mask]

# Default
X5clean = dl1[finiteY3mask]
r_rowdl1, p_valuedl1 = stats.pearsonr(X5clean,Y3clean)
slopedl1, interceptdl1, rdl1, pdl1, std_errdl1= stats.linregress(X5clean,Y3clean)

# InvOcean
X6clean = dl2[finiteY3mask]
r_rowdl2, p_valuedl2 = stats.pearsonr(X6clean,Y3clean)
slopedl2, interceptdl2, rdl2, pdl2, std_errdl2= stats.linregress(X6clean,Y3clean)
#------------------------------------------------------------------------------
# Plot the axis for each graph
fig = plt.figure()
plt.subplots_adjust(hspace=0.5)

# Graph 1
ax=plt.subplot(111) # options graph 2 (vertical no, horizontal no, graph no)

# SIMULATIONS
# Plot the monthly mean for the location
#line2, = plt.plot(wet1_Mavg,obs_Mavg,"o", color='red', label="Default:")    
line2, = plt.plot(wet1,y1, "o", color='red', label="Default:")
#line2, = plt.plot(dl_1,dl_Obs, "o", color='red', label="Default:")
#line3, = plt.plot(wet2_Mavg,obs_Mavg, "o", color='blue', label="InvOcean")
line3, = plt.plot(wet2,y1, "o", color='blue', label="InvOcean") 
#line3, = plt.plot(dl_2,dl_Obs, "o", color='blue', label="InvOcean") 

# Plot the regression line
#line4, = plt.plot(wet1_Mavg, interceptM1 + slopeM1*wet1_Mavg, color='red', label="Default:\n (slope: "+str("%7.4f"%(p_valueM1))+"$\pm$ "+str("%7.4f"%(std_errM1))+"%, r: "+str("%7.4f"%(r_rowM1))+")")
line4, = plt.plot(wet1, interceptD1 + slopeD1*wet1, color='red', label="Default:\n (slope: "+str("%7.4f"%(p_valueD1))+"$\pm$ "+str("%7.4f"%(std_errD1))+"%, r: "+str("%7.4f"%(r_rowD1))+")")
#line4, = plt.plot(dl1, interceptdl1 + slopedl1*dl1, color='red', label="Default:\n (slope: "+str("%7.4f"%(p_valuedl1))+"$\pm$ "+str("%7.4f"%(std_errdl1))+"%, r: "+str("%7.4f"%(r_rowdl1))+")")
#line5, = plt.plot(wet2_Mavg, interceptM2 + slopeM2*wet2_Mavg, color='blue',label="InvOcean:\n (slope: "+str("%7.4f"%(p_valueM2))+"$\pm$ "+str("%7.4f"%(std_errM2))+"%, r: "+str("%7.4f"%(r_rowM2))+")")
line5, = plt.plot(wet2, interceptD2 + slopeD2*wet2, color='blue',label="InvOcean:\n (slope: "+str("%7.4f"%(p_valueD2))+"$\pm$ "+str("%7.4f"%(std_errD2))+"%, r: "+str("%7.4f"%(r_rowD2))+")")
#line5, = plt.plot(dl2, interceptdl2 + slopedl2*dl2, color='blue',label="InvOcean:\n (slope: "+str("%7.4f"%(p_valuedl2))+"$\pm$ "+str("%7.4f"%(std_errdl2))+"%, r: "+str("%7.4f"%(r_rowdl2))+")")

ax.yaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.05))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.01))
#ax.set_xlim(64.75, 80.25)  # Monthly
#ax.set_ylim(0.772, 1.258)  # Monthly
#ax.set_xlim(49.75, 100.25) # Daily
#ax.set_ylim(0.542, 1.458)  # Daily
#ax.set_xlim(49.75, 100.25)  # Decimal logarithm
#ax.set_ylim(-0.208, 0.208)  # Decimal logarithm

# Plot the axis labels, legend and title
#plt.ylabel('Hg$^0$ (ng/m$^3$)', fontsize=10) # Daily and Monthly
plt.ylabel('Observation Hg$^0$ (ng/m$^3$)', fontsize=10)        # Decimal Logarithm
plt.xlabel('Simulation Hg$^0$ (ng/m$^3$)', fontsize=10)
plt.title("Relationship between simulated and observed values of Hg$^0$ (Cape Grim, Tas)", fontsize=15)
#plt.legend(handles=[line6,line4,line5], loc=1, handlelength=0, handletextpad=0)
#plt.annotate("Default:\n (slope: "+str("%7.4f"%(slopeM1))+"$\pm$ "+str("%7.4f"%(std_errM1))+" ng m$^-$$^3$ %$^-$$^1$, r: "+str("%7.4f"%(rM1))+")", xy=(76.1,1.21), color='red', fontweight='bold')
#plt.annotate("Default:\n (slope: "+str("%7.4f"%(slopeD1))+"$\pm$ "+str("%7.4f"%(std_errD1))+" ng m$^-$$^3$ %$^-$$^1$, r: "+str("%7.4f"%(rD1))+")", xy=(50.1,1.36), color='red', fontweight='bold')
plt.annotate("Default:\n (slope: "+str("%7.4f"%(slopedl1))+"$\pm$ "+str("%7.4f"%(std_errdl1))+" %$^-$$^1$, r: "+str("%7.4f"%(rdl1))+")", xy=(50.0,0.155), color='red', fontweight='bold')
#plt.annotate("InvOcean:\n (slope: "+str("%7.4f"%(slopeM2))+"$\pm$ "+str("%7.4f"%(std_errM2))+" ng m$^-$$^3$ %$^-$$^1$, r: "+str("%7.4f"%(rM2))+")", xy=(76.1,1.185), color='blue', fontweight='bold')
#plt.annotate("InvOcean:\n (slope: "+str("%7.4f"%(slopeD2))+"$\pm$ "+str("%7.4f"%(std_errD2))+" ng m$^-$$^3$ %$^-$$^1$, r: "+str("%7.4f"%(rD2))+")", xy=(50.1,1.31), color='blue', fontweight='bold')
plt.annotate("InvOcean:\n (slope: "+str("%7.4f"%(slopedl2))+"$\pm$ "+str("%7.4f"%(std_errdl2))+" %$^-$$^1$, r: "+str("%7.4f"%(rdl2))+")", xy=(50.0,0.130), color='blue', fontweight='bold')

# Graph 2
#ax=plt.subplot(212) # options graph 2 (vertical no, horizontal no, graph no)

