#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 10:50:58 2017

THIS SCRIPT CALCULATES THE MONTHLY MEAN AND STANDARD DEVIATION
NOTE: IT DOESN'T PRODUCE A GRAPH

@author: ncp532
"""
# File system packages
from netCDF4 import Dataset				# function used to open single netcdf file
from netCDF4 import MFDataset				# function used to open multiple netcdf files

# Drawing packages
import matplotlib.pyplot as plt             # import package as shorter nickname
import matplotlib.dates as mdates           # 

# Data handing packages
import numpy as np                          # import package as shorter nickname - Numpy as great at handling multidimensional data arrays.
import pandas as pd
from scipy import signal, stats

# Date and Time handling package
from datetime import datetime,timedelta		# functions to handle date and time

#------------------------------------------------------------------------------
# Define the dataset
dataset1 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Default_Files/coards.201*.nc')      # Default
dataset2 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Optimal_Files/coards.201*.nc')      # InvOcean
dataset3 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Anthro_Files/coards.201*.nc')       # Anthro
dataset4 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Biomass_Files/coards.201*.nc')      # Biomass
dataset5 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Geogenic_Files/coards.201*.nc')     # Geogenic
dataset6 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Soil_Files/coards.201*.nc')         # Soil
dataset7 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Horowitz/Offline_MITgcm_ocean*.nc') # NewChem/NewOcean                    
dataset8 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Horowitz/SlabOcean.201*.nc')        # NewChem/SlabOcean
dataset9 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Horowitz/AIRDENS.2010*.nc')         # AIRDENS

# OBSERVATION
# Retrieve the observation data
x1 = np.genfromtxt("/Users/ncp532/Desktop/Hg_data/Hg0SiteDaily_AMS_2013-2016_Daily_YearlyMean.csv", delimiter=",", skip_header=1, usecols=(0), invalid_raise=False) # day number
y1 = np.genfromtxt("/Users/ncp532/Desktop/Hg_data/Hg0SiteDaily_AMS_2013-2016_Daily_YearlyMean.csv", delimiter=",", skip_header=1, usecols=(1), invalid_raise=False) # Hg Concentration

# Default
# Get the GEOS-Chem data for the stations
airdena=dataset1.variables['TIME_SER__AIRDEN'][:]                           # dry air number density in (molecules air/m3)
wet1=(dataset1.variables['IJ_AVG_S__Hg0'][:])*(1e-3)*200.59/(6.0221*1e23)*airdena*1e6 # Hg0
# Convert the dataset from ppt) to ng/m3 = ppt*(10^-12)*(molecular mass Hg/molar volume Hg)*(molar volume Hg/Avogadro's number)*(10^9)*(air density in molecules/m3)
timea=dataset1.variables['time'][:]                           # in minutes
leva=dataset1.variables['lev'][:]                             # atmospheric level (0 to 46)
lata=dataset1.variables['lat'][:]                             # latitude 

# InvOcean
# Get the GEOS-Chem data for the stations
airdenb=dataset2.variables['TIME_SER__AIRDEN'][:]                           # dry air number density in (molecules air/m3)
wet2=(dataset2.variables['IJ_AVG_S__Hg0'][:])*(1e-3)*200.59/(6.0221*1e23)*airdenb*1e6
#(2.5*1e25)
timeb=dataset2.variables['time'][:]                           # in minutes
#levb=dataset2.variables['lev'][:]                             # atmospheric level (0 to 46)
latb=dataset2.variables['lat'][:]                             # latitude 

# Anthro
# Get the GEOS-Chem data for the stations
airdenc=dataset3.variables['TIME_SER__AIRDEN'][:]                           # dry air number density in (molecules air/m3)
wet3=(dataset3.variables['IJ_AVG_S__Hg0'][:])*(1e-3)*200.59/(6.0221*1e23)*airdenc*1e6
#(2.5*1e25)
timec=dataset3.variables['time'][:]                           # in minutes
#levc=dataset3.variables['lev'][:]                             # atmospheric level (0 to 46)
latc=dataset3.variables['lat'][:]                             # latitude 

# Biomass
# Get the GEOS-Chem data for the stations
airdend=dataset4.variables['TIME_SER__AIRDEN'][:]                           # dry air number density in (molecules air/m3)
wet4=(dataset4.variables['IJ_AVG_S__Hg0'][:])*(1e-3)*200.59/(6.0221*1e23)*airdend*1e6
#(2.5*1e25)
timed=dataset4.variables['time'][:]                           # in minutes
#levd=dataset4.variables['lev'][:]                             # atmospheric level (0 to 46)
latd=dataset4.variables['lat'][:]                             # latitude 

# Geogenic
# Get the GEOS-Chem data for the stations
airdene=dataset5.variables['TIME_SER__AIRDEN'][:]                           # dry air number density in (molecules air/m3)
wet5=(dataset5.variables['IJ_AVG_S__Hg0'][:])*(1e-3)*200.59/(6.0221*1e23)*airdene*1e6
#(2.5*1e25)
timee=dataset5.variables['time'][:]                           # in minutes
#leve=dataset5.variables['lev'][:]                             # atmospheric level (0 to 46)
late=dataset5.variables['lat'][:]                             # latitude 

# Soil
# Get the GEOS-Chem data for the stations
airdenf=dataset6.variables['TIME_SER__AIRDEN'][:]                           # dry air number density in (molecules air/m3)
wet6=(dataset6.variables['IJ_AVG_S__Hg0'][:])*(1e-3)*200.59/(6.0221*1e23)*airdenf*1e6
#(2.5*1e25)
timef=dataset6.variables['time'][:]                           # in minutes
#levf=dataset6.variables['lev'][:]                             # atmospheric level (0 to 46)
latf=dataset6.variables['lat'][:]                             # latitude 

# NewChem/NewOcean
# Get the GEOS-Chem data for the stations
airden3=dataset9.variables['BXHGHT_S__AIRNUMDE'][:,0,12,67]                           # dry air number density in (molecules air/m3)
wet7=(dataset7.variables['IJ_AVG_S__Hg0'][:,0,12,67])*(1e-3)*200.59/(6.0221*1e23)*airden3 #cm3 to m3 is a factor of 1e-6 #(2.5*1e25)
# Convert the dataset from ppt) to ng/m3 = ppt*(10^-12)*(molecular mass Hg/molar volume Hg)*(molar volume Hg/Avogadro's number)*(10^9)*(air density in molecules/m3)
timeg=dataset7.variables['time'][:]                           # in minutes
levg=dataset7.variables['lev'][:]                             # atmospheric level (0 to 46)
latg=dataset7.variables['lat'][:]                             # latitude 

# NewChem/SlabOcean
# Get the GEOS-Chem data for the stations
wet8=(dataset8.variables['IJ_AVG_S__Hg0'][:,0,12,67])*(1e-3)*200.59/(6.0221*1e23)*airden3 #cm3 to m3 is a factor of 1e-6 #(2.5*1e25)
# Convert the dataset from ppt) to ng/m3 = ppt*(10^-12)*(molecular mass Hg/molar volume Hg)*(molar volume Hg/Avogadro's number)*(10^9)*(air density in molecules/m3)
timeh=dataset8.variables['time'][:]                           # in minutes
levh=dataset8.variables['lev'][:]                             # atmospheric level (0 to 46)
lath=dataset8.variables['lat'][:]                             # latitude
#------------------------------------------------------------------------------
# Default
# Convert GEOS-Chem time into standard time (GEOS-Chem time is hours since 1985/01/01/0000)
d0=datetime(1985,1,1,0)
datea = []
for t in timea:
    hrs=timedelta(hours=int(t))
    datea.append(d0+hrs)
    
# Default
# Convert datea to monthsa
monthsa = [int(d.strftime('%Y%m')) for d in datea]
monthsa = set(monthsa)
monthsa = np.sort([str(m)+"01" for m in monthsa])
monthsa = [datetime.strptime(m,"%Y%m%d") for m in monthsa]

#------------------------------------------------------------------------------    
# Calculate the monthly mean
def monthlyM(x, date):
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
wet1_Mavg, date_avg=monthlyM(wet1[:,0,6,56],datea) # Default
wet2_Mavg, date_avg=monthlyM(wet2[:,0,6,56],datea) # InvOcean
wet3_Mavg, date_avg=monthlyM(wet3[:,0,6,56],datea) # Anthro
wet4_Mavg, date_avg=monthlyM(wet4[:,0,6,56],datea) # Biomass
wet5_Mavg, date_avg=monthlyM(wet5[:,0,6,56],datea) # Geogenic
wet6_Mavg, date_avg=monthlyM(wet6[:,0,6,56],datea) # Soil
wet7_Mavg, date_avg=monthlyM(y1,datea)             # Observations

#------------------------------------------------------------------------------
# Calculate the standard deviation
def monthlySTD(x, date):
    df = pd.DataFrame({'X':x}, index=date) 
    df = df.resample('MS').stdev()
    #Reset the index
    df =df.reset_index()
    #extract the values
    x=df['X']
    date=df['index']  
    #convert the pandas series date to list
    date = date.tolist()
    return x,date 
# Cape Grim
wet1_STDavg, date_avg=monthlySTD(wet1[:,0,6,56],datea) # Default
wet2_STDavg, date_avg=monthlySTD(wet2[:,0,6,56],datea) # InvOcean
wet3_STDavg, date_avg=monthlySTD(wet3[:,0,6,56],datea) # Anthro
wet4_STDavg, date_avg=monthlySTD(wet4[:,0,6,56],datea) # Biomass
wet5_STDavg, date_avg=monthlySTD(wet5[:,0,6,56],datea) # Geogenic
wet6_STDavg, date_avg=monthlySTD(wet6[:,0,6,56],datea) # Soil
wet7_STDavg, date_avg=monthlySTD(y1,datea)             # Observations

print "Done"
