#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May  2 10:50:22 2017

THIS SCRIPT PLOTS:
1) LINEGRAPHS OF DAILY/MONTHLY/SIM-MNB
2) LINEGRAPHS OF DAILY/VARIABILITY
3) SCATTERPLOT OF SIMULATED VS OBSERVED CONCENTRATIONS

The GEOS-Chem co-ordinates for Cape Grim at:
- 2x2.5 resolution are [:,0,6,56]
- 4x5 resolution are [:,0,12,67]

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
dataset3 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Anthro_Files/coards.201*.nc')   # Anthro
dataset4 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Biomass_Files/coards.201*.nc')  # Biomass
dataset5 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Geogenic_Files/coards.201*.nc') # Geogenic
dataset6 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Soil_Files/coards.201*.nc')     # Soil
#dataset7 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Bromine_x3/coards.201*.nc')    # Bromine_x3

# OBSERVATION
# Retrieve the observation data
x1 = np.genfromtxt("/Users/ncp532/Documents/Hg_data/Hg0SiteDaily_CG_2013-2016_Daily_YearlyMean.csv", delimiter=",", skip_header=1, usecols=(0), invalid_raise=False) # day number
y1 = np.genfromtxt("/Users/ncp532/Documents/Hg_data/Hg0SiteDaily_CG_2013-2016_Daily_YearlyMean.csv", delimiter=",", skip_header=1, usecols=(1), invalid_raise=False) # Hg Concentration

# Default
# Get the GEOS-Chem data for the stations
airdena=dataset1.variables['TIME_SER__AIRDEN'][:,0,6,56]                           # dry air number density in (molecules air/m3)
wet1=(dataset1.variables['IJ_AVG_S__Hg0'][:,0,6,56])*(1e-3)*200.59/(6.0221*1e23)*airdena*1e6 #cm3 to m3 is a factor of 1e-6 #(2.5*1e25)
# Convert the dataset from ppt) to ng/m3 = ppt*(10^-12)*(molecular mass Hg/molar volume Hg)*(molar volume Hg/Avogadro's number)*(10^9)*(air density in molecules/m3)
timea=dataset1.variables['time'][:]                           # in minutes

# InvOcean
# Get the GEOS-Chem data for the stations
airdenb=dataset2.variables['TIME_SER__AIRDEN'][:,0,6,56]                           # dry air number density in (molecules air/m3)
wet2=(dataset2.variables['IJ_AVG_S__Hg0'][:,0,6,56])*(1e-3)*200.59/(6.0221*1e23)*airdenb*1e6

# Anthro
# Get the GEOS-Chem data for the stations
airdenc=dataset3.variables['TIME_SER__AIRDEN'][:,0,6,56]                           # dry air number density in (molecules air/m3)
wet3=(dataset3.variables['IJ_AVG_S__Hg0'][:,0,6,56])*(1e-3)*200.59/(6.0221*1e23)*airdenc*1e6

# Biomass
# Get the GEOS-Chem data for the stations
airdend=dataset4.variables['TIME_SER__AIRDEN'][:,0,6,56]                           # dry air number density in (molecules air/m3)
wet4=(dataset4.variables['IJ_AVG_S__Hg0'][:,0,6,56])*(1e-3)*200.59/(6.0221*1e23)*airdend*1e6

# Geogenic
# Get the GEOS-Chem data for the stations
airdene=dataset5.variables['TIME_SER__AIRDEN'][:,0,6,56]                           # dry air number density in (molecules air/m3)
wet5=(dataset5.variables['IJ_AVG_S__Hg0'][:,0,6,56])*(1e-3)*200.59/(6.0221*1e23)*airdene*1e6

# Soil
# Get the GEOS-Chem data for the stations
airdenf=dataset6.variables['TIME_SER__AIRDEN'][:,0,6,56]                           # dry air number density in (molecules air/m3)
wet6=(dataset6.variables['IJ_AVG_S__Hg0'][:,0,6,56])*(1e-3)*200.59/(6.0221*1e23)*airdenf*1e6

# Bromine_x3
# Get the GEOS-Chem data for the stations
#airdeng=dataset7.variables['TIME_SER__AIRDEN'][:,0,6,56]                           # dry air number density in (molecules air/m3)
#wet7=(dataset7.variables['IJ_AVG_S__Hg0'][:,0,6,56])*(1e-3)*200.59/(6.0221*1e23)*airdeng*1e6
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
wet3_Mavg, date_avg=monthly(wet3,datea) # Anthro
wet4_Mavg, date_avg=monthly(wet4,datea) # Biomass
wet5_Mavg, date_avg=monthly(wet5,datea) # Geogenic
wet6_Mavg, date_avg=monthly(wet6,datea) # Soil
#wet7_Mavg, date_avg=monthly(wet7,datea) # Bromine_x3
obs_Mavg,  date_avg=monthly(y1,datea)             # observations
#------------------------------------------------------------------------------
# Format the monthly mean from (1 value per month) to (1 value per day)
# wet1_Mavg (Default)
df = pd.DataFrame((wet1_Mavg), index=date_avg) 
df = pd.DataFrame({'X':wet1}, index=datea) # Cape Grim
df = df.resample('MS').mean()
# Set the start and end dates
start_date = df.index.min() - pd.DateOffset(day=1)
end_date = df.index.max() + pd.DateOffset(day=31)
dates = pd.date_range(start_date, end_date, freq='D')
dates.name = 'date'
df = df.reindex(dates, method='ffill')
df = df.reset_index()
# Extract the values
nwetM1=df['X']
# Convert the pandas series date to list
new_wetM1 = nwetM1.tolist()

# wet2_Mavg (InvOcean)
df = pd.DataFrame((wet2_Mavg), index=date_avg) 
df = pd.DataFrame({'X':wet2}, index=datea) # Cape Grim
df = df.resample('MS').mean()
# Set the start and end dates
start_date = df.index.min() - pd.DateOffset(day=1)
end_date = df.index.max() + pd.DateOffset(day=31)
dates = pd.date_range(start_date, end_date, freq='D')
dates.name = 'date'
df = df.reindex(dates, method='ffill')
df = df.reset_index()
# Extract the values
nwetM2=df['X']
# Convert the pandas series date to list
new_wetM2 = nwetM2.tolist()

# wet3_Mavg (Anthro)
df = pd.DataFrame((wet3_Mavg), index=date_avg) 
df = pd.DataFrame({'X':wet3}, index=datea) # Cape Grim
df = df.resample('MS').mean()
# Set the start and end dates
start_date = df.index.min() - pd.DateOffset(day=1)
end_date = df.index.max() + pd.DateOffset(day=31)
dates = pd.date_range(start_date, end_date, freq='D')
dates.name = 'date'
df = df.reindex(dates, method='ffill')
df = df.reset_index()
# Extract the values
nwetM3=df['X']
# Convert the pandas series date to list
new_wetM3 = nwetM3.tolist()

# wet4_Mavg (Biomass)
df = pd.DataFrame((wet4_Mavg), index=date_avg) 
df = pd.DataFrame({'X':wet4}, index=datea) # Cape Grim
df = df.resample('MS').mean()
# Set the start and end dates
start_date = df.index.min() - pd.DateOffset(day=1)
end_date = df.index.max() + pd.DateOffset(day=31)
dates = pd.date_range(start_date, end_date, freq='D')
dates.name = 'date'
df = df.reindex(dates, method='ffill')
df = df.reset_index()
# Extract the values
nwetM4=df['X']
# Convert the pandas series date to list
new_wetM4 = nwetM4.tolist()

# wet5_Mavg (Geogenic)
df = pd.DataFrame((wet5_Mavg), index=date_avg) 
df = pd.DataFrame({'X':wet5}, index=datea) # Cape Grim
df = df.resample('MS').mean()
# Set the start and end dates
start_date = df.index.min() - pd.DateOffset(day=1)
end_date = df.index.max() + pd.DateOffset(day=31)
dates = pd.date_range(start_date, end_date, freq='D')
dates.name = 'date'
df = df.reindex(dates, method='ffill')
df = df.reset_index()
# Extract the values
nwetM5=df['X']
# Convert the pandas series date to list
new_wetM5 = nwetM5.tolist()

# wet6_Mavg (Soil)
df = pd.DataFrame((wet6_Mavg), index=date_avg) 
df = pd.DataFrame({'X':wet6}, index=datea) # Cape Grim
df = df.resample('MS').mean()
# Set the start and end dates
start_date = df.index.min() - pd.DateOffset(day=1)
end_date = df.index.max() + pd.DateOffset(day=31)
dates = pd.date_range(start_date, end_date, freq='D')
dates.name = 'date'
df = df.reindex(dates, method='ffill')
df = df.reset_index()
# Extract the values
nwetM6=df['X']
# Convert the pandas series date to list
new_wetM6 = nwetM6.tolist()

# wet7_Mavg (bromine_x3)
#df = pd.DataFrame((wet7_Mavg), index=date_avg) 
#df = pd.DataFrame({'X':wet7}, index=datea) # Cape Grim
#df = df.resample('MS').mean()
# Set the start and end dates
#start_date = df.index.min() - pd.DateOffset(day=1)
#end_date = df.index.max() + pd.DateOffset(day=31)
#dates = pd.date_range(start_date, end_date, freq='D')
#dates.name = 'date'
#df = df.reindex(dates, method='ffill')
#df = df.reset_index()
# Extract the values
#nwetM7=df['X']
# Convert the pandas series date to list
#new_wetM7 = nwetM7.tolist()

# Observations
df = pd.DataFrame((obs_Mavg), index=date_avg) 
df = pd.DataFrame({'X':y1}, index=datea) # Cape Grim
df = df.resample('MS').mean()
# Set the start and end dates
start_date = df.index.min() - pd.DateOffset(day=1)
end_date = df.index.max() + pd.DateOffset(day=31)
dates = pd.date_range(start_date, end_date, freq='D')
dates.name = 'date'
df = df.reindex(dates, method='ffill')
df = df.reset_index()
# Extract the values
nobsM=df['X']
# Convert the pandas series date to list
new_obsM = nobsM.tolist()
#------------------------------------------------------------------------------
# Calculate the daily variability (from Monthly Mean)
# (Daily Model Values) - (Model Monthly Mean) 

# Cape Grim
dvM1 = wet1 - new_wetM1 # Default
dvM2 = wet2 - new_wetM2 # InvOcean
dvM3 = wet3 - new_wetM3 # Anthro
dvM4 = wet4 - new_wetM4 # Biomass
dvM5 = wet5 - new_wetM5 # Geogenic
dvM6 = wet6 - new_wetM6 # Soil
#dvM7 = wet7 - new_wetM7 # bromine_x3
dvobsM = y1 - new_obsM            # Observations
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
wet1_Yavg, date_avg=yearly(wet1,datea) # default  
wet2_Yavg, date_avg=yearly(wet2,datea) # InvOcean
wet3_Yavg, date_avg=yearly(wet3,datea) # Anthro  
wet4_Yavg, date_avg=yearly(wet4,datea) # Biomass
wet5_Yavg, date_avg=yearly(wet5,datea) # Geogenic    
wet6_Yavg, date_avg=yearly(wet6,datea) # Soil
#wet7_Yavg, date_avg=yearly(wet7,datea) # bromine_x3
obs_Yavg,  date_avg=yearly(y1,datea)             # Observations

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
X3clean =[wet3[i] for i in goodvals]
X4clean =[wet4[i] for i in goodvals]
X5clean =[wet5[i] for i in goodvals]
X6clean =[wet6[i] for i in goodvals]
#X9clean =[wet9[i] for i in goodvals]
dateclean =[datea[i] for i in goodvals]

# 2nd filter
goodvals =np.where([i < (OM_OBS + 2*OSTD_OBS) for i in Y1_filtered])[0] #(x < mean + 2 * sd)
Y1clean =np.array([Y1_filtered[i] for i in goodvals])
X1clean =np.array([X1clean[i] for i in goodvals])
X2clean =np.array([X2clean[i] for i in goodvals])
X3clean =np.array([X3clean[i] for i in goodvals])
X4clean =np.array([X4clean[i] for i in goodvals])
X5clean =np.array([X5clean[i] for i in goodvals])
X6clean =np.array([X6clean[i] for i in goodvals])
#X9clean =np.array([X9clean[:,0,6,56][i] for i in goodvals])
dateclean =np.array([dateclean[i] for i in goodvals])

STD_DN=np.empty(len(y1))
STD_DN.fill(STDDN)
STD_UP=np.empty(len(y1))
STD_UP.fill(STDUP)
#------------------------------------------------------------------------------
# Calculate the overall mean
OM_1 = np.nanmean(wet1)   # Default
OM_2 = np.nanmean(wet2)   # InvOcean
OM_3 = np.nanmean(wet3)   # Anthro
OM_4 = np.nanmean(wet4)   # Biomass
OM_5 = np.nanmean(wet5)   # Geogenic
OM_6 = np.nanmean(wet6)   # Soil
#OM_7 = np.nanmean(wet7)   # bromine_x3
OM_OBS = np.nanmean(y1)             # Observations
#------------------------------------------------------------------------------
# Calculate the standard deviation
OSTD_1 = np.nanstd(wet1)   # Default
OSTD_2 = np.nanstd(wet2)   # InvOcean
OSTD_3 = np.nanstd(wet3)   # Anthro
OSTD_4 = np.nanstd(wet4)   # Biomass
OSTD_5 = np.nanstd(wet5)   # Geogenic
OSTD_6 = np.nanstd(wet6)   # Soil
#OSTD_7 = np.nanstd(wet7)   # bromine_x3
OSTD_OBS = np.nanstd(y1)             # Observations
#------------------------------------------------------------------------------
# Format the yearly mean from (1 value per year) to (1 value per day)
# wet1_Yavg (Default)
df = pd.DataFrame((wet1_Yavg), index=date_avg) 
df = pd.DataFrame({'X':wet1}, index=datea) # Cape Grim

df = df.resample('AS').mean()
# Set the start and end dates
start_date = df.index.min() - pd.DateOffset(day=1)
end_date = df.index.max() + pd.DateOffset(month=12, day=31)
dates = pd.date_range(start_date, end_date, freq='D')
dates.name = 'date'
df = df.reindex(dates, method='ffill')
df = df.reset_index()
# Extract the values
nwetY1=df['X']
# Convert the pandas series date to list
new_wetY1 = nwetY1.tolist()

# wet2_Yavg (InvOcean)
df = pd.DataFrame((wet2_Yavg), index=date_avg) 
df = pd.DataFrame({'X':wet2}, index=datea) # Cape Grim
df = df.resample('AS').mean()
# Set the start and end dates
start_date = df.index.min() - pd.DateOffset(day=1)
end_date = df.index.max() + pd.DateOffset(month=12, day=31)
dates = pd.date_range(start_date, end_date, freq='D')
dates.name = 'date'
df = df.reindex(dates, method='ffill')
df = df.reset_index()
# Extract the values
nwetY2=df['X']
# Convert the pandas series date to list
new_wetY2 = nwetY2.tolist()

# wet3_Yavg (Anthro)
df = pd.DataFrame((wet3_Yavg), index=date_avg) 
df = pd.DataFrame({'X':wet3}, index=datea) # Cape Grim

df = df.resample('AS').mean()
# Set the start and end dates
start_date = df.index.min() - pd.DateOffset(day=1)
end_date = df.index.max() + pd.DateOffset(month=12, day=31)
dates = pd.date_range(start_date, end_date, freq='D')
dates.name = 'date'
df = df.reindex(dates, method='ffill')
df = df.reset_index()
# Extract the values
nwetY3=df['X']
# Convert the pandas series date to list
new_wetY3 = nwetY3.tolist()

# wet4_Yavg (Biomass)
df = pd.DataFrame((wet4_Yavg), index=date_avg) 
df = pd.DataFrame({'X':wet4}, index=datea) # Cape Grim
df = df.resample('AS').mean()
# Set the start and end dates
start_date = df.index.min() - pd.DateOffset(day=1)
end_date = df.index.max() + pd.DateOffset(month=12, day=31)
dates = pd.date_range(start_date, end_date, freq='D')
dates.name = 'date'
df = df.reindex(dates, method='ffill')
df = df.reset_index()
# Extract the values
nwetY4=df['X']
# Convert the pandas series date to list
new_wetY4 = nwetY4.tolist()

# wet5_Yavg (Geogenic)
df = pd.DataFrame((wet5_Yavg), index=date_avg) 
df = pd.DataFrame({'X':wet5}, index=datea) # Cape Grim

df = df.resample('AS').mean()
# Set the start and end dates
start_date = df.index.min() - pd.DateOffset(day=1)
end_date = df.index.max() + pd.DateOffset(month=12, day=31)
dates = pd.date_range(start_date, end_date, freq='D')
dates.name = 'date'
df = df.reindex(dates, method='ffill')
df = df.reset_index()
# Extract the values
nwetY5=df['X']
# Convert the pandas series date to list
new_wetY5 = nwetY5.tolist()

# wet6_Yavg (Soil)
df = pd.DataFrame((wet6_Yavg), index=date_avg) 
df = pd.DataFrame({'X':wet6}, index=datea) # Cape Grim
df = df.resample('AS').mean()
# Set the start and end dates
start_date = df.index.min() - pd.DateOffset(day=1)
end_date = df.index.max() + pd.DateOffset(month=12, day=31)
dates = pd.date_range(start_date, end_date, freq='D')
dates.name = 'date'
df = df.reindex(dates, method='ffill')
df = df.reset_index()
# Extract the values
nwetY6=df['X']
# Convert the pandas series date to list
new_wetY6 = nwetY6.tolist()

# wet7_Yavg (bromine_x3)
#df = pd.DataFrame((wet7_Yavg), index=date_avg) 
#df = pd.DataFrame({'X':wet7}, index=datea) # Cape Grim
#df = df.resample('AS').mean()
# Set the start and end dates
#start_date = df.index.min() - pd.DateOffset(day=1)
#end_date = df.index.max() + pd.DateOffset(month=12, day=31)
#dates = pd.date_range(start_date, end_date, freq='D')
#dates.name = 'date'
#df = df.reindex(dates, method='ffill')
#df = df.reset_index()
# Extract the values
#nwetY7=df['X']
# Convert the pandas series date to list
#new_wetY7 = nwetY7.tolist()

# Observations
df = pd.DataFrame((obs_Yavg), index=date_avg) 
df = pd.DataFrame({'X':y1}, index=datea) # Cape Grim
df = df.resample('AS').mean()
# Set the start and end dates
start_date = df.index.min() - pd.DateOffset(day=1)
end_date = df.index.max() + pd.DateOffset(month=12, day=31)
dates = pd.date_range(start_date, end_date, freq='D')
dates.name = 'date'
df = df.reindex(dates, method='ffill')
df = df.reset_index()
# Extract the values
nobsY=df['X']
# Convert the pandas series date to list
new_obsY = nobsY.tolist()
#------------------------------------------------------------------------------
# Calculate the daily variability (from annual mean)
# (Daily Model Values) - (Model Annual Mean)

# Cape Grim
dvY1 = wet1 - new_wetY1 # Default
dvY2 = wet2 - new_wetY2 # InvOcean
dvY3 = wet3 - new_wetY3 # Anthro
dvY4 = wet4 - new_wetY4 # Biomass
dvY5 = wet5 - new_wetY5 # Geogenic
dvY6 = wet6 - new_wetY6 # Soil
#dvY7 = wet7 - new_wetY7 # bromine_x3
dvobsY = y1 - new_obsY            # Observations
#------------------------------------------------------------------------------
# Calculate the daily variability (from the sensitivity simulation)
# (InvOcean Values) - (Sensitivity Simulation Values)

# Cape Grim
#dvY1 = wet1 - new_wetY1 # Default
#dvS2 = wet2 - new_wetY2 # InvOcean
dvS3 = wet2 - wet3 # Anthro
dvS4 = wet2 - wet4 # Biomass
dvS5 = wet5 - wet2 # Geogenic
dvS6 = wet2 - wet6 # Soil
#dvS7 = wet2 - wet7 # bromine_x3
#------------------------------------------------------------------------------
# Calculate the mean bias from the difference (between simulation and observation)
MB1 = (np.mean((X1clean - Y1clean)/Y1clean)) # Default
MB2 = (np.mean((X2clean - Y1clean)/Y1clean)) # InvOcean
MB3 = (np.mean((X3clean - Y1clean)/Y1clean)) # Anthro
MB4 = (np.mean((X4clean - Y1clean)/Y1clean)) # Biomass
MB5 = (np.mean((X5clean - Y1clean)/Y1clean)) # Geogenic
MB6 = (np.mean((X6clean - Y1clean)/Y1clean)) # Soil
#MB9 = (np.mean((X9clean - Y1clean)/Y1clean) # Bromine_x3
#------------------------------------------------------------------------------
# Calculate the mean normalised bias (between simulation and observation)
MNB1 = (np.mean((X1clean - Y1clean)/Y1clean))*100 # Default
MNB2 = (np.mean((X2clean - Y1clean)/Y1clean))*100 # InvOcean
MNB3 = (np.mean((X3clean - Y1clean)/Y1clean))*100 # Anthro
MNB4 = (np.mean((X4clean - Y1clean)/Y1clean))*100 # Biomass
MNB5 = (np.mean((X5clean - Y1clean)/Y1clean))*100 # Geogenic
MNB6 = (np.mean((X6clean - Y1clean)/Y1clean))*100 # Soil
#MNB9 = (np.mean((X9clean - Y1clean)/Y1clean)*100 # Bromine_x3
#------------------------------------------------------------------------------
# Calculate the standard deviation (of the daily variability)
# Default
std1 = np.std(wet1) # Cape Grim

# InvOcean
std2 = np.std(wet2) # Cape Grim

# Anthro
std3 = np.std(wet3) # Cape Grim

# Biomass
std4 = np.std(wet4) # Cape Grim

# Geogenic
std5 = np.std(wet5) # Cape Grim

# Soil
std6 = np.std(wet6) # Cape Grim

# bromine_x3
#std7 = np.std(wet7) # Cape Grim

# OBSERVATIONS
stdOBS = np.nanstd(y1) # Standard deviation of Observations Daily Variability
#------------------------------------------------------------------------------
# Simulation - Overall Mean Bias

SMB1 = wet1 - MB1 # Default
SMB2 = wet2 - MB2 # InvOcean
SMB3 = wet3 - MB3 # Anthro
SMB4 = wet4 - MB4 # Biomass
SMB5 = wet5 - MB5 # Geogenic
SMB6 = wet6 - MB6 # Soil
#SMB7 = wet7 - MB7 # bromine_x3
#------------------------------------------------------------------------------
# Calculate the Coefficient of Correlation (r)

# Default
r_row1, p_value1 = stats.pearsonr(X1clean,Y1clean)
slope1, intercept1, r1, p1, std_err1= stats.linregress(X1clean,Y1clean)

# InvOcean
r_row2, p_value2 = stats.pearsonr(X2clean,Y1clean)
slope2, intercept2, r2, p2, std_err2= stats.linregress(X2clean,Y1clean)

# Anthro
r_row3, p_value3 = stats.pearsonr(X3clean,Y1clean)
slope3, intercept3, r3, p3, std_err3= stats.linregress(X3clean,Y1clean)

# Biomass
r_row4, p_value4 = stats.pearsonr(X4clean,Y1clean)
slope4, intercept4, r4, p4, std_err4= stats.linregress(X4clean,Y1clean)

# Geogenic
r_row5, p_value5 = stats.pearsonr(X5clean,Y1clean)
slope5, intercept5, r5, p5, std_err5= stats.linregress(X5clean,Y1clean)

# Soil
r_row6, p_value6 = stats.pearsonr(X6clean,Y1clean)
slope6, intercept6, r6, p6, std_err6= stats.linregress(X6clean,Y1clean)

# bromine_x3
#r_row7, p_value7 = stats.pearsonr(X7clean,Y1clean)
#slope7, intercept7, r7, p7, std_err7= stats.linregress(X7clean,Y1clean)

# Default (from monthly mean value)
finiteYmask = np.isfinite(obs_Mavg) # Scan for NaN values
YMclean = obs_Mavg[finiteYmask]     # Remove NaN values from the obs
XMclean = wet1_Mavg[finiteYmask]     # Remove NaN values from the obs

r_row8, p_value8 = stats.pearsonr(XMclean,YMclean)
slope8, intercept8, r8, p8, std_err8= stats.linregress(XMclean,YMclean)
#------------------------------------------------------------------------------
# Calculate the Coefficient of Determination (r2)

r2_1 = r_row1 * r_row1 # Default
r2_2 = r_row2 * r_row2 # InvOcean
r2_3 = r_row3 * r_row3 # Anthro
r2_4 = r_row4 * r_row4 # Biomass
r2_5 = r_row5 * r_row5 # Geogenic
r2_6 = r_row6 * r_row6 # Soil
#r2_7 = r_row7 * r_row7 # bromine_x3
r2_8 = r_row8 * r_row8 # Soil
#------------------------------------------------------------------------------
# Plot the axis for each graph
# Graph 1
fig = plt.figure()
plt.subplots_adjust(hspace=0.5)
ax=plt.subplot(311) # options graph 1 (vertical no, horizontal no, graph no)

# OBSERVATIONS
# Plot the observations for the location
plt.plot(datea,y1, 'k', linewidth=2, label="Observations")

# SIMULATIONS
# Plot the modelling simulations for each location
plt.plot(datea, wet1, ls='--', color='red', linewidth=2, label="Default \n   (r value: "+str("%7.4f"%(r_row1))+")\n   (r$^2$ value: "+str("%7.4f"%(r2_1))+")\n   (MNB: "+str("%7.4f"%(MNB1))+"%)")       # wet(time,altitude,latitude,longitude)
plt.plot(datea, wet2, ls='--', color='blue', linewidth=2, label="InvOcean \n   (r value: "+str("%7.4f"%(r_row2))+")\n   (r$^2$ value: "+str("%7.4f"%(r2_2))+")\n   (MNB: "+str("%7.4f"%(MNB2))+"%)")     # wet(time,altitude,latitude,longitude)
plt.plot(datea, wet3, ls='--', color='cyan', linewidth=2, label="Anthro \n   (r value: "+str("%7.4f"%(r_row3))+")\n   (r$^2$ value: "+str("%7.4f"%(r2_3))+")")       # wet(time,altitude,latitude,longitude)
plt.plot(datea, wet4, ls='--', color='magenta', linewidth=2, label="Biomass \n   (r value: "+str("%7.4f"%(r_row4))+")\n   (r$^2$ value: "+str("%7.4f"%(r2_4))+")")   # wet(time,altitude,latitude,longitude)
plt.plot(datea, wet5, ls='--', color='darkslategray', linewidth=2, label="Geogenic \n   (r value: "+str("%7.4f"%(r_row5))+")\n   (r$^2$ value: "+str("%7.4f"%(r2_5))+")")   # wet(time,altitude,latitude,longitude)
plt.plot(datea, wet6, ls='--', color='saddlebrown', linewidth=2, label="Soil \n   (r value: "+str("%7.4f"%(r_row6))+")\n   (r$^2$ value: "+str("%7.4f"%(r2_6))+")")  # wet(time,altitude,latitude,longitude)
#plt.plot(datea, wet7, ls='--', color='orange', linewidth=2, label="Bromine_x3 \n   (r value: "+str("%7.4f"%(r_row7))+")\n   (r$^2$ value: "+str("%7.4f"%(r2_7))+")") # wet(time,altitude,latitude,longitude)

# Plot the Monthly and Yearly mean for the location
#plt.plot(datea, wet1_Mavg[:], 'g--', linewidth=2, label="Sim 1 Monthly Avg, Cape Grim, Tas")   # wet(time,altitude,latitude,longitude)

# Plot the standard deviation for a location
#mu1=wet1+np.std(wet1)         # find the upper limit
#mu2=wet1-np.std(wet1)         # find the lower limit
#plt.plot(datea, mu1, 'r-', linewidth=2, alpha=0.2)   # plot the upper limit
#plt.plot(datea, mu2, 'r-', linewidth=2, alpha=0.2)   # plot the lower limit
#p2 = ax.fill_between(datea, wet1+np.std(wet1), wet1-np.std(wet1), facecolor='red', alpha=0.2) # fill the distribution
#ax.legend([(p2[0], p1[0]), ], ['Stuff'])

# Plot the standard deviation for the observations
#mu1=y1+np.nanstd(y1)         # find the upper limit
#mu2=y1-np.nanstd(y1)         # find the lower limit
#plt.plot(datea, mu1, 'k-', linewidth=2, alpha=0.2)   # plot the upper limit
#plt.plot(datea, mu2, 'k-', linewidth=2, alpha=0.2)   # plot the lower limit
#p2 = ax.fill_between(datea, y1+np.nanstd(y1), y1-np.nanstd(y1), facecolor='k', alpha=0.2) # fill the distribution

plt.xlim(datetime(2013,1,1),datetime(2016,12,31))   # set the date limits
plt.xticks(rotation=15)
date_formatter = mdates.DateFormatter('%b/%Y')      # format how the date is displayed
ax.xaxis.set_major_formatter(date_formatter)
ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3)) # set the interval between dispalyed dates
ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
ax.yaxis.set_ticks_position('both')

# Plot the axis labels, legend and title
plt.ylabel('Hg$^0$ (ng/m$^3$)', fontsize=10)
plt.title("Daily Average Hg$^0$ concentrations (Cape Grim, Tas)", fontsize=15)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#plt.plot(wet,lev , '--', linewidth=2)
#plt.plot(wet,lat , '--', linewidth=2)

# Graph 2
ax=plt.subplot(312) # options graph 2 (vertical no, horizontal no, graph no)

# OBSERVATIONS
# Plot the Monthly/Yearly mean for the location
plt.plot(monthsa,obs_Mavg, 'k', linewidth=5, label="Observations")

# SIMULATIONS
# Plot the monthly mean for the location
# Cape Grim
plt.plot(monthsa, wet1_Mavg, color='red', linewidth=3, label="Default")        # wet1_Mavg
plt.plot(monthsa, wet2_Mavg, color='blue', linewidth=3, label="InvOcean")     # wet2_Mavg
plt.plot(monthsa, wet3_Mavg, ls='--', color='cyan', linewidth=2, label="Anthro")       # wet3_Mavg
plt.plot(monthsa, wet4_Mavg, ls='--', color='magenta', linewidth=2, label="Biomass")   # wet4_Mavg
plt.plot(monthsa, wet5_Mavg, ls='--', color='darkslategray', linewidth=2, label="Geogenic")   # wet5_Mavg
plt.plot(monthsa, wet6_Mavg, ls='--', color='saddlebrown', linewidth=2, label="Soil")  # wet6_Mavg
#plt.plot(monthsa, wet7_Mavg, ls='--', color='orange', linewidth=2, label="Bromine_x3") # wet7_Mavg

# Cape Grim
#plt.plot(datea, new_wet1, 'r--', linewidth=2, label="Sim 1, Cape Grim, Tas")   # new_wet1
#plt.plot(datea, new_wet2, 'b--', linewidth=2, label="Sim 2, Cape Grim, Tas")   # new_wet2

plt.xlim(datetime(2013,1,1),datetime(2016,12,31))   # set the date limits
plt.xticks(rotation=15)
date_formatter = mdates.DateFormatter('%b/%Y')      # format how the date is displayed
ax.xaxis.set_major_formatter(date_formatter)
ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3)) # set the interval between dispalyed dates
ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=0.5))
ax.yaxis.set_ticks_position('both')

# Plot the axis labels, legend and title
#plt.ylabel('Hg$^0$ (ng/m$^3$)', fontsize=10)
#plt.title("Monthly Mean Hg$^0$ concentrations (Cape Grim, Tas)", fontsize=15)
#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#plt.plot(wet,lev , '--', linewidth=2)
#plt.plot(wet,lat , '--', linewidth=2)

# Graph 3
ax=plt.subplot(313) # options graph 2 (vertical no, horizontal no, graph no)

# OBSERVATIONS
# Plot the Monthly/Yearly mean for the location
plt.plot(datea,y1, 'k', linewidth=5, label="Observations")

# SIMULATIONS
# Plot the Simulation - MNB
# Cape Grim
plt.plot(datea, SMB1, color='red', linewidth=4, label="Default")        # Mean bias observations and Normal
plt.plot(datea, SMB2, color='blue', linewidth=4, label="InvOcean")     # Mean bias observations and InvOcean
#plt.plot(datea, SMB3, ls='--', color='cyan', linewidth=2, label="Anthro")       # Mean bias observations and Anthro
#plt.plot(datea, SMB4, ls='--', color='magenta', linewidth=2, label="Biomass")   # Mean bias observations and Biomass
#plt.plot(datea, SMB5, ls='--', color='darkslategray', linewidth=2, label="Geogenic")   # Mean bias observations and Geogenic
#plt.plot(datea, SMB6, ls='--', color='saddlebrown', linewidth=2, label="Soil")  # Mean bias observations and Soil
#plt.plot(datea, SMB7, ls='--', color='orange', linewidth=2, label="Bromine_x3") # Mean bias observations and Bromine_x3

plt.xlim(datetime(2013,1,1),datetime(2016,12,31))   # set the date limits
plt.xticks(rotation=15)
date_formatter = mdates.DateFormatter('%b/%Y')      # format how the date is displayed
ax.xaxis.set_major_formatter(date_formatter)
ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3)) # set the interval between dispalyed dates
ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax.yaxis.set_ticks_position('both')
ax.set_ylim(0.5, 1.4)  

# Plot the axis labels, legend and title
#plt.xlabel('Date', fontsize=10)
#plt.ylabel('Hg$^0$ (ng/m$^3$)', fontsize=10)
#plt.title("Simulation - MNB (Cape Grim, Tas)", fontsize=15)
#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#plt.axhline(0, linewidth=0.5, color='k')
#plt.plot(wet,lev , '--', linewidth=2)
#plt.plot(wet,lat , '--', linewidth=2)

# DAILY VARIABILITY
fig = plt.figure()
plt.subplots_adjust(hspace=0.5)
#Graph 1
ax=plt.subplot(311) # options graph 2 (vertical no, horizontal no, graph no)

# SIMULATIONS

# OBSERVATIONS

# Plot the Daily Variability from the Monthly/Yearly mean for the location
#plt.plot(datea,dvobsY, 'k', linewidth=2, label="Observations")
#
## Plot the Daily Variability
## Cape Grim
#plt.plot(datea, dvY1, ls='--', color='red', linewidth=2, label="Default \n(r value: "+str("%7.4f"%(r_row1))+")\n(r$^2$ value: "+str("%7.4f"%(r2_1))+")")        # dvM1
##plt.plot(datea, dvM2, ls='--', color='blue', linewidth=2, label="InvOcean \n(r value: "+str("%7.4f"%(r_row2))+")\n(r$^2$ value: "+str("%7.4f"%(r2_2))+")")     # dvM2
##plt.plot(datea, dvM3, ls='--', color='cyan', linewidth=2, label="Anthro \n(r value: "+str("%7.4f"%(r_row3))+")\n(r$^2$ value: "+str("%7.4f"%(r2_3))+")")       # dvM3
##plt.plot(datea, dvM4, ls='--', color='magenta', linewidth=2, label="Biomass \n(r value: "+str("%7.4f"%(r_row4))+")\n(r$^2$ value: "+str("%7.4f"%(r2_4))+")")   # dvM4
##plt.plot(datea, dvM5, ls='--', color='darkslategray', linewidth=2, label="Geogenic \n(r value: "+str("%7.4f"%(r_row5))+")\n(r$^2$ value: "+str("%7.4f"%(r2_5))+")")   # dvM5
##plt.plot(datea, dvM6, ls='--', color='saddlebrown', linewidth=2, label="Soil \n(r value: "+str("%7.4f"%(r_row6))+")\n(r$^2$ value: "+str("%7.4f"%(r2_6))+")")  # dvM6
##plt.plot(datea, dvM7, ls='--', color='orange', linewidth=2, label="Bromine_x3 \n(r value: "+str("%7.4f"%(r_row7))+")\n(r$^2$ value: "+str("%7.4f"%(r2_7))+")")  # dvM7
#
#plt.xlim(datetime(2013,1,1),datetime(2016,12,31))   # set the date limits
#plt.xticks(rotation=20)
#date_formatter = mdates.DateFormatter('%b/%Y')      # format how the date is displayed
#ax.xaxis.set_major_formatter(date_formatter)
#ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3)) # set the interval between dispalyed dates
#ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
#ax.yaxis.set_ticks_position('both')
#ax.set_ylim(-0.3, 0.3)
#
## Plot the axis labels, legend and title
#plt.ylabel('Hg$^0$ (ng/m$^3$)', fontsize=10)
#plt.title("Daily Variability of Default Simulation (Cape Grim, Tas)", fontsize=15)
##plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#plt.axhline(0, linewidth=0.5, color='k')
##plt.plot(wet,lev , '--', linewidth=2)
##plt.plot(wet,lat , '--', linewidth=2)
#
##Graph 2
#ax=plt.subplot(322) # options graph 2 (vertical no, horizontal no, graph no)
#
## SIMULATIONS
#
## OBSERVATIONS
#
## Plot the Daily Variability from the Monthly/Yearly mean for the location
#plt.plot(datea,dvobsY, 'k', linewidth=2, label="Observations")
#
## Plot the Daily Variability
## Cape Grim
##plt.plot(datea, dvM1, ls='--', color='red', linewidth=2, label="Default \n(r value: "+str("%7.4f"%(r_row1))+")\n(r$^2$ value: "+str("%7.4f"%(r2_1))+")")        # dvM1
#plt.plot(datea, dvY2, ls='--', color='blue', linewidth=2, label="InvOcean \n(r value: "+str("%7.4f"%(r_row2))+")\n(r$^2$ value: "+str("%7.4f"%(r2_2))+")")     # dvM2
##plt.plot(datea, dvM3, ls='--', color='cyan', linewidth=2, label="Anthro \n(r value: "+str("%7.4f"%(r_row3))+")\n(r$^2$ value: "+str("%7.4f"%(r2_3))+")")       # dvM3
##plt.plot(datea, dvM4, ls='--', color='magenta', linewidth=2, label="Biomass \n(r value: "+str("%7.4f"%(r_row4))+")\n(r$^2$ value: "+str("%7.4f"%(r2_4))+")")   # dvM4
##plt.plot(datea, dvM5, ls='--', color='darkslategray', linewidth=2, label="Geogenic \n(r value: "+str("%7.4f"%(r_row5))+")\n(r$^2$ value: "+str("%7.4f"%(r2_5))+")")   # dvM5
##plt.plot(datea, dvM6, ls='--', color='saddlebrown', linewidth=2, label="Soil \n(r value: "+str("%7.4f"%(r_row6))+")\n(r$^2$ value: "+str("%7.4f"%(r2_6))+")")  # dvM6
##plt.plot(datea, dvM7, ls='--', color='orange', linewidth=2, label="Bromine_x3 \n(r value: "+str("%7.4f"%(r_row7))+")\n(r$^2$ value: "+str("%7.4f"%(r2_7))+")")  # dvM7
#
#plt.xlim(datetime(2013,1,1),datetime(2016,12,31))   # set the date limits
#plt.xticks(rotation=20)
#date_formatter = mdates.DateFormatter('%b/%Y')      # format how the date is displayed
#ax.xaxis.set_major_formatter(date_formatter)
#ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3)) # set the interval between dispalyed dates
#ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
#ax.yaxis.set_ticks_position('both')
#ax.set_ylim(-0.3, 0.3)
#
## Plot the axis labels, legend and title
#plt.ylabel('Hg$^0$ (ng/m$^3$)', fontsize=10)
#plt.title("Daily Variability of InvOcean Simulation (Cape Grim, Tas)", fontsize=15)
##plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#plt.axhline(0, linewidth=0.5, color='k')
##plt.plot(wet,lev , '--', linewidth=2)
##plt.plot(wet,lat , '--', linewidth=2)

#Graph 3
ax=plt.subplot(311) # options graph 2 (vertical no, horizontal no, graph no)

# SIMULATIONS

# OBSERVATIONS

# Plot the Daily Variability from the Monthly/Yearly mean for the location
#plt.plot(datea,dvobsY, 'k', linewidth=2, label="Observations")

# Plot the Daily Variability
# Cape Grim
#plt.plot(datea, dvM1, ls='--', color='red', linewidth=2, label="Default \n(r value: "+str("%7.4f"%(r_row1))+")\n(r$^2$ value: "+str("%7.4f"%(r2_1))+")")        # dvM1
#plt.plot(datea, dvM2, ls='--', color='blue', linewidth=2, label="InvOcean \n(r value: "+str("%7.4f"%(r_row2))+")\n(r$^2$ value: "+str("%7.4f"%(r2_2))+")")     # dvM2
plt.plot(datea, dvS3, ls='--', color='cyan', linewidth=2, label="Anthro \n(r value: "+str("%7.4f"%(r_row3))+")\n(r$^2$ value: "+str("%7.4f"%(r2_3))+")")       # dvM3
#plt.plot(datea, dvM4, ls='--', color='magenta', linewidth=2, label="Biomass \n(r value: "+str("%7.4f"%(r_row4))+")\n(r$^2$ value: "+str("%7.4f"%(r2_4))+")")   # dvM4
#plt.plot(datea, dvM5, ls='--', color='darkslategray', linewidth=2, label="Geogenic \n(r value: "+str("%7.4f"%(r_row5))+")\n(r$^2$ value: "+str("%7.4f"%(r2_5))+")")   # dvM5
#plt.plot(datea, dvM6, ls='--', color='saddlebrown', linewidth=2, label="Soil \n(r value: "+str("%7.4f"%(r_row6))+")\n(r$^2$ value: "+str("%7.4f"%(r2_6))+")")  # dvM6
#plt.plot(datea, dvM7, ls='--', color='orange', linewidth=2, label="Bromine_x3 \n(r value: "+str("%7.4f"%(r_row7))+")\n(r$^2$ value: "+str("%7.4f"%(r2_7))+")")  # dvM7

#plt.xlim(datetime(2013,1,1),datetime(2016,12,31))   # set the date limits
#plt.xticks(rotation=20)
#date_formatter = mdates.DateFormatter('%b/%Y')      # format how the date is displayed
#ax.xaxis.set_major_formatter(date_formatter)
#ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3)) # set the interval between dispalyed dates
#ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
ax.yaxis.set_ticks_position('both')
#ax.set_ylim(-0.3, 0.3)
plt.yticks(fontsize=14)

# Plot the axis labels, legend and title
#plt.ylabel('Hg$^0$ (ng/m$^3$)', fontsize=10)
#plt.title("Daily Variability of Anthro Simulation (Cape Grim, Tas)", fontsize=15)
#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.axhline(0, linewidth=0.5, color='k')
#plt.plot(wet,lev , '--', linewidth=2)
#plt.plot(wet,lat , '--', linewidth=2)

#Graph 4
ax=plt.subplot(312) # options graph 2 (vertical no, horizontal no, graph no)

# SIMULATIONS

# OBSERVATIONS

# Plot the Daily Variability from the Monthly/Yearly mean for the location
#plt.plot(datea,dvobsY, 'k', linewidth=2, label="Observations")

# Plot the Daily Variability
# Cape Grim
#plt.plot(datea, dvM1, ls='--', color='red', linewidth=2, label="Default \n(r value: "+str("%7.4f"%(r_row1))+")\n(r$^2$ value: "+str("%7.4f"%(r2_1))+")")        # dvM1
#plt.plot(datea, dvM2, ls='--', color='blue', linewidth=2, label="InvOcean \n(r value: "+str("%7.4f"%(r_row2))+")\n(r$^2$ value: "+str("%7.4f"%(r2_2))+")")     # dvM2
#plt.plot(datea, dvM3, ls='--', color='cyan', linewidth=2, label="Anthro \n(r value: "+str("%7.4f"%(r_row3))+")\n(r$^2$ value: "+str("%7.4f"%(r2_3))+")")       # dvM3
plt.plot(datea, dvS4, ls='--', color='magenta', linewidth=2, label="Biomass \n(r value: "+str("%7.4f"%(r_row4))+")\n(r$^2$ value: "+str("%7.4f"%(r2_4))+")")   # dvM4
#plt.plot(datea, dvM5, ls='--', color='darkslategray', linewidth=2, label="Geogenic \n(r value: "+str("%7.4f"%(r_row5))+")\n(r$^2$ value: "+str("%7.4f"%(r2_5))+")")   # dvM5
#plt.plot(datea, dvM6, ls='--', color='saddlebrown', linewidth=2, label="Soil \n(r value: "+str("%7.4f"%(r_row6))+")\n(r$^2$ value: "+str("%7.4f"%(r2_6))+")")  # dvM6
#plt.plot(datea, dvM7, ls='--', color='orange', linewidth=2, label="Bromine_x3 \n(r value: "+str("%7.4f"%(r_row7))+")\n(r$^2$ value: "+str("%7.4f"%(r2_7))+")")  # dvM7

#plt.xlim(datetime(2013,1,1),datetime(2016,12,31))   # set the date limits
#plt.xticks(rotation=20)
#date_formatter = mdates.DateFormatter('%b/%Y')      # format how the date is displayed
#ax.xaxis.set_major_formatter(date_formatter)
#ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3)) # set the interval between dispalyed dates
#ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
ax.yaxis.set_ticks_position('both')
#ax.set_ylim(-0.3, 0.3)
plt.yticks(fontsize=14)

# Plot the axis labels, legend and title
#plt.ylabel('Hg$^0$ (ng/m$^3$)', fontsize=10)
#plt.title("Daily Variability of Biomass Simulation (Cape Grim, Tas)", fontsize=15)
#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.axhline(0, linewidth=0.5, color='k')
#plt.plot(wet,lev , '--', linewidth=2)
#plt.plot(wet,lat , '--', linewidth=2)

##Graph 5
#ax=plt.subplot(325) # options graph 2 (vertical no, horizontal no, graph no)
#
## SIMULATIONS
#
## OBSERVATIONS
#
## Plot the Daily Variability from the Monthly/Yearly mean for the location
##plt.plot(datea,dvobsY, 'k', linewidth=2, label="Observations")
#
## Plot the Daily Variability
## Cape Grim
##plt.plot(datea, dvM1, ls='--', color='red', linewidth=2, label="Default \n(r value: "+str("%7.4f"%(r_row1))+")\n(r$^2$ value: "+str("%7.4f"%(r2_1))+")")        # dvM1
##plt.plot(datea, dvM2, ls='--', color='blue', linewidth=2, label="InvOcean \n(r value: "+str("%7.4f"%(r_row2))+")\n(r$^2$ value: "+str("%7.4f"%(r2_2))+")")     # dvM2
##plt.plot(datea, dvM3, ls='--', color='cyan', linewidth=2, label="Anthro \n(r value: "+str("%7.4f"%(r_row3))+")\n(r$^2$ value: "+str("%7.4f"%(r2_3))+")")       # dvM3
##plt.plot(datea, dvM4, ls='--', color='magenta', linewidth=2, label="Biomass \n(r value: "+str("%7.4f"%(r_row4))+")\n(r$^2$ value: "+str("%7.4f"%(r2_4))+")")   # dvM4
#plt.plot(datea, dvS5, ls='--', color='darkslategray', linewidth=2, label="Geogenic \n(r value: "+str("%7.4f"%(r_row5))+")\n(r$^2$ value: "+str("%7.4f"%(r2_5))+")")   # dvM5
##plt.plot(datea, dvM6, ls='--', color='saddlebrown', linewidth=2, label="Soil \n(r value: "+str("%7.4f"%(r_row6))+")\n(r$^2$ value: "+str("%7.4f"%(r2_6))+")")  # dvM6
##plt.plot(datea, dvM7, ls='--', color='orange', linewidth=2, label="Bromine_x3 \n(r value: "+str("%7.4f"%(r_row7))+")\n(r$^2$ value: "+str("%7.4f"%(r2_7))+")")  # dvM7
#
#plt.xlim(datetime(2013,1,1),datetime(2016,12,31))   # set the date limits
#plt.xticks(rotation=20)
#date_formatter = mdates.DateFormatter('%b/%Y')      # format how the date is displayed
#ax.xaxis.set_major_formatter(date_formatter)
#ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3)) # set the interval between dispalyed dates
#ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
#ax.yaxis.set_ticks_position('both')
##ax.set_ylim(-0.3, 0.3)
#
## Plot the axis labels, legend and title
#plt.ylabel('Hg$^0$ (ng/m$^3$)', fontsize=10)
#plt.title("Daily Variability of Geogenic Simulation (Cape Grim, Tas)", fontsize=15)
##plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#plt.axhline(0, linewidth=0.5, color='k')
##plt.plot(wet,lev , '--', linewidth=2)
##plt.plot(wet,lat , '--', linewidth=2)

#Graph 6
ax=plt.subplot(313) # options graph 2 (vertical no, horizontal no, graph no)

# SIMULATIONS

# OBSERVATIONS

# Plot the Daily Variability from the Monthly/Yearly mean for the location
#plt.plot(datea,dvobsY, 'k', linewidth=2, label="Observations")

# Plot the Daily Variability
# Cape Grim
#plt.plot(datea, dvM1, ls='--', color='red', linewidth=2, label="Default \n(r value: "+str("%7.4f"%(r_row1))+")\n(r$^2$ value: "+str("%7.4f"%(r2_1))+")")        # dvM1
#plt.plot(datea, dvM2, ls='--', color='blue', linewidth=2, label="InvOcean \n(r value: "+str("%7.4f"%(r_row2))+")\n(r$^2$ value: "+str("%7.4f"%(r2_2))+")")     # dvM2
#plt.plot(datea, dvM3, ls='--', color='cyan', linewidth=2, label="Anthro \n(r value: "+str("%7.4f"%(r_row3))+")\n(r$^2$ value: "+str("%7.4f"%(r2_3))+")")       # dvM3
#plt.plot(datea, dvM4, ls='--', color='magenta', linewidth=2, label="Biomass \n(r value: "+str("%7.4f"%(r_row4))+")\n(r$^2$ value: "+str("%7.4f"%(r2_4))+")")   # dvM4
#plt.plot(datea, dvM5, ls='--', color='darkslategray', linewidth=2, label="Geogenic \n(r value: "+str("%7.4f"%(r_row5))+")\n(r$^2$ value: "+str("%7.4f"%(r2_5))+")")   # dvM5
plt.plot(datea, dvS6, ls='--', color='saddlebrown', linewidth=2, label="Soil \n(r value: "+str("%7.4f"%(r_row6))+")\n(r$^2$ value: "+str("%7.4f"%(r2_6))+")")  # dvM6
#plt.plot(datea, dvM7, ls='--', color='orange', linewidth=2, label="Bromine_x3 \n(r value: "+str("%7.4f"%(r_row7))+")\n(r$^2$ value: "+str("%7.4f"%(r2_7))+")")  # dvM7

plt.xlim(datetime(2013,1,1),datetime(2016,12,31))   # set the date limits
plt.xticks(rotation=40)
date_formatter = mdates.DateFormatter('%b/%Y')      # format how the date is displayed
ax.xaxis.set_major_formatter(date_formatter)
ax.xaxis.set_major_locator(mdates.MonthLocator(interval=6)) # set the interval between dispalyed dates
ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=1))
ax.yaxis.set_ticks_position('both')
#ax.set_ylim(-0.3, 0.3)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)

# Plot the axis labels, legend and title
#plt.ylabel('Hg$^0$ (ng/m$^3$)', fontsize=10)
#plt.title("Daily Variability of Soil Simulation  (Cape Grim, Tas)", fontsize=15)
#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.axhline(0, linewidth=0.5, color='k')
#plt.plot(wet,lev , '--', linewidth=2)
#plt.plot(wet,lat , '--', linewidth=2)

# MODELLED vs OBSERVED VALUES
# Plot the axis for each graph
fig = plt.figure()
plt.subplots_adjust(hspace=0.5)

# Graph 1
ax=plt.subplot(111) # options graph 2 (vertical no, horizontal no, graph no)

# SIMULATIONS
# Plot the monthly mean for the location
#line2, = plt.plot(wet1,y1, "o", color='red', label="Default:") # Daily Values
line2, = plt.plot(SMB1,y1, "o", color='red', label="Default:") # Default - MB
#line3, = plt.plot(wet2,y1, "o", color='blue', label="InvOcean") # Daily Values
line3, = plt.plot(SMB2,y1, "o", color='blue', label="InvOcean") # InvOcean - MB

# Plot the regression line
#line4, = plt.plot(wet1, intercept1 + slope1*wet1, color='red', label="Default:\n (slope: "+str("%7.4f"%(p_value1))+"$\pm$ "+str("%7.4f"%(std_err1))+"%, r: "+str("%7.4f"%(r_row1))+")")
line4, = plt.plot(SMB1, intercept1 + slope1*SMB1, color='red', label="Default:\n (slope: "+str("%7.4f"%(p_value1))+"$\pm$ "+str("%7.4f"%(std_err1))+", r: "+str("%7.4f"%(r_row1))+")")
#line5, = plt.plot(wet2, intercept2 + slope2*wet2, color='blue',label="InvOcean:\n (slope: "+str("%7.4f"%(p_value2))+"$\pm$ "+str("%7.4f"%(std_err2))+"%, r: "+str("%7.4f"%(r_row2))+")")
line5, = plt.plot(SMB2, intercept2 + slope2*SMB2, color='blue',label="InvOcean:\n (slope: "+str("%7.4f"%(p_value2))+"$\pm$ "+str("%7.4f"%(std_err2))+", r: "+str("%7.4f"%(r_row2))+")")
#line6 = plt.plot(wet1, STD_DN, color='black', alpha=0.5)
line6 = plt.plot(SMB1, STD_DN, color='black', alpha=0.5)
#line7 = plt.plot(wet1, STD_UP, color='black', alpha=0.5)
line7 = plt.plot(SMB1, STD_UP, color='black', alpha=0.5)

ax.yaxis.set_ticks_position('both')
ax.xaxis.set_major_locator(ticker.MultipleLocator(0.05))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.01))
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.10))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.01))
#ax.set_xlim(49.75, 100.25) # Daily
#ax.set_ylim(0.542, 1.458)  # Daily

# Plot the axis labels, legend and title
#plt.ylabel('Hg$^0$ (ng/m$^3$)', fontsize=10) # Daily and Monthly
plt.ylabel('Observation Hg$^0$ (ng/m$^3$)', fontsize=10)        
#plt.xlabel('Simulation Hg$^0$ (ng/m$^3$)', fontsize=10)
plt.xlabel('(Simulation - MNB) Hg$^0$ (ng/m$^3$)', fontsize=10)
#plt.title("Relationship between simulated and observed values of Hg$^0$ (Cape Grim, Tas)", fontsize=15)
plt.title("Relationship between (Simulation - MNB) and observed values of Hg$^0$ (Cape Grim, Tas)", fontsize=15)
#plt.legend(handles=[line6,line4,line5], loc=1, handlelength=0, handletextpad=0)
# ANNOTATE DAILY VALUES
#plt.annotate("Default:\n (slope: "+str("%7.4f"%(slope1))+"$\pm$ "+str("%7.4f"%(std_err1))+" %$^-$$^1$, r: "+str("%7.4f"%(r1))+")", xy=(1.18,1.85), color='red', fontweight='bold')
#plt.annotate("InvOcean:\n (slope: "+str("%7.4f"%(slope2))+"$\pm$ "+str("%7.4f"%(std_err2))+" %$^-$$^1$, r: "+str("%7.4f"%(r2))+")", xy=(1.18,1.75), color='blue', fontweight='bold')
#plt.annotate("Upper limit: "+str("%7.4f"%(STDUP)), xy=(1.18,1.70), color='black', fontweight='bold', alpha=0.5)
#plt.annotate("Lower limit: "+str("%7.4f"%(STDDN)), xy=(1.18,1.65), color='black', fontweight='bold', alpha=0.5)
# ANNOTATE SIM_MNB
plt.annotate("Default:\n (slope: "+str("%7.4f"%(slope1))+" $\pm$ "+str("%7.4f"%(std_err1))+", r: "+str("%7.4f"%(r1))+")", xy=(1.07,1.85), color='red', fontweight='bold')
plt.annotate("InvOcean:\n (slope: "+str("%7.4f"%(slope2))+" $\pm$ "+str("%7.4f"%(std_err2))+", r: "+str("%7.4f"%(r2))+")", xy=(1.07,1.75), color='blue', fontweight='bold')
plt.annotate("Upper limit: "+str("%7.4f"%(STDUP))+" ng/m$^3$", xy=(1.07,1.70), color='black', fontweight='bold', alpha=0.5)
plt.annotate("Lower limit: "+str("%7.4f"%(STDDN))+" ng/m$^3$", xy=(1.07,1.65), color='black', fontweight='bold', alpha=0.5)
