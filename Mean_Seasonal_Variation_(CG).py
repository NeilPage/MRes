#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May  2 10:50:22 2017

THIS SCRIPT PLOTS A LINEGRAPH OF THE MEAN SEASONAL VARIATION

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
dataset7 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Horowitz/Offline_MITgcm_ocean*.nc') # NewChem/NewOcean                    
dataset8 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Horowitz/SlabOcean.201*.nc') # NewChem/SlabOcean
dataset9 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Horowitz/AIRDENS.2010*.nc') # AIRDENS
#dataset10 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Bromine_x3/coards.201*.nc')        # Bromine_x3

# OBSERVATION
# Retrieve the observation data
x1 = np.genfromtxt("/Users/ncp532/Documents/Hg_data/Hg0SiteDaily_AMS_2013-2016_Daily_YearlyMean.csv", delimiter=",", skip_header=1, usecols=(0), invalid_raise=False) # day number
y1 = np.genfromtxt("/Users/ncp532/Documents/Hg_data/Hg0SiteDaily_AMS_2013-2016_Daily_YearlyMean.csv", delimiter=",", skip_header=1, usecols=(1), invalid_raise=False) # Hg Concentration

# Default
# Get the GEOS-Chem data for the stations
airdena=dataset1.variables['TIME_SER__AIRDEN'][:,0,6,56]                           # dry air number density in (molecules air/m3)
wet1=(dataset1.variables['IJ_AVG_S__Hg0'][:,0,6,56])*(1e-3)*200.59/(6.0221*1e23)*airdena*1e6 # Hg0
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

# NewChem/NewOcean
# Get the GEOS-Chem data for the stations
airden3=dataset9.variables['BXHGHT_S__AIRNUMDE'][:,0,12,67]                           # dry air number density in (molecules air/m3)
wet7=(dataset7.variables['IJ_AVG_S__Hg0'][:,0,12,67])*(1e-3)*200.59/(6.0221*1e23)*airden3 #cm3 to m3 is a factor of 1e-6 #(2.5*1e25)
# Convert the dataset from ppt) to ng/m3 = ppt*(10^-12)*(molecular mass Hg/molar volume Hg)*(molar volume Hg/Avogadro's number)*(10^9)*(air density in molecules/m3)
timeg=dataset7.variables['time'][:]                           # in minutes

# NewChem/SlabOcean
# Get the GEOS-Chem data for the stations
wet8=(dataset8.variables['IJ_AVG_S__Hg0'][:,0,12,67])*(1e-3)*200.59/(6.0221*1e23)*airden3 #cm3 to m3 is a factor of 1e-6 #(2.5*1e25)
# Convert the dataset from ppt) to ng/m3 = ppt*(10^-12)*(molecular mass Hg/molar volume Hg)*(molar volume Hg/Avogadro's number)*(10^9)*(air density in molecules/m3)

# Bromine_x3
# Get the GEOS-Chem data for the stations
#airdeni=dataset10.variables['TIME_SER__AIRDEN'][:,0,6,56]                           # dry air number density in (molecules air/m3)
#wet9=(dataset10.variables['IJ_AVG_S__Hg0'][:,0,6,56])*(1e-3)*200.59/(6.0221*1e23)*airdeni*1e6
#------------------------------------------------------------------------------
# Default
# Convert GEOS-Chem time into standard time (GEOS-Chem time is hours since 1985/01/01/0000)
d0=datetime(1985,1,1,0)
datea = []
for t in timea:
    hrs=timedelta(hours=int(t))
    datea.append(d0+hrs)
    
# NewChem/NewOcean
# Convert GEOS-Chem time into standard time (GEOS-Chem time is hours since 1985/01/01/0000)
d0=datetime(1989,1,1,0)
dateg = []
for t in timeg:
    hrs=timedelta(hours=int(t))
    dateg.append(d0+hrs)

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
#wet9_Mavg, date_avg=monthly(wet9,datea) # Bromine_x3
wetOBS_Mavg, date_avg=monthly(y1,datea)           # Observations
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
wet3_Yavg, date_avg=yearly(wet3,datea) # Anthro  
wet4_Yavg, date_avg=yearly(wet4,datea) # Biomass
wet5_Yavg, date_avg=yearly(wet5,datea) # Geogenic    
wet6_Yavg, date_avg=yearly(wet6,datea) # Soil
#wet9_Yavg, date_avg=yearly(wet9,datea) # Bromine_x3
wetOBS_Yavg, date_avg=yearly(y1,datea)           # Observations
#------------------------------------------------------------------------------    
# Clean up the data (keep outliers)
finiteY1mask = np.isfinite(y1) # Scan for NaN values
Y1cln = y1[finiteY1mask]     # Remove NaN values from the obs

X1cln = wet1[finiteY1mask]
X2cln = wet2[finiteY1mask]
X3cln = wet3[finiteY1mask]
X4cln = wet4[finiteY1mask]
X5cln = wet5[finiteY1mask]
X6cln = wet6[finiteY1mask]
#X9cln = wet9[finiteY1mask]
#------------------------------------------------------------------------------
# keep only the Good Values in the observations (remove outliers)
OM_OBS = np.nanmean(y1)  # mean of Observations
OSTD_OBS = np.nanstd(y1) # Std of Observations

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
#------------------------------------------------------------------------------
# Calculate the overall mean
OM_1 = np.nanmean(wet1)   # Default
OM_2 = np.nanmean(wet2)   # InvOcean
OM_3 = np.nanmean(wet3)   # Anthro
OM_4 = np.nanmean(wet4)   # Biomass
OM_5 = np.nanmean(wet5)   # Geogenic
OM_6 = np.nanmean(wet6)   # Soil
OM_7 = np.nanmean(wet7) # NewChem/NewOcean
OM_8 = np.nanmean(wet8) # NewChem/SlabOcean
#OM_9 = np.nanmean(wet9) # Bromine_x3
OM_OBS = np.nanmean(y1)            # Observations
#------------------------------------------------------------------------------
# Calculate the standard deviation (of the daily variability)
# Default
std1=np.zeros([12]) # Set the size
stdev1= [ dt.month for dt in datea ] # Extract the month values from the date 
for monthind in np.arange(12):
    dataind=stdev1 == monthind+1
    std1[monthind] = np.std(wet1[dataind])

# InvOcean
std2=np.zeros([12]) # Set the size
stdev2= [ dt.month for dt in datea ] # Extract the month values from the date 
for monthind in np.arange(12):
    dataind=stdev2 == monthind+1
    std2[monthind] = np.std(wet2[dataind])

# Anthro
std3=np.zeros([12]) # Set the size
stdev3= [ dt.month for dt in datea ] # Extract the month values from the date 
for monthind in np.arange(12):
    dataind=stdev3 == monthind+1
    std3[monthind] = np.std(wet3[dataind])

# Biomass
std4=np.zeros([12]) # Set the size
stdev4= [ dt.month for dt in datea ] # Extract the month values from the date 
for monthind in np.arange(12):
    dataind=stdev4 == monthind+1
    std4[monthind] = np.std(wet4[dataind])

# Geogenic
std5=np.zeros([12]) # Set the size
stdev5= [ dt.month for dt in datea ] # Extract the month values from the date 
for monthind in np.arange(12):
    dataind=stdev5 == monthind+1
    std5[monthind] = np.std(wet5[dataind])

# Soil
std6=np.zeros([12]) # Set the size
stdev6= [ dt.month for dt in datea ] # Extract the month values from the date 
for monthind in np.arange(12):
    dataind=stdev6 == monthind+1
    std6[monthind] = np.std(wet6[dataind])
    
# OBSERVATIONS
std7=np.zeros([12]) # Set the size
stdev7= [ dt.month for dt in datea ] # Extract the month values from the date 
for monthind in np.arange(12):
    dataind=stdev7 == monthind+1
    std7[monthind] = np.nanstd(y1[dataind]) # find the mean of months with same month value

# Bromine_x3
#std9=np.zeros([12]) # Set the size
#stdev9= [ dt.month for dt in datea ] # Extract the month values from the date 
#for monthind in np.arange(12):
#    dataind=stdev9 == monthind+1
#    std9[monthind] = np.std(wet9[dataind])    
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

# Anthro
monthavg3=np.zeros([12]) # Set the size
monthvals3= [ dt.month for dt in monthsa ] # Extract the month values from the date 
for monthind in np.arange(12):
    dataind=monthvals3 == monthind+1
    monthavg3[monthind] = np.mean(wet3_Mavg[dataind]) # find the mean of months with same month value

# Biomass
monthavg4=np.zeros([12]) # Set the size
monthvals4= [ dt.month for dt in monthsa ] # Extract the month values from the date 
for monthind in np.arange(12):
    dataind=monthvals4 == monthind+1
    monthavg4[monthind] = np.mean(wet4_Mavg[dataind]) # find the mean of months with same month value

# Geogenic
monthavg5=np.zeros([12]) # Set the size
monthvals5= [ dt.month for dt in monthsa ] # Extract the month values from the date 
for monthind in np.arange(12):
    dataind=monthvals5 == monthind+1
    monthavg5[monthind] = np.mean(wet5_Mavg[dataind]) # find the mean of months with same month value

# Soil
monthavg6=np.zeros([12]) # Set the size
monthvals6= [ dt.month for dt in monthsa ] # Extract the month values from the date 
for monthind in np.arange(12):
    dataind=monthvals6 == monthind+1
    monthavg6[monthind] = np.mean(wet6_Mavg[dataind]) # find the mean of months with same month value
    
# OBSERVATIONS
monthavg7=np.zeros([12]) # Set the size
monthvals7= [ dt.month for dt in monthsa ] # Extract the month values from the date 
for monthind in np.arange(12):
    dataind=monthvals7 == monthind+1
    monthavg7[monthind] = np.nanmean(wetOBS_Mavg[dataind]) # find the mean of months with same month value

# Bromine_x3
#monthavg9=np.zeros([12]) # Set the size
#monthvals9= [ dt.month for dt in monthsa ] # Extract the month values from the date 
#for monthind in np.arange(12):
#    dataind=monthvals9 == monthind+1
#    monthavg9[monthind] = np.mean(wet9_Mavg[dataind]) # find the mean of months with same month value
#------------------------------------------------------------------------------
# Calculate the mean bias (between the simulation and observation)
MB1 = np.mean(X1clean - Y1clean) # Default
MB2 = np.mean(X2clean - Y1clean) # InvOcean
MB3 = np.mean(X3clean - Y1clean) # Anthro
MB4 = np.mean(X4clean - Y1clean) # Biomass
MB5 = np.mean(X5clean - Y1clean) # Geogenic
MB6 = np.mean(X6clean - Y1clean) # Soil
MB7 = np.mean(wet7 - OM_OBS) # NewChem/NewOcean
MB8 = np.mean(wet8 - OM_OBS) # NewChem/SlabOcean
#MB9 = np.mean(X9clean - Y1clean) # Bromine_x3

#------------------------------------------------------------------------------
# Calculate the mean normalised bias (between simulation and observation)
# Without outliers
MNB1 = (np.mean((X1clean - Y1clean)/Y1clean))*100 # Default
MNB2 = (np.mean((X2clean - Y1clean)/Y1clean))*100 # InvOcean
MNB3 = (np.mean((X3clean - Y1clean)/Y1clean))*100 # Anthro
MNB4 = (np.mean((X4clean - Y1clean)/Y1clean))*100 # Biomass
MNB5 = (np.mean((X5clean - Y1clean)/Y1clean))*100 # Geogenic
MNB6 = (np.mean((X6clean - Y1clean)/Y1clean))*100 # Soil
#NMB9 = (np.mean((X9clean - Y1clean)/Y1clean)*100 # Bromine_x3

# Including outliers
MNB11 = (np.mean((X1cln - Y1cln)/Y1cln))*100 # Default
MNB22 = (np.mean((X2cln - Y1cln)/Y1cln))*100 # InvOcean
MNB33 = (np.mean((X3cln - Y1cln)/Y1cln))*100 # Anthro
MNB44 = (np.mean((X4cln - Y1cln)/Y1cln))*100 # Biomass
MNB55 = (np.mean((X5cln - Y1cln)/Y1cln))*100 # Geogenic
MNB66 = (np.mean((X6cln - Y1cln)/Y1cln))*100 # Soil
MNB77 = (np.mean((wet7 - OM_OBS)/OM_OBS))*100 # NewChem/NewOcean
MNB88 = (np.mean((wet8 - OM_OBS)/OM_OBS))*100 # NewChem/SlabOcean
#NMB99 = (np.mean((X9cln - Y1cln)/Y1cln)*100 # Bromine_x3         
#------------------------------------------------------------------------------
# Calculate the mean normalised mean error (between simulation and observation)
# Without outliers
MNE1 = (np.mean(abs(X1clean - Y1clean)/Y1clean))*100 # Default
MNE2 = (np.mean(abs(X2clean - Y1clean)/Y1clean))*100 # InvOcean
MNE3 = (np.mean(abs(X3clean - Y1clean)/Y1clean))*100 # Anthro
MNE4 = (np.mean(abs(X4clean - Y1clean)/Y1clean))*100 # Biomass
MNE5 = (np.mean(abs(X5clean - Y1clean)/Y1clean))*100 # Geogenic
MNE6 = (np.mean(abs(X6clean - Y1clean)/Y1clean))*100 # Soil
#NMB9 = (np.mean((X9clean - Y1clean)/Y1clean)*100 # Bromine_x3

# Including outliers
MNE11 = (np.mean(abs(X1cln - Y1cln)/Y1cln))*100 # Default
MNE22 = (np.mean(abs(X2cln - Y1cln)/Y1cln))*100 # InvOcean
MNE33 = (np.mean(abs(X3cln - Y1cln)/Y1cln))*100 # Anthro
MNE44 = (np.mean(abs(X4cln - Y1cln)/Y1cln))*100 # Biomass
MNE55 = (np.mean(abs(X5cln - Y1cln)/Y1cln))*100 # Geogenic
MNE66 = (np.mean(abs(X6cln - Y1cln)/Y1cln))*100 # Soil
MNE77 = (np.mean(abs(wet7 - OM_OBS)/OM_OBS))*100 # NewChem/NewOcean
MNE88 = (np.mean(abs(wet8 - OM_OBS)/OM_OBS))*100 # NewChem/SlabOcean
#NMB99 = (np.mean((X9cln - Y1cln)/Y1cln)*100 # Bromine_x3
#------------------------------------------------------------------------------
# Mean difference between the InvOcean simulation and sensitivity simulations

MD1 = np.mean(((wet1 - wet2) / wet2) * 100) # Default
#MD2 = np.mean(((wet2 - wet2) / wet2) * 100) # InvOcean
MD3 = np.mean(((wet3 - wet2) / wet2) * 100) # Anthro
MD4 = np.mean(((wet4 - wet2) / wet2) * 100) # Biomass
MD5 = np.mean(((wet5 - wet2) / wet2) * 100) # Geogenic
MD6 = np.mean(((wet6 - wet2) / wet2) * 100) # Soil
#MD7 = np.mean(((wet7 - wet2) / wet2) * 100) # NewChem/NewOcean
#MD8 = np.mean(((wet8 - wet2) / wet2) * 100) # NewChem/SlabOcean
#MD9 = np.mean(((wet9 - wet2) / wet2) * 100) # Bromine_x3
#------------------------------------------------------------------------------
# Mean contribution of the sensitivity simulations to InvOcean

MC1 = np.mean((wet1 / wet2) * 100) # Default
#MC2 = np.mean((wet2 / wet2) * 100) # InvOcean
MC3 = np.mean((wet3 / wet2) * 100) # Anthro
MC4 = np.mean((wet4 / wet2) * 100) # Biomass
MC5 = np.mean((wet5 / wet2) * 100) # Geogenic
MC6 = np.mean((wet6 / wet2) * 100) # Soil
#MC7 = np.mean((wet7 / wet2) * 100) # NewChem/NewOcean
#MC8 = np.mean((wet8 / wet2) * 100) # NewChem/SlabOcean
#MC9 = np.mean((wet9 / wet2) * 100) # Bromine_x3
#------------------------------------------------------------------------------
# Plot the axis for each graph
fig = plt.figure()
plt.subplots_adjust(hspace=0.5)

# Graph 1
ax=plt.subplot(211) # options graph 2 (vertical no, horizontal no, graph no)

# OBSERVATIONS
# Plot the Monthly/Yearly mean for the location
plt.plot(dateg,monthavg7, 'k', linewidth=5, label="Observations")
plt.errorbar(dateg, monthavg7, yerr=std7, fmt='k', capsize=3, linewidth=2)

# SIMULATIONS
# Plot the monthly mean for the location
plt.plot(dateg, monthavg1, color='red', linewidth=3, label="Normal \n  (MNB: "+str("%6.3f"%(MNB1))+"%)\n  (MNB*: "+str("%6.3f"%(MNB11))+"%)\n  (MNE: "+str("%6.3f"%(MNE1))+"%)\n  (MNE*: "+str("%6.3f"%(MNE11))+"%)")       # wet1_Mavg
plt.errorbar(dateg, monthavg1, yerr=std1, ls='--', color='red', capsize=3, linewidth=2)
plt.plot(dateg, monthavg2, color='blue', linewidth=3, label="InvOcean \n  (MNB: "+str("%6.3f"%(MNB2))+"%)\n  (MNB*: "+str("%6.3f"%(MNB22))+"%)\n  (MNE: "+str("%6.3f"%(MNE2))+"%)\n  (MNE*: "+str("%6.3f"%(MNE22))+"%)")     # wet2_Mavg
plt.errorbar(dateg, monthavg2, yerr=std2, ls='--', color='blue', capsize=3, linewidth=2)

#plt.plot(dateg, wet7, ls='--', color='green', linewidth=1, label="NewChem/NewOcean \n  (MNB: "+str("%6.3f"%(MNB77))+"%)\n  (MNE: "+str("%6.3f"%(MNE77))+"%)")   # wet7_Mavg
#plt.plot(dateh, wet8, ls='--', color='darkorange', linewidth=1, label="NewChem/SlabOcean \n  (MNB: "+str("%6.3f"%(MNB88))+"%)\n  (MNE: "+str("%6.3f"%(MNE88))+"%)")  # wet8_Mavg

#plt.xlim(datetime(2014,1,1),datetime(2014,12,31))   # set the date limits
plt.margins(0.05, 0.05)
plt.xticks(rotation=15)
date_formatter = mdates.DateFormatter('%b')      # format how the date is displayed
ax.xaxis.set_major_formatter(date_formatter)
ax.xaxis.set_major_locator(mdates.MonthLocator(interval=1)) # set the interval between dispalyed dates
ax.yaxis.set_ticks_position('both')


# Plot the axis labels, legend and title
plt.ylabel('Hg$^0$ (ng/m$^3$)', fontsize=10)
plt.title("Mean seasonal variation and temporal standard deviation of Hg$^0$ at Cape Grim, Tas", fontsize=15)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

# Graph 2
ax=plt.subplot(212) # options graph 2 (vertical no, horizontal no, graph no)

# OBSERVATIONS
# Plot the Monthly/Yearly mean for the location
plt.plot(dateg,monthavg7, 'k', linewidth=1, label="Observations")
plt.errorbar(dateg, monthavg7, yerr=std7, fmt='k', capsize=2, linewidth=1)

# SIMULATIONS
# Plot the monthly mean for the location
plt.plot(dateg, monthavg2, ls='--', color='blue', linewidth=1, label="InvOcean")                                                           # InvOcean
plt.errorbar(dateg, monthavg2, yerr=std2, ls='--', color='blue', capsize=2, linewidth=1)
plt.plot(dateg, monthavg3, ls='--', color='cyan', linewidth=1, label="Anthro \n  (Mean Difference: "+str("%6.3f"%(MD3))+"%)\n  (Mean Contribution: "+str("%6.3f"%(MC3))+"%)")          # Anthro
plt.errorbar(dateg, monthavg3, yerr=std3, ls='--', color='cyan', capsize=2, linewidth=1)
plt.plot(dateg, monthavg4, ls='--', color='magenta', linewidth=1, label="Biomass \n(Mean Difference: "+str("%6.3f"%(MD4))+"%)\n  (Mean Contribution: "+str("%6.3f"%(MC4))+"%)") # Biomass
plt.errorbar(dateg, monthavg4, yerr=std4, ls='--', color='magenta', capsize=2, linewidth=1)
plt.plot(dateg, monthavg5, ls='--', color='darkslategray', linewidth=1, label="Geogenic \n  (Mean Difference: "+str("%6.3f"%(MD5))+"%)\n  (Mean Contribution: n/a)")     # Geogenic
plt.errorbar(dateg, monthavg5, yerr=std5, ls='--', color='darkslategray', capsize=2, linewidth=1)
plt.plot(dateg, monthavg6, ls='--', color='saddlebrown', linewidth=1, label="Soil \n  (Mean Difference: "+str("%6.3f"%(MD6))+"%)\n  (Mean Contribution: "+str("%6.3f"%(MC6))+"%)")     # Soil
plt.errorbar(dateg, monthavg6, yerr=std6, ls='--', color='saddlebrown', capsize=2, linewidth=1)
#plt.plot(dateg, monthavg9*100, ls='--', color='orange', linewidth=1, label="Bromine_x3 \n(Mean Difference: "+str("%6.3f"%(MD9))+"%)")      # Bromine_x3
#plt.errorbar(dateg, monthavg9*100, yerr=std9, ls='--', color='orange', capsize=2, linewidth=1)

#plt.xlim(datetime(2014,1,1),datetime(2014,12,31))   # set the date limits
plt.margins(0.05, 0.05)
plt.xticks(rotation=15)
date_formatter = mdates.DateFormatter('%b')      # format how the date is displayed
ax.xaxis.set_major_formatter(date_formatter)
ax.xaxis.set_major_locator(mdates.MonthLocator(interval=1)) # set the interval between dispalyed dates
ax.yaxis.set_ticks_position('both')

# Plot the axis labels, legend and title
plt.ylabel('Hg$^0$ (ng/m$^3$)', fontsize=10)
plt.title("Senitivity simulations", fontsize=15)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
