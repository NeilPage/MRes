#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May  2 10:50:22 2017

THIS SCRIPT PLOTS THE HG SURFACE FLUXES AT THE SOUTHERN HEMISPHERE SCALE
FOR EACH MONTH

@author: ncp532
"""
# File system packages
from netCDF4 import Dataset				# function used to open single netcdf file
from netCDF4 import MFDataset				# function used to open multiple netcdf files
import netCDF4
from pylab import *

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
dataset1 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Geogenic_Files/trac_avg.201*.nc') # InvOcean

# Get the variables needed for calculations
airdena = dataset1.variables['BXHGHT_S__AIRNUMDE'][:] # dry air number density in (molecules air/m3)
area = dataset1.variables['DXYP__DXYP'][:,0:45,:]            # Grid box surface area (m2)
timea = dataset1.variables['time'][:]                 # in minutes
leva = dataset1.variables['lev'][:]                   # atmospheric level (0 to 46)
lata = dataset1.variables['lat'][:]                   # latitude 

# Get the source data
#wet0 = (dataset1.variables['IJ_AVG_S__Hg0'][:])*(1e-3)*200.59/(6.0221*1e23)*airdena*1e6
anthro1 = (dataset1.variables['HG_SRCE__Hg0_an'][:,0:45,:])  # Hg0 anthro emissions (kg/month)
anthro2 = (dataset1.variables['HG_SRCE__Hg2_an'][:,0:45,:])  # Hg2 anthro emissions (kg/month)
anthro3 = (dataset1.variables['HG_SRCE__HgP_an'][:,0:45,:])  # HgP anthro emissions (kg/month)
anthro = (anthro1 + anthro2 + anthro3)                # anthro emissions (kg/month)
biomass = (dataset1.variables['HG_SRCE__Hg_bb'][:,0:45,:])   # biomass emissions (kg/month)
soil = (dataset1.variables['HG_SRCE__Hg_so'][:,0:45,:])      # soil emissions (kg/month)
geogenic = (dataset1.variables['HG_SRCE__Hg0_na'][:,0:45,:]) # land natural emissions (kg/month)
land_leg = (dataset1.variables['HG_SRCE__Hg0_ln'][:,0:45,:]) # land legacy/re- emissions (kg/month)
ocean = (dataset1.variables['HG_SRCE__Hg0_oc'][:,0:45,:])    # ocean net emissions (kg/month)
hg0_dd=(dataset1.variables['DRYD_FLX__Hg0df'][:,0:45,:])     # Hg0 dry deposition flux (molec/cm2/s)
hgII_1=(dataset1.variables['WETDLS_S__Hg2'][:,:,0:45,:])       # Loss of Hg2 in precipitation (kg/s)
hgII_2=(dataset1.variables['WETDCV_S__Hg2'][:,:,0:45,:])       # Rainout loss of tracer in convective updrafts (kg/s)
hgII_wd = (hgII_1 + hgII_2)                           # HgII wet deposition flux (kg/s)
snow=(dataset1.variables['HG_SRCE__Hg0_snow'][:,0:45,:])     # Snow emission of Hg (kg)
# Convert the dataset from ppt to ng/m3 = ppt*(10^-12)*(molecular mass Hg/molar volume Hg)*(molar volume Hg/Avogadro's number)*(10^9)*(air density in molecules/m3)
#------------------------------------------------------------------------------
days=np.array([31,28,31,30,31,30,31,31,30,31,30,31,
               31,28,31,30,31,30,31,31,30,31,30,31,
               31,28,31,30,31,30,31,31,30,31,30,31,
               31,29,31,30,31,30,31,31,30,31,30,31])
#------------------------------------------------------------------------------
# Convert the source data from units (kg/month) into units (kg/m2/month)
# NOTE: if you want to calculate the Hg budget then don't divide by area
# you need to multiple by area in the case of HG0_DD
anthro = anthro/area
biomass = biomass/area
soil = soil/area
geogenic = geogenic/area
land_leg = land_leg/area
ocean = ocean/area

# Convert the source data from units (molec/cm2/s) into units (kg/m2/month) 
# kg/m2/month = (molecules/cm2/s) * (1e4 cm2 / m2) /  (6.0221e23 molecules/mol)
#               * (200.59 g/mol) * (1e-3 kg/g) * (86400 s/d) * (days_in_month d/month)
HgTimeD, HglatD, HglonD = np.shape(hg0_dd)
HG0_DD = np.copy(hg0_dd)
for i in range(HgTimeD):
    HG0_DD[i,:,:] = hg0_dd[i,:,:]*days[i] # Overall

HG0_DD = HG0_DD[:,:,:] * (1e4) / (6.0221*1e23) * 200.59 * (1e-3)*86400#*area

# Convert the source data from units (kg/s) into units (kg/m2/month)
# kg/m2/month = (kg/s) / (grid_box_area m2 ) * (86400 s/d) * (no of days_in_month)
hgII_wd = np.sum(hgII_wd, axis=1)
HgTime, Hglat, Hglon = np.shape(hgII_wd)
HGII_WD = np.copy(hgII_wd)
for i in range(HgTime):
    HGII_WD[i,:,:] = hgII_wd[i,:,:]*days[i] # Overall

HGII_WD = HGII_WD[:,:,:]*86400/area

#------------------------------------------------------------------------------
# DATETIME

# InvOcean
# Convert GEOS-Chem time into standard time (GEOS-Chem time is hours since 1985/01/01/0000)
d0=datetime(1985,1,1,0)
datea = []
for t in timea:
    hrs=timedelta(hours=int(t))
    datea.append(d0+hrs)
#------------------------------------------------------------------------------
# Calculate the seasonal variation

# Anthro
monthavg1=np.zeros([12]) # Set the size
monthvals1= [ dt.month for dt in datea ] # Extract the month values from the date 
for monthind in np.arange(12):
    dataind=monthvals1 == monthind+1
    monthavg1[monthind] = np.mean(anthro[:,:,:][dataind]) # find the mean of months with same month value

# Biomass
monthavg2=np.zeros([12]) # Set the size
monthvals2= [ dt.month for dt in datea ] # Extract the month values from the date 
for monthind in np.arange(12):
    dataind=monthvals2 == monthind+1
    monthavg2[monthind] = np.mean(biomass[:,:,:][dataind]) # find the mean of months with same month value

# Soil
monthavg3=np.zeros([12]) # Set the size
monthvals3= [ dt.month for dt in datea ] # Extract the month values from the date 
for monthind in np.arange(12):
    dataind=monthvals3 == monthind+1
    monthavg3[monthind] = np.mean(soil[:,:,:][dataind]) # find the mean of months with same month value

# Geogenic
monthavg4=np.zeros([12]) # Set the size
monthvals4= [ dt.month for dt in datea ] # Extract the month values from the date 
for monthind in np.arange(12):
    dataind=monthvals4 == monthind+1
    monthavg4[monthind] = np.mean(geogenic[:,:,:][dataind]) # find the mean of months with same month value

# Land Legacy
monthavg5=np.zeros([12]) # Set the size
monthvals5= [ dt.month for dt in datea ] # Extract the month values from the date 
for monthind in np.arange(12):
    dataind=monthvals5 == monthind+1
    monthavg5[monthind] = np.mean(land_leg[:,:,:][dataind]) # find the mean of months with same month value

# Ocean
monthavg6=np.zeros([12]) # Set the size
monthvals6= [ dt.month for dt in datea ] # Extract the month values from the date 
for monthind in np.arange(12):
    dataind=monthvals6 == monthind+1
    monthavg6[monthind] = np.mean(ocean[:,:,:][dataind]) # find the mean of months with same month value

# Hg0 Dry Deposition
monthavg7=np.zeros([12]) # Set the size
monthvals7= [ dt.month for dt in datea ] # Extract the month values from the date 
for monthind in np.arange(12):
    dataind=monthvals7 == monthind+1
    monthavg7[monthind] = np.mean(HG0_DD[:,:,:][dataind]) # find the mean of months with same month value

# HgII Wet Deposition
monthavg8=np.zeros([12]) # Set the size
monthvals8= [ dt.month for dt in datea ] # Extract the month values from the date 
for monthind in np.arange(12):
    dataind=monthvals8 == monthind+1
    monthavg8[monthind] = np.mean(HGII_WD[:,:,:][dataind]) # find the mean of months with same month value
#------------------------------------------------------------------------------
# Calculate the Net Emissions Flux
NEF = (monthavg1+monthavg2+monthavg3+monthavg4+monthavg5+monthavg6)-(monthavg7+monthavg8)
#------------------------------------------------------------------------------
# Multiply by 1e8 easier for plotting
monthavg1=monthavg1*1e10
monthavg2=monthavg2*1e10
monthavg3=monthavg3*1e10
monthavg4=monthavg4*1e10
monthavg5=monthavg5*1e10
monthavg6=monthavg6*1e10
monthavg7=monthavg7*1e10
monthavg8=monthavg8*1e10
NEF=NEF*1e10
#------------------------------------------------------------------------------
# Sum the monthly average for the Annual Total
Sum1 = np.sum(anthro)/4
Sum2 = np.sum(biomass)/4
Sum3 = np.sum(soil)/4
Sum4 = np.sum(geogenic)/4
Sum5 = np.sum(land_leg)/4
Sum6 = np.sum(ocean)/4
Sum7 = np.sum(HG0_DD)/4
Sum8 = np.sum(HGII_WD)/4
Sum9 = np.sum(snow)/4
#------------------------------------------------------------------------------
# Plot the axis for each graph
fig = plt.figure()
plt.subplots_adjust(hspace=0.5)

# Graph 1
ax=plt.subplot(111) # options graph 2 (vertical no, horizontal no, graph no)

# SIMULATIONS
# Plot the monthly mean for the location
lines1 = ax.bar(datea[0:12], monthavg1, width=22, color='red', label="Anthropogenic")
lines2 = ax.bar(datea[0:12], monthavg2, width=22, color='blue', bottom=monthavg1, label="Biomass")
lines3 = ax.bar(datea[0:12], monthavg3, width=22, bottom=monthavg1+monthavg2, color='yellow', label="Soil")
lines4 = ax.bar(datea[0:12], monthavg4, width=22, color='green', bottom=monthavg1+monthavg2+monthavg3, label="Geogenic")
lines5 = ax.bar(datea[0:12], monthavg5, width=22, color='cyan', bottom=monthavg1+monthavg2+monthavg3+monthavg4, label="Land Legacy")
lines6 = ax.bar(datea[0:12], monthavg6, width=22, color='magenta', bottom=monthavg1+monthavg2+monthavg3+monthavg4+monthavg5, label="Ocean (net)")
lines7 = ax.bar(datea[0:12], monthavg7*(-1), width=22, color='orange', label="Deposition Dry$^*$")
lines8 = ax.bar(datea[0:12], monthavg8*(-1), width=22, bottom=monthavg7*(-1), color='darkslategray', label="Deposition Wet")
lines9 = plt.plot(datea[0:12],NEF, "o--",color='black', linewidth=1, label="Net Emissions Flux")

#plt.xlim(datetime(2013,01,1),datetime(2013,12,31))   # set the date limits
date_formatter = mdates.DateFormatter('%b')      # format how the date is displayed
ax.xaxis.set_major_formatter(date_formatter)
ax.xaxis.set_major_locator(mdates.MonthLocator(interval=1)) # set the interval between dispalyed dates
plt.axhline(0, color='black')
ax.yaxis.set_ticks_position('both')
ax.set_ylim(-7.0, 8) 
#ax.set_ylim(-0.5e-9, 1.0e-9) 
#ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(1))

# Plot the axis labels, legend and title
#plt.ylabel('Hg surface fluxes (1e-9 Kg/m2/month)', fontsize=10) 
plt.ylabel('Hg surface fluxes ($\mu$g/m2/month)', fontsize=10)
plt.xlabel('Month', fontsize=10)
plt.title("Southern Hemisphere Hg surface fluxes", fontsize=15)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

