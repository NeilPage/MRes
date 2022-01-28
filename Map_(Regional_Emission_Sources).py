#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 16:19:13 2017

THIS SCIPT PLOTS A MAP OF THE MEAN REGIONAL EMISSIONS FROM DIFFERENT SOURCES

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
from matplotlib.ticker import LogFormatter, LogFormatterMathtext
import matplotlib.gridspec as gridspec
# Date and Time handling package
from datetime import datetime,timedelta		# functions to handle date and time

# Data handing packages
import math						         # package to do mathematical calculations (such as log, exponential)
import numpy as np					    # import package as shorter nickname - Numpy as great at handling multidimensional data arrays.

# Define the dataset
#dataset1 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_Default_Files/trac_avg.201*.nc')	#define the dataset
dataset1 = MFDataset('/Users/ncp532/Documents/Some_scripts_/nc_InvOcean_Files/trac_avg.201*.nc')	#define the dataset

# Get the GEOS-Chem data for the stations
#airden1=dataset1.variables['BXHGHT_S__AIRNUMDE'][:] # for Hg0
airden1=dataset1.variables['BXHGHT_S__AIRNUMDE'][:,0,:,:] # Air Density
#wet1=(dataset1.variables['IJ_AVG_S__Hg0'][:])*(1e-3)*200.59/(6.0221*1e23)*airden1*1e6 # Hg0
wet1=(dataset1.variables['HG_SRCE__Hg0_an'][:])  # Hg0 anthro emissions (kg)
wet2=(dataset1.variables['HG_SRCE__Hg2_an'][:])  # Hg2 anthro emissions (kg)
wet3=(dataset1.variables['HG_SRCE__HgP_an'][:])  # HgP anthro emissions (kg)
wet4=(dataset1.variables['HG_SRCE__Hg_bb'][:])   # biomass emissions (kg)
wet5=(dataset1.variables['HG_SRCE__Hg_so'][:])   # soil emissions (kg)
wet6=(dataset1.variables['HG_SRCE__Hg0_na'][:])  # land natural emissions (kg)
wet7=(dataset1.variables['HG_SRCE__Hg0_ln'][:])  # land legacy/re- emissions (kg)
wet8=(dataset1.variables['HG_SRCE__Hg0_oc'][:])  # ocean net emissions (kg)
wet9=(dataset1.variables['HG_SRCE__Hg_vg'][:])   # vegetation emissions (kg)
time1=dataset1.variables['time'][:]                           # in minutes
lat1=dataset1.variables['lat'][:]                             # in degrees
lon1=dataset1.variables['lon'][:]                             # in degrees

# CALCULATE THE AVERAGE EMISSIONS FROM DIFFERENT SOURCES
# Anthropogenic Emissions
avg1=wet1+wet2 
wet1_avg=np.sum(avg1,axis=0)/4

# Biomass Burning Emisssions
avg2=wet4
wet2_avg=np.sum(avg2,axis=0)/4

# Soil Emissions
avg3=wet5
wet3_avg=np.sum(avg3,axis=0)/4

# Land Natural Emissions
avg4=wet5+wet6
wet4_avg=np.sum(avg4,axis=0)/4

# Land Legacy Emissions
avg5=wet7
wet5_avg=np.sum(avg5,axis=0)/4

# Ocean Net Emissions
avg6=wet8
wet6_avg=np.sum(avg6,axis=0)/4

# Vegetation Emissions
avg7=wet9
wet7_avg=np.sum(avg7,axis=0)/4

# Fix lat and long (GMAO 2x25 GRID)
lons_e = np.array([ -181.25, -178.75, -176.25, -173.75, -171.25, -168.75,
                    -166.25, -163.75, -161.25, -158.75, -156.25, -153.75,
                    -151.25, -148.75, -146.25, -143.75, -141.25, -138.75,
                    -136.25, -133.75, -131.25, -128.75, -126.25, -123.75,
                    -121.25, -118.75, -116.25, -113.75, -111.25, -108.75,
                    -106.25, -103.75, -101.25,  -98.75,  -96.25,  -93.75,
                    -91.25,  -88.75,  -86.25,  -83.75, -81.25,  -78.75,
                    -76.25,  -73.75,  -71.25,  -68.75,  -66.25,  -63.75,
                    -61.25,  -58.75,  -56.25,  -53.75,  -51.25,  -48.75,
                    -46.25,  -43.75,  -41.25,  -38.75,  -36.25,  -33.75,
                    -31.25,  -28.75,  -26.25,  -23.75, -21.25,  -18.75,
                    -16.25,  -13.75,  -11.25,   -8.75,   -6.25,   -3.75,
                    -1.25,    1.25,    3.75,    6.25,    8.75,   11.25,
                    13.75,   16.25,   18.75,   21.25,   23.75,   26.25,
                      28.75,   31.25,   33.75,   36.25,  38.75,   41.25,
                      43.75,   46.25,   48.75,   51.25,   53.75,   56.25,
                      58.75,   61.25,   63.75,   66.25,   68.75,   71.25,
                      73.75,   76.25,   78.75,   81.25,   83.75,   86.25,
                      88.75,   91.25,   93.75,   96.25,   98.75,  101.25,
                      103.75,  106.25,  108.75,  111.25,  113.75,  116.25,
                      118.75,  121.25,  123.75,  126.25,  128.75,  131.25,
                      133.75,  136.25,  138.75,  141.25,  143.75,  146.25,
                      148.75,  151.25,  153.75,  156.25,  158.75,  161.25,
                      163.75,  166.25,  168.75,  171.25,  173.75,  176.25,
                      178.75,])
lons_m= np.array([-180.0, -177.5, -175.0, -172.5, -170.0, -167.5, -165.0, -162.5,
                       -160.0, -157.5, -155.0, -152.5, -150.0, -147.5, -145.0, -142.5,
                       -140.0, -137.5, -135.0, -132.5, -130.0, -127.5, -125.0, -122.5,
                       -120.0, -117.5, -115.0, -112.5, -110.0, -107.5, -105.0, -102.5,
                       -100.0,  -97.5,  -95.0,  -92.5,  -90.0,  -87.5,  -85.0,  -82.5,
                       -80.0,  -77.5,  -75.0,  -72.5,  -70.0,  -67.5,  -65.0,  -62.5,
                       -60.0,  -57.5,  -55.0,  -52.5,  -50.0,  -47.5,  -45.0,  -42.5,
                       -40.0,  -37.5,  -35.0,  -32.5,  -30.0,  -27.5,  -25.0,  -22.5,
                       -20.0,  -17.5,  -15.0,  -12.5,  -10.0,   -7.5,   -5.0,   -2.5,
                       0.0,    2.5,    5.0,    7.5,   10.0,   12.5,   15.0,   17.5,
                       20.0,   22.5,   25.0,   27.5,   30.0,   32.5,   35.0,   37.5,
                       40.0,   42.5,   45.0,   47.5,   50.0,   52.5,   55.0,   57.5,
                       60.0,   62.5,   65.0,   67.5,   70.0,   72.5,   75.0,   77.5,
                       80.0,   82.5,   85.0,   87.5,   90.0,   92.5,   95.0,   97.5,
                       100.0,  102.5,  105.0,  107.5,  110.0,  112.5,  115.0,  117.5,
                       120.0,  122.5,  125.0,  127.5,  130.0,  132.5,  135.0,  137.5,
                       140.0,  142.5,  145.0,  147.5,  150.0,  152.5,  155.0,  157.5,
                       160.0,  162.5,  165.0,  167.5,  170.0,  172.5,  175.0,  177.5,
                       ])

lats_e=np.array([  -90.,  -89.,  -87.,  -85.,  -83.,  -81.,  -79.,  -77.,
                 -75.,  -73.,  -71.,  -69.,  -67.,  -65.,  -63.,  -61.,
                 -59.,  -57.,  -55.,  -53.,  -51.,  -49.,  -47.,  -45.,
                 -43.,  -41.,  -39.,  -37.,  -35.,  -33.,  -31.,  -29.,
                 -27.,  -25.,  -23.,  -21.,  -19.,  -17.,  -15.,  -13.,
                 -11.,   -9.,   -7.,   -5.,   -3.,   -1.,    1.,    3.,
                 5.,    7.,    9.,   11.,   13.,   15.,   17.,   19.,
                 21.,   23.,   25.,   27.,   29.,   31.,   33.,   35.,
                 37.,   39.,   41.,   43.,   45.,   47.,   49.,   51.,
                 53.,   55.,   57.,   59.,   61.,   63.,   65.,   67.,
                 69.,   71.,   73.,   75.,   77.,   79.,   81.,   83.,
                 85.,   87.,   89.,   90., ])

lats_m=np.array([  -89.5, -88.,  -86.,  -84.,  -82.,  -80.,  -78.,  -76.,
                -74.,  -72.,  -70.,  -68.,  -66.,  -64.,  -62.,  -60.,
                -58.,  -56.,  -54.,  -52.,  -50.,  -48.,  -46.,  -44.,
                -42.,  -40.,  -38.,  -36.,  -34.,  -32.,  -30.,  -28.,
                -26.,  -24.,  -22.,  -20.,  -18.,  -16.,  -14.,  -12.,
                -10.,   -8.,   -6.,   -4.,   -2.,    0.,    2.,    4.,
                6.,    8.,   10.,   12.,   14.,   16.,   18.,   20.,
                22.,   24.,   26.,   28.,   30.,   32.,   34.,   36.,
                38.,   40.,   42.,   44.,   46.,   48.,   50.,   52.,
                54.,   56.,   58.,   60.,   62.,   64.,   66.,   68.,
                70.,   72.,   74.,   76.,   78.,   80.,   82.,   84.,
                86.,   88.,   89.5, ])


# Convert GEOS-Chem time into standard time (GEOS-Chem time is hours since 1985/01/01/0000)
d0=datetime(1985,1,1,0)
date1 = []
for t in time1:
    hrs=timedelta(hours=int(t))
    date1.append(d0+hrs)

def plotkg(x, blah, tit):
#def plotkg(x, tit):    
    lons, lats = np.meshgrid(lons_e,lats_e)
    m=Basemap(llcrnrlat=-47,  urcrnrlat=-2,        # set the latitude limits
              llcrnrlon=10, urcrnrlon=160,         # set the longitude limits
              resolution='c',projection='merc')
        # np.transpose because it mixes the lat lon order
    
    mlons,mlats = m(lons, lats) # convert lons/lat from middle of grid-box to to corner of grid-box
    
    cs=m.pcolormesh(mlons, mlats, x, latlon=False, norm=colors.LogNorm(),cmap=cm.plasma)   # define the colorbar attributes (colors used, add "_r" to flip the color scale)
    #cs=m.pcolormesh(mlons, mlats, x, latlon=False, vmin=1, cmap=cm.plasma)   # define the colorbar attributes (colors used, add "_r" to flip the color scale)    
    formatter = LogFormatterMathtext(10, labelOnlyBase=False) 
    m.drawcoastlines()
    
    if blah=='Y':
        cb=m.colorbar(cs,"bottom",size="10%", pad="30%",ticks=[0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000],format=formatter)
        cb.set_label(r"$kg$ $year$$^-$$^1$", fontsize =10)
        cb.ax.tick_params(labelsize=10)
        cb.ax.minorticks_on()
    
    #cb=m.colorbar(cs,"right",size="5%", pad="2%")
    #cb.set_label(r"$kg$", fontsize =10)
    #cb.ax.tick_params(labelsize=10)
    parallels = np.arange(-45,-2, 5.)
    # Plot parallels and meridians
    m.drawparallels(parallels,labels=[True,False,True,False], fontsize=10)
    meridians = np.arange(15,160,15.)
    m.drawmeridians(meridians,labels=[True,False,False,True], fontsize=10)
    
    plt.xlabel('Longitude [$^\circ$]', fontsize=10, labelpad=15)
    plt.ylabel('Latitude [$^\circ$]', fontsize=10, labelpad=40)
    plt.title(tit, fontsize=15)
    #ax1.xaxis.set_label_coords(0.5, -0.1)
    #ax1.yaxis.set_label_coords(-0.07, 0.5)
    return m

fig = plt.figure()
gs = gridspec.GridSpec(3,2)
#height_ratios=[3, 3, 6])
                       #height_ratios=[1.55, 1.55, 1.55, 2])
                       #height_ratios=[3, 3, 3])

#MAP 1
ax1=plt.subplot(gs[0,0]) # options graph 1 (vertical no, horizontal no, graph no)
ax1 = plotkg(wet1_avg[:,:], 'Y', "Anthropogenic Emissions")
#ax1 = plotkg(wet1_avg[:,:], "TGM Anthropogenic Emissions")

#MAP 2
ax2=plt.subplot(gs[0,1]) # options graph 2 (vertical no, horizontal no, graph no)
ax2 = plotkg(wet2_avg[:,:], 'Y', "Biomass Burning Emissions")
#ax2 = plotkg(wet2_avg[:,:],"TGM Biomass Burning Emissions")

#MAP 3
ax3=plt.subplot(gs[1,0]) # options graph 3 (vertical no, horizontal no, graph no)
ax3 = plotkg(wet3_avg[:,:], 'Y', "Soil Emissions")
#ax3 = plotkg(wet3_avg[:,:],"TGM Soil Emissions")

#MAP 4
#ax4=plt.subplot(gs[1,1]) # options graph 4 (vertical no, horizontal no, graph no)
#ax4 = plotkg(wet4_avg[:,:], 'Y', "Land Natural Emissions")
#ax4 = plotkg(wet4_avg[:,:], "TGM Land Natural Emissions")

#MAP 5
ax4=plt.subplot(gs[1,1]) # options graph 5 (vertical no, horizontal no, graph no)
ax4 = plotkg(wet5_avg[:,:], 'Y', "Land Legacy Emissions")
#ax5 = plotkg(wet5_avg[:,:],"TGM Land Legacy Emissions")

#MAP 6
ax5=plt.subplot(gs[2,:]) # options graph 6 (vertical no, horizontal no, graph no)
ax5 = plotkg(wet6_avg[:,:], 'Y', "Ocean Net Emissions")
#ax6 = plotkg(wet6_avg[:,:],"TGM Ocean Net Emissions")

#MAP 7
#ax7=plt.subplot(gs[3,0]) # options graph 7 (vertical no, horizontal no, graph no)
#ax7 = plotkg(wet7_avg[:,:], 'Y', "TGM Vegetation Emissions")
#ax7 = plotkg(wet7_avg[:,:],"TGM Vegetation Emissions")

plt.subplots_adjust(hspace=0.4)