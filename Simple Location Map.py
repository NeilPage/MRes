#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 16:19:59 2017

THIS SCRIPT CREATES A SIMPLE LOCATION MAP 

@author: ncp532
"""

import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from datetime import datetime
# miller projection
map = Basemap(projection='mill')
# plot coastlines, draw label meridians and parallels.
map.drawcoastlines()
map.drawparallels(np.arange(-90,90,15),labels=[1,0,0,0], fontsize=15)
map.drawmeridians(np.arange(map.lonmin,map.lonmax+30,60),labels=[0,0,0,1], fontsize=15)
# fill continents 'coral' (with zorder=0), color wet areas 'aqua'
map.drawmapboundary(fill_color='aqua')
map.fillcontinents(color='coral',lake_color='aqua', zorder=0)
lons = [144.6832, 131.0447, 151.1018, 77.3423, 18.4897]
lats = [-40.6832, -12.2491, -32.4777, -37.4748, -34.3535]
x,y = map(lons, lats)
map.scatter(x, y, c='white', edgecolors='black', s=200)
    
# Plot location labels
labels = ['Cape Grim', 'Gunn Point', 'Glenville', 'Amsterdam Island', 'Cape Point']
x_offsets = [-800000, -5500000, -750000, -2000000, -900000]
y_offsets = [-950000, -350000, 350000, -900000, -900000]
for label, xpt, ypt, x_offset, y_offset in zip(labels, x, y, x_offsets, y_offsets):
    plt.text(xpt+x_offset, ypt+y_offset, label, fontsize=15, color='black', fontweight='bold')
plt.xlabel('Longitude [$^\circ$]', fontsize=15, labelpad=30)
plt.ylabel('Latitude [$^\circ$]', fontsize=15, labelpad=50)
plt.title('Locations of the ground based monitoring sites', fontsize=25)
plt.show()