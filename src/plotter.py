#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 13:14:40 2021

@author: Maureen Cohen
"""
import iris
import iris.plot as iplt
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
from matplotlib.colors import TwoSlopeNorm


reds = mpl_cm.get_cmap('brewer_Reds_09')

def plotter(cubes, time=-1, level=47):
    
    for cube in cubes:
        if cube.long_name == 'deep convection indicator':
            data = cube.copy()            
    
    run_length = np.arange(0,data.shape[0])
    longitudes = data.shape[2]
    latitudes = data.shape[1]
    
    plt.figure(figsize=(10,5))
    plt.contourf(np.arange(-longitudes/2, longitudes/2), np.arange(-latitudes/2, latitudes/2), np.roll(data[time,:,:].data, 72, axis=1), reds.N, cmap=reds)
    plt.title('Deep Convection Indicator, t=270.0 to 300.0 days')
    plt.xlabel('Longitude [degrees]')
    plt.ylabel('Latitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.yticks((-45,-30,-15,0,15,30,45),('90S','60S','30S','0','30N','60N','90N'))  
    plt.colorbar(pad=0.1)
    plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/laso/epsfigs/deepconvection.eps', format='eps')   
    plt.show()