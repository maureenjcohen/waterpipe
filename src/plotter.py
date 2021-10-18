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
redblu = mpl_cm.get_cmap('brewer_RdBu_11')


def plotter(cubes, time=-1, level=47):
    
    for cube in cubes:
        if cube.long_name == 'change_over_time_in_air_temperature_due_to_gravity_wave_drag':
            data = cube.copy()            

    run_length = np.arange(0,data.shape[0])
    longitudes = data.shape[3]
    latitudes = data.shape[2]
    heights = np.round(data.coord('level_height').points*1e-03,0)
    
    plt.figure(figsize=(10,5))
    plt.contourf(np.arange(-longitudes/2, longitudes/2), np.arange(-latitudes/2, latitudes/2), np.roll(data[time,level,:,:].data, 72, axis=1), redblu.N, cmap=redblu, norm=TwoSlopeNorm(0))
    plt.title('Temperature change due to gravity wave drag, h=%s km' %(heights[level]))
    plt.xlabel('Longitude [degrees]')
    plt.ylabel('Latitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.yticks((-45,-30,-15,0,15,30,45),('90S','60S','30S','0','30N','60N','90N'))  
    cbar = plt.colorbar(pad=0.1)
    cbar.ax.set_title('K')
    # plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/laso/epsfigs/deepconvection.eps', format='eps')   
    plt.show()
    
    
def line_plot(cubes, time_slice=-1, lat=45, level=74):
    
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube[time_slice,:,lat,:].copy()
    
    heights = np.round(x_wind.coord('level_height').points*1e-03, 0)
    longitudes = x_wind.shape[1]/2
    
    plt.figure(figsize=(10,5))
    plt.plot(np.arange(-longitudes, longitudes), np.roll(x_wind[level,:].data, 72, axis=0), color='b')
    plt.title('U at Equator, t=%s months, h=%s km' %(time_slice, heights[level]))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.ylabel('Wind speed [m s-1]')
    plt.show()