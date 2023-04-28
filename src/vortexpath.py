#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 13:13:12 2023

@author: Mo Cohen
"""
import warnings, iris
import matplotlib.pyplot as plt
import numpy as np
from numpy import unravel_index

warnings.filterwarnings('ignore')

def vortexpath(cubes, start=0, end=-1, level=8, lat=73, meaning=3):
    
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube[start:end, :, :, :].copy()
            longterm_x_wind = cube.copy()
        if cube.standard_name == 'y_wind':
            y_wind = cube[start:end, :, :, :].copy()
            
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
    # Regrid y-wind onto coordinates of x-wind cube
    km_heights = np.round(x_wind.coord('level_height').points*1e-03, 2)
    # Extract heights in km for labelling plots
    latitudes = x_wind.coord('latitude').points
    longitudes = x_wind.coord('longitude').points
    # Extract latitudes and longitudes
    lat_deg = int(latitudes[lat])
    print('Calculating at ' + str(lat_deg))

    time_axis = np.arange(0, y_wind.shape[0]-meaning+1)
    
    v = y_wind[:, level, lat, 36:108].data
    plt.contourf(longitudes, latitudes, y_wind[end,level,:,:].data, cmap='RdBu_r')
    plt.show()
    print('Plotting just the meridional wind shows a sharp transition...')

    lon_deg = []
    lon_ind = []
    for ytime in range(0, v.shape[0]):
        y_ind = np.where(np.diff(np.sign(v[ytime, :])) == 2.)[0] + 36
        # Find where the sign of y_wind changes. We add 36 because we're only
        # looking at the nightside, so our 0 index is the 36 index in the
        # longitudes array.
        print(y_ind, longitudes[y_ind])
        # Prints out the longitude index where the y_wind changes sign, as well
        # as the longitude in degrees at that index

        if len(y_ind) == 1:
            lon_deg.append(longitudes[y_ind][0])
            lon_ind.append(y_ind[0])
        else:
            lon_deg.append(0)
            lon_ind.append(0)

    lon_deg = np.convolve(np.array(lon_deg).flatten(),
                                 np.ones(meaning), 'valid')/meaning
    # Applies a rolling mean with the period (in days) specified in the args
    longterm_mean = np.mean(lon_deg)

    fig, ax1 = plt.subplots(figsize=(6,4))
    ax1.set_xlabel('Time [days]', fontsize=14)
    ax1.set_ylabel('Longitude [deg E]', fontsize=14)
    ax1.plot(time_axis, np.ones_like(lon_deg)*longterm_mean, 
             color='r', label='Long-term mean')
    ax1.plot(time_axis, lon_deg, color='b', label='Gyre longitude')

    ax1.set_ylim(0,360)
    plt.title('Longitude of gyre centre at %s N' %lat_deg, fontsize=14)
    plt.legend()
    fig.tight_layout()
    plt.show()
    
    
    
    