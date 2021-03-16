#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 10 16:56:00 2021

@author: Mo Cohen

Plotting functions for studying atmospheric variability

"""

import os, iris
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import iris.plot as iplt
import iris.coords
from iris.coords import DimCoord
import numpy as np
import scipy as sp
import windspharm
# Import packages


brewer_red = mpl_cm.get_cmap('brewer_YlOrRd_09')
brewer_bg = mpl_cm.get_cmap('brewer_PuBu_09')
brewer_reds = mpl_cm.get_cmap('brewer_Reds_09')
brewer_redblu = mpl_cm.get_cmap('brewer_RdBu_11')
# Load colormaps to use in plots

def get_heights(cubes):
    
    """ Extract altitude information from air pressure cube
        Can be used with matplotlib contourf when Iris isn't playing nice """ 
    
    for cube in cubes:
        if cube.standard_name == 'air_pressure':
            pressure = cube.copy()
            
    heights = pressure.coord('level_height').points
    
    return heights
    

def cloudswirls(cubes):
    
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube.copy()
        if cube.long_name == 'cloud_area_fraction_assuming_maximum_random_overlap':
            cloud_cover = cube.copy()
    
    height = x_wind.shape[1]
    run_length = x_wind.shape[0]
    lats = x_wind.coord('latitude')
    longs = x_wind.coord('longitude')

    if lats.bounds == None:
        x_wind.coord('latitude').guess_bounds()
    if longs.bounds == None:
        x_wind.coord('longitude').guess_bounds()
        
    strat = x_wind.extract(iris.Constraint(longitude=lambda v: 355 < v <= 359 or 0 <= v <= 4, latitude=lambda v: -4 <= v <= 4, model_level_number=40))
    trop = x_wind.extract(iris.Constraint(longitude=lambda v: 355 < v <= 359 or 0 <= v <= 4, latitude=lambda v: -4 <= v <= 4, model_level_number=24))
#    high = x_wind.extract(iris.Constraint(longitude=lambda v: 355 < v <= 359 or 0 <= v <= 4, latitude=lambda v: -4 <= v <= 4, model_level_number=50))
#    low = x_wind.extract(iris.Constraint(longitude=lambda v: 355 < v <= 359 or 0 <= v <= 4, latitude=lambda v: -4 <= v <= 4, model_level_number=29))
    strat_grid = iris.analysis.cartography.area_weights(strat)
    trop_grid = iris.analysis.cartography.area_weights(trop)
    # high_grid = iris.analysis.cartography.area_weights(high)
    # low_grid = iris.analysis.cartography.area_weights(low)
    
    strat_mean = strat.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=strat_grid)
    trop_mean = trop.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=trop_grid)
    
    # high_mean = high.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=high_grid)
    # low_mean = low.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=low_grid)
    
    clear_list = []
    for day in range(0,cloud_cover.shape[0]):
        data = (cloud_cover[day,:,36].data + cloud_cover[day,:,108].data)/2
        clear = data[np.where( data <= 0.2)]
        count = len(clear.tolist())
        clear_list.append(count)
        
    time_axis = np.arange(0,cloud_cover.shape[0])
    plt.plot(time_axis,(np.array(clear_list)/90)*100)
#    plt.plot(np.arange(0,run_length), high_mean.data, linestyle='--', color='r', label='50 km')
    plt.plot(np.arange(0,run_length), strat_mean.data, linestyle='-', color='r', label='40 km')
#    plt.plot(np.arange(0,run_length), low_mean.data, linestyle='--', color='b', label='30 km')
    plt.plot(np.arange(0,run_length), trop_mean.data, linestyle='-', color='b', label='25 km')
    plt.title('Cloud Variability')
    plt.xlabel('Time [months]') 
    plt.ylabel('Velocity [m s-1]')
    plt.legend()
    plt.show()

