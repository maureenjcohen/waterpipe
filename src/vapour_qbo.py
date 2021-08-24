#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 10:38:54 2021

@author: Maureen Cohen
"""

import iris
import iris.plot as iplt
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
from matplotlib.colors import TwoSlopeNorm
from iris.coord_systems import GeogCS


bg = mpl_cm.get_cmap('brewer_PuBu_09')


def vapour_series(cubes, radius=7160000, level=47, x=(0,4), y=(43,47)):
    
    for cube in cubes:
        if cube.standard_name == 'specific_humidity':
            vapour = cube.copy()
            
    vapour.coord('latitude').coord_system = GeogCS(radius)
    vapour.coord('longitude').coord_system = GeogCS(radius)
    # Sets planet radius in m for area-weighted average. Default is radius of Proxima Centauri b
    
    run_length, height = vapour.shape[0], vapour.shape[1]
    # Gives time in (Earth) days if the sampling rate is in samples per day
    heights = np.round(vapour.coord('Hybrid height').points*1e-03,0)
    lats, lat_points = vapour.coord('latitude'), vapour.coord('latitude').points
    longs, long_points = vapour.coord('longitude'), vapour.coord('longitude').points
    
    if lats.bounds == None:
        vapour.coord('latitude').guess_bounds()
    if longs.bounds == None:
        vapour.coord('longitude').guess_bounds()       
    # Set grid box bounds if there are none
        
    patch = vapour[:,level,y[0]:y[1],x[0]:x[1]].copy()    
    patch_grid = iris.analysis.cartography.area_weights(patch)
    patch_mean = patch.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=patch_grid)
    data = patch_mean.data
    
    lat_index, long_index = int((y[0]+y[1])/2), int((x[0]+x[1])/2)

    plt.plot(np.arange(0,run_length), data, linestyle='-', color='b')
    plt.title('Specific humidity at lat=%s, long=%s, h=%s km' %(lat_points[lat_index], long_points[long_index], heights[level]))
    plt.xlabel('Time [days]') 
    plt.ylabel('Water vapour [kg kg-1]')
    plt.show()
    


def vapour_anomaly(cubes, lat=45, long=0, time_slice=-1):
    
    for cube in cubes:
        if cube.standard_name == 'specific_humidity':
            vapour = cube.copy()
            
    mean_vapour = vapour.collapsed('t', iris.analysis.MEAN)
    anomaly = vapour - mean_vapour
    data = vapour.data
    
    run_length = vapour.shape[0]
    
    x_axis = vapour.coord('longitude').points
    y_axis = vapour.coord('Hybrid height').points
    
    plt.contourf(np.arange(0,run_length)*0.25, y_axis, data[:,:,lat,long].T, bg.N, cmap=bg, norm=TwoSlopeNorm(0))
    plt.title('Water vapour abundance at equator, long=%s' %x_axis[long])
    plt.xlabel('Time [days]')
    plt.ylabel('Height [m]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    