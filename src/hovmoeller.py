#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 17:52:54 2021

@author: Mo Cohen

Hovmoeller plotting function

"""
import iris
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import iris.plot as iplt
import iris.coords
from iris.coords import DimCoord
import numpy as np

# Import packages


brewer_red = mpl_cm.get_cmap('brewer_YlOrRd_09')
brewer_bg = mpl_cm.get_cmap('brewer_PuBu_09')
brewer_reds = mpl_cm.get_cmap('brewer_Reds_09')
brewer_redblu = mpl_cm.get_cmap('brewer_RdBu_11')


def plot_hovmoellerx(cubes):
    
    """ Extracts horizontal wind cube from CubeList and makes Hovmoeller plot
        Input: Iris CubeList
        Outputs: Hovmoeller plots of dayside and nightside zonal equatorial winds"""
    
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube.copy()
            
            
    run_length = x_wind.shape[0]
    
    lats = x_wind.coord('latitude')
    longs = x_wind.coord('longitude')

    if lats.bounds == None:
        x_wind.coord('latitude').guess_bounds()
    if longs.bounds == None:
        x_wind.coord('longitude').guess_bounds()
        
    dayside = x_wind.extract(iris.Constraint(longitude=lambda v: 270 < v <= 359 or 0 <= v <= 90, latitude=lambda v: -10 <= v <= 10))
    nightside = x_wind.extract(iris.Constraint(longitude=lambda v: 90 < v <= 270, latitude=lambda v: -10 <= v <= 10))
    # Extract latitudes -10 to +10 on dayside and nightside
    
    day_grid = iris.analysis.cartography.area_weights(dayside)
    night_grid = iris.analysis.cartography.area_weights(nightside)
    # Get areas of extracted gridboxes
    
    dayside_mean = dayside.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=day_grid)
    nightside_mean = nightside.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=night_grid)
    # Calculate dayside and nightside area-weighted averages
    
  #  months = DimCoord((np.arange(0,run_length)), standard_name='time', units='months')
     # Iris automatically plots this in days since 1 January 1970, with "now" being May 2004
     # Not got around to fixing the time coordinates yet
     
    iplt.contourf(dayside_mean, brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Dayside Mean Zonal Equatorial Wind [m s-1]')
    plt.xlabel('Time')
    plt.ylabel('Height [m]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    iplt.contourf(nightside_mean, brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Nightside Mean Zonal Equatorial Wind [m s-1]')
    plt.xlabel('Time')
    plt.ylabel('Height [m]')
    plt.colorbar(pad=0.1)
    plt.show()