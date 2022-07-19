#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 16:30:25 2022

@author: Maureen Cohen
"""
import iris
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
from matplotlib.colors import TwoSlopeNorm
from iris.coord_systems import GeogCS
import numpy as np
# Import packages

redblu = mpl_cm.get_cmap('RdBu')

def vert_profile(cubes,start=500,end=600,level=8):

    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube[start:end,:,:,:].copy()

    heights = np.round(x_wind.coord('level_height').points*1e-03,0)
    time_axis = np.arange(0,x_wind.shape[0])
    lats = x_wind.coord('latitude')
    lons = x_wind.coord('longitude')

    if lats.bounds == None:
        x_wind.coord('latitude').guess_bounds()
    if lons.bounds == None:
        x_wind.coord('longitude').guess_bounds()

    grid_areas = iris.analysis.cartography.area_weights(x_wind)
    global_mean = x_wind.collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=grid_areas)

    fig1, ax1 = plt.subplots()
    ax1.set_xlabel('Time [days]')
    ax1.set_ylabel('Wind [m/s]')
    ax1.set_title('Global mean zonal wind at h=%s km' %heights[level])
    plt.plot(global_mean[:,level].data)
    plt.show()

    fig2, ax2 = plt.subplots()
    ax2.set_xlabel('Time [days]')
    ax2.set_ylabel('Height [km]')
    ax2.set_title('Global mean zonal wind')
    plt.contourf(time_axis, heights[:30], global_mean[:,:30].data.T, redblu.N, cmap=redblu, norm=TwoSlopeNorm(0))
    plt.colorbar(pad=0.1)
    plt.show()

