#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 15:21:07 2021

@author: s1144983
"""
import iris
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import iris.plot as iplt
from matplotlib.colors import TwoSlopeNorm

# Import packages


brewer_redblu = mpl_cm.get_cmap('brewer_RdBu_11')

def meridional_wind(cubes, n=3, time_slice=-1, level=38):
    
    for cube in cubes:
        if cube.standard_name == 'y_wind':
            y_wind = cube.copy()
        if cube.standard_name == 'air_pressure':
            pressure = cube.copy()            
    
    height = [('level_height', y_wind.coord('level_height').points)]
    pressure = pressure.interpolate(height, iris.analysis.Linear())
    # Regrid so that all three cubes are on the same x, y, z grid. Uses the x_wind as reference for the others
    
    p_heights = np.round(pressure.data*1e-05,4)
    km_heights = np.round(pressure.coord('level_height').points*1e-03,1)
    # For labelling the quiver plot

    mean_y = y_wind[time_slice,:,:,:].collapsed('longitude', iris.analysis.MEAN)
    y_prime = y_wind - mean_y
    # Zonal mean meridional winds
    
    y_wind = y_wind[time_slice,level,:,:].data    
    x_wind = np.zeros_like(y_wind) 

    X,Y = np.meshgrid(np.arange(-72,72), np.arange(-45,46))
    fig1, ax1 = plt.subplots(figsize = (10,5)) 
    q1 = ax1.quiver(X[::n,::n], Y[::n,::n], np.roll(x_wind[::n,::n], 72, axis=1), np.roll(y_wind[::n,::n], 72, axis=1), angles='xy', scale_units='xy', scale=3)
    # Create a quiver plot. The np.roll function moves the cube data so the substellar point is centred.
    ax1.quiverkey(q1, X=0.9, Y=1.05, U=3, label='3 m/s', labelpos='E', coordinates='axes')
    # This creates the key with the arrow size. The value of U, the label string, and the value of scale in the previous line should all match.
    plt.title('Meridional wind vectors [m s-1], h=%s bar, %s km' %(p_heights[time_slice,level,0,0], km_heights[level]))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.yticks((-45,-30,-15,0,15,30,45),('90S','60S','30S','0','30N','60N','90N')) 
    # Labelled the longs and lats manually 
    plt.show()
            
    CS = iplt.contourf(mean_y, brewer_redblu.N, cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    plt.title('Zonal Mean Meridional Wind [m s-1]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.clabel(CS, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    plt.show()