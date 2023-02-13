#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 15:54:49 2023

@author: Mo Cohen
"""

import iris, windspharm
from windspharm.iris import VectorWind
import iris.plot as iplt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
from matplotlib.colors import TwoSlopeNorm
from iris.analysis import calculus
from numpy import unravel_index


redblu = mpl_cm.get_cmap('RdBu_r')

def plot_vorts(cubes, time_slice=-1, level=8, omega=1.19e-05, g=9.12, lat=50):
    
    """ Calculate and plot the potential vorticity
    Inputs:
        - Iris DataCube
        - Time index - default is the last time coordinate
        - Level - Default is 2.96 km
        - Rotation rate - Default is TRAPPIST-1e
        - Gravitational constant - Default is TRAPPIST-1e
        - Latitude - Default is 50 degrees north
        
    Outputs:
        - Contour fill plot of potential vorticity for the chosen
        level and time
        - Iris cube of potential vorticity

    Optional:
    Can also plot the relative and absolute vorticity by uncommenting
    sections of code.                                                  """
    

    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube[time_slice,:,:,:].copy()
        if cube.standard_name == 'y_wind':
            y_wind = cube[time_slice,:,:,:].copy()
        if cube.standard_name =='air_potential_temperature':
            theta = cube[time_slice,:,:,:].copy()
        if cube.standard_name == 'air_pressure':
            pressure = cube[time_slice,:,:,:].copy()
            
    heights = np.round(x_wind.coord('level_height').points*1e-03,2) # in km
    longitudes = x_wind.shape[2]
    latitudes = x_wind.shape[1]
    # Extract some data for labelling/creating plots

    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())    
    theta = theta.regrid(x_wind, iris.analysis.Linear())
    pressure = pressure.regrid(x_wind, iris.analysis.Linear())
    
    vertical = [('level_height', x_wind.coord('level_height').points)]

    theta = theta.interpolate(vertical, iris.analysis.Linear())
    pressure = pressure.interpolate(vertical, iris.analysis.Linear())
    # Regrid/interpolate so all data cubes are on the same 3-D grid
    
    dtheta = iris.analysis.calculus.differentiate(theta, 'level_height')
    dpressure = iris.analysis.calculus.differentiate(pressure, 'level_height')
    # Calculate theta and pressure deltas
    dth_dp = dtheta/dpressure
    # Change of potential temperature with change in pressure (approx. derivative)
    
    lat_rads = lat*(np.pi/180) # Convert input latitude to radians
    f_constant = 2*omega*np.sin(lat_rads) # Calculate Coriolis parameter

    winds = VectorWind(x_wind, y_wind) # Windspharm VectorWind object
    
    # winds2 = windspharm.iris.VectorWind(x_wind[level,:,:],y_wind[level,:,:])
    # uchi,vchi,upsi,vpsi = winds2.helmholtz(truncation=21)
   
    # zonal_upsi = upsi.collapsed('longitude', iris.analysis.MEAN)
    # zonal_vpsi = vpsi.collapsed('longitude', iris.analysis.MEAN)
    # eddy_upsi = upsi - zonal_upsi
    # eddy_vpsi = vpsi - zonal_vpsi
    # eddies = windspharm.iris.VectorWind(eddy_upsi, eddy_vpsi)
    # rel_vort = eddies.vorticity()

    rel_vort = winds.vorticity()
#    abs_vort = winds.absolutevorticity()
    # Calculate relative and absolute vorticity
    
    f = np.full_like(rel_vort, f_constant) 
    # Make an array of f to add to relative vorticity
    pot_vort = -g*(rel_vort[0:-1,:,:].data + f)*dth_dp.data
    # Calculate potential vorticity
    
    # Can quickly look at relative and absolute vorticity below if you want
    rcolor_levs = np.arange(-10, 10, 1)
    r_constant = 1e05
    rexponent = -5
    plt.figure(figsize=(8,5))
    plt.contourf(np.arange(-longitudes/2, longitudes/2), np.arange(-latitudes/2, latitudes/2),
                 np.roll(rel_vort[level,:,:,].data*r_constant, 72, axis=1),levels=rcolor_levs, cmap=redblu, norm=TwoSlopeNorm(0))
    plt.title('Relative Vorticity, h = %s km' %(heights[level]))
    plt.xlabel('Longitude [degrees]')
    plt.ylabel('Latitude [degrees]')
    plt.xticks((-72, -60, -48, -36, -24, -12, 0, 12, 24, 36, 48, 60, 72), ('180W', '150W',
                '120W', '90W', '60W', '30W', '0', '30E', '60E', '90E', '120E', '150E', '180E'))
    plt.yticks((-45, -30, -15, 0, 15, 30, 45),
                ('90S', '60S', '30S', '0', '30N', '60N', '90N'))
    cbar = plt.colorbar(pad=0.2, orientation='horizontal')
    cbar.set_ticks(rcolor_levs)
    cbar.ax.set_title('$10^{%s}$ s$^{-1}$' %rexponent)
    plt.show()
    
    
    rv_min = unravel_index(np.argmax(rel_vort[level,:20,0:72].data, axis=None), rel_vort[level,:20,0:72].shape)
    print(rv_min[0], rv_min[1])
    print(x_wind.coord('latitude').points[rv_min[0]], x_wind.coord('longitude').points[rv_min[1]])
    print(rel_vort[level,rv_min[0],rv_min[1]].data)
    
    # iplt.contourf(abs_vort[level,:,:], redblu.N, cmap=redblu, norm=TwoSlopeNorm(0))
    # plt.title('Absolute Vorticity, h = %s km' %(heights[level]), y=1.20)
    # plt.ylabel('Latitude [degrees]')
    # plt.xlabel('Longitude [degrees]')
    # ax = plt.gca()
    # ax.gridlines(draw_labels=True)
    # plt.colorbar(orientation='horizontal')
    # plt.show()
    
    # The values below need to be adjusted to find the right colorbar settings
    # color_levs = np.arange(-10,10,1) # Levels for the colorbar
    # multiply_constant = 1e07 # Vorticity values are typically 10^-5, 10^-6, 10^-7, etc.
    # # Multiply your data cube by this constant to display as 10^0
    # exponent = -7 # Use the exponent for the colorbar title
    
    # plt.figure(figsize=(8,5))
    # plt.contourf(np.arange(-longitudes/2, longitudes/2), np.arange(-latitudes/2, latitudes/2),
    #              np.roll(pot_vort[level,:,:].data*multiply_constant, 72, axis=1),levels=color_levs, cmap=redblu, norm=TwoSlopeNorm(0))
    # plt.title('Potential Vorticity, h = %s km' %(heights[level]))
    # plt.xlabel('Longitude [degrees]')
    # plt.ylabel('Latitude [degrees]')
    # plt.xticks((-72, -60, -48, -36, -24, -12, 0, 12, 24, 36, 48, 60, 72), ('180W', '150W',
    #             '120W', '90W', '60W', '30W', '0', '30E', '60E', '90E', '120E', '150E', '180E'))
    # plt.yticks((-45, -30, -15, 0, 15, 30, 45),
    #             ('90S', '60S', '30S', '0', '30N', '60N', '90N'))
    # cbar = plt.colorbar(pad=0.2, orientation='horizontal')
    # cbar.set_ticks(color_levs)
    # cbar.ax.set_title('$10^{%s}$ s$^{-1}$' %exponent)
    # plt.show()
    
   
    

    
    
    
    