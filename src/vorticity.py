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


redblu = mpl_cm.get_cmap('brewer_RdBu_11')
plasma = mpl_cm.get_cmap('plasma')


def plot_vorts(cubes, time_slice=-1, level=8, omega=1.19e-05, g=9.12, lat=50):
    
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube[time_slice,:,:].copy()
        if cube.standard_name == 'y_wind':
            y_wind = cube[time_slice,:,:,:].copy()
        if cube.standard_name =='air_potential_temperature':
            theta = cube[time_slice,:,:,:].copy()
        if cube.standard_name == 'air_pressure':
            pressure = cube[time_slice,:,:,:].copy()
            
    heights = np.round(x_wind.coord('level_height').points*1e-03,2) # in km
    longitudes = x_wind.shape[2]
    latitudes = x_wind.shape[1]

    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())    
    theta = theta.regrid(x_wind, iris.analysis.Linear())
    pressure = pressure.regrid(x_wind, iris.analysis.Linear())
    
    vertical = [('level_height', x_wind.coord('level_height').points)]

    theta = theta.interpolate(vertical, iris.analysis.Linear())
    pressure = pressure.interpolate(vertical, iris.analysis.Linear())
    
    dtheta = iris.analysis.calculus.differentiate(theta, 'level_height')
    dpressure = iris.analysis.calculus.differentiate(pressure, 'level_height')
    dth_dp = dtheta/dpressure
    
    lat_rads = lat*(np.pi/180)
    f_constant = 2*omega*np.sin(lat_rads)

    winds = VectorWind(x_wind, y_wind)
    
    rel_vort = winds.vorticity()
    abs_vort = winds.absolutevorticity()
    
    f = np.full_like(rel_vort, f_constant)   
    pot_vort = -g*(rel_vort[0:-1,:,:].data + f)*dth_dp.data
    
    # iplt.contourf(rel_vort[level,:,:], redblu.N, cmap=redblu, norm=TwoSlopeNorm(0))
    # plt.title('Relative Vorticity, h = %s km' %(heights[level]), y=1.20)
    # plt.ylabel('Latitude [degrees]')
    # plt.xlabel('Longitude [degrees]')
    # ax = plt.gca()
    # ax.gridlines(draw_labels=True)
    # plt.colorbar(orientation='horizontal')
    # plt.show()
    
    # iplt.contourf(abs_vort[level,:,:], redblu.N, cmap=redblu, norm=TwoSlopeNorm(0))
    # plt.title('Absolute Vorticity, h = %s km' %(heights[level]), y=1.20)
    # plt.ylabel('Latitude [degrees]')
    # plt.xlabel('Longitude [degrees]')
    # ax = plt.gca()
    # ax.gridlines(draw_labels=True)
    # plt.colorbar(orientation='horizontal')
    # plt.show()
    
    color_levs = np.arange(-10,10,1)
    plt.figure(figsize=(8,5))
    plt.contourf(np.arange(-longitudes/2, longitudes/2), np.arange(-latitudes/2, latitudes/2),
                 np.roll(rel_vort[level,:,:].data, 72, axis=1), cmap=redblu)
    plt.title('Relative Vorticity, h = %s km' %(heights[level]))
    plt.xlabel('Longitude [degrees]')
    plt.ylabel('Latitude [degrees]')
    plt.xticks((-72, -60, -48, -36, -24, -12, 0, 12, 24, 36, 48, 60, 72), ('180W', '150W',
                '120W', '90W', '60W', '30W', '0', '30E', '60E', '90E', '120E', '150E', '180E'))
    plt.yticks((-45, -30, -15, 0, 15, 30, 45),
                ('90S', '60S', '30S', '0', '30N', '60N', '90N'))
    cbar = plt.colorbar(pad=0.2, orientation='horizontal')
#    cbar.set_ticks(color_levs)
    cbar.ax.set_title('$10^{-5}$ s$^{-1}$')
    # plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/laso/epsfigs/deepconvection.eps', format='eps')
    plt.show()
    
    

    
    
    
    