#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 11:07:06 2021

@author: s1144983
"""


import iris
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


def vertical_wind(cubes, level=47, long=0, lat=45):
    
    for cube in cubes:
        if cube.standard_name == 'upward_air_velocity':
            z_wind = cube[:,:,lat,long].copy()    
    
    heights = np.round(z_wind.coord('Hybrid height').points*1e-03,0)
    run_length = z_wind.shape[0]
    data = z_wind.data

    plt.plot(np.arange(0,run_length), data[:,level])
    plt.title('Vertical wind at substellar point, h= %s km' %heights[level])
    plt.xlabel('Time [6-hours]') 
    plt.ylabel('Velocity [m s-1]')
    plt.show()
    
def qbo_cell(cubes, long=0, level=47, time1=0, time2=720):
    
    for cube in cubes:
        if cube.standard_name == 'air_potential_temperature':
            theta = cube[:,:,:,long].copy()
        if cube.standard_name == 'air_pressure':
            pressure = cube[:,:,:,long].copy()
    
    p0 = iris.coords.AuxCoord(100000.0, long_name='reference_pressure', units='Pa')
    p0.convert_units(pressure.units)
    temp = theta[:,level,:]*((pressure[:,level,:]/p0)**(287.05/1005)) # R and cp in J/kgK for 300K
    temperature = temp.data
    
    x_axis = theta.coord('latitude').points
    heights = np.round(theta.coord('Hybrid height').points*1e-03,0)
    plt.plot(x_axis, temperature[time1,:], color='b', label='%s days' %(time1/4))
    plt.plot(x_axis, temperature[time2,:], color='r', label='%s days' %(time2/4))
    plt.title('Temperature structure at %s km' %heights[level])
    plt.xlabel('Latitude [degrees]')
    plt.ylabel('Temperature [K]')
    plt.legend()
    plt.show()

    
    