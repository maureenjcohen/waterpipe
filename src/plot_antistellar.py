#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  1 14:01:01 2021

@author: Maureen Cohen
"""

import os, iris
import matplotlib.pyplot as plt
import iris.plot as iplt
import numpy as np


def plot_antiwind(cubes, lat=45, long=72, level=47):
    
    for cube in cubes:
        if cube.standard_name == 'eastward_wind':
            x_wind = cube[:,:,lat,long].copy()
            
    data = x_wind.data
    time = x_wind.shape[0]
    heights = np.round(x_wind.coord('Hybrid height').points*1e-03,0)
    
    plt.plot(np.arange(0,time)*0.25, data[:,level])
    plt.title('Zonal wind at h=%s km' %heights[level])
    plt.xlabel('Time [days]')
    plt.ylabel('Wind velocity [m s-1]')
    plt.show()
    


def plot_antitemp(cubes, lat=45, long=72, level=47):
        
    for cube in cubes:
        if cube.standard_name == 'air_potential_temperature':
            theta = cube[:,:,lat,long].copy()
        if cube.standard_name == 'air_pressure':
            pressure = cube[:,:,lat,long].copy()
        
    p0 = iris.coords.AuxCoord(100000.0, long_name='reference_pressure', units='Pa')
    p0.convert_units(pressure.units)
    temperature = theta*((pressure/p0)**(287.05/1005)) # R and cp in J/kgK for 300K        
   
    data = temperature.data
    time = temperature.shape[0]
    heights = np.round(temperature.coord('level_height').points*1e-03,0)
    
    avg = np.mean(data[:,level])
    smallest = np.min(data[:,level])
    biggest = np.max(data[:,level])
    print('The average temperature is %s K' %avg)
    print('The min and max are %s and %s' %(smallest, biggest))

    
    plt.plot(np.arange(0,time)*0.25, data[:,level])
    plt.title('Temperature at h=%s km' %heights[level])
    plt.xlabel('Time [days]')
    plt.ylabel('Temperature [K]')
    plt.show()