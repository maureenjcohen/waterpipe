#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 16:03:25 2021

@author: Maureen Cohen
"""
import iris, cartopy, pywt
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
from matplotlib.colors import TwoSlopeNorm
from iris.coord_systems import GeogCS
import numpy as np
import scipy as sp
# Import packages


magma = mpl_cm.get_cmap('magma')


def wavelets(cubes, radius=7160000, scales=512, wavelet='mexh', x=(106,110), y=(43,47), level=47, sampling=0.25):
    """ Performs wavelet transform for 1-D data
    Inputs: Iris CubeList, radius of planet, max scale for wavelet, mother wavelet shape,
    range of longitude columns to average over, range of latitude columns to average over,
    atmospheric level, number of samples per day
    
    Outputs: Plot of data, plot of scaleogram"""
    
    for cube in cubes:
        if cube.standard_name == 'eastward_wind':
            x_wind = cube.copy()    
    # Extract zonal wind data
    
    x_wind.coord('latitude').coord_system = GeogCS(radius)
    x_wind.coord('longitude').coord_system = GeogCS(radius)
    # Sets planet radius in m for area-weighted average. Default is radius of Proxima Centauri b
    
    run_length, height = x_wind.shape[0], x_wind.shape[1]
    time = np.arange(0,run_length)*sampling
    # Gives time in (Earth) days if the sampling rate is in samples per day
    heights = np.round(x_wind.coord('Hybrid height').points*1e-03,0)
    lats, lat_points = x_wind.coord('latitude'), x_wind.coord('latitude').points
    longs, long_points = x_wind.coord('longitude'), x_wind.coord('longitude').points
    # Extract info for labelling plots

    if lats.bounds == None:
        x_wind.coord('latitude').guess_bounds()
    if longs.bounds == None:
        x_wind.coord('longitude').guess_bounds()
        
    # Set grid box bounds if there are none
        
    patch = x_wind[:,level,y[0]:y[1],x[0]:x[1]].copy()    
    patch_grid = iris.analysis.cartography.area_weights(patch)
    patch_mean = patch.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=patch_grid)
    # Extract data from a small patch of the atmosphere and find the mean

    data = patch_mean.data
    coeffs, freqs = pywt.cwt(data, np.arange(1,scales), wavelet, 1./run_length)

    power = (abs(coeffs))**2
    # Find power from coefficients
    period = (1./freqs)*run_length*sampling
    # Gives period in days
    
    lat_index, long_index = int((y[0]+y[1])/2), int((x[0]+x[1])/2)
    # Indices for central lat and long coordinates, for the plot titles
    
    plt.plot(np.arange(0,run_length)*sampling, patch_mean.data, linestyle='-', color='b')
    plt.title('Zonal Wind at lat=%s, long=%s, h=%s km' %(lat_points[lat_index], long_points[long_index], heights[level]))
    plt.xlabel('Time [days]') 
    plt.ylabel('Velocity [m s-1]')
    plt.show()
    
    plt.pcolormesh(time, period, power, cmap=magma)
    plt.title('Scaleogram for lat=%s, long=%s, h=%s km' %(lat_points[lat_index], long_points[long_index], heights[level]))
    plt.xlabel('Time [days]')
    plt.ylabel('Period [days]')
    plt.show()
    


def wavelets2d(cubes, scales=512, wavelet='mexh', start=0, end=480, x=108, y=45, level=47, sampling=0.25):
    
    for cube in cubes:
        if cube.standard_name == 'eastward_wind':
            x_wind = cube[start:end,:,:,:].copy()    
    # Extract zonal wind data
    
    run_length, height, longitudes = x_wind.shape[0], x_wind.shape[1], x_wind.shape[3]
    time = np.arange(0,run_length)*sampling
    # Gives time in (Earth) days if the sampling rate is in samples per day
    heights = np.round(x_wind.coord('Hybrid height').points*1e-03,0)
    lat_points = x_wind.coord('latitude').points
    long_points = x_wind.coord('longitude').points
    # Extract info for labelling plots
    
    x_data = x_wind.data
    
    u_prime_list = []
    
    for timex in range(0,x_wind.shape[0]):
        time_list = []
        
        for levelx in range(0,x_wind.shape[1]):
            level_list = []
            
            for latx in range(0,x_wind.shape[2]):
            
                u_fft = sp.fftpack.fft(x_data[timex,levelx,latx,:])
                u_psd = np.abs(u_fft)**2
                u_freq = sp.fftpack.fftfreq(len(u_psd), 1./144)
                
                highpass = u_fft.copy()
                highpass[np.abs(u_freq) < 5.1] = 0
                
                u_cleaned = np.real(sp.fftpack.ifft(highpass))
                u_bar = np.mean(u_cleaned)
                u_prime = u_cleaned - u_bar
                level_list.append(u_cleaned)
                
            time_list.append(level_list)
        u_prime_list.append(time_list)
        
    u_prime = np.array(u_prime_list)
            
    intime = u_prime[:,level,y,x].data
    inspace = u_prime[end-1,level,y,:].data
    timecoeffs, freqs = pywt.cwt(intime, np.arange(1,scales), wavelet, 1./run_length)
    spacecoeffs, wavenums = pywt.cwt(inspace, np.arange(1,scales), wavelet, 1./longitudes)
    
    timepower = (abs(timecoeffs))**2
    timescale = (1./freqs)*run_length*sampling
    spacepower = (abs(spacecoeffs))**2
    spacescale = (1./wavenums)*longitudes
    
    
    plt.plot(np.arange(0,run_length)*sampling, intime, linestyle='-', color='b')
    plt.title('Zonal Wind at lat=%s, long=%s, h=%s km' %(lat_points[y], long_points[x], heights[level]))
    plt.xlabel('Time [days]') 
    plt.ylabel('Velocity [m s-1]')
    plt.show()
    
    plt.pcolormesh(time, freqs, timepower, cmap=magma)
    plt.title('Scaleogram for lat=%s, long=%s, h=%s km' %(lat_points[y], long_points[x], heights[level]))
    plt.xlabel('Time [days]')
    plt.ylabel('Frequency')
    plt.show()
    
    plt.plot(np.arange(-longitudes/2,longitudes/2), inspace, linestyle='-', color='b')
    plt.title('Zonal Wind at lat=%s, time=%s days, h=%s km' %(lat_points[y], time[end-1], heights[level]))
    plt.xlabel('Longitude') 
    plt.ylabel('Velocity [m s-1]')
    plt.show()
    
    plt.pcolormesh(np.arange(-longitudes/2,longitudes/2), wavenums, spacepower, cmap=magma)
    plt.title('Scaleogram for lat=%s, time=%s days, h=%s km' %(lat_points[y], time[end-1], heights[level]))
    plt.xlabel('Longitude')
    plt.ylabel('Wavenumber')
    plt.show()
