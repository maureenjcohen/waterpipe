#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 13:15:21 2021

@author: Mo Cohen

Plots u'w' 

"""

import iris
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
from matplotlib.colors import TwoSlopeNorm

brewer_red = mpl_cm.get_cmap('brewer_YlOrRd_09')
brewer_bg = mpl_cm.get_cmap('brewer_PuBu_09')
brewer_reds = mpl_cm.get_cmap('brewer_Reds_09')
brewer_redblu = mpl_cm.get_cmap('brewer_RdBu_11')


def plot_uwprime(cubes, time_slice=1780, lat=45, level=47, plots=True):
    
    for cube in cubes:
        if cube.standard_name == 'eastward_wind' or cube.standard_name=='x_wind':
            x_wind = cube[time_slice,:,:,:].copy()
        if cube.standard_name == 'upward_air_velocity':
            z_wind = cube[time_slice,:,:,:].copy()
            
    longitudes = x_wind.shape[2]/2
    
    z_wind = z_wind.regrid(x_wind, iris.analysis.Linear())
    
    for coord in x_wind.coords():
        if coord.long_name == 'Hybrid height':
            heights = np.round(x_wind.coord('Hybrid height').points*1e-03,0)
            vertical = [('Hybrid height', x_wind.coord('Hybrid height').points)]
        elif coord.long_name == 'level_height':
            heights = np.round(x_wind.coord('level_height').points*1e-03,0)
            vertical = [('level_height', x_wind.coord('level_height').points)]
            
    z_wind = z_wind.interpolate(vertical, iris.analysis.Linear()) 
    z_bar = z_wind.collapsed('longitude', iris.analysis.MEAN)
    z_prime = z_wind - z_bar           

    u_fft = sp.fftpack.fft(x_wind[level,lat,:].data)
    u_psd = np.abs(u_fft)**2
    u_freq = sp.fftpack.fftfreq(len(u_psd), 1./144)
    pos = u_freq > 0
    
    highpass = u_fft.copy()
    highpass[np.abs(u_freq) < 1.1] = 0
    highpass_psd = np.abs(highpass)**2
    
    if plots == True:
        
        plt.figure(figsize=(10,5))
        plt.plot(np.arange(-longitudes, longitudes), np.roll(x_wind[level,lat,:].data, 72, axis=0))
        plt.title('U at Equator, t=%s days, h=%s km' %(time_slice/4, heights[level]))
        plt.xlabel('Longitude [degrees]')
        plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
        plt.ylabel('Wind speed [m s-1]')
        plt.show()
        
        plt.figure(figsize=(10,5))
        plt.plot(np.arange(-longitudes, longitudes), np.roll(z_prime[level,lat,:].data, 72, axis=0))
        plt.title('$W^{\prime}$ at Equator, t=%s days, h=%s km' %(time_slice/4, heights[level]))
        plt.xlabel('Longitude [degrees]')
        plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
        plt.ylabel('Wind speed [m s-1]')
        plt.show()
    
        fig, ax = plt.subplots(1,1,figsize=(8,4))
        ax.plot(u_freq[pos], highpass_psd[pos])
        ax.set_xlabel('Wavenumber')
        ax.set_ylabel('PSD')
        ax.set_title('U wavenumbers with k=1 removed')
        plt.show()
    
    u_cleaned = np.real(sp.fftpack.ifft(highpass))
    u_bar = np.mean(u_cleaned)
    u_prime = u_cleaned - u_bar
    cross_term = u_prime*z_prime[level,lat,:].data
    
    if plots == True: 
        
        plt.figure(figsize=(10,5))    
        plt.plot(np.arange(-longitudes, longitudes), np.roll(u_prime, 72))
        plt.title('$U^{\prime}$ at Equator, t=%s days, h=%s km' %(time_slice/4, heights[level]))
        plt.xlabel('Longitude [degrees]')
        plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
        plt.ylabel('Wind speed [m s-1]')
        plt.show()
    
        plt.figure(figsize=(10,5))    
        plt.plot(np.arange(-longitudes, longitudes), np.roll(cross_term, 72))
        plt.title('$U^{\prime} \cdot W^{\prime}$ at Equator, t=%s days, h=%s km' %(time_slice/4, heights[level]))
        plt.xlabel('Longitude [degrees]')
        plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
        plt.ylabel('Vertical Eddy Momentum [m2 s-2]')
        plt.show()
        
    return u_prime, cross_term


def create_xsection(cubes, time_slice=1780, lat=45):
    
    for cube in cubes:
        if cube.standard_name == 'eastward_wind' or cube.standard_name=='x_wind':
            x_wind = cube[time_slice,:,:,:].copy()
            
    longitudes, heights = x_wind.shape[2]/2, np.round(x_wind.coord('Hybrid height').points*1e-03,0)    
    
    u_prime_array = []
    xterm_array = []
    
    for level in range(0,60):
    
        u_prime, cross_term = plot_uwprime(cubes=cubes, time_slice=time_slice, lat=lat, level=level, plots=False)
        
        u_prime_array.append(u_prime)
        xterm_array.append(cross_term)
    
    u_prime_array = np.array(u_prime_array)
    xterm_array = np.array(xterm_array)
    
    plt.figure(figsize=(10,5))    
    plt.contourf(np.arange(-longitudes, longitudes), np.array(heights), np.roll(u_prime_array, 72, axis=1), brewer_redblu.N, cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    plt.title('$U^{\prime}$ at Equator, t=%s days' %(time_slice/4))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.ylabel('Height [km]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    plt.figure(figsize=(10,5))    
    plt.contourf(np.arange(-longitudes, longitudes), np.array(heights), np.roll(xterm_array, 72, axis=1), brewer_redblu.N, cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    plt.title('$U^{\prime} \cdot W^{\prime}$ at Equator, t=%s days' %(time_slice/4))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.ylabel('Height [km]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    plt.figure(figsize=(10,5))    
    plt.contourf(np.arange(-longitudes, longitudes), np.array(heights[33:42]), np.roll(u_prime_array[33:42], 72, axis=1), np.linspace(-5,5,80), cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    plt.title('$U^{\prime}$ at Equator, t=%s days' %(time_slice/4))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.ylabel('Height [km]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    plt.figure(figsize=(10,5))    
    plt.contourf(np.arange(-longitudes, longitudes), np.array(heights[35:58]), np.roll(xterm_array[35:58], 72, axis=1), brewer_redblu.N, cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    plt.title('$U^{\prime} \cdot W^{\prime}$ at Equator, t=%s days' %(time_slice/4))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.ylabel('Height [km]')
    plt.colorbar(pad=0.1)
    plt.show()
    