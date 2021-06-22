#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 14:37:36 2021

@author: Maureen Cohen
"""
import iris
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
from matplotlib.colors import TwoSlopeNorm


brewer_redblu = mpl_cm.get_cmap('RdBu_r')

def wave_acceleration(cubes, lat=45, long=0, start=0, end=240, plot=False):
    
    for cube in cubes:
        if cube.standard_name == 'eastward_wind':
            x_wind = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'upward_air_velocity':
            z_wind = cube[start:end,:,:,:].copy()
            
    z_wind = z_wind.regrid(x_wind, iris.analysis.Linear())
    
    vertical = [('Hybrid height', x_wind.coord('Hybrid height').points)]
    z_wind = z_wind.interpolate(vertical, iris.analysis.Linear())

    longitudes = x_wind.shape[3]/2
    heights = x_wind.coord('Hybrid height').points
    
    """ Filter zonal wind, remove longer wavelengths, calculate u-prime"""
    
    x_data = x_wind.data
    
    u_prime_list = []
    
    for time in range(0,x_wind.shape[0]):
        time_list = []
        
        for level in range(0,x_wind.shape[1]):
            level_list = []
            
            for latx in range(0,x_wind.shape[2]):
            
                u_fft = sp.fftpack.fft(x_data[time,level,latx,:])
                u_psd = np.abs(u_fft)**2
                u_freq = sp.fftpack.fftfreq(len(u_psd), 1./144)
                
                highpass = u_fft.copy()
                highpass[np.abs(u_freq) < 5.1] = 0
                
                u_cleaned = np.real(sp.fftpack.ifft(highpass))
                u_bar = np.mean(u_cleaned)
                u_prime = u_cleaned - u_bar
                level_list.append(u_prime)
                
            time_list.append(level_list)
        u_prime_list.append(time_list)
        
    u_prime = np.array(u_prime_list)
    print(u_prime.shape)

    if plot == True: 
        
        plt.figure(figsize=(10,5))    
        plt.contourf(np.arange(-longitudes, longitudes), np.array(heights), np.roll(u_prime[-1,:,lat,:], 72, axis=1), np.linspace(-20,20,40), cmap=brewer_redblu, norm=TwoSlopeNorm(0))
        plt.title('$U^{\prime}$ at Equator, t=%s days' %(end/4))
        plt.xlabel('Longitude [degrees]')
        plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
        plt.ylabel('Height [m]')
        plt.colorbar(pad=0.1)
        plt.show()    

    """ Calculate w-prime"""

    z_mean = z_wind.collapsed('t', iris.analysis.MEAN)    
    z_anomaly = z_wind - z_mean    
    
    if plot == True:
        
        plt.figure(figsize=(10,5))
        plt.contourf(np.arange(-longitudes,longitudes), np.array(heights), np.roll(z_anomaly[-1,:,lat,:].data, 72, axis=1), np.linspace(-0.09,0.09,20), cmap=brewer_redblu, norm=TwoSlopeNorm(0))
        plt.title('Vertical Wind Anomaly at Equator [m s-1], t = %s days' %(end/4))
        plt.xlabel('Longitude [degrees]')
        plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
        plt.ylabel('Height [m]')
        plt.colorbar(pad=0.1)
        plt.show()
    
    """ Differentiate wrt height"""

    array = u_prime*z_anomaly.data

    h = np.array([heights])
    h = np.repeat(h[:,np.newaxis], (end-start), axis=0)
    h = np.reshape(h,((end-start),60))
    h = np.repeat(h[:,:,np.newaxis], 90, axis=2)
    h = np.repeat(h[:,:,:,np.newaxis], 144, axis=3)
    
    print(array.shape, h.shape)
    
    acceleration = array.copy()
    acceleration[:,0,...] = (array[:,1,...]-array[:,0,...])/(h[:,1,...]-h[:,0,...])
    acceleration[:,-1,...] = (array[:,-1,...]-array[:,-2,...])/(h[:,-1,...]-h[:,-2,...])
    acceleration[:,1:-1,...] = (array[:,2:,...]-array[:,0:-2,...])/(h[:,2:,...]-h[:,0:-2,...])
    
    """ Plot """
    net_acc = np.mean(acceleration, axis=3) # + x_wind[start,:,31:59,long].data
    net_acc = np.sum(net_acc, axis=0)
    
    x_mean = x_wind.collapsed('longitude', iris.analysis.MEAN)
    mean_acc = np.mean(acceleration, axis=0)    

    plt.figure(figsize=(10,5))
    plt.contourf(np.arange(-longitudes,longitudes), np.array(heights), np.roll(-mean_acc[:,lat,:], 72, axis=1), np.linspace(-5e-06,5e-06,20), cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    plt.title('Wave-Induced Acceleration at Equator [m s-2], t = %s days' %(end/4))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.ylabel('Height [m]')
    plt.colorbar(pad=0.1)
    plt.show()
                
    x_axis = x_wind.coord('latitude').points
    y_axis = x_wind.coord('Hybrid height').points
    plt.figure(figsize=(8,10))
    plt.contourf(x_axis, y_axis, -net_acc, brewer_redblu.N, cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    plt.colorbar(pad=0.1)
    CS = plt.contour(x_axis, y_axis, (x_mean[-1,:,:].data - x_mean[0,:,:].data), colors='black', linewidths=1.5)
    plt.title('Mean wave-induced acceleration from day=%s to %s' %(start/4, end/4))
    plt.xlabel('Latitude [degrees]')
    plt.ylabel('Height [m]')
    plt.clabel(CS, inline=False, colors='k', fmt='%1.1f')
    plt.show()
    
    