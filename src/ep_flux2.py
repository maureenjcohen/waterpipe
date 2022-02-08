#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 10:45:21 2022

@author: Mo Cohen
"""

import iris
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
from matplotlib.colors import TwoSlopeNorm
import iris.plot as iplt
import iris.coords
from iris.coord_systems import GeogCS
from iris.analysis import calculus
import numpy as np
import scipy as sp
# Import packages


brewer_red = mpl_cm.get_cmap('brewer_YlOrRd_09')
brewer_bg = mpl_cm.get_cmap('brewer_PuBu_09')
brewer_reds = mpl_cm.get_cmap('brewer_Reds_09')
brewer_redblu = mpl_cm.get_cmap('brewer_RdBu_11')


def ep_zderivative(array,h):
    work = array.copy()
    work[:,0,...] = (array[:,1,...]-array[:,0,...])/(h[:,1,...]-h[:,0,...])
    work[:,-1,...] = (array[:,-1,...]-array[:,-2,...])/(h[:,-1,...]-h[:,-2,...])
    work[:,1:-1,...] = (array[:,2:,...]-array[:,1:-1,...])/(h[:,2:,...]-h[:,1:-1,...])
    
    return work

def ep_latderivative(array,l):
    
    work = array.copy()
    work[:,:,0] = (array[:,:,1]-array[:,:,0])/(l[1]-l[0])
    work[:,:,-1] = (array[:,:,-1]-array[:,:,-2])/(l[-1]-l[-2])
    work[:,:,1:-1] = (array[:,:,2:]-array[:,:,1:-1])/(l[2:]-l[1:-1])
    
    return work

def ep_flux2(cubes, omega=0.64617667e-05, long_slice=(0,-1), start=2880, end=3240):
    
    for cube in cubes:
        if cube.standard_name == 'air_potential_temperature':
            theta = cube[start:end,:,:,long_slice[0]:long_slice[1]].copy()
        if cube.standard_name == 'eastward_wind':
            x_wind = cube[start:end,:,:,long_slice[0]:long_slice[1]].copy()
        if cube.standard_name == 'northward_wind':
            y_wind = cube[start:end,:,:,long_slice[0]:long_slice[1]].copy()
        if cube.standard_name == 'upward_air_velocity':
            z_wind = cube[start:end,:,:,long_slice[0]:long_slice[1]].copy()
        if cube.standard_name == 'air_pressure':
            pressure = cube[start:end,:,:,long_slice[0]:long_slice[1]].copy()
    
    heights, longitudes, longs = np.round(pressure.coord('Hybrid height').points*1e-03,0), pressure.coord('longitude').points, pressure.shape[3]/2
    
    y_wind, x_wind = y_wind.regrid(pressure, iris.analysis.Linear()), x_wind.regrid(pressure, iris.analysis.Linear())
    theta = theta.regrid(pressure, iris.analysis.Linear())
    vertical = [('Hybrid height', pressure.coord('Hybrid height').points)]
    x_wind, y_wind = x_wind.interpolate(vertical, iris.analysis.Linear()), y_wind.interpolate(vertical, iris.analysis.Linear())
    theta = theta.interpolate(vertical, iris.analysis.Linear())
    
    x_wind_data = x_wind.data
    y_wind_data = y_wind.data
    pressure_data = pressure.data
    logpressure = np.log(pressure_data)
    
    theta_zonal_mean = theta.collapsed('longitude', iris.analysis.MEAN)
    theta_zonal_mean_data = theta_zonal_mean.data
    logpressure_zonal_mean = np.mean(logpressure, axis=3)
    pressure_zonal_mean = np.mean(pressure_data, axis=3)
    
    THETAp = ep_zderivative(theta_zonal_mean_data, pressure_zonal_mean)
    THETAp /= pressure_zonal_mean
    
    u_prime_list = []

    for timex in range(0, pressure.shape[0]):
        timex_list = []

        for levelx in range(0, pressure.shape[1]):
            levelx_list = []

            for latx in range(0, pressure.shape[2]):

                u_fft = sp.fftpack.fft(x_wind_data[timex, levelx, latx, :])
                u_psd = np.abs(u_fft)**2
                u_freq = sp.fftpack.fftfreq(len(u_psd), 1./143)

                highpass = u_fft.copy()
                highpass[np.abs(u_freq) < 1.1] = 0

                u_cleaned = np.real(sp.fftpack.ifft(highpass))
                u_bar = np.mean(u_cleaned)
                u_prime = u_cleaned - u_bar
                levelx_list.append(u_prime)

            timex_list.append(levelx_list)
            
        u_prime_list.append(timex_list)
            
    v_prime_list = []
    
    for timey in range(0, pressure.shape[0]):
        v_time_list = []
    
        for levely in range(0,pressure.shape[1]):
            v_level_list = []
            
            for laty in range(0,pressure.shape[2]):
            
                v_fft = sp.fftpack.fft(y_wind_data[timey,levely,laty,:])
                v_psd = np.abs(v_fft)**2
                v_freq = sp.fftpack.fftfreq(len(v_psd), 1./90)
                
                vhighpass = v_fft.copy()
                vhighpass[np.abs(v_freq) < 1.1] = 0
                
                v_cleaned = np.real(sp.fftpack.ifft(vhighpass))
                v_bar = np.mean(v_cleaned)
                v_prime = v_cleaned - v_bar
            
                v_level_list.append(v_prime)
                
            v_time_list.append(v_level_list)
                
        v_prime_list.append(v_time_list)

    x_anomaly_data = np.array(u_prime_list)    
    y_anomaly_data = np.array(v_prime_list)
    theta_anomaly = theta - theta_zonal_mean    
    theta_anomaly_data = theta_anomaly.data    
    
    uvprimes = x_anomaly_data*y_anomaly_data
    vthprimes = y_anomaly_data*theta_anomaly_data
    uvprimes_zonal_mean = np.mean(uvprimes, axis=3)
    vthprimes_zonal_mean = np.mean(vthprimes, axis=3)
    vertical_eddy = vthprimes_zonal_mean/THETAp

    lats = (pressure.coord('latitude').points)*(np.pi/180) # Latitudes in radians
    coriolis = 2*omega*np.sin(lats) # Coriolis force
    radius = 7160000 # radius of ProxB in meters
    
    u_bar = x_wind.collapsed('longitude', iris.analysis.MEAN)
    u_bar_data = u_bar.data
    
    EP_phi = np.cos(lats)*ep_latderivative(-uvprimes_zonal_mean,radius*np.sin(lats)) + 2*np.sin(lats)*uvprimes_zonal_mean
    EP_phi = EP_phi/(radius*np.cos(lats))
    EP_up = ep_zderivative(coriolis*vertical_eddy, pressure_zonal_mean)
    
    print(EP_phi[-1,-1,20], EP_up[-1,-1,20])
    print(coriolis[20],uvprimes_zonal_mean[-1,-1,20], THETAp[-1,-1,20])
    
    EP_phi_mean = np.mean(EP_phi, axis=0)
    EP_up_mean = np.mean(EP_up, axis=0)
    EP_div = EP_phi_mean + EP_up_mean
    
    x_mean = np.mean(u_bar_data, axis=0)
        
    x_axis = pressure.coord('latitude').points
    # y_axis = pressure_zonal_mean
    y_axis = pressure.coord('Hybrid height').points
    plt.contourf(x_axis[15:75], y_axis, EP_div[:,15:75], brewer_redblu.N, cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    cbar = plt.colorbar(pad=0.1)
    cbar.ax.set_title('units')
    CS = plt.contour(x_axis[15:75], y_axis, x_mean[:,15:75], colors='black', linewidths=0.5)
    plt.title('Mean EP flux divergence, day=%s to %s' %(start/4, end/4))
    plt.xlabel('Latitude [degrees]')
    plt.ylabel('Height [m]')
    plt.clabel(CS, inline=False, colors='k', fmt='%1.1f')
    plt.show()
    


    