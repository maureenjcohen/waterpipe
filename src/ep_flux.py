#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 17 14:15:57 2021

@author: s1144983
"""

import os, iris, cartopy, windspharm, pywt
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
    work[0,...] = (array[1,...]-array[0,...])/(h[1,...]-h[0,...])
    work[-1,...] = (array[-1,...]-array[-2,...])/(h[-1,...]-h[-2,...])
    work[1:-1,...] = (array[2:,...]-array[0:-2,...])/(h[2:,...]-h[0:-2,...])
    
    return work

def ep_latderivative(array,l):
    
    work = array.copy()
    work[:,0,...] = (array[:,1,...]-array[:,0,...])/(l[1,...]-l[0,...])
    work[:,-1,...] = (array[:,-1,...]-array[:,-2,...])/(l[-1,...]-l[-2,...])
    work[:,1:-1,...] = (array[:,2:,...]-array[:,0:-2,...])/(l[2:,...]-l[0:-2,...])
    
    return work

def ep_flux(cubes, omega=0.64617667e-05, long_slice=(0,-1), time_slice=1780):
    
    for cube in cubes:
        if cube.standard_name == 'air_potential_temperature':
            theta = cube[time_slice,:,:,long_slice[0]:long_slice[1]].copy()
        if cube.standard_name == 'eastward_wind':
            x_wind = cube[time_slice,:,:,long_slice[0]:long_slice[1]].copy()
        if cube.standard_name == 'northward_wind':
            y_wind = cube[time_slice,:,:,long_slice[0]:long_slice[1]].copy()
        if cube.standard_name == 'air_pressure':
            pressure = cube[time_slice,:,:,long_slice[0]:long_slice[1]].copy()
    
    heights, longitudes, longs = np.round(pressure.coord('Hybrid height').points*1e-03,0), pressure.coord('longitude').points, pressure.shape[2]/2
    
    y_wind, x_wind = y_wind.regrid(pressure, iris.analysis.Linear()), x_wind.regrid(pressure, iris.analysis.Linear())
    
    vertical = [('Hybrid height', pressure.coord('Hybrid height').points)]
    x_wind, y_wind = x_wind.interpolate(vertical, iris.analysis.Linear()), y_wind.interpolate(vertical, iris.analysis.Linear())
    
    # for cube in (theta, x_wind, y_wind, pressure):
    #     cube = cube.intersection(longitude=long_slice)
    
    theta_data = theta.data
    x_wind_data = x_wind.data
    y_wind_data = y_wind.data
    pressure_data = pressure.data
    logpressure = np.log(pressure_data)
    
    theta_zonal_mean = theta.collapsed('longitude', iris.analysis.MEAN)
    theta_zonal_mean_data = theta_zonal_mean.data
    logpressure_zonal_mean = np.mean(logpressure, axis=2)
    pressure_zonal_mean = np.mean(pressure_data, axis=2)
    
    THETAp = ep_zderivative(theta_zonal_mean_data, logpressure_zonal_mean)
    THETAp /= pressure_zonal_mean
    
    u_prime_lats = []
    
    for level in range(0,61):
        u_prime_heights = []
        
        for lat in range(0,90):
        
            u_fft = sp.fftpack.fft(x_wind_data[level,lat,:])
            u_psd = np.abs(u_fft)**2
            u_freq = sp.fftpack.fftfreq(len(u_psd), 1./144)
            
            highpass = u_fft.copy()
            highpass[np.abs(u_freq) < 1.1] = 0
            
            u_cleaned = np.real(sp.fftpack.ifft(highpass))
            u_bar = np.mean(u_cleaned)
            u_prime = u_cleaned - u_bar
        
            u_prime_heights.append(u_prime)
            
        u_prime_lats.append(u_prime_heights)

    y_zonal_mean = y_wind.collapsed('longitude', iris.analysis.MEAN)
    x_anomaly_data = np.array(u_prime_lats)
    y_anomaly = y_wind - y_zonal_mean 
    theta_anomaly = theta - theta_zonal_mean
    
    y_anomaly_data = y_anomaly.data
    theta_anomaly_data = theta_anomaly.data
    
    XY_anomaly = x_anomaly_data*y_anomaly_data
    YTH_anomaly = y_anomaly_data*theta_anomaly_data
    XY_anomaly_zonal_mean = np.mean(XY_anomaly, axis=2)
    YTH_anomaly_zonal_mean = np.mean(YTH_anomaly, axis=2)

    lats = (pressure.coord('latitude').points)*(np.pi/180) # Latitudes in radians
    coriolis = 2*omega*np.sin(lats) # Coriolis force
    radius = 7160000 # radius of ProxB in meters
    rcosphi = radius*np.cos(lats)
    rsinphi = radius*np.sin(lats)
    latfac = rcosphi*np.cos(lats)
    
    EP_phi = -XY_anomaly_zonal_mean*latfac
    EP_up = (coriolis*rcosphi*YTH_anomaly_zonal_mean)/THETAp
    
    EP_phi_div = ep_latderivative(EP_phi, rsinphi)
    EP_up_div = ep_zderivative(EP_up, pressure_zonal_mean)
    EP_div = (EP_phi_div + EP_up_div)
    
    # x_zonal_mean = x_wind.collapsed('longitude', iris.analysis.MEAN)
    x_mean = x_wind.collapsed('longitude', iris.analysis.MEAN)
        
    x_axis = pressure.coord('latitude').points
    # y_axis = pressure_zonal_mean
    y_axis = pressure.coord('Hybrid height').points
    plt.contourf(x_axis, y_axis[35:59], EP_div[35:59,:], brewer_redblu.N, cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    plt.colorbar(pad=0.1)
    CS = plt.contour(x_axis, y_axis[35:59], x_mean[35:59,:].data, colors='black', linewidths=0.5)
    plt.title('EP flux divergence for day=%s, long=%s to %s' %(time_slice/4, longitudes[0], longitudes[-1]))
    plt.xlabel('Latitude [degrees]')
    plt.ylabel('Height [m]')
    plt.clabel(CS, inline=False, colors='k', fmt='%1.1f')
    plt.show()
    
    # THETAp_non = ep_zderivative(theta_data, logpressure)
    # THETAp_non /= logpressure

    # EP_phi_non = -XY_anomaly*np.repeat(latfac[:,np.newaxis], 144,1)
    # EP_up_non = (np.repeat(coriolis[:,np.newaxis], 144,1)*np.repeat(rcosphi[:,np.newaxis], 144,1)*YTH_anomaly)/THETAp_non
    # EP_phi_non_div = ep_latderivative(EP_phi_non, np.repeat(rsinphi[:,np.newaxis], 144,1))
    # EP_up_non_div = ep_zderivative(EP_up_non, logpressure)
    # EP_div_non = EP_phi_non_div + EP_up_non_div
    # EP_div_non_mean = np.mean(EP_div_non, axis=0)
    
    # plt.contourf(np.roll(pressure.coord('longitude').points, 72), pressure.coord('latitude').points, EP_div_non_mean[level,:,:], brewer_redblu.N, cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    # plt.title('EP flux divergence for day=%s to %s, h=%s km' %(start/4, end/4, heights[level]))
    # plt.xlabel('Longitude [degrees]')
    # plt.ylabel('Latitude [degrees]')
    # plt.colorbar(pad=0.1)
    # plt.show()


    