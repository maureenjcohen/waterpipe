#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 17 14:15:57 2021

@author: Mo Cohen
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

def calculate_density(cubes):
    
    for cube in cubes:
        if cube.standard_name == 'air_potential_temperature':
            theta = cube.copy()
        if cube.standard_name == 'specific_humidity':
            q = cube.copy()
        if cube.standard_name == 'air_pressure':
            pressure = cube.copy()  
            
    R = 287.05
    cp = 1005        
    p0 = iris.coords.AuxCoord(100000.0, long_name='reference_pressure', units='Pa')
    p0.convert_units(pressure.units)
    absolute_temp = theta*((pressure/p0)**(R/cp)) # R and cp in J/kgK for 300K
    
    density = pressure/(R*absolute_temp*(1+0.61*q))
    
    return density.data


def ep_zderivative(array,h):
    work = array.copy()
    work[:,0,...] = (array[:,1,...]-array[:,0,...])/(h[:,1,:]-h[:,0,:])
    work[:,-1,...] = (array[:,-1,...]-array[:,-2,...])/(h[:,-1,:]-h[:,-2,:])
    work[:,1:-1,...] = (array[:,2:,...]-array[:,1:-1,...])/(h[:,2:,:]-h[:,1:-1,:])
    
    return work

def ep_latderivative(array,l):
    
    work = array.copy()
    work[:,:,0] = (array[:,:,1]-array[:,:,0])/(l[1]-l[0])
    work[:,:,-1] = (array[:,:,-1]-array[:,:,-2])/(l[-1]-l[-2])
    work[:,:,1:-1] = (array[:,:,2:]-array[:,:,1:-1])/(l[2:]-l[1:-1])
    
    return work

def ep_flux(cubes, density, g=10.9, omega=0.64617667e-05, long_slice=(0,-1), start=3235, end=3240):
    
    for cube in cubes:
        if cube.standard_name == 'air_potential_temperature':
            theta = cube[start:end,:,:,long_slice[0]:long_slice[1]].copy()
        if cube.standard_name == 'eastward_wind':
            x_wind = cube[start:end,:,:,long_slice[0]:long_slice[1]].copy()
        if cube.standard_name == 'northward_wind':
            y_wind = cube[start:end,:,:,long_slice[0]:long_slice[1]].copy()
        if cube.standard_name == 'air_pressure':
            pressure = cube[start:end,:,:,long_slice[0]:long_slice[1]].copy()
    
    heights, longitudes, longs = np.round(pressure.coord('Hybrid height').points*1e-03,0), pressure.coord('longitude').points, pressure.shape[3]/2
    
    y_wind, x_wind = y_wind.regrid(pressure, iris.analysis.Linear()), x_wind.regrid(pressure, iris.analysis.Linear())
    theta = theta.regrid(pressure, iris.analysis.Linear())
    vertical = [('Hybrid height', pressure.coord('Hybrid height').points)]
    x_wind, y_wind = x_wind.interpolate(vertical, iris.analysis.Linear()), y_wind.interpolate(vertical, iris.analysis.Linear())
    theta = theta.interpolate(vertical, iris.analysis.Linear())

    theta_data = theta.data
    x_wind_data = x_wind.data
    y_wind_data = y_wind.data
    pressure_data = pressure.data
    logpressure = np.log(pressure_data)
    
    theta_zonal_mean = theta.collapsed('longitude', iris.analysis.MEAN)
    dth_dz = iris.analysis.calculus.differentiate(theta_zonal_mean, 'Hybrid height')
    dth_dz = dth_dz.data
    density = np.mean(density[start:end,:-1,:,:], axis=3)
    THETAp = (-1/(g*density))*dth_dz
    
    # theta_zonal_mean_data = theta_zonal_mean.data
    logpressure_zonal_mean = np.mean(logpressure, axis=3)
    # pressure_zonal_mean = np.mean(pressure_data, axis=3)
    
    # THETAp = ep_zderivative(theta_zonal_mean_data, logpressure_zonal_mean)
    # THETAp /= pressure_zonal_mean
    
    u_prime_list = []

    for time in range(0, pressure.shape[0]):
        time_list = []

        for level in range(0, pressure.shape[1]):
            level_list = []

            for latx in range(0, pressure.shape[2]):

                u_fft = sp.fftpack.fft(x_wind_data[time, level, latx, :])
                u_psd = np.abs(u_fft)**2
                u_freq = sp.fftpack.fftfreq(len(u_psd), 1./143)

                highpass = u_fft.copy()
                highpass[np.abs(u_freq) < 5.1] = 0

                u_cleaned = np.real(sp.fftpack.ifft(highpass))
                u_bar = np.mean(u_cleaned)
                u_prime = u_cleaned - u_bar
                level_list.append(u_prime)

            time_list.append(level_list)
            
        u_prime_list.append(time_list)

            
    v_prime_list = []
    
    for time in range(0, pressure.shape[0]):
        v_time_list = []
    
        for level in range(0,pressure.shape[1]):
            v_level_list = []
            
            for laty in range(0,pressure.shape[2]):
            
                v_fft = sp.fftpack.fft(y_wind_data[time,level,laty,:])
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

    # y_zonal_mean = y_wind.collapsed('longitude', iris.analysis.MEAN)
    x_anomaly_data = np.array(u_prime_list)    
    # y_anomaly = y_wind - y_zonal_mean
    y_anomaly_data = np.array(v_prime_list)
    print(x_anomaly_data.shape, y_anomaly_data.shape)
    theta_anomaly = theta - theta_zonal_mean
    
    # y_anomaly_data = y_anomaly.data
    theta_anomaly_data = theta_anomaly.data
    
    XY_anomaly = x_anomaly_data*y_anomaly_data
    YTH_anomaly = y_anomaly_data*theta_anomaly_data
    XY_anomaly_zonal_mean = np.mean(XY_anomaly, axis=3)
    YTH_anomaly_zonal_mean = np.mean(YTH_anomaly, axis=3)

    lats = (pressure.coord('latitude').points)*(np.pi/180) # Latitudes in radians
    coriolis = 2*omega*np.sin(lats) # Coriolis force
    radius = 7160000 # radius of ProxB in meters
    rcosphi = radius*np.cos(lats)
    rsinphi = radius*np.sin(lats)
    latfac = rcosphi*np.cos(lats)
    
    # EP_phi = -XY_anomaly_zonal_mean*np.cos(lats)*np.cos(lats)
    # EP_up = (coriolis*latfac*YTH_anomaly_zonal_mean[:,:-1,:])/THETAp
    EP_phi = -XY_anomaly_zonal_mean*latfac
    EP_up = (coriolis*rcosphi*YTH_anomaly_zonal_mean[:,:-1,:])/THETAp
    
    h = pressure.coord('Hybrid height').points
    h = np.array([h])
    h = np.repeat(h[:, np.newaxis], (end-start), axis=0)
    h = np.reshape(h, ((end-start), 61))
    h = np.repeat(h[:, :, np.newaxis], 90, axis=2)
    
    EP_phi_div = ep_latderivative(EP_phi, rsinphi)
    EP_up_div = (-1/(g*density))*ep_zderivative(EP_up, h[:,:-1,:])

    
    EP_div = np.mean((EP_phi_div[:,:-1,:] + EP_up_div), axis=0)*(1/rcosphi)
    # EP_div = np.mean((EP_phi_div[:,:-1,:]), axis=0)*(1/rcosphi)

    
    # x_zonal_mean = x_wind.collapsed('longitude', iris.analysis.MEAN)
    # x_zonal_mean = x_zonal_mean.data
    x_mean = x_wind.collapsed(['longitude'], iris.analysis.MEAN)
        
    x_axis = pressure.coord('latitude').points
    # y_axis = pressure_zonal_mean
    y_axis = pressure.coord('Hybrid height').points
    plt.figure(figsize=(8, 6))

    plt.contourf(x_axis, y_axis[35:59], EP_div[35:59,:], np.arange(-2e-05,2e-05,0.5e-06), cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    cbar = plt.colorbar(pad=0.1)
    cbar.set_label('m/s$^2$', size=15)
    CS = plt.contour(x_axis, y_axis[35:59], x_mean[-1,35:59,:].data - x_mean[0,35:59,:].data, colors='black', linewidths=1.5)
    plt.title('Mean EP flux divergence, day=%s to %s' %(start/4, end/4))
    plt.xlabel('Latitude [degrees]')
    plt.ylabel('Height [m]')
    plt.clabel(CS, inline=False, colors='k', fmt='%1.1f')
    plt.show()
    


def plot_ywind(cubes, level=47, time_slice=1780):

    for cube in cubes:
        if cube.standard_name =='northward_wind':
            y_wind = cube.copy()
    
    heights = np.round(y_wind.coord('Hybrid height').points*1e-03,0)
    
    iplt.contourf(y_wind[time_slice,level,:,:], cmap=brewer_redblu)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Meridional wind, h=%s, day=%s' %(heights[level],time_slice*0.25))
    cbar = plt.colorbar(orientation='horizontal')
    cbar.ax.set_title('m/s')
    plt.show()
    