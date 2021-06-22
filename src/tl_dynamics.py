#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 07:46:08 2021

@author: Maureen Cohen
"""
import os, iris
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import iris.plot as iplt
import iris.coords
from iris.coords import DimCoord
import numpy as np
import scipy as sp
import windspharm
from matplotlib.colors import TwoSlopeNorm

# Import packages


brewer_red = mpl_cm.get_cmap('brewer_YlOrRd_09')
brewer_bg = mpl_cm.get_cmap('brewer_PuBu_09')
brewer_reds = mpl_cm.get_cmap('brewer_Reds_09')
brewer_redblu = mpl_cm.get_cmap('brewer_RdBu_11')


def plot_uw(cubes, start=0, end=120, lat=45, levels=(37,44,47,53)):
    
    for cube in cubes:
        if cube.standard_name == 'eastward_wind' or cube.standard_name == 'x_wind':
            x_wind = cube[start:end,:,:,:].copy() 
        if cube.standard_name == 'northward_wind' or cube.standard_name == 'y_wind':
            y_wind = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'upward_air_velocity':
            z_wind = cube[start:end,:,:,:].copy()
    
    run_length, longitudes, latitudes = x_wind.shape[0], x_wind.shape[3]/2, x_wind.shape[2]/2
            
    y_wind, z_wind = y_wind.regrid(x_wind, iris.analysis.Linear()), z_wind.regrid(x_wind, iris.analysis.Linear())
    
    for coord in x_wind.coords():
        if coord.long_name == 'Hybrid height':
            heights = np.round(x_wind.coord('Hybrid height').points*1e-03,0)
            vertical = [('Hybrid height', x_wind.coord('Hybrid height').points)]
        elif coord.long_name == 'level_height':
            heights = np.round(x_wind.coord('level_height').points*1e-03,0)
            vertical = [('level_height', x_wind.coord('level_height').points)]
            
    z_wind = z_wind.interpolate(vertical, iris.analysis.Linear())
    
    zonal_mean_u = x_wind.collapsed('t', iris.analysis.MEAN)
    u_prime = x_wind - zonal_mean_u
    zonal_mean_w = z_wind.collapsed('t', iris.analysis.MEAN)
    w_prime = z_wind - zonal_mean_w
    cross_term = u_prime*w_prime
    
    
    for level in levels:
        
        plt.figure(figsize=(10,5))
        plt.contourf(np.arange(-longitudes,longitudes), np.arange(0, run_length), np.roll(cross_term[:,level,lat,:].data, 72, axis=1),brewer_redblu.N, cmap=brewer_redblu, norm=TwoSlopeNorm(0))
        plt.title('$U^{\prime} \cdot W^{\prime}$ at Equator [m2 s-2], h=%s km, %s to %s days' %(heights[level], start/4, end/4))
        plt.xlabel('Longitude [degrees]')
        plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
        plt.ylabel('Time [6-hours]')
        plt.colorbar(pad=0.1)
        plt.show()
        
        plt.figure(figsize=(10,5))
        plt.contourf(np.arange(-longitudes,longitudes), np.arange(-latitudes,latitudes), np.roll(cross_term[-1,level,:,:].data, 72, axis=1),brewer_redblu.N, cmap=brewer_redblu, norm=TwoSlopeNorm(0))
        plt.title('$U^{\prime} \cdot W^{\prime}$ [m2 s-2], h=%s km, day=%s' %(heights[level], end/4))
        plt.xlabel('Longitude [degrees]')
        plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
        plt.ylabel('Latitude [degrees]')
        plt.colorbar(pad=0.1)
        plt.show()
        
def plot_uv(cubes, start=0, end=240, levels=(37,44,47,53)):
    
    for cube in cubes:
        if cube.standard_name == 'eastward_wind' or cube.standard_name == 'x_wind':
            x_wind = cube[start:end,:,:,:].copy() 
        if cube.standard_name == 'northward_wind' or cube.standard_name == 'y_wind':
            y_wind = cube[start:end,:,:,:].copy()
    
    longitudes, latitudes = x_wind.shape[3]/2, x_wind.shape[2]/2
            
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
    
    for coord in x_wind.coords():
        if coord.long_name == 'Hybrid height':
            heights = np.round(x_wind.coord('Hybrid height').points*1e-03,0)
        elif coord.long_name == 'level_height':
            heights = np.round(x_wind.coord('level_height').points*1e-03,0)            
    
    zonal_mean_u = x_wind.collapsed('t', iris.analysis.MEAN)
    u_prime = x_wind - zonal_mean_u
    zonal_mean_v = y_wind.collapsed('t', iris.analysis.MEAN)
    v_prime = y_wind - zonal_mean_v
    cross_term = u_prime*v_prime
        
    for level in levels:
        
        plt.figure(figsize=(10,5))
        plt.contourf(np.arange(-longitudes,longitudes), np.arange(-latitudes,latitudes), np.roll(cross_term[end-1,level,:,:].data, 72, axis=1),brewer_redblu.N, cmap=brewer_redblu, norm=TwoSlopeNorm(0))
        plt.title('$U^{\prime} \cdot V^{\prime}$ [m2 s-2], h=%s km, day=%s' %(heights[level], end/4))
        plt.xlabel('Longitude [degrees]')
        plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
        plt.ylabel('Latitude [degrees]')
        plt.colorbar(pad=0.1)
        plt.show() 
        

def wave_acceleration(cubes, lat=45, long=0, start=0, end=240):
    
    for cube in cubes:
        if cube.standard_name == 'eastward_wind':
            x_wind = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'northward_wind':
            y_wind = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'upward_air_velocity':
            z_wind = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'air_pressure':
            pressure = cube[start:end,:,:,:].copy()
            
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
    z_wind = z_wind.regrid(x_wind, iris.analysis.Linear())
    z_wind = -1*z_wind
    pressure = pressure.regrid(x_wind, iris.analysis.Linear())
    
    heights = pressure.coord('Hybrid height').points
    
    height = [('Hybrid height', pressure.coord('Hybrid height').points)]
    x_wind, y_wind = x_wind.interpolate(height, iris.analysis.Linear()), y_wind.interpolate(height, iris.analysis.Linear())

    longitudes = x_wind.shape[3]/2
    # x_zonal_mean = x_wind.collapsed('longitude', iris.analysis.MEAN)
    # avg_zonal_mean = x_zonal_mean.collapsed('t', iris.analysis.MEAN)

    x_mean = x_wind.collapsed('t', iris.analysis.MEAN)
    z_mean = z_wind.collapsed('t', iris.analysis.MEAN)
    
    x_anomaly = x_wind - x_mean
    z_anomaly = z_wind - z_mean
    
    result = x_wind.copy()
    result.data = x_anomaly.data*z_anomaly.data
    array = result.data
    # h = np.log(pressure.data)
    h = np.array([heights])
    h = np.repeat(h[:,np.newaxis], (end-start), axis=0)
    h = np.reshape(h, ((end-start),61))
    h = np.repeat(h[:,:,np.newaxis], 90, axis=2)
    h = np.repeat(h[:,:,:,np.newaxis], 144, axis=3)

    print(h.shape)
    
    acceleration = array.copy()
    acceleration[:,0,...] = (array[:,1,...]-array[:,0,...])/(h[:,1,...]-h[:,0,...])
    acceleration[:,-1,...] = (array[:,-1,...]-array[:,-2,...])/(h[:,-1,...]-h[:,-2,...])
    acceleration[:,1:-1,...] = (array[:,2:,...]-array[:,0:-2,...])/(h[:,2:,...]-h[:,0:-2,...])
    
    # acceleration = array.copy()
    # acceleration[:,0,...] = (array[:,1,...]-array[:,0,...])/(h[1]-h[0])
    # acceleration[:,-1,...] = (array[:,-1,...]-array[:,-2,...])/(h[-1]-h[-2])
    # acceleration[:,1:-1,...] = (array[:,2:,...]-array[:,0:-2,...])/(h[2:]-h[0:-2])
    
    # zonal_acc = np.mean(acceleration, axis=3)
    # net_acc = np.sum(acceleration, axis=3)
    # net_acc_total = np.mean(zonal_acc, axis=0)
    
    substellar = acceleration[:,:,31:59:,long]
    net_acc = np.mean(substellar, axis=0) # + x_wind[start,:,31:59,long].data

    plt.figure(figsize=(10,5))
    plt.contourf(np.arange(-longitudes,longitudes), np.array(heights), np.roll(x_anomaly[end-start-1,:,lat,:].data, 72, axis=1), brewer_redblu.N, cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    plt.title('Zonal Wind Anomaly at Equator [m s-1], t = %s days' %(end/4))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.ylabel('Height [m]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    plt.figure(figsize=(10,5))
    plt.contourf(np.arange(-longitudes,longitudes), np.array(heights), np.roll(z_anomaly[end-start-1,:,lat,:].data, 72, axis=1), np.linspace(-0.05,0.05,20), cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    plt.title('Vertical Wind Anomaly at Equator [m s-1], t = %s days' %(end/4))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.ylabel('Height [m]')
    plt.colorbar(pad=0.1)
    plt.show()

    plt.figure(figsize=(10,5))
    plt.contourf(np.arange(-longitudes,longitudes), np.array(heights), np.roll(acceleration[end-start-1,:,lat,:], 72, axis=1), np.linspace(-0.0024,0.0024,20), cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    plt.title('Wave-Induced Acceleration at Equator [m s-2], t = %s days' %(end/4))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.ylabel('Height [m]')
    plt.colorbar(pad=0.1)
    plt.show()
                
    x_axis = pressure.coord('latitude').points
    y_axis = pressure.coord('Hybrid height').points
    plt.figure(figsize=(6,10))
    plt.contourf(x_axis[31:59], y_axis, net_acc, brewer_redblu.N, cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    plt.colorbar(pad=0.1)
    CS = plt.contour(x_axis[31:59], y_axis, x_anomaly[end-start-1,:,31:59,long].data, colors='black', linewidths=0.5)
    plt.title('Mean acceleration at long %s from day=%s to %s' %(long*2.5, start/4, end/4))
    plt.xlabel('Latitude [degrees]')
    plt.ylabel('Height [m]')
    plt.clabel(CS, inline=False, colors='k', fmt='%1.1f')
    plt.show()
    
    
    # acc_cube = x_wind.copy()
    # acc_cube.data = acceleration
    # substellar = acc_cube.intersection(longitude=(-90,90))
    # plt.figure(figsize=(5,5))
    # plt.contourf(np.arange(-36,37), np.array(heights[:58]), substellar[end-1,:58,:].data, np.linspace(-0.5,0.5,20), cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    # plt.title('Wave-Induced Acceleration at Equator [m s-2], t = %s days' %(end/4))
    # plt.xlabel('Longitude [degrees]')
    # plt.xticks((-36,-24,-12,0,12,24,36),('90W','60W','30W','0','30E','60E','90E'))
    # plt.ylabel('Height [m]')
    # plt.colorbar(pad=0.1)
    # plt.show()  
    
    
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

def ep_flux(cubes, omega=0.64617667e-05, start=0, end=240):
    
    for cube in cubes:
        if cube.standard_name == 'air_potential_temperature':
            theta = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'eastward_wind':
            x_wind = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'northward_wind':
            y_wind = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'air_pressure':
            pressure = cube[start:end,:,:,:].copy()
    
    heights, longitudes, longs = np.round(pressure.coord('Hybrid height').points*1e-03,0), pressure.coord('longitude').points, pressure.shape[3]/2
    
    y_wind, x_wind = y_wind.regrid(pressure, iris.analysis.Linear()), x_wind.regrid(pressure, iris.analysis.Linear())
    
    vertical = [('Hybrid height', pressure.coord('Hybrid height').points)]
    x_wind, y_wind = x_wind.interpolate(vertical, iris.analysis.Linear()), y_wind.interpolate(vertical, iris.analysis.Linear())
    
    
    theta_data = theta.data
    x_wind_data = x_wind.data
    y_wind_data = y_wind.data
    pressure_data = pressure.data
    logpressure = np.log(pressure_data)
    
    theta_zonal_mean = theta.collapsed('t', iris.analysis.MEAN)
    theta_zonal_mean_data = theta_zonal_mean.data
    logpressure_zonal_mean = np.mean(logpressure, axis=2)
    pressure_zonal_mean = np.mean(pressure_data, axis=2)
    
    THETAp = ep_zderivative(theta_zonal_mean_data, logpressure_zonal_mean)
    THETAp /= pressure_zonal_mean
    
    # u_prime_lats = []
    
    # for level in range(0,61):
    #     u_prime_heights = []
        
    #     for lat in range(0,90):
        
    #         u_fft = sp.fftpack.fft(x_wind_data[level,lat,:])
    #         u_psd = np.abs(u_fft)**2
    #         u_freq = sp.fftpack.fftfreq(len(u_psd), 1./144)
            
    #         highpass = u_fft.copy()
    #         highpass[np.abs(u_freq) < 1.1] = 0
            
    #         u_cleaned = np.real(sp.fftpack.ifft(highpass))
    #         u_bar = np.mean(u_cleaned)
    #         u_prime = u_cleaned - u_bar
        
    #         u_prime_heights.append(u_prime)
            
    #     u_prime_lats.append(u_prime_heights)

    y_zonal_mean = y_wind.collapsed('t', iris.analysis.MEAN)
    x_zonal_mean = x_wind.collapsed('t', iris.analysis.MEAN)
    # x_anomaly_data = np.array(u_prime_lats)
    x_anomaly = x_wind - x_zonal_mean
    y_anomaly = y_wind - y_zonal_mean 
    theta_anomaly = theta - theta_zonal_mean
    
    y_anomaly_data = y_anomaly.data
    theta_anomaly_data = theta_anomaly.data
    x_anomaly_data = x_anomaly.data
    
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