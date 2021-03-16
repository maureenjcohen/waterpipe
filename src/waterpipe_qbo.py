#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 15:24:49 2021

@author: Mo Cohen

Functions relating to QBO-like phenomena

"""

import os, iris
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import iris.plot as iplt
import iris.coords
from iris.coords import DimCoord
from iris.analysis import calculus
import numpy as np
import scipy as sp
import windspharm
# Import packages


brewer_red = mpl_cm.get_cmap('brewer_YlOrRd_09')
brewer_bg = mpl_cm.get_cmap('brewer_PuBu_09')
brewer_reds = mpl_cm.get_cmap('brewer_Reds_09')
brewer_redblu = mpl_cm.get_cmap('brewer_RdBu_11')

def get_heights(cubes):
    
    """ Extract altitude information from air pressure cube
        Can be used with matplotlib contourf when Iris isn't playing nice """ 
    
    for cube in cubes:
        if cube.standard_name == 'air_pressure':
            pressure = cube.copy()
            
    heights = pressure.coord('level_height').points
    
    return heights


def plot_hovmoellerx(cubes):
    
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube.copy()
    
    heights = x_wind.coord('level_height').points
    run_length = x_wind.shape[0]
    lats = x_wind.coord('latitude')
    longs = x_wind.coord('longitude')

    if lats.bounds == None:
        x_wind.coord('latitude').guess_bounds()
    if longs.bounds == None:
        x_wind.coord('longitude').guess_bounds()
        
    dayside = x_wind.extract(iris.Constraint(longitude=lambda v: 270 < v <= 359 or 0 <= v <= 90, latitude=lambda v: -10 <= v <= 10))
    nightside = x_wind.extract(iris.Constraint(longitude=lambda v: 90 < v <= 270, latitude=lambda v: -10 <= v <= 10))
    dayside_full = x_wind.extract(iris.Constraint(longitude=lambda v: 270 < v <= 359 or 0 <= v <= 90, latitude=lambda v: -90 <= v <= 90))
    nightside_full = x_wind.extract(iris.Constraint(longitude=lambda v: 90 < v <= 270, latitude=lambda v: -90 <= v <= 90))
    day_grid = iris.analysis.cartography.area_weights(dayside)
    night_grid = iris.analysis.cartography.area_weights(nightside)
    # fullday_grid = iris.analysis.cartography.area_weights(dayside_full)
    # fullnight_grid = iris.analysis.cartography.area_weights(nightside_full)
    
    dayside_mean = dayside.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=day_grid)
    nightside_mean = nightside.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=night_grid)
    
    dayside_zonal_mean = dayside_full.collapsed('longitude', iris.analysis.MEAN)
    nightside_zonal_mean = nightside_full.collapsed('longitude', iris.analysis.MEAN)  
    print(dayside_zonal_mean.shape)
    

    plt.contourf(np.arange(0,run_length), np.array(heights), dayside_mean.data.T, levels=np.linspace(-80,130,20), cmap=brewer_redblu)
    plt.title('Dayside Mean Zonal Equatorial Wind [m s-1]')
    plt.xlabel('Time [months]')
    plt.ylabel('Height [m]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    plt.contourf(np.arange(0,run_length), np.array(heights), nightside_mean.data.T, levels=np.linspace(-80,130,20), cmap=brewer_redblu)
    plt.title('Nightside Mean Zonal Equatorial Wind [m s-1]')
    plt.xlabel('Time [months]')
    plt.ylabel('Height [m]')
    plt.colorbar(pad=0.1)
    plt.show()    
    
    plt.contourf(np.arange(0,run_length), np.arange(-45,45), dayside_zonal_mean[:,47].data.T, levels=np.linspace(-80,130,20), cmap=brewer_redblu)
    plt.title('Dayside Mean Zonal Wind at 40 km [m s-1]')
    plt.xlabel('Time [months]')
    plt.ylabel('Latitude')
    plt.yticks((-45,0,45),('90S', '0', '90N'))
    plt.colorbar(pad=0.1)
    plt.show()
    
    plt.contourf(np.arange(0,run_length), np.arange(-45,45), nightside_zonal_mean[:,47].data.T, levels=np.linspace(-80,130,20), cmap=brewer_redblu)
    plt.title('Nightside Mean Zonal Wind at 40 km [m s-1]')
    plt.xlabel('Time [months]')
    plt.ylabel('Latitude')
    plt.yticks((-45,0,45),('90S', '0', '90N'))
    plt.colorbar(pad=0.1)
    plt.show()


def plot_temp_anomaly(cubes):
    
    for cube in cubes:
        if cube.standard_name == 'air_potential_temperature':
            theta = cube.copy()
    # Extract potential temperature cube
        
    run_length = theta.shape[0]
    longitudes = theta.shape[3]/2
    theta_time_mean = theta.collapsed('time', iris.analysis.MEAN)
    theta_zonal_mean = theta.collapsed('longitude',iris.analysis.MEAN)
    
    anomaly = theta - theta_time_mean
    spatial_anomaly = theta - theta_zonal_mean
            
    plt.contourf(np.arange(-longitudes,longitudes), np.arange(0, run_length), np.roll(anomaly[:,35,45,:].data, 72, axis=1),brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Temperature Anomaly at Equator [K], h=25 km')
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-52,-32,-12,0,12,32,52,72),('180W','140W','100W','60W','0','60E','100E','140E','180E'))
    plt.ylabel('Time [months]')
    plt.colorbar(pad=0.1)
    plt.show()  
    
    plt.contourf(np.arange(-longitudes,longitudes), np.arange(0, run_length), np.roll(anomaly[:,44,45,:].data, 72, axis=1), brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Temperature Anomaly at Equator [K], h=35 km')
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-52,-32,-12,0,12,32,52,72),('180W','140W','100W','60W','0','60E','100E','140E','180E'))
    plt.ylabel('Time [months]')
    plt.colorbar(pad=0.1)
    plt.show()  
    
    plt.contourf(np.arange(-longitudes,longitudes), np.arange(0, run_length), np.roll(anomaly[:,47,45,:].data, 72, axis=1), brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Temperature Anomaly at Equator [K], h=40 km')
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-52,-32,-12,0,12,32,52,72),('180W','140W','100W','60W','0','60E','100E','140E','180E'))
    plt.ylabel('Time [months]')
    plt.colorbar(pad=0.1)
    plt.show() 
        
    plt.contourf(np.arange(-longitudes,longitudes), np.arange(0, run_length), np.roll(spatial_anomaly[:,35,45,:].data, 72, axis=1), brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Temperature Spatial Anomaly at Equator [K],, h=25 km')
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-52,-32,-12,0,12,32,52,72),('180W','140W','100W','60W','0','60E','100E','140E','180E'))
    plt.ylabel('Time [months]')
    plt.colorbar(pad=0.1)
    plt.show()  
    
    plt.contourf(np.arange(-longitudes,longitudes), np.arange(0, run_length), np.roll(spatial_anomaly[:,44,45,:].data, 72, axis=1), brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Temperature Spatial Anomaly at Equator [K], h=35 km')
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-52,-32,-12,0,12,32,52,72),('180W','140W','100W','60W','0','60E','100E','140E','180E'))
    plt.ylabel('Time [months]')
    plt.colorbar(pad=0.1)
    plt.show()  
    
    plt.contourf(np.arange(-longitudes,longitudes), np.arange(0, run_length), np.roll(spatial_anomaly[:,47,45,:].data, 72, axis=1), brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Temperature Spatial Anomaly at Equator [K],, h=40 km')
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-52,-32,-12,0,12,32,52,72),('180W','140W','100W','60W','0','60E','100E','140E','180E'))
    plt.ylabel('Time [months]')
    plt.colorbar(pad=0.1)
    plt.show() 
    

def plot_wind_anomaly(cubes):
    
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube.copy()
    # Extract potential temperature cube
    
    heights = x_wind.coord('level_height').points
    run_length = x_wind.shape[0]
    longitudes = x_wind.shape[3]/2

    x_wind_time_mean = x_wind.collapsed('time', iris.analysis.MEAN)
    x_wind_zonal_mean = x_wind.collapsed('longitude',iris.analysis.MEAN)
    
    anomaly = x_wind - x_wind_time_mean
    spatial_anomaly = x_wind - x_wind_zonal_mean
    
#   anomaly = theta.copy()
    # for t in range(0,run_length):
    #     anomaly = theta[t,:,:,:] - theta_time_mean
        
    
    plt.contourf(np.arange(-longitudes,longitudes), np.array(heights), np.roll(anomaly[39,:,45,:].data, 72, axis=1), brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Zonal Wind Time Anomaly at Equator, t = 40 months')
    plt.xlabel('Longitude [degrees]')
    plt.ylabel('Height [m]')
    plt.xticks((-72,-52,-32,-12,0,12,32,52,72),('180W','140W','100W','60W','0','60E','100E','140E','180E'))
    plt.colorbar(pad=0.1)
    plt.show()
    
    plt.contourf(np.arange(-longitudes,longitudes), np.array(heights), np.roll(anomaly[47,:,45,:].data, 72, axis=1), brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Zonal Wind Time Anomaly at Equator, t = 48 months')
    plt.xlabel('Longitude [degrees]')
    plt.ylabel('Height [m]')
    plt.xticks((-72,-52,-32,-12,0,12,32,52,72),('180W','140W','100W','60W','0','60E','100E','140E','180E'))
    plt.colorbar(pad=0.1)
    plt.show()  
    
    plt.contourf(np.arange(-longitudes,longitudes), np.array(heights), np.roll(spatial_anomaly[39,:,45,:].data, 72, axis=1), brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Zonal Wind Spatial Anomaly at Equator, t = 40 months')
    plt.xlabel('Longitude [degrees]')
    plt.ylabel('Height [m]')
    plt.xticks((-72,-52,-32,-12,0,12,32,52,72),('180W','140W','100W','60W','0','60E','100E','140E','180E'))
    plt.colorbar(pad=0.1)
    plt.show()
    
    plt.contourf(np.arange(-longitudes,longitudes), np.array(heights), np.roll(spatial_anomaly[47,:,45,:].data, 72, axis=1), brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Zonal Wind Spatial Anomaly at Equator, t = 48 months')
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-52,-32,-12,0,12,32,52,72),('180W','140W','100W','60W','0','60E','100E','140E','180E'))
    plt.ylabel('Height [m]')
    plt.colorbar(pad=0.1)
    plt.show()  
 

def qbo_period(cubes, periodicity=False):
    
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube.copy()    
    
    height = x_wind.shape[1]
    run_length = x_wind.shape[0]
    lats = x_wind.coord('latitude')
    longs = x_wind.coord('longitude')

    if lats.bounds == None:
        x_wind.coord('latitude').guess_bounds()
    if longs.bounds == None:
        x_wind.coord('longitude').guess_bounds()
        
    # strat = x_wind.extract(iris.Constraint(longitude=lambda v: 176 < v <= 184, latitude=lambda v: -4 <= v <= 4, model_level_number=47))
    # trop = x_wind.extract(iris.Constraint(longitude=lambda v: 176 < v <= 184, latitude=lambda v: -4 <= v <= 4, model_level_number=35))
    # high = x_wind.extract(iris.Constraint(longitude=lambda v: 176 < v <= 184, latitude=lambda v: -4 <= v <= 4, model_level_number=51))
    # low = x_wind.extract(iris.Constraint(longitude=lambda v: 176 < v <= 184, latitude=lambda v: -4 <= v <= 4, model_level_number=40))
        
    strat = x_wind.extract(iris.Constraint(longitude=lambda v: 355 < v <= 359 or 0 <= v <= 4, latitude=lambda v: -4 <= v <= 4, model_level_number=47))
    trop = x_wind.extract(iris.Constraint(longitude=lambda v: 355 < v <= 359 or 0 <= v <= 4, latitude=lambda v: -4 <= v <= 4, model_level_number=35))
    high = x_wind.extract(iris.Constraint(longitude=lambda v: 355 < v <= 359 or 0 <= v <= 4, latitude=lambda v: -4 <= v <= 4, model_level_number=51))
    low = x_wind.extract(iris.Constraint(longitude=lambda v: 355 < v <= 359 or 0 <= v <= 4, latitude=lambda v: -4 <= v <= 4, model_level_number=40))
    strat_grid = iris.analysis.cartography.area_weights(strat)
    trop_grid = iris.analysis.cartography.area_weights(trop)
    high_grid = iris.analysis.cartography.area_weights(high)
    low_grid = iris.analysis.cartography.area_weights(low)
    
    strat_mean = strat.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=strat_grid)
    trop_mean = trop.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=trop_grid)
    
    high_mean = high.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=high_grid)
    low_mean = low.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=low_grid)
    
    plt.plot(np.arange(0,run_length), high_mean.data, linestyle='--', color='r', label='50 km')
    plt.plot(np.arange(0,run_length), strat_mean.data, linestyle='-', color='r', label='40 km')
    plt.plot(np.arange(0,run_length), low_mean.data, linestyle='--', color='b', label='30 km')
    plt.plot(np.arange(0,run_length), trop_mean.data, linestyle='-', color='b', label='25 km')
    plt.title('Mean Substellar Zonal Wind')
    plt.xlabel('Time [months]') 
    plt.ylabel('Velocity [m s-1]')
    plt.legend()
    plt.show()
    
    if periodicity == True:
        
        data = np.array(high_mean.data)
        fft = sp.fftpack.fft(data)
        psd = np.abs(fft)**2
        freq = sp.fftpack.fftfreq(len(psd), 1./run_length)
        i = freq > 0
        
        data2 = np.array(strat_mean.data)
        fft2 = sp.fftpack.fft(data2)
        psd2 = np.abs(fft2)**2
        freq2 = sp.fftpack.fftfreq(len(psd2), 1./run_length)
        i2 = freq2 > 0
        
        periods_50km = np.round(run_length/freq[i], 2)
        fig, ax = plt.subplots(1,1,figsize=(8,4))
        ax.plot(periods_50km, psd[i])
        ax.set_xlabel('Period [months]')
        ax.set_ylabel('PSD')
        ax.set_title('Periodicity of mean substellar zonal wind at 50 km')
        psd_sorted = -np.sort(-psd[i])
        top_three = psd_sorted[:3]
        for i,j in zip(periods_50km,psd[i]):
            if j in top_three:
                ax.annotate('%s' %i, xy=(i,j), xytext=(2, 5), textcoords='offset points')
        # plt.annotate(str(maxpsd), location)
        plt.show()
        
        periods_40km = np.round(run_length/freq2[i2], 2)
        fig, ax = plt.subplots(1,1,figsize=(8,4))
        ax.plot(periods_40km, psd2[i2])
        ax.set_xlabel('Period [months]')
        ax.set_ylabel('PSD')
        ax.set_title('Periodicity of mean substellar zonal wind at 40 km')
        psd_sorted2 = -np.sort(-psd2[i2])
        top_three2 = psd_sorted2[:3]
        for i,j in zip(periods_40km,psd2[i2]):
            if j in top_three2:
                ax.annotate('%s' %i, xy=(i,j), xytext=(2, 5), textcoords='offset points')
        # plt.annotate(str(maxpsd), location)
        plt.show()


def wave_acceleration(cubes, time_slice=47):
    
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube.copy()
        if cube.standard_name == 'y_wind':
            y_wind = cube.copy()
        if cube.standard_name == 'upward_air_velocity':
            z_wind = cube[:,1:,:,:].copy()
        if cube.standard_name == 'air_pressure':
            pressure = cube.copy()
            
    longitudes = x_wind.shape[3]/2
    heights = pressure.coord('level_height').points
    x_zonal_mean = x_wind.collapsed('longitude', iris.analysis.MEAN)
    z_zonal_mean = z_wind.collapsed('longitude', iris.analysis.MEAN)
    x_anomaly = x_wind - x_zonal_mean
    z_anomaly = z_wind - z_zonal_mean
    print(x_anomaly.shape,z_anomaly.shape)
    
    result = x_wind.copy()
    result.data = x_anomaly.data*z_anomaly.data
    acceleration = iris.analysis.calculus.differentiate(result, 'level_height')
    
    plt.contourf(np.arange(-longitudes,longitudes), np.array(heights[2:]), np.roll(acceleration[time_slice,:,45,:].data, 72, axis=1), brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Wave-Induced Acceleration at Equator, t = %s months' %(time_slice+1))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-52,-32,-12,0,12,32,52,72),('180W','140W','100W','60W','0','60E','100E','140E','180E'))
    plt.ylabel('Height [m]')
    plt.colorbar(pad=0.1)
    plt.show()  
        


# def ep_flux(cubes):
#     for cube in cubes:
#         if cube.standard_name == 'air_ potential_temperature':
#             theta = cube.copy()
#         if cube.standard_name == 'x_wind':
#             x_wind = cube.copy()
#         if cube.standard_name == 'y_wind':
#             y_wind = cube.copy()
#         if cube.standard_name == 'air_pressure':
#             pressure = cube.copy()
            
    
    
    