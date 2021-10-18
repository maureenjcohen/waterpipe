#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 15:24:49 2021

@author: Mo Cohen

Functions relating to QBO-like phenomena

"""

import os, iris, cartopy, windspharm, pywt
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
from matplotlib import ticker
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

def get_heights(cubes):
    
    """ Extract altitude information from air pressure cube
        Can be used with matplotlib contourf when Iris isn't playing nice """ 
    
    for cube in cubes:
        if cube.standard_name == 'air_pressure':
            pressure = cube.copy()
            
    heights = pressure.coords()[1].points
    
    return heights


def plot_hovmoellerx(cubes, radius=7160000, time='days'):
    
    """ Plot Hovmoeller diagrams of zonal wind 
        Arguments: CubeList, time (can be '6-hours' or 'days')
        Outputs: 1) Hovmoeller plot of mean equatorial wind between -10 and 10 lat, axes time vs. height
                 2) Hovmoeller plot of zonal mean wind at 40 km, axes time vs. latitude (top view) """
    
    for cube in cubes:
        if cube.standard_name == 'eastward_wind' or cube.standard_name == 'x_wind':
            x_wind = cube.copy()
            
    x_wind.coord('latitude').coord_system = GeogCS(radius)
    x_wind.coord('longitude').coord_system = GeogCS(radius)
    
    run_length, lats, longs = x_wind.shape[0], x_wind.coord('latitude'), x_wind.coord('longitude')
    
    for coord in x_wind.coords():
        if coord.long_name == 'Hybrid height':
            heights = np.round(x_wind.coord('Hybrid height').points*1e-03,0)
        elif coord.long_name == 'level_height':
            heights = np.round(x_wind.coord('level_height').points*1e-03,0)
        
    if lats.bounds == None:
        x_wind.coord('latitude').guess_bounds()
    if longs.bounds == None:
        x_wind.coord('longitude').guess_bounds()
        
    dayside, nightside = x_wind.intersection(longitude=(-90,89), latitude=(-10,10)), x_wind.intersection(longitude=(90,270), latitude=(-10,10))
    dayside_full, nightside_full = x_wind.intersection(longitude=(-90,89)), x_wind.intersection(longitude=(90,270))
    day_grid, night_grid = iris.analysis.cartography.area_weights(dayside), iris.analysis.cartography.area_weights(nightside)
    
    dayside_mean = dayside.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=day_grid)
    nightside_mean = nightside.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=night_grid)
    
    dayside_zonal_mean = dayside_full.collapsed('longitude', iris.analysis.MEAN)
    nightside_zonal_mean = nightside_full.collapsed('longitude', iris.analysis.MEAN)    

    for cube in (dayside_mean, nightside_mean): 
        if cube == dayside_mean:
            side = 'Dayside'
        else:
            side = 'Nightside'
        plt.contourf(np.arange(0,run_length), np.array(heights), cube.data.T, brewer_redblu.N, cmap=brewer_redblu, norm=TwoSlopeNorm(0))
        plt.title('%s Mean Zonal Equatorial Wind' %(side))
        plt.xlabel('Time [%s]' %(time))
        plt.ylabel('Height [km]')
        cbar = plt.colorbar(pad=0.1)
        cbar.ax.set_title('m/s')
        # plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/laso/epsfigs/hovmoeller_%s.eps' %(side), format='eps')   
        plt.show()
    
    for cube in (dayside_zonal_mean, nightside_zonal_mean):
        if cube == dayside_zonal_mean:
            side = 'Dayside'
        else:
            side = 'Nightside'  
        plt.contourf(np.arange(0,run_length), np.arange(-45,45), cube[:,47,:].data.T, brewer_redblu.N, cmap=brewer_redblu, norm=TwoSlopeNorm(0))
        plt.title('%s Mean Zonal Wind at 41 km' %(side))
        plt.xlabel('Time [%s]' %(time))
        plt.ylabel('Latitude')
        plt.yticks((-45,0,45),('90S', '0', '90N'))
        mbar = plt.colorbar(pad=0.1)
        mbar.ax.set_title('m/s')
        # plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/laso/epsfigs/topview_%s.eps' %(side), format='eps')   
        plt.show()
    


def plot_temp_anomaly(cubes, period=(0,220), lat=45, levels=(35,44,47,53)):
    
    """ Plot temperature anomaly, temperature minus the local time-averaged temperature
        Default latitude is the equator
        Arguments: CubeList, latitude, and atmospheric level"""
    
    for cube in cubes:
        if cube.standard_name == 'air_potential_temperature':
            theta = cube[period[0]:period[1],:,:,:].copy()
        if cube.standard_name == 'air_pressure':
            pressure = cube[period[0]:period[1],:,:,:].copy()
        
    run_length, longitudes = theta.shape[0], theta.shape[3]/2
    
    for coord in theta.coords():
        if coord.long_name == 'Hybrid height':
            heights = np.round(theta.coord('Hybrid height').points*1e-03,0)
        elif coord.long_name == 'level_height':
            heights = np.round(theta.coord('level_height').points*1e-03,0)
    
    p0 = iris.coords.AuxCoord(100000.0, long_name='reference_pressure', units='Pa')
    p0.convert_units(pressure.units)
    temperature = theta*((pressure/p0)**(287.05/1005)) # R and cp in J/kgK for 300K
    
    temp_time_mean = temperature.collapsed('time', iris.analysis.MEAN)
    anomaly = temperature - temp_time_mean
    
    for level in levels:
            
        plt.figure(figsize=(10,5))
        plt.contourf(np.arange(-longitudes,longitudes), np.arange(0, run_length), np.roll(anomaly[:,level,lat,:].data, 72, axis=1),brewer_redblu.N, cmap=brewer_redblu, norm=TwoSlopeNorm(0))
        plt.title('Temperature Anomaly at Equator, h=%s km' %(heights[level]))
        plt.xlabel('Longitude [degrees]')
        plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
        plt.ylabel('Time [months]')
        cbar = plt.colorbar(pad=0.1)
        cbar.ax.set_title('K')
        # plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/laso/epsfigs/abs_temp_anomaly_%s.eps' %(heights[level]), format='eps')
#        plt.show()  
    
 

def plot_eddy_temp(cubes, lat=45, start=0, end=160, levels=(38,44,47,53)):
    
    """ Plot the eddy (varying, non-mean) component of air temperature. 
        The eddy temperature is the temperature minus the zonal mean temperature.
        Arguments: CubeList, start and end time for period, levels of atmosphere
        Outputs: Plots of eddy temperature, axes longitude vs. time"""
    
    for cube in cubes:
        if cube.standard_name == 'air_potential_temperature':
            theta = cube[start:end,:,lat,:].copy()
        if cube.standard_name == 'air_pressure':
            pressure = cube[start:end,:,lat,:].copy()

    run_length, heights, longitudes = pressure.shape[0], np.round(pressure.coords()[1].points*1e-03, 0), pressure.shape[2]
    
    p0 = iris.coords.AuxCoord(100000.0, long_name='reference_pressure', units='Pa')
    p0.convert_units(pressure.units)
    temperature = theta*((pressure/p0)**(287.05/1005)) # R and cp in J/kgK for 300K

    zonal_mean = temperature.collapsed('longitude', iris.analysis.MEAN)
    temp_prime = temperature - zonal_mean 
    
    daily = np.mean(temp_prime.data.reshape(-1,4,61,144),axis=1)
                    
    plt.figure(figsize=(10,5))
    plt.contourf(np.arange(-longitudes/2,longitudes/2), np.array(heights[47:53]), np.roll(temp_prime[-1,47:53,:].data, 72, axis=1), brewer_redblu.N, cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    plt.title('$T^{\prime}$ at Equator [K], day %s' %(end/4))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.ylabel('Height [km]')
    plt.colorbar(pad=0.1)
    plt.show() 
    
    for level in levels:       
                
        plt.figure(figsize=(10,5))
        plt.contourf(np.arange(-longitudes/2,longitudes/2), np.arange(0, run_length), np.roll(temp_prime[:,level,:].data, 72, axis=1), np.linspace(-10,10,100), cmap=brewer_redblu, norm=TwoSlopeNorm(0))
        plt.title('$T^{\prime}$ at Equator [K], h=%s km, %s to %s days' %(heights[level], start/4, end/4))
        plt.xlabel('Longitude [degrees]')
        plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
        plt.ylabel('Time [6-hours]')
        plt.colorbar(pad=0.1)
        plt.show() 


    
def plot_water_anomaly(cubes):
    
    for cube in cubes: 
        if cube.standard_name == 'specific_humidity':
            water = cube.copy()
            
    run_length = water.shape[0]
    longitudes = water.shape[3]/2
    
    water = water.collapsed('model_level_number', iris.analysis.SUM)
    
    time_mean = water.collapsed('time', iris.analysis.MEAN)
    anomaly = water - time_mean
    anomaly /= time_mean    
    
    plt.figure(figsize=(10,5))
    plt.contourf(np.arange(-longitudes,longitudes), np.arange(0, run_length), np.roll(anomaly[:,45,:].data, 72, axis=1), levels=np.linspace(-2,2),  cmap=brewer_redblu)
    plt.title('Water Column Anomaly at Equator [fractional]')
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.ylabel('Time [days]')
    plt.colorbar(pad=0.1)
    plt.show() 

        
def plot_wind_anomaly(cubes, lat=45, times=(0,1500,1780,2060)):
    
    for cube in cubes:
        if cube.standard_name == 'upward_air_velocity' or cube.standard_name == 'x_wind':
            x_wind = cube[:,:,lat,:].copy()
    
    longitudes = x_wind.shape[2]/2
    
    for coord in x_wind.coords():
        if coord.long_name == 'Hybrid height':
            heights = np.round(x_wind.coord('Hybrid height').points*1e-03,0)
            x_wind_time_mean = x_wind.collapsed('t', iris.analysis.MEAN)
            time_unit = 'days'
        elif coord.long_name == 'level_height':
            heights = np.round(x_wind.coord('level_height').points*1e-03,0)
            x_wind_time_mean = x_wind.collapsed('time', iris.analysis.MEAN)
            time_unit = 'months'
    
    anomaly = x_wind - x_wind_time_mean
    
    zonal_mean = anomaly.collapsed('longitude', iris.analysis.MEAN)
    zonal_anomaly = anomaly - zonal_mean
    
    for time in times:
      
        plt.figure(figsize=(10,5))    
        plt.contourf(np.arange(-longitudes,longitudes), np.array(heights), np.roll(zonal_anomaly[time,:,:].data, 72, axis=1), brewer_redblu.N, cmap=brewer_redblu, norm=TwoSlopeNorm(0))
        plt.title('Zonal Wind Anomaly at Equator [m s-1], t = %s %s' %(time/4, time_unit))
        plt.xlabel('Longitude [degrees]')
        plt.ylabel('Height [m]')
        plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
        plt.colorbar(pad=0.1)
        plt.show()

    

def plot_u(cubes, start=640, end=800,lat=45, level=47):
    
    """ Plots u vs longitude and time
        Arguments: CubeList, latitude, and atmospheric level"""
    
    for cube in cubes:
        if cube.standard_name == 'eastward_wind' or cube.standard_name == 'x_wind':
            x_wind = cube[start:end,:,lat,:].copy()
    
    run_length, longitudes = x_wind.shape[0], x_wind.shape[2]/2
    
    for coord in x_wind.coords():
        if coord.long_name == 'Hybrid height':
            heights = np.round(x_wind.coord('Hybrid height').points*1e-03,0)
        elif coord.long_name == 'level_height':
            heights = np.round(x_wind.coord('level_height').points*1e-03,0)
    
    plt.figure(figsize=(10,5))
    plt.contourf(np.arange(-longitudes,longitudes), np.arange(0, run_length), np.roll(x_wind[:,level,:].data, 72, axis=1),brewer_redblu.N, cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    plt.title('Eastward Wind at Equator, h=%s km, t=%s to %s days' %(heights[level], start, end))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.ylabel('Time [days]')
    cbar = plt.colorbar(pad=0.1)
    cbar.ax.set_title('m/s')
    plt.show()  
    
    print(np.min(x_wind[:,level,:].data), np.max(x_wind[:,level,:].data))
    plt.figure(figsize=(10,5))
    plt.contourf(np.arange(-longitudes,longitudes), np.arange(0, run_length), x_wind[:,level,:].data,brewer_redblu.N, cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    plt.title('Eastward Wind at Equator, h=%s km, t=%s to %s days' %(heights[level], start/4, end/4))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('0','30E','60E','90E','120E','150E','180E/W','150W','120W','90W','60W','30W','0'))
    plt.ylabel('Time [days]')
    mbar = plt.colorbar(pad=0.1) 
    mbar.locator = ticker.AutoLocator()
    mbar.update_ticks()
    # mbar.ax.set_yticklabels(['140','120','100','80','60','40','20','0','-20','-40','-60','-80', '-100'])
    mbar.ax.set_title('m/s')
    # plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/laso/epsfigs/zonalwind_%s_to_%s_ticks.eps' %(start, end), format='eps')  
    plt.show() 
    
    
def u_prime(cubes, time_slice=1780, lat=45, level=47, bandpass=False):
    
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube[time_slice,:,lat,:].copy()
    
    heights = np.round(x_wind.coord('level_height').points*1e-03, 0)
    longitudes = x_wind.shape[1]/2
    
    dayside = x_wind.intersection(longitude=(-90,89))
    nightside = x_wind.intersection(longitude=(90,269))
    
    dayside_zonal_mean = dayside.collapsed('longitude',iris.analysis.MEAN)
    nightside_zonal_mean = nightside.collapsed('longitude',iris.analysis.MEAN)
    global_zonal_mean = x_wind.collapsed('longitude', iris.analysis.MEAN)

    dayside_spatial_anomaly = dayside - dayside_zonal_mean 
    nightside_spatial_anomaly = nightside - nightside_zonal_mean
    global_spatial_anomaly = x_wind - global_zonal_mean
    
    plt.figure(figsize=(10,5))
    plt.plot(np.arange(-longitudes, longitudes), np.roll(x_wind[level,:].data, 72, axis=0), color='k', label='U')
    plt.plot(np.arange(-longitudes, longitudes), np.roll(global_spatial_anomaly[level,:].data, 72, axis=0), color='r', label='U$^{\prime}$')
    plt.title('U at Equator, t=%s days, h=%s km' %(time_slice/4, heights[level]))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.ylabel('Wind speed [m s-1]')
    plt.axhline(y=global_zonal_mean[level].data, color='b', linestyle='--', label='$\overline{U}$')
    plt.legend()
    plt.show()
    
    plt.figure(figsize=(10,5))    
    plt.contourf(np.arange(-longitudes,longitudes), np.array(heights[35:58]), np.roll(global_spatial_anomaly[35:58].data, 72,axis=0), brewer_redblu.N, cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    plt.title('$U^{\prime}$ at Equator [m s-1], t = %s days' %(time_slice/4))
    plt.xlabel('Longitude [degrees]')
    plt.ylabel('Height [km]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.colorbar(pad=0.1)
    plt.show()
   
    plt.figure(figsize=(10,5))    
    plt.contourf(np.arange(-longitudes/2,longitudes/2), np.array(heights[35:58]), nightside_spatial_anomaly[35:58,:].data, brewer_redblu.N, cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    plt.title('Nightside $U^{\prime}$ at Equator [m s-1], t = %s days' %(time_slice/4))
    plt.xlabel('Longitude [degrees]')
    plt.ylabel('Height [km]')
    plt.xticks((-36,-18,0,18,36), ('90E', '135E', '180', '135W','90W'))
    plt.colorbar(pad=0.1)
    plt.show()
    
    plt.figure(figsize=(10,5))    
    plt.contourf(np.arange(-longitudes/2,longitudes/2), np.array(heights[35:58]), dayside_spatial_anomaly[35:58,:].data, brewer_redblu.N, cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    plt.title('Dayside $U^{\prime}$ at Equator [m s-1], t = %s days' %(time_slice/4))
    plt.xlabel('Longitude [degrees]')
    plt.ylabel('Height [km]')
    plt.xticks((-36,-18,0,18,36), ('90W', '45W', '0', '45E','90E'))
    plt.colorbar(pad=0.1)
    plt.show()
    
        
    if bandpass==True:
        
        data = np.array(global_spatial_anomaly[level,45,:].data)
        fft = sp.fftpack.rfft(data)
        psd = np.abs(fft)**2
        freq = sp.fftpack.rfftfreq(len(psd), 1./144)
        psd[freq<=3] = 0
        # i = freq > 0
        
        fig, ax = plt.subplots(1,1,figsize=(8,4))
        ax.plot(freq, psd)
        ax.set_xlabel('Wavenumber')
        ax.set_ylabel('PSD')
        ax.set_title('Wavenumbers with k<=3 removed')
        plt.show()
        
        back = sp.fftpack.irfft(psd)
        plt.figure(figsize=(10,5))    
        plt.plot(np.arange(-longitudes, longitudes), back)
        plt.title('Residual U$^{\prime}$')
        plt.xlabel('Longitude [degrees]')
        plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
        plt.ylabel('Wind speed [m s-1]')
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
        
    # strat = x_wind.extract(iris.Constraint(longitude=lambda v: 178 < v <= 182, latitude=lambda v: -2 <= v <= 2, model_level_number=47))
    # trop = x_wind.extract(iris.Constraint(longitude=lambda v: 178 < v <= 182, latitude=lambda v: -2 <= v <= 2, model_level_number=35))
    # high = x_wind.extract(iris.Constraint(longitude=lambda v: 178 < v <= 182, latitude=lambda v: -2 <= v <= 2, model_level_number=51))
    # low = x_wind.extract(iris.Constraint(longitude=lambda v: 178 < v <= 182, latitude=lambda v: -2 <= v <= 2, model_level_number=40))
        
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


def wave_acceleration(cubes, time_slice=2000):
    
    for cube in cubes:
        if cube.standard_name == 'eastward_wind':
            x_wind = cube[time_slice,:,:,:].copy()
        if cube.standard_name == 'northward_wind':
            y_wind = cube[time_slice,:,:,:].copy()
        if cube.standard_name == 'upward_air_velocity':
            z_wind = cube[time_slice,:,:,:].copy()
        if cube.standard_name == 'air_pressure':
            pressure = cube[time_slice,:,:,:].copy()
            
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
    z_wind = z_wind.regrid(x_wind, iris.analysis.Linear())
    pressure = pressure.regrid(x_wind, iris.analysis.Linear())
    
    heights = pressure.coord('Hybrid height').points
    
    height = [('Hybrid height', pressure.coord('Hybrid height').points)]
    x_wind, y_wind = x_wind.interpolate(height, iris.analysis.Linear()), y_wind.interpolate(height, iris.analysis.Linear())
    
    x_wind = x_wind.extract(iris.Constraint(latitude=45))
    z_wind = z_wind.extract(iris.Constraint(latitude=45))
    pressure = pressure.extract(iris.Constraint(latitude=45))

    longitudes = x_wind.shape[1]/2
    x_zonal_mean = x_wind.collapsed('longitude', iris.analysis.MEAN)
    z_zonal_mean = z_wind.collapsed('longitude', iris.analysis.MEAN)
    
    x_anomaly = x_wind - x_zonal_mean
    z_anomaly = z_wind - z_zonal_mean
    
    result = x_wind.copy()
    result.data = x_anomaly.data*z_anomaly.data
    array = result.data
    h = np.log(pressure.data)

#    acceleration = iris.analysis.calculus.differentiate(result, 'Hybrid height')
    
    acceleration = array.copy()
    acceleration[0,:] = (array[1,:]-array[0,:])/(h[1,:]-h[0,:])
    acceleration[-1,:] = (array[-1,:]-array[-2,:])/(h[-1,:]-h[-2,:])
    acceleration[1:-1,:] = (array[2:,:]-array[0:-2,:])/(h[2:,:]-h[0:-2,:])

    plt.figure(figsize=(10,5))
    plt.contourf(np.arange(-longitudes,longitudes), np.array(heights[:57]), np.roll(acceleration[:57,:], 72, axis=1), levels=np.linspace(-5,5,50), cmap=brewer_redblu)
    plt.title('Wave-Induced Acceleration at Equator [m s-2], t = %s days' %(time_slice/4))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.ylabel('Height [m]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    acc_cube = x_wind.copy()
    acc_cube.data = acceleration
    substellar = acc_cube.intersection(longitude=(-90,90))
    plt.figure(figsize=(5,5))
    plt.contourf(np.arange(-36,37), np.array(heights[:58]), substellar[:58,:].data, np.linspace(-3,3,24), cmap=brewer_redblu)
    plt.title('Wave-Induced Acceleration at Equator [m s-2], t = %s days' %(time_slice/4))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-36,-24,-12,0,12,24,36),('90W','60W','30W','0','30E','60E','90E'))
    plt.ylabel('Height [m]')
    plt.colorbar(pad=0.1)
    plt.show()  
    

def equatorial_wind(cubes, time_slice=1780, low=47, high=53, bandpass=False):
    
    for cube in cubes:
        if cube.standard_name == 'eastward_wind':
            x_wind = cube.copy()
    
#    acceleration = iris.analysis.calculus.differentiate(x_wind, 'longitude')
    longitudes = x_wind.shape[3]/2     
    heights = np.round(x_wind.coord('Hybrid height').points*1e-03,0)
    trop = np.roll(x_wind[time_slice,37,45,:].data, 72)
    westerly = np.roll(x_wind[time_slice,low,45,:].data, 72)
    easterly = np.roll(x_wind[time_slice,high,45,:].data, 72)
    # westerly_acc = np.roll(acceleration[time_slice,low,45,:].data, 72)
    # easterly_acc = np.roll(acceleration[time_slice,high,45,:].data, 72)
    
    plt.figure(figsize=(10,5))
    plt.plot(np.arange(-longitudes, longitudes), trop)
    plt.title('Tropospheric Zonal Jet at Equator, t=%s months, h=%s km' %(time_slice+1, heights[37]))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.axvline(x=-36, color='r', linestyle='--')
    plt.axvline(x=36, color='r', linestyle='--')
    plt.ylabel('Wind speed [m s-1]')
    plt.show()

    plt.figure(figsize=(10,5))
    plt.plot(np.arange(-longitudes, longitudes), westerly)
    plt.title('Zonal Wind at Equator, t=%s days, h=%s km' %(time_slice+1, heights[low]))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.axvline(x=-36, color='r', linestyle='--')
    plt.axvline(x=36, color='r', linestyle='--')
    plt.ylabel('Wind speed [m s-1]')
    plt.show()
    
    plt.figure(figsize=(10,5))
    plt.plot(np.arange(-longitudes, longitudes), easterly)
#    plt.plot(np.arange(-longitudes, longitudes), np.zeros_like(easterly), color='g', linestyle='--')
    plt.title('Zonal Wind at Equator, t=%s days, h=%s km' %(time_slice+1, heights[high]))
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.axvline(x=-36, color='r', linestyle='--')
    plt.axvline(x=36, color='r', linestyle='--')    
    plt.ylabel('Wind speed [m s-1]')
    plt.show()
    
    if bandpass==True:
        
        data = np.array(trop)
        fft = sp.fftpack.fft(data)
        psd = np.abs(fft)**2
        freq = sp.fftpack.fftfreq(len(psd), 1./144)
        psd[psd>1e5] = 0
        # i = freq > 0
        
        fig, ax = plt.subplots(1,1,figsize=(8,4))
        ax.plot(freq, psd)
        ax.set_xlabel('Wavenumber')
        ax.set_ylabel('PSD')
        ax.set_title('Frequencies')
        plt.show()
        
        back = sp.fftpack.ifft(psd)
        plt.plot(np.arange(-longitudes, longitudes), back)
        plt.show()
        
    

def equatorial_pressure(cubes, time_slice=10, low=47, high=53):
    
    for cube in cubes:
        if cube.standard_name == 'air_pressure':
            pressure = cube.copy()
    
    longitudes = pressure.shape[3]/2     
    heights = np.round(pressure.coord('level_height').points*1e-03,0)
    westerly = np.roll(pressure[time_slice,low,45,:].data, 72)
    easterly = np.roll(pressure[time_slice,high,45,:].data, 72)
    
    plt.figure(figsize=(10,5))
    plt.plot(np.arange(-longitudes, longitudes), np.roll(pressure[time_slice,38,45,:].data, 72))
    plt.title('Pressure at Equator, t=%s months, h=%s km' %(time_slice+1, heights[38]))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.axvline(x=-36, color='r', linestyle='--')
    plt.axvline(x=36, color='r', linestyle='--')
    plt.ylabel('Pressure [Pa]')
    plt.show()
    
    iplt.contourf(pressure[time_slice,38,:,:], brewer_red.N, cmap=brewer_red)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Pressure at h=%s km, t=%s months' %(heights[38], time_slice+1), y=1.20)
    plt.colorbar(pad=0.1)
    plt.show()

    plt.figure(figsize=(10,5))
    plt.plot(np.arange(-longitudes, longitudes), westerly)
    plt.title('Pressure at Equator, t=%s months, h=%s km' %(time_slice+1, heights[low]))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.axvline(x=-36, color='r', linestyle='--')
    plt.axvline(x=36, color='r', linestyle='--')
    plt.ylabel('Pressure [Pa]')
    plt.show()
    
    iplt.contourf(pressure[time_slice,low,:,:], brewer_red.N, cmap=brewer_red)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Pressure at h=%s km, t=%s months' %(heights[low], time_slice+1), y=1.20)
    plt.colorbar(pad=0.1)
    plt.show()
    
    plt.figure(figsize=(10,5))
    plt.plot(np.arange(-longitudes, longitudes), easterly)
    plt.title('Pressure at Equator, t=%s months, h=%s km' %(time_slice+1, heights[high]))
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.axvline(x=-36, color='r', linestyle='--')
    plt.axvline(x=36, color='r', linestyle='--')    
    plt.ylabel('Pressure [Pa]')
    plt.show()
    
    iplt.contourf(pressure[time_slice,high,:,:], brewer_red.N, cmap=brewer_red)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Pressure at h=%s km, t=%s months' %(heights[high], time_slice+1), y=1.20)
    plt.colorbar(pad=0.1)
    plt.show()
    

def equatorial_temp(cubes, time_slice=10, low=47, high=53):
    
    for cube in cubes:
        if cube.standard_name == 'air_pressure':
            pressure = cube.copy()
        if cube.standard_name == 'air_potential_temperature':
            theta = cube.copy()
            
    p0 = iris.coords.AuxCoord(100000.0, long_name='reference_pressure', units='Pa')
    p0.convert_units(pressure.units)
    absolute_temp = theta*((pressure/p0)**(287.05/1005)) # R and cp in J/kgK for 300K
    
    longitudes = pressure.shape[3]/2     
    heights = np.round(pressure.coord('Hybrid height').points*1e-03,0)
    westerly = np.roll(absolute_temp[time_slice,low,45,:].data, 72)
    easterly = np.roll(absolute_temp[time_slice,high,45,:].data, 72)
    
    plt.figure(figsize=(10,5))
    plt.plot(np.arange(-longitudes, longitudes), np.roll(absolute_temp[time_slice,38,45,:].data, 72))
    plt.title('Temperature at Equator, t=%s months, h=%s km' %(time_slice+1, heights[38]))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.axvline(x=-36, color='r', linestyle='--')
    plt.axvline(x=36, color='r', linestyle='--')
    plt.ylabel('Temperature [K]]')
    plt.show()

    plt.figure(figsize=(10,5))
    plt.plot(np.arange(-longitudes, longitudes), westerly)
    plt.title('Temperature at Equator, t=%s months, h=%s km' %(time_slice+1, heights[low]))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.axvline(x=-36, color='r', linestyle='--')
    plt.axvline(x=36, color='r', linestyle='--')
    plt.ylabel('Temperature [K]]')
    plt.show()
    
    plt.figure(figsize=(10,5))
    plt.plot(np.arange(-longitudes, longitudes), easterly)
    plt.title('Temperature at Equator, t=%s months, h=%s km' %(time_slice+1, heights[high]))
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.axvline(x=-36, color='r', linestyle='--')
    plt.axvline(x=36, color='r', linestyle='--')    
    plt.ylabel('Temperature [K]')
    plt.show()


def ep_zderivative(array,h):
    work = array.copy()
    work[:,0,...] = (array[:,1,...]-array[:,0,...])/(h[:,1,...]-h[:,0,...])
    work[:,-1,...] = (array[:,-1,...]-array[:,-2,...])/(h[:,-1,...]-h[:,-2,...])
    work[:,1:-1,...] = (array[:,2:,...]-array[:,0:-2,...])/(h[:,2:,...]-h[:,0:-2,...])
    
    return work

def ep_latderivative(array,l):
    
    work = array.copy()
    work[:,:,0,...] = (array[:,:,1,...]-array[:,:,0,...])/(l[1,...]-l[0,...])
    work[:,:,-1,...] = (array[:,:,-1,...]-array[:,:,-2,...])/(l[-1,...]-l[-2,...])
    work[:,:,1:-1,...] = (array[:,:,2:,...]-array[:,:,0:-2,...])/(l[2:,...]-l[0:-2,...])
    
    return work

def ep_flux(cubes, omega=0.64617667e-05, level=47, long_slice=(120,132), start=1780, end=1784):
    
    for cube in cubes:
        if cube.standard_name == 'air_potential_temperature':
            theta = cube[start:end,:,:,long_slice[0]:long_slice[-1]].copy()
        if cube.standard_name == 'eastward_wind':
            x_wind = cube[start:end,:,:,long_slice[0]:long_slice[-1]].copy()
        if cube.standard_name == 'northward_wind':
            y_wind = cube[start:end,:,:,long_slice[0]:long_slice[-1]].copy()
        if cube.standard_name == 'air_pressure':
            pressure = cube[start:end,:,:,long_slice[0]:long_slice[-1]].copy()
    
    heights, longitudes, longs = np.round(pressure.coord('Hybrid height').points*1e-03,0), pressure.coord('longitude').points, pressure.shape[3]/2
    
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
    logpressure_zonal_mean = np.mean(logpressure, axis=3)
    pressure_zonal_mean = np.mean(pressure_data, axis=3)
    
    THETAp = ep_zderivative(theta_zonal_mean_data, logpressure_zonal_mean)
    THETAp /= logpressure_zonal_mean
    
    x_zonal_mean = x_wind.collapsed('longitude', iris.analysis.MEAN)
    y_zonal_mean = y_wind.collapsed('longitude', iris.analysis.MEAN)
    x_anomaly = x_wind - x_zonal_mean
    y_anomaly = y_wind - y_zonal_mean 
    theta_anomaly = theta - theta_zonal_mean
    
    x_anomaly_data = x_anomaly.data
    y_anomaly_data = y_anomaly.data
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
    
    EP_phi = -XY_anomaly_zonal_mean*latfac
    EP_up = (coriolis*rcosphi*YTH_anomaly_zonal_mean)/THETAp
    
    EP_phi_div = ep_latderivative(EP_phi, rsinphi)
    EP_up_div = ep_zderivative(EP_up, logpressure_zonal_mean)
    EP_div = EP_phi_div + EP_up_div
    EP_div_mean = np.mean(EP_div, axis=0)  
    x_mean = x_zonal_mean.collapsed('t', iris.analysis.MEAN)
        
    x_axis = pressure.coord('latitude').points
    y_axis = pressure.coord('Hybrid height').points
    plt.contourf(x_axis, y_axis[35:59], EP_div_mean[35:59,:], brewer_redblu.N, cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    plt.colorbar(pad=0.1)
    CS = plt.contour(x_axis, y_axis[35:59], x_mean[35:59].data, colors='black', linewidths=0.5)
    plt.title('EP flux divergence for day=%s, long=%s to %s' %(start/4, longitudes[0], longitudes[-1]))
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


def momentum_transport(cubes, start=1500, end=1780):
    for cube in cubes:
        if cube.standard_name == 'eastward_wind':
            x_wind = cube[start:end,:,:,:].copy() 
        if cube.standard_name == 'northward_wind':
            y_wind = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'upward_air_velocity':
            z_wind = cube[start:end,1:,:,:].copy()

    heights = [('Hybrid height', x_wind.coord('Hybrid height').points)]
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
    z_wind = z_wind.regrid(x_wind, iris.analysis.Linear())
    z_wind = z_wind.interpolate(heights, iris.analysis.Linear())
            
    time_x = x_wind.collapsed('t', iris.analysis.MEAN)
    time_y = y_wind.collapsed('t', iris.analysis.MEAN)
    time_z = z_wind.collapsed('t', iris.analysis.MEAN)
    product_xy, product_xz = x_wind*y_wind, x_wind*z_wind
    time_product_xy, time_product_xz = product_xy.collapsed('t', iris.analysis.MEAN), product_xz.collapsed('t', iris.analysis.MEAN)
    
    # zonal_x = x_wind.collapsed('longitude', iris.analysis.MEAN)
    # zonal_y = y_wind.collapsed('longitude', iris.analysis.MEAN)
    # eddy_x = x_wind - zonal_x
    # eddy_y = y_wind - zonal_y
    # eddy_product = eddy_x*eddy_y
    # time_eddy_product = eddy_product.collapsed('t', iris.analysis.MEAN)
    
    zonal_uv = time_product_xy.collapsed('longitude', iris.analysis.MEAN)
    zonal_uw = time_product_xz.collapsed('longitude', iris.analysis.MEAN)
    zonal_u = time_x.collapsed('longitude', iris.analysis.MEAN)
    zonal_v = time_y.collapsed('longitude', iris.analysis.MEAN)
    zonal_w = time_z.collapsed('longitude', iris.analysis.MEAN)
    
    momentum_uv = zonal_uv - zonal_u*zonal_v
    momentum_uw = zonal_uw - zonal_u*zonal_w
    
    iplt.contourf(momentum_uv[:55], np.linspace(-240,240,240), cmap=brewer_redblu)
    plt.title('Horizontal Eddy Momentum Flux [m2 s-2]')
    plt.xlabel('Latitude [degrees]')
    plt.ylabel('Height [m]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    iplt.contourf(momentum_uw[:55], np.linspace(-0.25,0.25,20), cmap=brewer_redblu)
    plt.title('Vertical Eddy Momentum Flux [m2 s-2]')
    plt.xlabel('Latitude [degrees]')
    plt.ylabel('Height [m]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    
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
    
    zonal_mean_u = x_wind.collapsed('longitude', iris.analysis.MEAN)
    u_prime = x_wind - zonal_mean_u
    zonal_mean_w = z_wind.collapsed('longitude', iris.analysis.MEAN)
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
        
    
def plot_uv(cubes, time_slice=1500, levels=(37,44,47,53)):
    
    for cube in cubes:
        if cube.standard_name == 'eastward_wind' or cube.standard_name == 'x_wind':
            x_wind = cube[time_slice,:,:,:].copy() 
        if cube.standard_name == 'northward_wind' or cube.standard_name == 'y_wind':
            y_wind = cube[time_slice,:,:,:].copy()
    
    longitudes, latitudes = x_wind.shape[2]/2, x_wind.shape[1]/2
            
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
    
    for coord in x_wind.coords():
        if coord.long_name == 'Hybrid height':
            heights = np.round(x_wind.coord('Hybrid height').points*1e-03,0)
        elif coord.long_name == 'level_height':
            heights = np.round(x_wind.coord('level_height').points*1e-03,0)            
    
    zonal_mean_u = x_wind.collapsed('longitude', iris.analysis.MEAN)
    u_prime = x_wind - zonal_mean_u
    zonal_mean_v = y_wind.collapsed('longitude', iris.analysis.MEAN)
    v_prime = y_wind - zonal_mean_v
    cross_term = u_prime*v_prime
        
    for level in levels:
        
        plt.figure(figsize=(10,5))
        plt.contourf(np.arange(-longitudes,longitudes), np.arange(-latitudes,latitudes), np.roll(cross_term[level,:,:].data, 72, axis=1),brewer_redblu.N, cmap=brewer_redblu, norm=TwoSlopeNorm(0))
        plt.title('$U^{\prime} \cdot V^{\prime}$ [m2 s-2], h=%s km, day=%s' %(heights[level], time_slice/4))
        plt.xlabel('Longitude [degrees]')
        plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
        plt.ylabel('Latitude [degrees]')
        plt.colorbar(pad=0.1)
        plt.show() 
        
        

def wavelets(cubes, plot=False):
    
    for cube in cubes:
        if cube.standard_name == 'eastward_wind':
            x_wind = cube.copy()    
    
    height = x_wind.shape[1]
    heights = x_wind.coord('Hybrid height').points
    run_length = x_wind.shape[0]
    lats = x_wind.coord('latitude')
    longs = x_wind.coord('longitude')

    if lats.bounds == None:
        x_wind.coord('latitude').guess_bounds()
    if longs.bounds == None:
        x_wind.coord('longitude').guess_bounds()
        
    strat = x_wind.extract(iris.Constraint(longitude=lambda v: 178 < v <= 182, latitude=lambda v: -2 <= v <= 2, model_level_number=47))
    trop = x_wind.extract(iris.Constraint(longitude=lambda v: 178 < v <= 182, latitude=lambda v: -2 <= v <= 2, model_level_number=35))
    high = x_wind.extract(iris.Constraint(longitude=lambda v: 178 < v <= 182, latitude=lambda v: -2 <= v <= 2, model_level_number=51))
    low = x_wind.extract(iris.Constraint(longitude=lambda v: 178 < v <= 182, latitude=lambda v: -2 <= v <= 2, model_level_number=40))
        
    # strat = x_wind.extract(iris.Constraint(longitude=lambda v: 355 < v <= 359 or 0 <= v <= 4, latitude=lambda v: -4 <= v <= 4))
    # print(strat.shape)
    # trop = x_wind.extract(iris.Constraint(longitude=lambda v: 355 < v <= 359 or 0 <= v <= 4, latitude=lambda v: -4 <= v <= 4, Hybrid height=heights[35]))
    # high = x_wind.extract(iris.Constraint(longitude=lambda v: 355 < v <= 359 or 0 <= v <= 4, latitude=lambda v: -4 <= v <= 4, Hybrid height=heights[51]))
    # low = x_wind.extract(iris.Constraint(longitude=lambda v: 355 < v <= 359 or 0 <= v <= 4, latitude=lambda v: -4 <= v <= 4, Hybrid height=heights[40]))
    strat_grid = iris.analysis.cartography.area_weights(strat)
    trop_grid = iris.analysis.cartography.area_weights(trop)
    high_grid = iris.analysis.cartography.area_weights(high)
    low_grid = iris.analysis.cartography.area_weights(low)
    
    strat_mean = strat.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=strat_grid)
    trop_mean = trop.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=trop_grid)
    
    high_mean = high.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=high_grid)
    low_mean = low.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=low_grid)
    
    if plot==True: 
        plt.plot(np.arange(0,run_length), high_mean.data, linestyle='--', color='r', label='50 km')
        plt.plot(np.arange(0,run_length), strat_mean.data, linestyle='-', color='r', label='40 km')
        plt.plot(np.arange(0,run_length), low_mean.data, linestyle='--', color='b', label='30 km')
        plt.plot(np.arange(0,run_length), trop_mean.data, linestyle='-', color='b', label='25 km')
        plt.title('Mean Substellar Zonal Wind')
        plt.xlabel('Time [6-hours]') 
        plt.ylabel('Velocity [m s-1]')
        plt.legend()
        plt.show()
    
    data = strat_mean.data
    coeffs, freqs = pywt.cwt(data, np.arange(1,65), 'gaus1')
    plt.contourf(coeffs, cmap=brewer_redblu)
    plt.title('Wavelet transform of h = 40km')
    plt.xlabel('Time')
    plt.ylabel('Scale')
    plt.show()
    
    
def wk_filter(arg):
    
    """ 1-2-1 filter for smoothing Fourier coefficients in the Wheeler-Kiladis diagram
        Three neighbouring values are weighted 1*, 2*, 1*, added together, and divided by 4 """
        
    out = np.zeros_like(arg)
    out[0] = (3*arg[0] + arg[1])/4
    out[-1] = (arg[-2] + 3*arg[-1])/4
    out[1:-1] = (2*arg[1:-1]+ arg[0:-2] + arg[2:])/4
    return out

def wheeler_kiladis(cubes, omega=0.64617667e-05, radius=7160000, period=(0,160), level=38, lat=45, smooth=400):
    
    for cube in cubes:
        if cube.standard_name == 'air_potential_temperature':
            theta = cube[period[0]:period[1],:,:,:].copy()
        if cube.standard_name == 'air_pressure':
            pressure = cube[period[0]:period[1],:,:,:].copy()
            
    # theta.coord('latitude').coord_system = GeogCS(radius)
    # theta.coord('longitude').coord_system = GeogCS(radius)
    run_length, heights, longs = theta.shape[0], np.round(theta.coord('Hybrid height').points*1e-03,0), theta.shape[3]

    p0 = iris.coords.AuxCoord(100000.0, long_name='reference_pressure', units='Pa')
    p0.convert_units(pressure.units)
    temperature = theta*((pressure/p0)**(287.05/1005)) # R and cp in J/kgK for 300K
    # Calculate absolute temperature from potential temperature
    
    zonal_mean_temperature = temperature.collapsed('t', iris.analysis.MEAN)
    eddy_temperature = temperature - zonal_mean_temperature
    print(eddy_temperature.shape)
    # Calculate eddy temperatures
    
    # daily = np.mean(eddy_temperature.data.reshape(-1,4,61,90,144),axis=1)
    # print(daily.shape)    
    data = eddy_temperature[:,level,lat,:].data
    # Extract eddy temperature data for requested atmospheric level and latitude
    # Default latitude is the equator
    
    """ Here you do the Fourier stuff"""
    
    coeffs = np.abs(np.fft.fftshift(np.fft.fft2(data)))**2
    print(coeffs.shape)
    # Compute power spectral density 
    
    wavenumbers = -np.fft.fftshift(np.fft.fftfreq(coeffs.shape[1])*coeffs.shape[1])
    freqs = np.fft.fftshift(np.fft.fftfreq(coeffs.shape[0], d=1./float(4)))
    # Find temporal and spatial frequencies for dataset. Sampling frequency is 4x/day. 
    # Sample size for time is the run length, sample size for space is 144 (number of columns for longitude in the UM)
    pos = freqs > 0  
    # Boolean mask for positive frequencies
    
    background = coeffs.copy()
    for i in range(0,smooth):
        background = wk_filter(background)
        # Smooth a copy of the powers x number of times, where x is given by the arg 'smooth'
     
    signal = wk_filter(coeffs) # Smooth the powers 1 time
    signal = signal/background # Ratio of signal to background
    
    """ Calculate dispersion relation """
    g = 10.9
    lats = pressure.coord('latitude').points*(np.pi/180)
    omega = omega*86400/(2*np.pi) # Convert rotation rate from rad/s to cycles/day
    f0 = 2*omega*np.cos(abs(lats[lat]))
    beta = f0/radius
    
    c = np.sqrt(g*1)

    # freq_scale = (2*np.pi/86400)*freqs[pos]/np.sqrt(beta*c)
    # wavenumber_scale = (wavenumbers/radius)*np.sqrt(c/beta)
    
    # grav_0 = (np.sqrt(wavenumber_scale**2) + 1)*(86400/(2*np.pi))*np.sqrt(beta*c)
    # grav_1 = (np.sqrt(wavenumber_scale**2) + 3)
    # plan_1 = -wavenumber_scale/(3 + wavenumber_scale**2)
    # plan_2 = -wavenumber_scale/(5 + wavenumber_scale**2)
    # kelvin = c*wavenumber_scale
        
    plt.contourf(wavenumbers, freqs[pos], signal[pos,:], brewer_reds.N, cmap=brewer_reds)
    # plt.plot(wavenumbers, grav_0, color='k', linewidth=0.5)
    # plt.plot(wavenumbers, grav_1, color='k', linewidth=0.5)
    # plt.plot(wavenumbers, plan_1, color='k', linewidth=0.5)
    # plt.plot(wavenumbers, plan_2, color='k', linewidth=0.5)
    # plt.plot(wavenumbers, kelvin, color='k', linewidth=0.5)

    plt.title('Wavenumber-frequency spectrum, h=%s km' %(heights[level]))
    plt.xlabel('Zonal wavenumber')
    plt.ylabel('Frequency [CPD]')
    # plt.ylim([0,2])
    plt.colorbar(pad=0.1)
    plt.show()
    
    return coeffs

    