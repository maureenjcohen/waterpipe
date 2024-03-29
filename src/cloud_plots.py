#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 11:49:52 2021

@author: Mo Cohen
"""

import os, iris
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import iris.plot as iplt
import iris.coords
from iris.coords import DimCoord
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import scipy as sp
import windspharm
# import pyvista as pv

brewer_bg = mpl_cm.get_cmap('brewer_PuBu_09')
brewer_red = mpl_cm.get_cmap('brewer_YlOrRd_09')
brewer_redblu = mpl_cm.get_cmap('brewer_RdBu_11')



def plot_clouds(cubes, time_slice=-1, periodicity=False):
    
    """ Plot cloud area fraction for all longs/lats on final output data dump
        Plot cloud volume fraction at longitude 0, final output
        Plot cloud volume fraction at terminators (averaged), final output 
        Plot percentage of latitudes with clear skies at terminators over time 
        If periodicity=True, plot periodicity of clear sky days      
        
        Returns: numpy array of percentage of latitudes with clear skies at terminators """

        
    for cube in cubes:
        if cube.long_name == 'cloud_area_fraction_assuming_maximum_random_overlap':
            cloud_cover = cube.copy()
        if cube.long_name == 'cloud_volume_fraction_in_atmosphere_layer':
            cloud_volume = cube.copy()
        if cube.long_name == 'liquid_cloud_volume_fraction_in_atmosphere_layer':
            liquid_cloud = cube.copy()
        if cube.long_name == 'ice_cloud_volume_fraction_in_atmosphere_layer':
            ice_cloud = cube.copy()
        if cube.standard_name == 'mass_fraction_of_cloud_ice_in_air':
            ice_condensate = cube.copy()
        if cube.standard_name == 'mass_fraction_of_cloud_liquid_water_in_air':
            liquid_condensate = cube.copy()
    
    iplt.contourf(cloud_cover[time_slice,:,:], brewer_bg.N, cmap=brewer_bg)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Cloud Area Fraction', y=1.20)
    plt.colorbar(pad=0.1)
    plt.show()

    iplt.contourf(cloud_volume[time_slice,:,:,0], brewer_bg.N, cmap=brewer_bg)
    plt.title('Cloud Volume Fraction at Longitude 0', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    iplt.contourf((cloud_volume[time_slice,:,:,36]+cloud_volume[time_slice,:,:,108])/2, brewer_bg.N, cmap=brewer_bg)
    plt.title('Cloud Volume Fraction at Terminators [kg kg-1]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    clear_list = []
    for day in range(0,cloud_cover.shape[0]):
        data = (cloud_cover[day,:,36].data + cloud_cover[day,:,108].data)/2
        clear = data[np.where( data <= 0.2)]
        count = len(clear.tolist())
        clear_list.append(count)
    
    average = round(np.mean(np.array(clear_list)), 3)
    std_dev = round(np.std(np.array(clear_list)), 3)
    
    time_axis = np.arange(0,cloud_cover.shape[0])
    plt.plot(time_axis,(np.array(clear_list)/90)*100)
    plt.title('Percent of latitudes with clear skies (cloud frac <= 0.2)')
    plt.xlabel('Days')
    plt.ylabel('Percent')
    plt.text(x=0,y=0, s="Mean = %s, std = %s"%(average, std_dev))
    plt.show()
    # Get an overview of what percentage of days might allow observation
    # of atmospheric spectrum in transmission spectroscopy
    
    if periodicity == True:
        
        daily_cloud = np.array(clear_list)
        run_length = cloud_cover.shape[0]
        cloud_fft = sp.fftpack.fft(daily_cloud)
        cloud_psd = np.abs(cloud_fft)**2
        cloudfreq = sp.fftpack.fftfreq(len(cloud_psd), 1./run_length)
        i = cloudfreq > 0
        
        fig, ax = plt.subplots(1,1,figsize=(8,4))
        ax.plot(cloudfreq[i], cloud_psd[i])
        ax.set_xlim(0,15)
        ax.set_xlabel('Frequency [1/run_length]')
        ax.set_ylabel('PSD')
        ax.set_title('Periodicity of clear sky')
        # Check if there is periodicity in cloud cover over time
        # Uses scipy discrete Fast Fourier Transform. Units of 1/run_length as
        # there is no meaningful time period in a tidally-locked planet without
        # eccentricity or obliquity. The plot tells you the amplitude of a cycle
        # and how many times it repeated during the model run time.
    
    return np.array(clear_list)


def cloud_lats(cubes, time_slice=-1, long=0):
    
    for cube in cubes:

        # if cube.long_name == 'cloud_volume_fraction_in_atmosphere_layer':
        #     cloud_volume = cube.copy()
        # if cube.long_name == 'liquid_cloud_volume_fraction_in_atmosphere_layer':
        #     liquid_cloud = cube.copy()
        # if cube.long_name == 'ice_cloud_volume_fraction_in_atmosphere_layer':
        #     ice_cloud = cube.copy()
        if cube.standard_name == 'mass_fraction_of_cloud_ice_in_air':
            ice_condensate = cube.copy()
        if cube.standard_name == 'mass_fraction_of_cloud_liquid_water_in_air':
            liquid_condensate = cube.copy()
    
    longitudes = liquid_condensate.coord('longitude').points
    for cube in (
                 ice_condensate, liquid_condensate):
        
        if cube.standard_name == None:
            title = cube.long_name
        else:
            title = cube.standard_name
            

        iplt.contourf(cube[time_slice,:,:,long], brewer_bg.N, cmap=brewer_bg)
        plt.title('%s at %s' %(title,longitudes[long]), y=1.05)
        plt.ylabel('Height [m]')
        plt.xlabel('Latitude [degrees]')
        cbar = plt.colorbar(pad=0.1)
        cbar.ax.set_title('kg/kg')
        plt.show()
        
    
def cloud_type(cubes, time_slice=-1, long1=36, long2=108, filtering=True, 
               ice_threshold = 0.02, liq_threshold = 0.02):
    
    for cube in cubes:

        # if cube.long_name == 'cloud_volume_fraction_in_atmosphere_layer':
        #     cloud_volume_raw = cube.copy()
        # if cube.long_name == 'liquid_cloud_volume_fraction_in_atmosphere_layer':
        #     liquid_cloud_raw = cube.copy()
        # if cube.long_name == 'ice_cloud_volume_fraction_in_atmosphere_layer':
        #     ice_cloud_raw = cube.copy()
        if cube.standard_name == 'mass_fraction_of_cloud_ice_in_air':
            ice_condensate_raw = cube.copy()
        if cube.standard_name == 'mass_fraction_of_cloud_liquid_water_in_air':
            liquid_condensate_raw = cube.copy()
    
    if long1==36 and long2==108:
        titleloc = 'Terminators'
    elif long1==143 and long2==0:
        titleloc = 'Substellar Point'
    else:
        titleloc = '20W'
    
    # liquid_cloud = liquid_cloud_raw.data
    # ice_cloud = ice_cloud_raw.data
    ice_condensate = ice_condensate_raw.data
    liquid_condensate = liquid_condensate_raw.data
    
    y_axis = np.round(liquid_condensate_raw.coord('level_height').points*1e-03,0)
    x_axis = liquid_condensate_raw.coord('latitude').points
    time_axis = np.arange(0,liquid_condensate.shape[0])
    
    # plt.figure(figsize=(10,5))
    # plt.contourf(x_axis, y_axis, (liquid_cloud[time_slice,:,:,long1]+liquid_cloud[time_slice,:,:,long2])/2, brewer_bg.N, cmap=brewer_bg)
    # cbarl = plt.colorbar(pad=0.1)
    # cbarl.ax.set_title('liq')
    # plt.title('Liquid Cloud Volume Fraction at Terminators, day=%s' %time_slice, y=1.05)
    # plt.ylabel('Height [km]')
    # plt.xlabel('Latitude [degrees]')
    # plt.show()
 
    # plt.figure(figsize=(10,5))
    # plt.contourf(x_axis, y_axis, (ice_cloud[time_slice,:,:,long1]+ice_cloud[time_slice,:,:,long2])/2, brewer_bg.N, cmap=brewer_bg)
    # cbari = plt.colorbar(pad=0.1)
    # cbari.ax.set_title('ice')
    # plt.title('Ice Cloud Volume Fraction at Terminators, day=%s' %time_slice, y=1.05)
    # plt.ylabel('Height [km]')
    # plt.xlabel('Latitude [degrees]')
    # plt.show()
    
    ice_east = np.mean(ice_condensate[:,:,:,long1], axis=(1,2))
    ice_west = np.mean(ice_condensate[:,:,:,long2], axis=(1,2))
    total_ice = (ice_east + ice_west)/2
    
    liq_east = np.mean(liquid_condensate[:,:,:,long1], axis=(1,2))
    liq_west = np.mean(liquid_condensate[:,:,:,long2], axis=(1,2))
    total_liq = (liq_east + liq_west)/2
    
    if filtering == True:
        fft_ice = sp.fftpack.fft(total_ice)
        psd_ice = np.abs(fft_ice)**2
        freqs_ice = sp.fftpack.fftfreq(len(psd_ice), 1./1)
        lowpass_ice = fft_ice.copy()
        lowpass_ice[np.abs(freqs_ice) > ice_threshold] = 0
        total_ice = np.real(sp.fftpack.ifft(lowpass_ice))
        
        fft_liq = sp.fftpack.fft(total_liq)
        psd_liq = np.abs(fft_liq)**2
        freqs_liq = sp.fftpack.fftfreq(len(psd_liq), 1./1)
        lowpass_liq = fft_liq.copy()
        lowpass_liq[np.abs(freqs_liq) > liq_threshold] = 0
        total_liq = np.real(sp.fftpack.ifft(lowpass_liq))
    
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('Time [days]')
    ax1.set_ylabel('Ice condensate [kg/kg]')
    ax1.plot(time_axis, total_ice, color='b', label='Ice')
    ax1.tick_params(axis='y', labelcolor='b')
    
    ax2 = ax1.twinx()
    ax2.set_ylabel('Liquid condensate [kg/kg]')
    ax2.plot(time_axis, total_liq, color='r', label='Liquid')
    ax2.tick_params(axis='y', labelcolor='r')
    
    plt.title('Mean Ice and Liquid Condensate at %s' %titleloc)
    fig.tight_layout()
    plt.show()
    
    # all_ice = np.sum(ice_condensate, axis=(1,2,3))
    # all_liq = np.sum(liquid_condensate, axis=(1,2,3))
    # total_all = all_ice + all_liq

    # plt.plot(time_axis, total_all)
    # plt.title('Total Planetary Condensates')
    # plt.xlabel('Time [days]')
    # plt.ylabel('Cloud condensate [kg/kg]')
    # plt.show()
    plt.plot(time_axis, total_ice, color='b', label='Ice')
    plt.plot(time_axis, total_liq, color='r', label='Liquid')
    plt.title('Mean Ice and Liquid Condensate at %s' %titleloc)
    plt.xlabel('Time [days]')
    plt.ylabel('Cloud condensate [kg/kg]')
    plt.legend()
    plt.show()
    

def wind_speed(cubes, level=25, lat=45):
    
    for cube in cubes:
        if cube.standard_name == 'x_wind': 
            x_wind = cube.copy()
            y_axis = np.round(x_wind.coord('level_height').points*1e-03,0)
            time_unit = 'days'

        elif cube.standard_name == 'eastward_wind':
            x_wind = cube.copy()
            y_axis = np.round(x_wind.coord('Hybrid height').points*1e-03,0)
            time_unit = 'days'
            
    x_axis = x_wind.coord('latitude').points
    time_axis = np.arange(0,x_wind.shape[0])

    zonal_mean = x_wind.collapsed('longitude', iris.analysis.MEAN) 
    data = zonal_mean[:,level,lat].data
    
    plt.plot(time_axis, data)
    plt.title('Zonal mean zonal wind at equator, h = %s km [m/s]' %y_axis[level])
    plt.xlabel('Time [%s]' %time_unit)
    plt.ylabel('Wind speed [m/s]')
    plt.show()
    


def plot_vorticity(cubes, level=15, time_slice=-1, omega=0.64617667):
    
    """ Uses windspharm package to plot the vorticity       """
    
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube.copy()
        if cube.standard_name == 'y_wind':
            y_wind = cube.copy()
   
    
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
    heights = np.round(x_wind.coord('level_height').points*1e-03,0)

    wind = windspharm.iris.VectorWind(x_wind, y_wind)

    planet_vort = wind.planetaryvorticity(omega=omega)
    relative_vort = wind.vorticity()
    absolute_vort = wind.absolutevorticity()
    
    iplt.contourf(relative_vort[time_slice,level,:,:], brewer_redblu.N, cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    plt.title('Relative Vorticity, h = %s km' %(heights[level]), y=1.20)
    plt.ylabel('Latitude [degrees]')
    plt.xlabel('Longitude [degrees]')
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.colorbar(orientation='horizontal')
    plt.show()
    

def vorticity_series(cubes, level=25, lats=(60,76), longs=(72,109), omega=0.64617667):    
    
    for cube in cubes:
        if cube.standard_name == 'x_wind' or cube.standard_name == 'eastward_wind':
            x_wind = cube.copy()
        if cube.standard_name == 'y_wind' or cube.standard_name == 'northward_wind':
            y_wind = cube.copy()
   
    
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
    # heights = np.round(x_wind.coord('level_height').points*1e-03,0)

    wind = windspharm.iris.VectorWind(x_wind, y_wind)

    planet_vort = wind.planetaryvorticity(omega=omega)
    relative_vort = wind.vorticity()
    
    rlats, rlongs =  relative_vort.coord('latitude'), relative_vort.coord('longitude') 
    
    # rlats.bounds = None
    # relative_vort.coord('latitude').guess_bounds()
    # rlongs.bounds = None
    # relative_vort.coord('latitude').guess_bounds()
    
    patch = relative_vort[:,level,lats[0]:lats[1], longs[0]:longs[1]]
    # patch_weights = iris.analysis.cartography.area_weights(patch)
    patch_mean = patch.collapsed(['latitude', 'longitude'], iris.analysis.MEAN)
    data = patch_mean.data
    
    plt.plot(np.arange(0,x_wind.shape[0]), data)
    plt.title('Mean relative vorticity')
    plt.xlabel('Time [days]')
    plt.show()
    
    
def plot_energy(cubes, lats=(30,61), longs=(0,13)):

    for cube in cubes:
        if cube.standard_name == 'surface_net_downward_shortwave_flux':
            shortwave = cube.copy()
        if cube.standard_name == 'surface_upward_latent_heat_flux':
            latent = cube.copy()
        if cube.standard_name == 'surface_upward_sensible_heat_flux':
            sensible = cube.copy()
            
    shlats, lalats, selats = shortwave.coord('latitude'), latent.coord('latitude'), sensible.coord('latitude')
    shlongs, lalongs, selongs = shortwave.coord('longitude'), latent.coord('longitude'), sensible.coord('longitude')
    time_axis = np.arange(0,shortwave.shape[0])

    if shlats.bounds == None:
        shortwave.coord('latitude').guess_bounds()
    if lalats.bounds == None:
        latent.coord('latitude').guess_bounds()
    if selats.bounds == None:
        sensible.coord('latitude').guess_bounds()
    if shlongs.bounds == None:
        shortwave.coord('longitude').guess_bounds()
    if lalongs.bounds == None:
        latent.coord('longitude').guess_bounds()
    if selongs.bounds == None:
        sensible.coord('longitude').guess_bounds()
            
    for cube in (shortwave, latent, sensible):  
        
        patch = cube[:,lats[0]:lats[1],longs[0]:longs[1]]
        patch_weights = iris.analysis.cartography.area_weights(patch)
        patch_mean = patch.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=patch_weights)

        plt.plot(time_axis, patch_mean.data)
        plt.title('mean %s' %cube.standard_name)
        plt.xlabel('Time [days]')
        plt.show()
        
def composites2D(cubes, start=0, end=10, nlev=38, level=25, n=3, meaning=7):
    
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube.copy()
        if cube.standard_name == 'y_wind':
            y_wind = cube.copy()
        if cube.long_name == 'cloud_area_fraction_assuming_maximum_random_overlap':
            cloud_cover = cube.copy()

    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())            
    meaned_x = np.mean(x_wind.data.reshape(-1,meaning,nlev,90,144),axis=1)
    meaned_y = np.mean(y_wind.data.reshape(-1,meaning,nlev,90,144),axis=1)
    meaned_cloud = np.mean(cloud_cover.data.reshape(-1,meaning,90,144), axis=1)
    
    # X,Y = np.meshgrid(np.arange(-meaned_x.shape[3]/2,meaned_x.shape[3]/2), np.arange(-meaned_x.shape[2]/2,meaned_x.shape[2]/2))
    X,Y = np.meshgrid(np.arange(0,144), np.arange(0,90))

    
    for time_slice in range(start, end+1):
        
        fig, ax = plt.subplots(figsize=(10,5))
        plt.imshow(np.roll(meaned_cloud[time_slice,:,:],72,axis=1), cmap=brewer_bg)
        plt.colorbar()
    
        
        plt.quiver(X[::n,::n],Y[::n,::n], np.roll(meaned_x[time_slice,level,::n,::n],72,axis=1), 
                    np.roll(-meaned_y[time_slice, level,::n,::n],72,axis=1),scale_units='xy',scale=5)

        plt.title('Cloud cover and horizontal wind at fortnight %s' %(time_slice+1))
        # plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
        # plt.yticks((-45,-30,-15,0,15,30,45),('90S','60S','30S','0','30N','60N','90N'))         
        # plt.xlabel('Longitude')
        # plt.ylabel('Latitude')
        
        # ax2.quiverkey(q1, X=0.9, Y=1.05, U=10, label='10 m/s', labelpos='E', coordinates='axes')

    
        plt.show()
        
        
def composites_levels(cubes, start=1,end=20, nlat=90, nlon=144, nlev=38, level=25, meaning=1, n=4, cloudtype='ice', fractype='mass'):

    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube.copy()
        if cube.standard_name == 'y_wind':
            y_wind = cube.copy() 
        if cube.long_name == 'ice_cloud_volume_fraction_in_atmosphere_layer' and fractype=='volume':
            ice = cube.copy()
        if cube.long_name == 'liquid_cloud_volume_fraction_in_atmosphere_layer' and fractype=='volume':
            liq = cube.copy()
        if cube.standard_name == 'mass_fraction_of_cloud_ice_in_air' and fractype=='mass':
            ice = cube.copy()
        if cube.standard_name == 'mass_fraction_of_cloud_liquid_water_in_air' and fractype=='mass':
            liq = cube.copy()
            
        
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
    heights = np.round(x_wind.coord('level_height').points*1e-03,2)
            
    meaned_x = np.mean(x_wind.data.reshape(-1,meaning,nlev,nlat,nlon),axis=1)
    meaned_y = np.mean(y_wind.data.reshape(-1,meaning,nlev,nlat,nlon),axis=1)
    
    if cloudtype=='ice':
        meaned_cloud = np.mean(ice.data.reshape(-1,meaning,nlev+1,nlat,nlon), axis=1)
        titleterm = 'Ice cloud'
    elif cloudtype=='liq':
        meaned_cloud = np.mean(liq.data.reshape(-1,meaning,nlev+1,nlat,nlon), axis=1)
        titleterm = 'Liquid cloud'
    else:
        print('Argument cloudtype must be ice or liq. Default is ice.')

    X,Y = np.meshgrid(np.arange(0,nlon), np.arange(0,nlat))   

    for time_slice in range(start, end+1):
        
        fig, ax = plt.subplots(figsize=(10,5))
        # plt.imshow(np.roll(meaned_cloud[time_slice,level, :,:],int(nlon/2),axis=1), cmap=brewer_bg)
        # cbar = plt.colorbar()
        
        # if fractype=='mass':
        #     cbar.ax.set_title('$10^{-4}$kg/kg')
    
        
        plt.quiver(X[::n,::n],Y[::n,::n], np.roll(meaned_x[time_slice,level,::n,::n],int(nlon/(2*n)),axis=1), 
                    np.roll(meaned_y[time_slice, level,::n,::n],int(nlon/(2*n)),axis=1),scale_units='xy',scale=5)

        plt.title('%s and horizontal wind, days %s to %s, h=%s km' %(titleterm, time_slice*meaning-meaning,time_slice*meaning, heights[level]))
        plt.xticks((0,12,24,36,48,60,72,84,96,108,120,132,144),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
        plt.yticks((90,75,60,45,30,15,0),('90S','60S','30S','0','30N','60N','90N'))   
        # plt.xticks((0,128,256,384,512),('180W','90W','0','90E','180E'))
        # plt.yticks((256,128,0),('90S','0','90N')) 
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        
        # ax2.quiverkey(q1, X=0.9, Y=1.05, U=10, label='10 m/s', labelpos='E', coordinates='axes')
    
        plt.show()
        

def cloud_profile(cubes, time_slice=340, lon=36, cloudtype='ice', fractype='volume'):
    
    for cube in cubes:
        if cube.long_name == 'ice_cloud_volume_fraction_in_atmosphere_layer' and fractype=='volume':
            ice = cube[time_slice,:,:,lon].copy()
        if cube.long_name == 'liquid_cloud_volume_fraction_in_atmosphere_layer' and fractype=='volume':
            liq = cube[time_slice,:,:,lon].copy()
        if cube.standard_name == 'mass_fraction_of_cloud_ice_in_air' and fractype=='mass':
            ice = cube[time_slice,:,:,lon].copy()
        if cube.standard_name == 'mass_fraction_of_cloud_liquid_water_in_air' and fractype=='mass':
            liq = cube[time_slice,:,:,lon].copy()
                            
    heights = np.round(ice.coord('level_height').points*1e-03,0)
    lats = ice.coord('latitude').points
    
    if cloudtype=='ice':
        cloud = ice.data
        titleterm = 'Ice cloud'
    elif cloudtype=='liq':
        cloud = liq.data
        titleterm = 'Liquid cloud'
    else:
        print('Argument cloudtype must be ice or liq. Default is ice.')
    
    fig, ax = plt.subplots(figsize=(10,5))
    plt.contourf(lats, heights, cloud, cmap=brewer_bg)
    cbar = plt.colorbar()
    
    # if fractype=='mass':
    #     cbar.ax.set_title('kg/kg')
    
    plt.title('%s at lon %s, day %s' %(titleterm, lon*2.5, time_slice))
    plt.ylabel('Height [km]')
    plt.xlabel('Latitude [degrees]')
    plt.show()


def cloud_structure(cubes, lat=45, time_slice=-1, fractype='mass'):
    
    for cube in cubes:
        if cube.long_name == 'ice_cloud_volume_fraction_in_atmosphere_layer' and fractype=='volume':
            ice = cube[time_slice,:,lat,:].copy()
        if cube.long_name == 'liquid_cloud_volume_fraction_in_atmosphere_layer' and fractype=='volume':
            liq = cube[time_slice,:,lat,:].copy()
        if cube.standard_name == 'mass_fraction_of_cloud_ice_in_air' and fractype=='mass':
            ice = cube[time_slice,:,lat,:].copy()
        if cube.standard_name == 'mass_fraction_of_cloud_liquid_water_in_air' and fractype=='mass':
            liq = cube[time_slice,:,lat,:].copy()
            
    heights = np.round(ice.coord('level_height').points*1e-03,0)
    lons = ice.coord('longitude').points
    
    if len(heights) > 39:
        final_height = 43
    else:
        final_height = -1
    
    fig, ax = plt.subplots(figsize=(10,5))
    
    plota = ax.contourf(np.roll(lons,72), heights[:final_height], ice[:final_height,:].data*10**4, cmap='Blues')
    cba = plt.colorbar(plota)
    cba.ax.set_title('$10^{-4}$ kg/kg', size=10)
    cba.set_label('Ice', size=15)
    
    plotb = ax.contourf(np.roll(lons,72), heights[:final_height], liq[:final_height,:].data*10**4, cmap='Reds', alpha=0.5)
    cbb = plt.colorbar(plotb)
    cbb.ax.set_title('$10^{-4}$ kg/kg', size=10)
    cbb.set_label('Liquid', size=15)
    
    plt.xticks([0,30,60,90,120,150,180,210,240,270,300,330,360], ['180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'])

    ax.set_xlabel('Longitude [degrees]')
    ax.set_ylabel('Height [km]')
    
    ax.set_title('Cloud structure at lat=%s' %((lat-45)*2))
    plt.show()
    

def y_gradient(cubes,start=0,end=-1,lon1=0,lon2=-1):
    
    for cube in cubes:
        if cube.standard_name == 'air_potential_temperature':
            theta = cube[start:end,:,:,lon1:lon2].copy()
        if cube.standard_name == 'air_pressure':
            pressure = cube[start:end,:,:,lon1:lon2].copy()
            
    p0 = iris.coords.AuxCoord(100000.0, long_name='reference_pressure', units='Pa')
    p0.convert_units(pressure.units)
    temperature = theta*((pressure/p0)**(287.05/1005)) # R and cp in J/kgK for 300K
    # Calculate absolute temperature from potential temperature
    
    zm_temp = temperature.collapsed('longitude',iris.analysis.MEAN)
    tm_temp = temperature.collapsed('time',iris.analysis.MEAN)
    zmtm = temperature.collapsed(['longitude','time'],iris.analysis.MEAN)
    
    zmtm = zmtm.data
    
    heights = np.round(theta.coord('level_height').points*1e-03,0)
    lats = np.round(theta.coord('latitude').points,0)
    
    plt.contourf(lats,heights,zmtm,cmap=brewer_redblu)
    plt.title('Air temperature [K], zonal mean')
    plt.xlabel('Latitude [degrees]')
    plt.ylabel('Height [km]')
    plt.colorbar()
    plt.show()
    
        