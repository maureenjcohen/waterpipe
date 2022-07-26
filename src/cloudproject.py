#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 09:42:00 2022

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
from numpy import unravel_index
import scipy as sp
import windspharm
from iris.analysis import calculus


redblu = mpl_cm.get_cmap('RdBu')
plasma = mpl_cm.get_cmap('plasma')
brewer_bg = mpl_cm.get_cmap('brewer_PuBu_09')

def composite(cubes, time_slice=500, nlat=90, nlon=144, nlev=38, level=8, n=4, cloudtype='both', fractype='mass', save='no'):

    """ Plot composites of the cloud cover and horizontal wind vectors
    
    Inputs: Iris CubeList, start and end times of data, number of latitudes, number of longitudes,
    number of levels, level to be plotted, meaning period (1 = no meaning), sparsity of quiver arrows,
    cloud type (options: ice, liq, both, none), fraction type (options: mass or volume), whether to save

    Outputs: Plot of horizontal wind vectors, with or without the cloud cover superimposed"""

    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube[time_slice,:,:,:].copy()
        if cube.standard_name == 'y_wind':
            y_wind = cube[time_slice,:,:,:].copy() 
        if cube.long_name == 'ice_cloud_volume_fraction_in_atmosphere_layer' and fractype=='volume':
            ice = cube[time_slice,:,:,:].copy()
        if cube.long_name == 'liquid_cloud_volume_fraction_in_atmosphere_layer' and fractype=='volume':
            liq = cube[time_slice,:,:,:].copy()
        if cube.standard_name == 'mass_fraction_of_cloud_ice_in_air' and fractype=='mass':
            ice = cube[time_slice,:,:,:].copy()
        if cube.standard_name == 'mass_fraction_of_cloud_liquid_water_in_air' and fractype=='mass':
            liq = cube[time_slice,:,:,:].copy()
            
        
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
    heights = np.round(x_wind.coord('level_height').points*1e-03,2)
            
    
    if cloudtype=='ice':
        cloud = ice
        titleterm = 'Ice cloud'
    elif cloudtype=='liq':
        cloud = liq
        titleterm = 'Liquid cloud'
    elif cloudtype=='both':
        cloud = ice + liq
        titleterm = 'Total cloud'
    elif cloudtype=='none':
        titleterm = 'Horizontal wind'
    else:
        print('Argument cloudtype must be ice, liq, both, or none. Default is both.')

    X,Y = np.meshgrid(np.arange(0,nlon), np.arange(0,nlat))   

    if cloudtype=='none':
        fig, ax = plt.subplots(figsize=(8.5,5))
        q1 = ax.quiver(X[::n,::n],Y[::n,::n], np.roll(x_wind[level,::n,::n].data,int(nlon/(2*n)),axis=1), 
                   np.roll(y_wind[level,::n,::n].data,int(nlon/(2*n)),axis=1),scale_units='xy',scale=5)
        ax.quiverkey(q1, X=0.9, Y=1.05, U=20, label='20 m/s', labelpos='E', coordinates='axes')
        plt.title('%s, day %s, h=%s km' %(titleterm, time_slice, heights[level]))
        plt.xticks((0,12,24,36,48,60,72,84,96,108,120,132,144),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
        plt.yticks((0,15,30,45,60,75,90),('90S','60S','30S','0','30N','60N','90N'))    
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
  
        if save == 'yes':
            plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/cloudproject/epsfigs/quiver_nocloud_%s.eps' %(time_slice), format='eps')
        else:
            pass
        plt.show()

    else:
        fig, ax = plt.subplots(figsize=(10,5))
        plt.imshow(np.roll(cloud[level,:,:].data*1e4,int(nlon/2),axis=1), cmap=brewer_bg)
        cbar = plt.colorbar()

        if fractype=='mass':
            cbar.set_label('$10^{-4}$ kg/kg', loc='center')    

        q1 = ax.quiver(X[::n,::n],Y[::n,::n], np.roll(x_wind[level,::n,::n].data,int(nlon/(2*n)),axis=1), 
                np.roll(-y_wind[level,::n,::n].data,int(nlon/(2*n)),axis=1),scale_units='xy',scale=5)
        ax.quiverkey(q1, X=0.95, Y=1.05, U=20, label='20 m/s', labelpos='E', coordinates='axes')

        plt.title('%s and horizontal wind, day %s, h=%s km' %(titleterm, time_slice, heights[level]))
        plt.xticks((0,12,24,36,48,60,72,84,96,108,120,132,144),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
        plt.yticks((90,75,60,45,30,15,0),('90S','60S','30S','0','30N','60N','90N'))    
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        if save == 'yes':
            plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/cloudproject/epsfigs/quiver_withcloud_%s.eps' %(time_slice), format='eps')
        else:
            pass            
        plt.show()

def rwave_velocity(cubes,start=500,end=600,nlat=90,nlon=144,level=8,omega=1.19e-05,radius=5797818,lat=55,save='no'):

    """ This function calculates the Rossby wave phase speed (including zonal wind) over time.

    The code performs a Helmholtz decomposition of the horizontal wind field at the input level.
    The wind vectors are converted to a wind speed (magnitude) and this 2-D data is input into a 
    2-D Fourier transform. The Fourier transform finds the spatial frequencies in the x- and y-directions.
    The function finds the maximum intensity wave and extracts the zonal and meridional wavenumbers.
    It then constructs a time series of the Rossby wave phase speed using the varying wavenumbers and 
    varying background zonal wind over time as inputs into the phase speed formula.

    At certain latitudes, the phase velocity alternates between positive and negative. This is the region
    where the cyclonic structure appears to travel back and forth periodically. """

    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'y_wind':
            y_wind = cube[start:end,:,:,:].copy()
        
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
    km_heights = np.round(x_wind.coord('level_height').points*1e-03,2)
    
    winds = windspharm.iris.VectorWind(x_wind[:,level,:,:],y_wind[:,level,:,:])
    uchi,vchi,upsi,vpsi = winds.helmholtz(truncation=21)

    
    zonal_upsi = upsi.collapsed('longitude', iris.analysis.MEAN)
    zonal_vpsi = vpsi.collapsed('longitude', iris.analysis.MEAN)
    eddy_upsi = upsi - zonal_upsi
    eddy_vpsi = vpsi - zonal_vpsi
    magnitude = np.sqrt(eddy_upsi.data**2 + eddy_vpsi.data**2)
    
    fft2 = sp.fft.fftshift(sp.fftpack.fft2(sp.fft.ifftshift(magnitude)))
    yfreqs = sp.fft.fftshift(sp.fft.fftfreq(fft2.shape[1],d=1./nlat))
    xfreqs = sp.fft.fftshift(sp.fft.fftfreq(fft2.shape[2],d=1./nlon))
    psd = np.abs(fft2)**2
    
    quadrant = psd[:,45:51,72:78] # Select positive zonal and meridional wavenumbers <= 5 
    xfreqs = xfreqs[72:78] # Select corresponding x- and y-frequencies
    yfreqs = yfreqs[45:51]

    x_wavenumbers = []
    y_wavenumbers = []
    for time in range(0,psd.shape[0]):
        peak = unravel_index(quadrant[time,:,:].argmax(), quadrant[time,:,:].shape)
        x_wavenumbers.append(peak[1])
        y_wavenumbers.append(peak[0])
        
    x_wavenumbers, y_wavenumbers = np.array(x_wavenumbers)+1, np.array(y_wavenumbers)+1
    
    zmzw_lat = lat 
    lat = np.abs(lat-45)*2 # Convert input row number to latitude in degrees north
    zmzw = x_wind[:,level,zmzw_lat].collapsed('longitude',iris.analysis.MEAN)
    zmzw = zmzw.data
        
    lat_rad = lat*(np.pi/180) # Convert input latitude to radians
    beta = 2*omega*np.cos(lat_rad)/radius # Beta factor    
    circum = 2*np.pi*radius*np.cos(lat_rad) # Circumference in meters at input latitude    
    lambda_x = circum/x_wavenumbers # Wavelength in x-direction at input latitude
    x_num = 2*np.pi/lambda_x
    lambda_y = 2*np.pi*radius/y_wavenumbers # Wavelength in y-direction at input latitude
    y_num = 2*np.pi/lambda_y
    
    c_phase = (zmzw - (beta/(x_num**2+y_num**2)))
    
    plt.plot(c_phase,color='b')
    plt.title('Rossby wave phase velocity at %s N' %(lat))
    plt.xlabel('Time [days]')
    plt.ylabel('Velocity [m/s]')
    if save == 'yes':
        plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/cloudproject/epsfigs/rwave_vel.eps', format='eps')
    else:
        pass  
    plt.show()


def vertical_profile(cubes,start=500,end=600,level=8,top_level=30,select='absolute',save='no'):

    """ This function plots the vertical profile of the chosen data input over time.
     Possible inputs: air temperature ('absolute'), potential temperature ('potential'),
     air pressure ('pressure'), vertical wind ('z_wind'), cloud ('cloud') """

    for cube in cubes:
        if cube.standard_name == 'air_potential_temperature':
            theta = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'air_pressure':
            pressure = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'upward_air_velocity':
            z_wind = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'mass_fraction_of_cloud_ice_in_air':
            ice = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'mass_fraction_of_cloud_liquid_water_in_air':
            liq = cube[start:end,:,:,:].copy()

    heat = mpl_cm.get_cmap('gist_heat')
    blues = mpl_cm.get_cmap('Blues')
    redblu = mpl_cm.get_cmap('RdBu')

    p0 = iris.coords.AuxCoord(
        100000.0, long_name='reference_pressure', units='Pa')
    p0.convert_units(pressure.units)
    temperature = theta*((pressure/p0)**(287.05/1005)) # R and cp in J/kgK for 300K

    if select == 'absolute':
        datacube = temperature.copy()
        titleterm = 'absolute temperature'
        y_axis = 'Temperature [K]'
        colors = heat
        unit = 'K'
        norm = TwoSlopeNorm(0)
    elif select == 'potential':
        datacube = theta.copy()
        titleterm = 'potential temperature'
        y_axis = 'Temperature [K]'
        colors = heat
        unit = 'K'
        norm = None
    elif select == 'pressure':
        datacube = pressure.copy()
        titleterm = 'air pressure'
        y_axis = 'Pressure [Pa]'
        colors = blues
        unit = 'Pa'
        norm = None
    elif select == 'z_wind':
        datacube = (z_wind.copy())*1e4
        titleterm = 'vertical wind'
        y_axis = 'Wind speed [m/s]'
        colors = redblu
        unit = '$10^{-4}$ m/s'
        norm = TwoSlopeNorm(0)
    elif select == 'cloud':
        datacube = (ice.copy() + liq.copy())*1e6
        titleterm = 'total cloud (ice and liquid)'
        y_axis = 'Cloud mass [kg/kg]'
        colors = blues
        unit = '$10^{-6}$ kg/kg'
        norm = None

    heights = np.round(datacube.coord('level_height').points*1e-03,0)
    time_axis = np.arange(0,datacube.shape[0])
    lats = datacube.coord('latitude')
    lons = datacube.coord('longitude')

    if lats.bounds == None:
        datacube.coord('latitude').guess_bounds()
    if lons.bounds == None:
        datacube.coord('longitude').guess_bounds()

    datacube = datacube.extract(iris.Constraint(longitude=lambda v: 270 < v <= 359 or 0 <= v <= 90, latitude=lambda v: -90 <= v <= 90))
    grid_areas = iris.analysis.cartography.area_weights(datacube)
    dayside_mean = datacube.collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=grid_areas)
    z_axis = heights[:top_level]
    dayside_time = dayside_mean[:,:top_level].data

    fig1, ax1 = plt.subplots(figsize=(8,5))
    ax1.set_xlabel('Time [days]')
    ax1.set_ylabel('%s' %y_axis)
    ax1.set_title('Dayside mean %s at h=%s km' %(titleterm, heights[level]))
    plt.plot(dayside_mean[:,level].data)
    plt.show()

    fig2, ax2 = plt.subplots(figsize=(8,5))
    ax2.set_xlabel('Time [days]')
    ax2.set_ylabel('Height [km]')
    ax2.set_title('Dayside mean %s' %titleterm)
    plt.contourf(time_axis, z_axis, dayside_time.T, np.arange(np.round(np.min(dayside_time),0),np.round(np.max(dayside_time),0)), cmap=colors, norm=norm)
    cbar = plt.colorbar(pad=0.1)
    cbar.ax.set_title('%s' %unit)

    if save == 'yes':
        plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/cloudproject/epsfigs/%s_%s_%s.eps' %(select,start,end), format='eps')
    else:
        pass

    plt.show() 