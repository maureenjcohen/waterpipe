#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 20:51:32 2022

@author: Mo Cohen
"""
import iris, windspharm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import scipy as sp
from numpy import unravel_index
from iris.coord_systems import GeogCS
from matplotlib.colors import TwoSlopeNorm



redblu = mpl_cm.get_cmap('RdBu')
plasma = mpl_cm.get_cmap('plasma')



def model_vel(cubes,start=500,end=600,level=8,omega=7.93e-06,radius=7160000,lat=55):
    
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'y_wind':
            y_wind = cube[start:end,:,:,:].copy()
        
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
    km_heights = np.round(x_wind.coord('level_height').points*1e-03,2)
    
    winds = windspharm.iris.VectorWind(x_wind[:,level,:,:],y_wind[:,level,:,:])
    # Create a VectorWind data object from the x and y wind cubes
    uchi,vchi,upsi,vpsi = winds.helmholtz(truncation=21)
    # Calculate the Helmholtz decomposition. Truncation is set to 21 because 
    # this is what Hammond and Lewis 2021 used.
    
    zonal_upsi = upsi.collapsed('longitude', iris.analysis.MEAN)
    zonal_vpsi = vpsi.collapsed('longitude', iris.analysis.MEAN)
    # Calculate zonal means of the x and y components of the rotational component
    eddy_upsi = upsi - zonal_upsi
    eddy_vpsi = vpsi - zonal_vpsi
    magnitude = np.sqrt(eddy_upsi.data**2 + eddy_vpsi.data**2)
    
    fft2 = sp.fft.fftshift(sp.fftpack.fft2(sp.fft.ifftshift(magnitude)))
    yfreqs = sp.fft.fftshift(sp.fft.fftfreq(fft2.shape[1],d=1./90))
    xfreqs = sp.fft.fftshift(sp.fft.fftfreq(fft2.shape[2],d=1./144))
    psd = np.abs(fft2)**2
    
    quadrant = psd[:,45:51,72:78]
    xfreqs = xfreqs[72:78]
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
    # c_group = (zmzw + (beta*(x_num**2-y_num**2))/((x_num**2+y_num**2)**2))
    
    # distance = np.cumsum(c_phase*60*60*24)*10e-03 # Convert m/s to m/day = distance travelled in a day
    # deg_dist = (distance/circum)
    
    
    fig, ax1 = plt.subplots()
    ax1.plot(x_wavenumbers,color='r',label='Zonal')
    # ax1.plot(y_wavenumbers,color='b',label='Meridional')

    ax1.set_xlabel('Time [days]')
    ax1.set_ylabel('Wavenumber')
    plt.legend()
    
    ax2 = ax1.twinx()
    ax2.plot(zmzw,color='k',label='ZMZW')
    ax2.set_ylabel('m/s')
    
    plt.legend()
    plt.show()
    
    
    plt.plot(c_phase,color='b',label='Phase vel')
    plt.title('Rossby wave phase velocity at %sN'%(lat))
    plt.xlabel('Time [days]')
    plt.ylabel('Velocity [m/s]')
    plt.legend()
    plt.show()
    
    # plt.plot(distance)
    # plt.title('Path travelled by cyclone at %sN'%(lat))
    # plt.xlabel('Time [days]')
    # plt.ylabel('Cumulative distance [km]')
    # plt.show()
    

def hbarotropic(cubes,start=500,end=600,level=8,omega=7.93e-06,radius=716000):
    
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'y_wind':
            y_wind = cube[start:end,:,:,:].copy()
            
            
    x_wind.coord('latitude').coord_system = GeogCS(radius)
    x_wind.coord('longitude').coord_system = GeogCS(radius)
        
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
    lats = np.round(x_wind.coord('latitude').points,0)
    km_heights = np.round(x_wind.coord('level_height').points*1e-03,2)
    
    winds = windspharm.iris.VectorWind(x_wind[:,level,:,:],y_wind[:,level,:,:])
    # Create a VectorWind data object from the x and y wind cubes
    uchi,vchi,upsi,vpsi = winds.helmholtz(truncation=21)
    # Calculate the Helmholtz decomposition. Truncation is set to 21 because 
    # this is what Hammond and Lewis 2021 used.
    
    zonal_upsi = upsi.collapsed('longitude', iris.analysis.MEAN)
    zonal_vpsi = vpsi.collapsed('longitude', iris.analysis.MEAN)
    # Calculate zonal means of the x and y components of the rotational component
    eddy_upsi = upsi - zonal_upsi
    eddy_vpsi = vpsi - zonal_vpsi
    magnitude = np.sqrt(eddy_upsi.data**2 + eddy_vpsi.data**2)
    
    fft2 = sp.fft.fftshift(sp.fftpack.fft2(sp.fft.ifftshift(magnitude)))
    yfreqs = sp.fft.fftshift(sp.fft.fftfreq(fft2.shape[1],d=1./90))
    xfreqs = sp.fft.fftshift(sp.fft.fftfreq(fft2.shape[2],d=1./144))
    psd = np.abs(fft2)**2
    
    quadrant = psd[:,46:51,73:78]
    xfreqs = xfreqs[73:78]
    yfreqs = yfreqs[46:51]

    x_wavenumbers = []
    y_wavenumbers = []
    for time in range(0,psd.shape[0]):
        peak = unravel_index(quadrant[time,:,:].argmax(), quadrant[time,:,:].shape)
        x_wavenumbers.append(peak[1])
        y_wavenumbers.append(peak[0])
        
    x_wavenumbers, y_wavenumbers = np.array(x_wavenumbers)+1, np.array(y_wavenumbers)+1
    
    zmzw = x_wind[:,level,:].collapsed('longitude',iris.analysis.MEAN)
    zmzw = zmzw.data    
    
    lats_rad = lats*(np.pi/180) # Convert input latitude to radians
    beta = 2*omega*np.cos(lats_rad)/radius # Beta factor
    
    circum = 2*np.pi*radius*np.cos(lats_rad) # Circumference in meters at input latitude
    
    x_set = []
    for day in range(0,psd.shape[0]):
        lambda_x = circum/x_wavenumbers[day] # Wavelength in x-direction at input latitude
        x_set.append(lambda_x)
    
    print(np.array(x_set).shape,zmzw.shape)
    x_num = 2*np.pi/np.array(x_set)
    
    lambda_y = 2*np.pi*radius/y_wavenumbers # Wavelength in y-direction at input latitude
    y_num = 2*np.pi/lambda_y

    cphases = []    
    for latitude in range(0,lats.shape[0]):
        c_phase = (zmzw[:,latitude] - (float(beta[latitude])/(x_num[:,latitude]**2+y_num**2)))
        cphases.append(c_phase)
    cphases = np.array(cphases)
    
    plt.contourf(np.arange(start,end),np.array(lats),np.array(cphases),levels=redblu.N,cmap=redblu,norm=TwoSlopeNorm(0))
    plt.title('Barotropic Rossby wave phase speed, h=%s km' %km_heights[level])
    plt.xlabel('Time [days]')
    plt.ylabel('Latitude [degrees]')
    plt.colorbar()
    plt.show()
    
    totalpath = []    
    for ylat in range(0,lats.shape[0]):
        distance = np.cumsum(cphases[ylat,:]*60*60*24)*10e-03 # Convert m/s to m/day = distance travelled in a day
        deg_dist = (distance/circum[ylat])
        totalpath.append(deg_dist)
        
    plt.contourf(np.arange(start,end),np.array(lats),np.array(totalpath),levels=redblu.N,cmap=redblu,norm=TwoSlopeNorm(0))
    plt.title('Total distance travelled by strongest Rossby wave, h=%s km' %km_heights[level])
    plt.xlabel('Time [days]')
    plt.ylabel('Latitude [degrees]')
    plt.colorbar()
    plt.show()    
    
    
def resonance(cubes,start=0,end=-1,level=8):

    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'y_wind':
            y_wind = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'surface_net_downward_shortwave_flux':
            heat = cube[start:end,:,:,:].copy()

    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
    km_heights = np.round(x_wind.coord('level_height').points*1e-03,2)
    
    winds = windspharm.iris.VectorWind(x_wind[:,level,:,:],y_wind[:,level,:,:])
    # Create a VectorWind data object from the x and y wind cubes
    uchi,vchi,upsi,vpsi = winds.helmholtz(truncation=21)
    # Calculate the Helmholtz decomposition. Truncation is set to 21 because 
    # this is what Hammond and Lewis 2021 used.
    
    zonal_upsi = upsi.collapsed('longitude', iris.analysis.MEAN)
    zonal_vpsi = vpsi.collapsed('longitude', iris.analysis.MEAN)
    # Calculate zonal means of the x and y components of the rotational component
    eddy_upsi = upsi - zonal_upsi
    eddy_vpsi = vpsi - zonal_vpsi
    magnitude = np.sqrt(eddy_upsi.data**2 + eddy_vpsi.data**2)
    
    fft2 = sp.fft.fftshift(sp.fftpack.fft2(sp.fft.ifftshift(magnitude)))
    yfreqs = sp.fft.fftshift(sp.fft.fftfreq(fft2.shape[1],d=1./90))
    xfreqs = sp.fft.fftshift(sp.fft.fftfreq(fft2.shape[2],d=1./144))
    psd = np.abs(fft2)**2

    one_zero = psd[:,45,73]
    print(xfreqs[73],yfreqs[45])
    two_one = psd[:,46,74]
    two_two = psd[:,47,74]
    three_two = psd[:,47,75]
    one_one = psd[:,46,73]
    wave_sum = one_zero + two_one + two_two + three_two + one_one

    lats = x_wind.coord('latitude')
    lons = x_wind.coord('longitude')

    if lats.bounds == None:
        x_wind.coord('latitude').guess_bounds()
    if lons.bounds == None:
        x_wind.coord('longitude').guess_bounds()

    grid_areas = iris.analysis.cartography.area_weights(x_wind[:,level,:,:])
    global_mean = x_wind[:,level,:,:].collapsed(['latitude','longitude'],iris.analysis.MEAN, weights=grid_areas)

    fig, ax1 = plt.subplots()
    ax1.set_xlabel('Time [days]')
    ax1.set_ylabel('Wind [m/s]')
    ax1.plot(global_mean.data,color='k',label='Wind')

    ax2 = ax1.twinx()
    ax2.set_ylabel('PSD')
    ax2.plot(one_zero,color='r',label='1-0 wave')
    ax2.plot(wave_sum,color='b',label='Wave sum')
    plt.legend()

    plt.title('Global mean zonal wind and PSD of 1-0 Rossby wave at h=%s km' %km_heights[level])
    plt.show()
