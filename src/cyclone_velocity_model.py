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
from iris.analysis import calculus
from matplotlib.colors import TwoSlopeNorm



redblu = mpl_cm.get_cmap('RdBu')
plasma = mpl_cm.get_cmap('plasma')



def model_vel(cubes,startlon=30,start=218,end=250,nlat=90,nlon=144,level=8,omega=1.19e-05,g=9.12,radius=5797818,lat=80,save='no'):
    
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube[start:end,:,:,:].copy()
            longterm_x_wind = cube.copy()
        if cube.standard_name == 'y_wind':
            y_wind = cube[start:end,:,:,:].copy()
        if cube.standard_name =='air_potential_temperature':
            theta = cube[start:end,:,:,:].copy()
        if cube.standard_name =='air_pressure':
            pressure = cube[start:end,:,:,:].copy()
            
    cphase_list = []
    cphase_hlist = []
    time_axis = np.arange(start, end)
    
    longterm_zmzw = longterm_x_wind[:,level,lat,0:72].collapsed(['longitude','time'], iris.analysis.MEAN)
    print(longterm_zmzw.data)
    
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
    km_heights = np.round(x_wind.coord('level_height').points*1e-03,2)
    latitudes = x_wind.coord('latitude').points
    
    winds = windspharm.iris.VectorWind(x_wind[:,level,:,:],y_wind[:,level,:,:])
    # Create a VectorWind data object from the x and y wind cubes
    
    lat_deg = int(latitudes[lat]) # Convert input row number to latitude in degrees
    print(lat_deg)
    shortterm_zmzw = x_wind[:,level,lat,0:72].collapsed('longitude',iris.analysis.MEAN)
    zmzw = shortterm_zmzw.data
    
    lat_rad = lat_deg*(np.pi/180) # Convert input latitude to radians
    beta = 2*omega*np.cos(lat_rad)/radius # Beta factor    
    circum = 2*np.pi*radius*np.cos(lat_rad) # Circumference in meters at input latitude   
    x_num = 2*np.pi/circum
    
    d_theta = iris.analysis.calculus.differentiate(theta, 'level_height')
    bv_freq = np.mean(np.sqrt(np.abs((g/theta[:,:-1,:,:].data)*d_theta.data)), axis=-1)
    Ld = (bv_freq[:,level,lat]*6800)/(2*omega*np.sin(lat_rad))
    
    wave_num = beta + zmzw*((1./Ld)**2)
    wave_denom = x_num**2 + (1./Ld)**2
    wave_component = wave_num/wave_denom
    
    c_phase = zmzw - wave_component                               
    cphase_list.append(c_phase)
    c_phase_h = (zmzw - (beta/(x_num**2))) 
    cphase_hlist.append(c_phase_h)

    distance = np.cumsum(c_phase*60*60*24) # Convert m/s to m/day = distance travelled in a day
    deg_dist = distance*(360/circum)+startlon
    stationary = beta/(x_num**2)
    print(deg_dist)
    
    markers_on = [0, 5, 10, 15]
    plt.plot(time_axis,c_phase, color='b', label='Phase vel')
    plt.plot(time_axis, c_phase, 'o', color='r', markevery=markers_on)
    plt.plot(time_axis, zmzw, color='r', label='ZMZW')
    plt.plot(time_axis, wave_component, color='r', linestyle='dashed', label='Wave comp')
    plt.plot(time_axis, longterm_zmzw.data*np.ones_like(zmzw), color='g', label='Longterm ZMZW')
#    plt.plot(np.array([0,5,10,15]), np.array([c_phase[0], c_phase[5], c_phase[10], c_phase[15]]), 'r' )
#    plt.plot(np.array([0, 5, 10, 15]),np.array([c_phase[0], c_phase[5], c_phase[10], c_phase[15]]), marker='o', mfc='r')
    plt.title('Rossby wave phase velocity at %sN'%(lat_deg))
    plt.xlabel('Time [days]')
    plt.ylabel('Velocity [m/s]')
    plt.legend(fontsize='small')
    plt.show()
    
    plt.plot(time_axis, deg_dist)
    plt.plot(time_axis, deg_dist,'o', color='r', markevery=markers_on)
    plt.title('Path travelled by cyclone at %sN'%(lat_deg))
    plt.xlabel('Time [days]')
    plt.ylabel('Cumulative distance [deg lon]')
    plt.show()
    

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
