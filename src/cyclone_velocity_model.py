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
    
    quadrant = psd[:,46:51,73:78]
    xfreqs = xfreqs[73:78]
    yfreqs = yfreqs[46:51]

    x_wavenumbers = []
    y_wavenumbers = []
    for time in range(0,psd.shape[0]):
        peak = unravel_index(quadrant[time,:,:].argmax(), quadrant[time,:,:].shape)
        x_wavenumbers.append(peak[1])
        y_wavenumbers.append(peak[0])
        
    x_wavenumbers, y_wavenumbers = np.array(x_wavenumbers), np.array(y_wavenumbers)
    
    zmzw_lat = np.abs(lat-45)*2
    zmzw = x_wind[:,level,zmzw_lat].collapsed('longitude',iris.analysis.MEAN)
    zmzw = zmzw.data
        
    lat_rad = lat*(np.pi/180) # Convert latitude to radians
    beta = 2*omega*np.cos(lat_rad)/radius # Beta factor
    
    circum = 2*np.pi*radius*np.cos(lat_rad) # Circumference in meters at input latitude
    
    lambda_x = circum/x_wavenumbers # Wavelength in x-direction at input latitude
    x_num = 2*np.pi/lambda_x
    
    lambda_y = 2*np.pi*radius/y_wavenumbers # Wavelength in y-direction at input latitude
    y_num = 2*np.pi/lambda_y
    
    c_phase = (zmzw - (beta/(x_num**2+y_num**2)))
    c_group = (zmzw + (beta*(x_num**2-y_num**2))/((x_num**2+y_num**2)**2))
    
    distance = c_phase*60*60*24
    
    plt.plot(c_phase,color='b')
    plt.plot(c_phase+zmzw,color='r')
    plt.show()