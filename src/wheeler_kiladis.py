#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 16:56:43 2021

@author: Maureen Cohen

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

brewer_redblu = mpl_cm.get_cmap('RdBu_r')
brewer_reds = mpl_cm.get_cmap('brewer_Reds_09')



def wk_filter(arg):
    
    """ 1-2-1 filter for smoothing Fourier coefficients in the Wheeler-Kiladis diagram
        Three neighbouring values are weighted 1*, 2*, 1*, added together, and divided by 4 """
        
    out = np.zeros_like(arg)
    out[0] = (3*arg[0] + arg[1])/4
    out[-1] = (arg[-2] + 3*arg[-1])/4
    out[1:-1] = (2*arg[1:-1]+ arg[0:-2] + arg[2:])/4
    return out

def wheeler_kiladis(cubes, omega=0.64617667e-05, radius=7160000, start=0, end=160, level=38, lat=45, smooth=400):
    
    for cube in cubes:
        if cube.standard_name == 'eastward_wind':
            x_wind = cube[start:end,:,:,:].copy()

    longitudes = x_wind.shape[3]/2
    heights = np.round(x_wind.coord('Hybrid height').points*1e-3,0)
    time = np.arange(0,x_wind.shape[0])*0.25
    
    """ Filter zonal wind, remove longer wavelengths, calculate u-prime"""
    
    # u_mean = x_wind.collapsed('longitude', iris.analysis.MEAN)
    # u_prime = x_wind - u_mean
    
    x_data = x_wind.data
    
    u_prime_list = []
    
    for timex in range(0,x_wind.shape[0]):
        time_list = []
        
        for levelx in range(0,x_wind.shape[1]):
            level_list = []
            
            for latx in range(0,x_wind.shape[2]):
            
                u_fft = sp.fftpack.fft(x_data[timex,levelx,latx,:])
                u_psd = np.abs(u_fft)**2
                u_freq = sp.fftpack.fftfreq(len(u_psd), 1./144)
                
                highpass = u_fft.copy()
                highpass[np.abs(u_freq) < 5.1] = 0
                
                u_cleaned = np.real(sp.fftpack.ifft(highpass))
                # u_bar = np.mean(u_cleaned)
                # u_prime = u_cleaned - u_bar
                level_list.append(u_cleaned)
                
            time_list.append(level_list)
        u_prime_list.append(time_list)
        
    u_prime = np.array(u_prime_list)
    
    # nu_prime_list = []
    
    # for timey in range(0,x_wind.shape[0]):
    #     timey_list = []
        
    #     for levely in range(0,x_wind.shape[1]):
    #         levely_list = []
            
    #         for laty in range(0,x_wind.shape[2]):
    #             laty_list = []
                
    #             for longy in range(0,x_wind.shape[3]):
                    
    #                 nu_fft = sp.fftpack.fft(u_prime[:,levely,laty,longy])
    #                 nu_psd = np.abs(nu_fft)**2
    #                 nu_freq = sp.fftpack.fftfreq(len(nu_psd), 1./int(4))
                    
    #                 nu_highpass = nu_fft.copy()
    #                 nu_highpass[np.abs(nu_freq) < 0.25] = 0
                    
    #                 returned = np.real(sp.fftpack.ifft(nu_highpass))
                    
    #                 laty_list.append(returned)
    #             levely_list.append(laty_list)
    #         timey_list.append(levely_list)
    #     nu_prime_list.append(timey_list)
        
    # nu_prime = np.array(nu_prime_list)
    # print(nu_prime.shape)
        
    plt.figure(figsize=(10,5))    
    plt.contourf(np.arange(-longitudes, longitudes), time, np.roll(u_prime[:,level,lat,:].data, 72, axis=1), brewer_redblu.N, cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    plt.title('$U^{\prime}$ at Equator, h=%s km' %(heights[level]))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.ylabel('Time [days]')
    plt.colorbar(pad=0.1)
    plt.show()            

    
    """ Here you do the Fourier stuff"""
    
    data = u_prime[:,level,lat,:].data
    coeffs = np.abs(np.fft.fftshift(np.fft.fft2(data)))**2
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
    lats = x_wind.coord('latitude').points*(np.pi/180)
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
        
    plt.contourf(wavenumbers, freqs[pos], signal[pos,:], np.linspace(0,-3,18), cmap=brewer_reds)
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