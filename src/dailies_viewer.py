#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 17:35:42 2021

@author: Mo Cohen
"""
import iris
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm



def plot_all(cubes, lat=45, long=36, level=25):
    
    for cube in cubes:
        # if cube.shape[0] == 1800:
        if cube.standard_name == None:
            title = cube.long_name
        else:
            title = cube.standard_name
            
        
        time_axis = cube.shape[0]
        units = cube.units
        
        if len(cube.shape) == 4:
            data = cube[:,level,lat,long].data
            heights = np.round(cube.coord('level_height').points*1e-03,2)
            print('at %s' %heights[level])
        else:
            data = cube[:,lat,long].data
        
        plt.plot(np.arange(0,time_axis), data)
        plt.title('%s' %title)
        plt.ylabel(units)
        plt.xlabel('days')
        plt.show()
        

def fourier(cubes, long1=36, long2=108, periodicity=False):
    
    for cube in cubes:
        if cube.long_name == 'cloud_area_fraction_assuming_maximum_random_overlap':
            cloud_cover = cube.copy()

    cloud_array = cloud_cover.data    
    clear_list = []
    for day in range(0,cloud_cover.shape[0]):
        data = (cloud_array[day,:,long1] + cloud_array[day,:,long2])/2
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