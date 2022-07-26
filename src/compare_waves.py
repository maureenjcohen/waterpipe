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
    
def compare_waves(cubes,start=0,end=-1,level=8):

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

    one_zero = psd[:,45,73]
    one_one = psd[:,46,73]
    two_zero = psd[:,45,74]
    two_one = psd[:,46,74]
    two_two = psd[:,47,74]
    three_one = psd[:,46,75]
    three_two = psd[:,47,75]

    lats = x_wind.coord('latitude')
    lons = x_wind.coord('longitude')

    if lats.bounds == None:
        x_wind.coord('latitude').guess_bounds()
    if lons.bounds == None:
        x_wind.coord('longitude').guess_bounds()

    grid_areas = iris.analysis.cartography.area_weights(x_wind[:,level,:,:])
    global_mean = x_wind[:,level,:,:].collapsed(['latitude','longitude'],iris.analysis.MEAN, weights=grid_areas)

    fig, ax = plt.subplots()
    ax.plot(one_zero,color='r',label='1-0 wave')
    ax.plot(one_one,color='b',label='1-1 wave')
    ax.plot(two_one,color='g',label='2-1 wave')
    ax.plot(two_two,color='k',label='2-2 wave')
    ax.set_ylabel('PSD')
    ax.set_xlabel('Time [days]')
    plt.legend()

    plt.title('PSD of Rossby waves at h=%s km' %km_heights[level])
    plt.show()
