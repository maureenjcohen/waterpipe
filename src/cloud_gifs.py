#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 14 20:24:57 2022

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
import gif

brewer_bg = mpl_cm.get_cmap('brewer_PuBu_09')
brewer_red = mpl_cm.get_cmap('brewer_YlOrRd_09')
brewer_redblu = mpl_cm.get_cmap('brewer_RdBu_11')


@gif.frame
def cloud_frame(cubes, time_slice=-1, nlat=256, nlon=512, nlev=38, level=25, meaning=10, n=8, cloudtype='ice', fractype='mass'):

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
        
    fig, ax = plt.subplots(figsize=(10,5))
    plt.imshow(np.roll(meaned_cloud[time_slice,level, :,:]*(10**4),int(nlon/2),axis=1), cmap=brewer_bg)
    cbar = plt.colorbar()
    
    if fractype=='mass':
        cbar.ax.set_title('$10^{-4}$kg/kg')

    
    plt.quiver(X[::n,::n],Y[::n,::n], np.roll(meaned_x[time_slice,level,::n,::n],int(nlon/(2*n)),axis=1), 
                np.roll(-meaned_y[time_slice, level,::n,::n],int(nlon/(2*n)),axis=1),scale_units='xy',scale=5)

    plt.title('%s and horizontal wind, days %s to %s, h=%s km' %(titleterm, time_slice*meaning-meaning,time_slice*meaning, heights[level]))
    # plt.xticks((0,12,24,36,48,60,72,84,96,108,120,132,144),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    # plt.yticks((90,75,60,45,30,15,0),('90S','60S','30S','0','30N','60N','90N'))   
    plt.xticks((0,128,256,384,512),('180W','90W','0','90E','180E'))
    plt.yticks((256,128,0),('90S','0','90N')) 
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
        
    # ax2.quiverkey(q1, X=0.9, Y=1.05, U=10, label='10 m/s', labelpos='E', coordinates='axes')
    
frames = []
for i in range(29):
    frame = cloud_frame(quad,time_slice=i)
    frames.append(frame)

gif.save(frames,'/exports/csce/datastore/geos/users/s1144983/um_data/cloudproject/gifs_quad/level.gif', duration = 30, unit = 's', between='startend')
        