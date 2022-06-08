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
import scipy as sp
import windspharm


brewer_bg = mpl_cm.get_cmap('brewer_PuBu_09')
brewer_red = mpl_cm.get_cmap('brewer_YlOrRd_09')
brewer_redblu = mpl_cm.get_cmap('brewer_RdBu_11')
plasma = mpl_cm.get_cmap('plasma')



@gif.frame
def cloud_frame(cubes, time_slice=-1, nlat=90, nlon=144, nlev=38, level=8, meaning=5, n=4, cloudtype='ice', fractype='mass'):

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
    plt.imshow(np.roll(meaned_cloud[time_slice,level, :,:],int(nlon/2),axis=1), cmap=brewer_bg, vmax=0.0002,vmin=0.0)
    cbar = plt.colorbar()
    
    # if fractype=='mass':
    #     cbar.ax.set_title('$10^{-4}$kg/kg')

    
    plt.quiver(X[::n,::n],Y[::n,::n], np.roll(meaned_x[time_slice,level,::n,::n],int(nlon/(2*n)),axis=1), 
                np.roll(-meaned_y[time_slice, level,::n,::n],int(nlon/(2*n)),axis=1),scale_units='xy',scale=5)

    plt.title('%s and horizontal wind, days %s to %s, h=%s km' %(titleterm, time_slice*meaning-meaning,time_slice*meaning, heights[level]))
    plt.xticks((0,12,24,36,48,60,72,84,96,108,120,132,144),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.yticks((90,75,60,45,30,15,0),('90S','60S','30S','0','30N','60N','90N'))   
    # plt.xticks((0,128,256,384,512),('180W','90W','0','90E','180E'))
    # plt.yticks((256,128,0),('90S','0','90N')) 
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
        
    # ax2.quiverkey(q1, X=0.9, Y=1.05, U=10, label='10 m/s', labelpos='E', coordinates='axes')
    
frames = []
for i in range(100,130):
    frame = cloud_frame(fartrap,time_slice=i)
    frames.append(frame)

gif.save(frames,'/exports/csce/datastore/geos/users/s1144983/um_data/cloudproject/gifs_trapfar_80km/level8ice.gif', duration = 29, unit = 's', between='startend')


# @gif.frame
# def ps_gif(cubes, time_slice=0, level=8):  
    
#     for cube in cubes:
#         if cube.standard_name == 'x_wind':
#             x_wind = cube[time_slice,:,:,:].copy()
#         if cube.standard_name == 'y_wind':
#             y_wind = cube[time_slice,:,:,:].copy()

#     # Select the cubes we want from the cube list
 
#     y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
    
#     height = [('level_height', x_wind.coord('level_height').points)]
#     # Regrid so that all three cubes are on the same x, y, z grid. Uses the x_wind as reference for the others
    
#     km_heights = np.round(x_wind.coord('level_height').points*1e-03,2)
#     # Extract pressure and kilometer values for labeling the plots
    
#     winds = windspharm.iris.VectorWind(x_wind,y_wind)
#     # Create a VectorWind data object from the x and y wind cubes
#     uchi,vchi,upsi,vpsi = winds.helmholtz(truncation=21)
#     # Calculate the Helmholtz decomposition. Truncation is set to 21 because 
#     # this is what Hammond and Lewis 2021 used.
    
#     zonal_upsi = upsi.collapsed('longitude', iris.analysis.MEAN)
#     zonal_vpsi = vpsi.collapsed('longitude', iris.analysis.MEAN)
#     # Calculate zonal means of the x and y components of the rotational component
#     eddy_upsi = upsi - zonal_upsi
#     eddy_vpsi = vpsi - zonal_vpsi
#     magnitude = np.sqrt(eddy_upsi.data**2 + eddy_vpsi.data**2)
#     fft2 = sp.fft.fftshift(sp.fftpack.fft2(sp.fft.ifftshift(magnitude[level,:,:])))
#     yfreqs = sp.fft.fftshift(sp.fft.fftfreq(fft2.shape[0],d=1./90))
#     xfreqs = sp.fft.fftshift(sp.fft.fftfreq(fft2.shape[1],d=1./144))
#     psd = np.abs(fft2)**2
    
#     fig2,ax2 = plt.subplots(figsize = (10,5))
#     im = ax2.contourf(xfreqs[73:78],yfreqs[46:51],psd[46:51,73:78], cmap=plasma)
#     plt.title('Power spectrum of eddy rotational wind magnitude, day %s' %(time_slice))
#     plt.xlabel('Zonal wavenumber')
#     plt.ylabel('Meridional wavenumber')
#     fig2.colorbar(im)

# ps_frames = []
# for i in range(280,350):
#     ps_frame = ps_gif(trap,time_slice=i)
#     ps_frames.append(ps_frame)

# gif.save(ps_frames,'/exports/csce/datastore/geos/users/s1144983/um_data/cloudproject/gifs_trapcontrol_80km/pslevel8_zoom.gif', duration = 35, unit = 's', between='startend')
