#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 19 15:16:33 2021

@author: Mo Cohen
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.cm as mpl_cm
import iris, cartopy
import cartopy.crs as ccrs
from matplotlib.colors import TwoSlopeNorm

import windspharm
import gif


hot = mpl_cm.get_cmap('hot')
blues = mpl_cm.get_cmap('Blues_r')
redblu = mpl_cm.get_cmap('bwr')

savepath = '/exports/csce/datastore/geos/users/s1144983/um_data/control/gifs'


@gif.frame
def sphere_qbo(cubes, level=47, time_slice=-1, long=0, lat=0):
    
    for cube in cubes:
        if cube.standard_name == 'x_wind' or cube.standard_name == 'eastward_wind':
            x_wind = cube.copy()
            
    longs = np.linspace(-180,180,144)
    lats = np.linspace(-90,90,90)
    # x_wind = x_wind[:,:,25:65,:]
    
    ortho = ccrs.Orthographic(central_longitude=long,central_latitude=lat)
    
    fig = plt.figure()
    # fig, ax0 = plt.subplots(figsize=(5,5))
    # ax0.imshow(image)
    fig.patch.set_facecolor('black')
    ax1 = plt.axes(projection=ortho)
    ax1.set_global()
    
    # ax.gridlines()
    ax1.contourf(longs, lats, np.roll(x_wind[time_slice,level, :,:].data, 72, axis=1), transform=ccrs.PlateCarree(), cmap=redblu, norm=TwoSlopeNorm(0))

    # plt.show()


temp_frames = []
for i in range(0,800):
    temp_frame = sphere_qbo(qbo, time_slice=i)
    temp_frames.append(temp_frame)    

gif.save(temp_frames, str(savepath) + '/qbo_200days.gif', duration = 200, unit = 's', between='startend')

    

# def sphere_streamlines(cubes, level=47, time_slice=-1, long=36, lat=0):
    
#     """ Inputs: Iris CubeList and model level you want to look at
#         wpharm=True allows windspharm library to be used on Linux machines
#         omega is the rotation rate of the planet in 10^-5 rad/s. Default is Proxima Centauri b.
        

#    """
        
#     for cube in cubes:
#         if cube.standard_name == 'x_wind' or cube.standard_name == 'eastward_wind':
    #         x_wind = cube.copy()
    #     if cube.standard_name == 'y_wind' or cube.standard_name == 'northward_wind':
    #         y_wind = cube.copy()
 
    # y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
   
    # speed = iris.analysis.maths.apply_ufunc(np.sqrt, (x_wind**2 + y_wind**2))
    # ortho = ccrs.Orthographic(central_longitude=long,central_latitude=lat)

    # X,Y = np.meshgrid(np.linspace(-180,180,144), np.linspace(-90,90,90))
    # x = x_wind[time_slice, level, :, :].data
    # y = y_wind[time_slice, level, :, :].data
    
    # fig, ax = plt.subplots(subplot_kw={'projection': ortho})
    # fig.set_size_inches([8,8])
    # ax.set_global()
    # ax.streamplot(X,Y,x,y, transform=ccrs.PlateCarree())
    # strm = plt.streamplot(X, Y, np.roll(x_wind[time_slice,level,:,:].data, 72, axis=1), np.roll(y_wind[time_slice,level,:,:].data, 72, axis=1), density = 0.5, color=np.roll(speed[time_slice,level,:,:].data, 72, axis=1), cmap=hot)
    # ax.contour(long, lat, strm, transform=ccrs.PlateCarree(), cmap=hot)

    # plt.show()
    


def sphere_rossby(cubes, time_slice=-1, long=36, lat=85):
    
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube.copy()
        if cube.standard_name == 'y_wind':
            y_wind = cube.copy()
 
    ortho = ccrs.Orthographic(central_longitude=long,central_latitude=lat)
 
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
    
    winds = windspharm.iris.VectorWind(x_wind,y_wind)
    
    eta = winds.absolutevorticity()
    div = winds.divergence()
    uchi, vchi = winds.irrotationalcomponent()
    etax, etay = winds.gradient(eta)
    etax.units = 'm**-1 s**-1'
    etay.units = 'm**-1 s**-1'
    
    S = eta*-1.*div-(uchi*etax+vchi*etay)
    S.coord('longitude').attributes['circular'] = True
    
    fig = plt.figure()

    fig.patch.set_facecolor('black')
    ax1 = plt.axes(projection=ortho)
    ax1.set_global()
    
    ax1.contourf(long, lat, np.roll(S[-1,:,:].data, 72, axis=1), transform=ccrs.PlateCarree(), cmap=redblu)


    clevs = [-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30]
    iplt.contourf(S[time_slice,level,:,:]*1e11,clevs,cmap=brewer_redblu, extend='both')
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.colorbar(orientation='horizontal')
    plt.title('Rossby Wave Source ($10^{-11}$s$^{-1}$)')
    plt.show()