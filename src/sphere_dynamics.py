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
import windspharm


hot = mpl_cm.get_cmap('hot')
blues = mpl_cm.get_cmap('Blues_r')
redblu = mpl_cm.get_cmap('RdBu')



def sphere_streamlines(cubes, level=55, time_slice=-1, long=36, lat=0):
    
    """ Inputs: Iris CubeList and model level you want to look at
        wpharm=True allows windspharm library to be used on Linux machines
        omega is the rotation rate of the planet in 10^-5 rad/s. Default is Proxima Centauri b.
        
        Plot zonal wind speed at level
        Plot streamplot at level
        Plot zonal mean zonal winds for planet
        
        If wpharm=True, plot stream function and relative vorticity    """
        
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube.copy()
        if cube.standard_name == 'y_wind':
            y_wind = cube.copy()
 
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
   
    speed = iris.analysis.maths.apply_ufunc(np.sqrt, (x_wind**2 + y_wind**2))
    ortho = ccrs.Orthographic(central_longitude=long,central_latitude=lat)

    X,Y = np.meshgrid(np.arange(-72,72), np.arange(-44,44))
    x = x_wind[time_slice, level, 0:88, :].data
    y = y_wind[time_slice, level, 0:88, :].data
    print(X.shape, Y.shape, x.shape, y.shape)
    
    fig, ax = plt.subplots(subplot_kw={'projection': ortho})
    fig.set_size_inches([8,8])
    ax.set_global()
    ax.streamplot(X,Y,x,y, transform=ccrs.PlateCarree())
    # strm = plt.streamplot(X, Y, np.roll(x_wind[time_slice,level,:,:].data, 72, axis=1), np.roll(y_wind[time_slice,level,:,:].data, 72, axis=1), density = 0.5, color=np.roll(speed[time_slice,level,:,:].data, 72, axis=1), cmap=hot)
    # ax.contour(long, lat, strm, transform=ccrs.PlateCarree(), cmap=hot)

    plt.show()
    


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


    # clevs = [-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30]
    # iplt.contourf(S[time_slice,level,:,:]*1e11,clevs,cmap=brewer_redblu, extend='both')
    # ax = plt.gca()
    # ax.gridlines(draw_labels=True)
    # plt.colorbar(orientation='horizontal')
    # plt.title('Rossby Wave Source ($10^{-11}$s$^{-1}$)')
    # plt.show()