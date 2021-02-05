# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 14:34:14 2021

@author: Mo Cohen

Pipeline for post-processing UM output data.
- Analyses model dynamics

Model: University of Exeter Stand Alone model of Proxima Centauri b, UM vn11.8

"""
import os, iris
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import iris.plot as iplt
import iris.coords
from iris.coords import DimCoord
import numpy as np
import scipy as sp
import windspharm
# Import packages


brewer_red = mpl_cm.get_cmap('brewer_YlOrRd_09')
brewer_bg = mpl_cm.get_cmap('brewer_PuBu_09')
brewer_reds = mpl_cm.get_cmap('brewer_Reds_09')
brewer_redblu = mpl_cm.get_cmap('brewer_RdBu_11')
# Load colormaps to use in plots

def plot_zonal_wind(cubes, time_slice=-1):
    
    """ Inputs: Iris CubeList and model level you want to look at
        Plots dayside, nightside, and global zonal mean winds     """
        
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube.copy()
    
    dayside = x_wind.extract(iris.Constraint(longitude=lambda v: 270 < v <= 359 or 0 <= v <= 90, latitude=lambda v: -90 <= v <= 90))
    nightside = x_wind.extract(iris.Constraint(longitude=lambda v: 90 < v <= 270, latitude=lambda v: -90 <= v <= 90))
    
    dayside_zonal_mean = dayside.collapsed('longitude', iris.analysis.MEAN)
    nightside_zonal_mean = nightside.collapsed('longitude', iris.analysis.MEAN)
    
    CS_day = iplt.contourf(dayside_zonal_mean[time_slice,:,:], brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Dayside Zonal Mean Zonal Wind [m s-1]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.clabel(CS_day, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    plt.show()
    
    CS_night = iplt.contourf(nightside_zonal_mean[time_slice,:,:], brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Nightside Zonal Mean Zonal Wind [m s-1]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.clabel(CS_night, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    plt.show()


def plot_streamfunction(cubes, level=14, time_slice=-1, omega=0.64617667):
    
    """ Uses windspharm package to plot the streamfunction       """
    
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube.copy()
        if cube.standard_name == 'y_wind':
            y_wind = cube.copy()
        if cube.standard_name == 'upward_air_velocity':
            z_wind = cube.copy()    
    
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())

    wind = windspharm.iris.VectorWind(x_wind, y_wind)
    streamfunction, velpotential = wind.sfvp()
    clevs = [-200, -180, -160, -120, -100, -80, -60, -40, -20, 0, 40, 80, 120, 160, 200]
    iplt.contourf(streamfunction[time_slice,level,:,:]*1e-06, clevs, cmap=brewer_redblu, extend='both')
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.colorbar(orientation='horizontal')
    plt.title('Streamfunction [$10^6$ m2 s-1], h = %s km' %(level+1), y=1.20)
    plt.show()

    planet_vort = wind.planetaryvorticity(omega=omega)
    relative_vort = wind.vorticity()
    absolute_vort = wind.absolutevorticity()
    
    iplt.contourf(relative_vort[time_slice,level,:,:], brewer_bg.N, cmap=brewer_bg)
    plt.title('Relative Vorticity, h = %s km' %(level+1), y=1.20)
    plt.ylabel('Latitude [degrees]')
    plt.xlabel('Longitude [degrees]')
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.colorbar(orientation='horizontal')
    plt.show()

        

def plot_streamlines(cubes, level=14, time_slice=-1):
    
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
    
    iplt.contourf(x_wind[time_slice,level,:,:], brewer_reds.N, cmap=brewer_reds)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Zonal Wind Speed, h=%s km [m s-1]' %(level+1), y=1.20)
    plt.colorbar(pad=0.1)
    plt.show()

    X,Y = np.meshgrid(np.arange(-72,72), np.arange(-45,45))
    fig = plt.figure(figsize = (12, 7)) 
    strm = plt.streamplot(X, Y, np.roll(x_wind[time_slice,level,:,:].data, 72, axis=1), np.roll(y_wind[time_slice,level,:,:].data, 72, axis=1), density = 0.5, color=np.roll(speed[time_slice,level,:,:].data, 72, axis=1), cmap=brewer_reds)
    # Since .data method extracts the numpy array and strips the metadata, the longitude/latitude information is lost.
    # To align plot so that (0,0) is at the center as in the Iris plots, use numpy.roll to shift columns (longitude) 180 degrees (72 places)
    fig.colorbar(strm.lines)
    plt.title('Wind speed and direction [m s-1], h=%s km' %(level+1))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.xticks((-72,-52,-32,-12,0,12,32,52,72),('180W','140W','100W','60W','0','60E','100E','140E','180E'))
    plt.yticks((-45,-30,-15,0,15,30,45),('90S','60S','30S','0','30N','60N','90N'))    
    plt.show() 


def plot_zwind(cubes, time_slice=-1):
    
    for cube in cubes:
        if cube.standard_name == 'upward_air_velocity':
            z_wind = cube.copy()
            
    dayside = z_wind.extract(iris.Constraint(longitude=lambda v: 270 < 359 or 0 <= v <= 90, latitude=lambda v: -90 <= v <= 90))
    nightside = z_wind.extract(iris.Constraint(longitude=lambda v: 90 < v <= 270, latitude=lambda v: -90 <= v <= 90))
    
    dayside_zonal_mean = dayside.collapsed('longitude', iris.analysis.MEAN)
    nightside_zonal_mean = nightside.collapsed('longitude', iris.analysis.MEAN)
    
    CS_day = iplt.contourf(dayside_zonal_mean[time_slice,:,:], brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Dayside Zonal Mean Z Wind [m s-1]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
#    plt.clabel(CS_day, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    plt.show()
    
    CS_night = iplt.contourf(nightside_zonal_mean[time_slice,:,:], brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Nightside Zonal Mean Z Wind [m s-1]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
#   plt.clabel(CS_night, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    plt.show()

    dayside_meridional_mean = dayside.collapsed('latitude', iris.analysis.MEAN)
    nightside_meridional_mean = nightside.collapsed('latitude', iris.analysis.MEAN)
    
    CS_day = iplt.contourf(dayside_meridional_mean[time_slice,:,:], brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Dayside Meridional Mean Z Wind [m s-1]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Longitude [degrees]')
#    plt.clabel(CS_day, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    plt.show()
    
    CS_night = iplt.contourf(nightside_meridional_mean[time_slice,:,:], brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Nightside Meridional Mean Z Wind [m s-1]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Longitude [degrees]')
#    plt.clabel(CS_night, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    plt.show()