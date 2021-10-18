#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 15:17:12 2021

@author: Maureen Cohen
"""

import iris
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
from matplotlib import ticker
from matplotlib.colors import TwoSlopeNorm
from iris.analysis import calculus

reds = mpl_cm.get_cmap('Reds')


def plot_number(cubes, radius=7160000, omega=0.64617667e-05, g=10.9, start=1080, end=1200, time_slice=-1, lat_slice=45, level=47):
    """ Level 47 is 41 km, level 20 is 8 km"""
    
    for cube in cubes:
        if cube.standard_name == 'eastward_wind' or cube.standard_name == 'x_wind':
            x_wind = cube[start:end,:-5,25:66,:].copy()
        if cube.standard_name == 'northward_wind' or cube.standard_name == 'y_wind':
            y_wind = cube[start:end,:-5,25:66,:].copy()
        # if cube.standard_name =='air_potential_temperature':
        #     theta = cube[start:end,:,:,:].copy()
    
    vertical, latitudes, longitudes = x_wind.shape[1], x_wind.shape[2], x_wind.shape[3]
    x_axis, y_axis = x_wind.coord('longitude').points, x_wind.coord('latitude').points
    
    for coord in x_wind.coords():
        if coord.long_name == 'Hybrid height':
            heights = np.round(x_wind.coord('Hybrid height').points*1e-03,0)
            tcoord = 't'
            hcoord = 'Hybrid height'
        elif coord.long_name == 'level_height':
            heights = np.round(x_wind.coord('level_height').points*1e-03,0)
            tcoord = 'time'
            hcoord = 'level_height'
    
    # d_theta = iris.analysis.calculus.differentiate(theta, 'Hybrid height')
    # bv_freq = np.sqrt(np.abs((g/theta[:,:-1,:,:].data)*d_theta.data))
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
    
    # x_wind = x_wind.collapsed(tcoord, iris.analysis.MEAN) 
    # y_wind = y_wind.collapsed(tcoord, iris.analysis.MEAN)
     
    x_diff = iris.analysis.calculus.differentiate(x_wind, tcoord)
    y_diff = iris.analysis.calculus.differentiate(y_wind, tcoord)
    x_grad = iris.analysis.calculus.differentiate(x_wind, 'longitude')
    y_grad = iris.analysis.calculus.differentiate(y_wind, 'latitude')
    print(np.min(y_diff.data), np.max(y_diff.data))

    
    mag = np.sqrt(x_wind.data**2 + y_wind.data**2)
    lats = x_wind.coord('latitude').points*(np.pi/180)
    # u_bar = x_wind.collapsed(['longitude', hcoord, tcoord], iris.analysis.MEAN)
    # f = 2*np.sin(lats)*(u_bar.data/(radius*np.cos(lats)))
    f10 = 2*omega*np.sin(lats[25])
    # f10 = 2*np.sin(lats[49])*(u_bar[49].data/(radius*np.cos(lats[49])))
    f = 2*omega*np.sin(lats)
    for i in range(15,25):
        f[i] = f10
    # print(f)

    f = np.repeat(f[:,np.newaxis], (end-start), axis=0)
    f = np.reshape(f,((end-start),latitudes))    
    f = np.repeat(f[:,np.newaxis,:], vertical, axis=1)
    f = np.repeat(f[:,:,:,np.newaxis], longitudes, axis=3)
    # f = 2*np.sin(lats)*(u_bar.data/(radius*np.cos(lats)))
    # f = np.repeat(f[:,:,:,np.newaxis], 144, axis=3)
    print(f.shape)
    
    # numerator = np.sqrt(x_diff[:,:,:-1,:].data**2 + y_diff[:,:,:-1,:].data**2 + (x_wind[:-1,:,:-1,:].data*x_grad[:-1,:,:-1,:].data)**2 + (y_wind[:-1,:,:-1,:].data*y_grad[:-1,:,:,:].data)**2)
    numerator = np.sqrt((x_wind[:-1,:,:-1,:].data*x_grad[:-1,:,:-1,:].data)**2 + (y_wind[:-1,:,:-1,:].data*y_grad[:-1,:,:,:].data)**2)
    print(np.min(numerator), np.max(numerator))
    # numerator = np.sqrt(x_diff.data**2 + y_diff.data**2)
    denominator = np.abs(f[:-1,:,:-1,:]*mag[:-1,:,:-1,:])
    print(np.min(denominator), np.max(denominator))
    rol = numerator/denominator
    # rol = np.log(numerator/denominator)

    # rol = rol/np.mean(rol)
    # print(rol.shape)
    rol = np.mean(rol, axis=0)
    rol = rol/np.min(rol)
    # rol = (rol - np.min(rol))/(np.max(rol) - np.min(rol))
    # print(np.where(rol==1.0))
    print(np.min(rol), np.max(rol))
        
    plt.figure(figsize=(10,5))
    plt.contourf(np.arange(-longitudes/2, longitudes/2), np.array(heights), np.roll(rol[:,lat_slice,:], 72, axis=1), np.linspace(0,60,30), cmap=reds)
    plt.title('Mean Lagrangian Rossby Number at Equator, t=%s to %s days' %(start/4, end/4))
    # plt.title('Mean Lagrangian Rossby Number at Equator')
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.ylabel('Height [km]')
    cbar = plt.colorbar(pad=0.1)
    cbar.locator = ticker.AutoLocator()
    cbar.update_ticks()
    # cbar.set_ticklabels(['0','10','20','30','40','50','60'])
    # plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/laso/epsfigs/lrn_equator_normed.eps', format='eps')  
    plt.show()
    
    
    plt.figure(figsize=(10,5))
    plt.contourf(np.arange(-longitudes/2, longitudes/2),  y_axis[:-1], np.roll(rol[level,:,:], 72, axis=1), np.linspace(0,60,30), cmap=reds)
    plt.title('Mean Lagrangian Rossby Number at h=%s km, t=%s to %s days' %(heights[level], start/4, end/4))
    # plt.title('Mean Lagrangian Rossby Number at h=%s km' %(heights[level]))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.ylabel('Latitude [degrees]')
    plt.yticks((-40,-20,0,20,40), ('40S', '20S', '0', '20N','40N'))
    mbar = plt.colorbar(pad=0.1)
    mbar.locator = ticker.AutoLocator()
    mbar.update_ticks()
    # mbar.set_ticks([1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0])
    # mbar.set_ticklabels(['1','0.9', '0.8', '0.7', '0.6', '0.5', '0.4', '0.3', '0.2', '0.1', '0.0'])
    plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/laso/epsfigs/lrn_%s_normed.eps' %(heights[level]), format='eps')  
    plt.show()
    
    

