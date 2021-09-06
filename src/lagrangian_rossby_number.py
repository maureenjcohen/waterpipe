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
from matplotlib.colors import TwoSlopeNorm
from iris.analysis import calculus

reds = mpl_cm.get_cmap('Reds')


def plot_number(cubes, radius=7160000, omega=0.64617667e-05, g=10.9, start=1080, end=1200, time_slice=-1, lat_slice=45, level=47):
    for cube in cubes:
        if cube.standard_name == 'eastward_wind':
            x_wind = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'northward_wind':
            y_wind = cube[start:end,:,:,:].copy()
        if cube.standard_name =='air_potential_temperature':
            theta = cube[start:end,:,:,:].copy()
    
    longitudes, latitudes = x_wind.shape[3], x_wind.shape[2]
    heights = x_wind.coord('Hybrid height').points
    h_label = np.round(heights*1e-03,0)
    
    # d_theta = iris.analysis.calculus.differentiate(theta, 'Hybrid height')
    # bv_freq = np.sqrt(np.abs((g/theta[:,:-1,:,:].data)*d_theta.data))
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())      
    x_diff = iris.analysis.calculus.differentiate(x_wind, 't')
    y_diff = iris.analysis.calculus.differentiate(y_wind, 't')
    x_grad = iris.analysis.calculus.differentiate(x_wind, 'longitude')
    y_grad = iris.analysis.calculus.differentiate(y_wind, 'latitude')
    
    mag = np.abs(np.sqrt(x_wind.data**2 + y_wind.data**2))
    lats = x_wind.coord('latitude').points*(np.pi/180)
    # u_bar = x_wind.collapsed(['longitude'], iris.analysis.MEAN)
    f = 2*omega*np.sin(lats)
    f = np.repeat(f[:,np.newaxis], (end-start), axis=0)
    f = np.reshape(f,((end-start),90))
    f = np.repeat(f[:,np.newaxis,:], 60, axis=1)
    f = np.repeat(f[:,:,:,np.newaxis], 144, axis=3)
    # f = 2*np.sin(lats)*(u_bar.data/(radius*np.cos(lats)))
    # f = np.repeat(f[:,:,:,np.newaxis], 144, axis=3)

    
    numerator = np.sqrt(x_diff.data[:,:,:-1,:]**2 + y_diff[:,:,:-1,:].data**2 + (x_wind[:-1,:,:-1,:].data*x_grad[:-1,:,:-1,:].data)**2 + (y_wind[:-1,:,:-1,:].data*y_grad[:-1,:,:,:].data)**2)
    # numerator = np.sqrt(x_diff.data**2 + y_diff.data**2)
    denominator = f[:-1,:,:-1,:]*mag[:-1,:,:-1,:]
    
    rol = np.log(numerator/denominator)
    rol = np.mean(rol, axis=0)
    
    plt.figure(figsize=(10,5))
    plt.contourf(np.arange(-longitudes/2,longitudes/2), np.array(h_label[:57]), np.roll(rol[:57,lat_slice,:], 72, axis=1), reds.N, cmap=reds)
    plt.title('Mean Lagrangian Rossby Number at Equator, t=%s to %s days' %(start/4, end/4))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.ylabel('Height [km]')
    plt.colorbar(pad=0.1)
    plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/laso/epsfigs/lrn_equator.eps', format='eps')  
    plt.show()
    
    plt.figure(figsize=(10,5))
    plt.contourf(np.arange(-longitudes/2,longitudes/2), np.arange(-latitudes/2, latitudes/2-1), np.roll(rol[level,:,:], 72, axis=1), reds.N, cmap=reds)
    plt.title('Mean Lagrangian Rossby Number at h=%s km, t=%s to %s days' %(h_label[level], start/4, end/4))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.ylabel('Latitude [degrees]')
    plt.yticks((-36,-18,0,18,36), ('90W', '45W', '0', '45E','90E'))
    plt.colorbar(pad=0.1)
    plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/laso/epsfigs/lrn_%s.eps' %(h_label[level]), format='eps')  
    plt.show()
    
    

