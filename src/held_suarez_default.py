#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 15:19:27 2021

@author: Mo Cohen
"""

import netCDF4 as nc
from netCDF4 import MFDataset
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
from matplotlib.colors import TwoSlopeNorm

reds = mpl_cm.get_cmap('Reds')
redblu = mpl_cm.get_cmap('RdBu')

path = ''
data = nc.Dataset(path)

u = data['ucomp'][:]
v = data['vcomp'][:]


def plot_streamlines(u, v, level=14, time_slice=-1):
    
   
    speed = np.sqrt(u**2 + v**2)

    X,Y = np.meshgrid(np.arange(-64,64), np.arange(-32,32))
    fig = plt.figure(figsize = (12, 7)) 
    strm = plt.streamplot(X, Y, u[time_slice,level,:,:], v[time_slice,level,:,:], density = 0.5, color=speed[time_slice,level,:,:], cmap=reds)
    fig.colorbar(strm.lines)
    plt.title('Wind speed and direction [m s-1], level %s' %(level))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.xticks((-64,-48,-32,-16,0,16,32,48,64),('180W','135W','90W','45W','0','45E','90E','135E','180E'))
    plt.yticks((-32,-24,-12,0,12,24,32),('90S','60S','30S','0','30N','60N','90N'))    
    plt.show() 
    
    
def plot_hovmoeller(u, bk, lat=32):

    x_axis = np.arange(0,u.shape[0])
    y_axis = np.arange(0,u.shape[1])
    
    zonal_mean = np.mean(u,axis=3)
    
    plt.contourf(x_axis, y_axis[:], zonal_mean[:,:,lat].T, redblu.N, cmap=redblu, norm=TwoSlopeNorm(0))
    plt.title('Zonal Mean Equatorial Wind [m s-1]')
    plt.xlabel('Time [months]')
    plt.ylabel('Level')
    plt.colorbar(pad=0.1)
    plt.show()
    