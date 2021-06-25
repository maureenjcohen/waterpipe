#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 11:43:11 2021

@author: Maureen Cohen
"""
import iris
import iris.plot as iplt
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import gif
import numpy as np
from matplotlib.colors import TwoSlopeNorm


brewer_redblu = mpl_cm.get_cmap('RdBu_r')

x_wind = control[69].copy()
data = x_wind.data
heights = x_wind.coord('level_height').points
longs = x_wind.coord('longitude').points
run_length = x_wind.shape[0]


@gif.frame
def plot_1d(cube, time):
    plt.plot(data[time,44:,45,0], heights[44:])
    plt.title('Month %s' %time)
    plt.xlabel('Wind velocity [m s-1]')
    plt.xlim(-140,140)
    plt.ylabel('Height [m]')
    
frames = []
for i in range(0,run_length):
    frame = plot_1d(data,i)
    frames.append(frame)

gif.save(frames, str('/exports/csce/datastore/geos/users/s1144983/um_data/control/gifs') + '/um_qbo_1dstrat.gif', duration = 30, unit = 's', between='startend')


""" 2D version """ 

@gif.frame
def plot_2d(cube, time):
    plt.contourf(longs, heights[44:], data[time,44:,45,:], np.linspace(-140,140,120), cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    plt.title('Month %s' %time)
    plt.xlabel('Longitude [degrees]')
    plt.ylabel('Height [m]')
    plt.colorbar(pad=0.1)
    
frames = []
for i in range(0,run_length):
    frame = plot_2d(data,i)
    frames.append(frame)
    
gif.save(frames, str('/exports/csce/datastore/geos/users/s1144983/um_data/control/gifs') + '/um_qbo_2d_strat.gif', duration = 30, unit = 's', between='startend')
