#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 21:13:33 2021

@author: Maureen Cohen
"""
import iris
import iris.plot as iplt
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
from matplotlib.colors import TwoSlopeNorm


brewer_redblu = mpl_cm.get_cmap('RdBu_r')

def gravity_wave_drag(cubes, time=-1, level=47):
    
    for cube in cubes:
        if cube.long_name == 'change_over_time_in_air_temperature_due_to_gravity_wave_drag':
            drag = cube.copy()            
    
    heights = np.round(drag.coord('level_height').points*1e-03,0)
    run_length = np.arange(0,drag.shape[0])
    # longitudes = drag.shape[3]
    # latitudes = drag.coord('latitude').points
    
    plt.figure(figsize=(10,5))
    iplt.contour(drag[time,level,:,:], brewer_redblu.N, cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Temp change due to gravity wave drag [K], month %s, h=%s km' %(run_length[time],heights[level]))
    plt.colorbar(pad=0.1)
    plt.show()
       
    # plt.figure(figsize=(10,5))
    # plt.contourf(np.arange(-longitudes/2, longitudes/2), latitudes, np.roll(drag[time,level,:,:].data, 72, axis=1), brewer_redblu.N, cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    # plt.title('Temp change due to gravity wave drag [K], month %s, h=%s km' %(run_length[time],heights[level]))
    # plt.xlabel('Longitude [degrees]')
    # plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    # plt.ylabel('Latitude [degrees]')
    # plt.colorbar(pad=0.1)
    # plt.show()