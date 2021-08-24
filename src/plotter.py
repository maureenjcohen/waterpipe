#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 25 13:14:40 2021

@author: Maureen Cohen
"""
import iris
import iris.plot as iplt
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
from matplotlib.colors import TwoSlopeNorm


reds = mpl_cm.get_cmap('brewer_Reds_09')

def plotter(cubes, time=-1, level=47):
    
    for cube in cubes:
        if cube.standard_name == 'specific_humidity':
            data = cube.copy()            
    
    run_length = np.arange(0,data.shape[0])
    # longitudes = drag.shape[3]
    # latitudes = drag.coord('latitude').points
    
    plt.figure(figsize=(12,6))
    iplt.contourf(data[time,level,:,:], reds.N, cmap=reds)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Specific humidity, month %s' %(run_length[time]+1))
    plt.colorbar(pad=0.1)
    plt.show()