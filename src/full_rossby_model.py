#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 16:32:45 2023

@author: Mo Cohen
"""
import iris, windspharm
from windspharm.iris import VectorWind
import iris.plot as iplt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
from matplotlib.colors import TwoSlopeNorm
from iris.analysis import calculus
from numpy import unravel_index

period = (505, 525)

def extract_core(cubes, start=period[0], end=period[1], level=8):
        
    
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'y_wind':
            y_wind = cube[start:end,:,:,:].copy()
            
    heights = np.round(x_wind.coord('level_height').points*1e-03,2) # in km
    longitudes = x_wind.shape[2]
    latitudes = x_wind.shape[1]
    # Extract some data for labelling/creating plots

    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())    

    winds = VectorWind(x_wind, y_wind) # Windspharm VectorWind object   
    rel_vort = winds.vorticity()   

    rv_min = unravel_index(np.argmax(rel_vort[level,:20,0:72].data, axis=None), rel_vort[level,:20,0:72].shape)
    print(rv_min[0], rv_min[1])
    print(x_wind.coord('latitude').points[rv_min[0]], x_wind.coord('longitude').points[rv_min[1]])
    print(rel_vort[level,rv_min[0],rv_min[1]].data)      
            
    return rv_min


def model_rwave(cubes,startlon=30,start=218,end=250,nlat=90,nlon=144,level=8,omega=1.19e-05,g=9.12,radius=5797818,lat=80,save='no'):
    
    return