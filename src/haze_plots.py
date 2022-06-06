#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 18 16:21:44 2022

@author: Mo Cohen
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import scipy as sp


plasma = mpl_cm.get_cmap('plasma')



def plot_aero(data,time_slice=-1,level=5, scroll=False):
    
    """ This function plots the global haze distribution.
    
    Inputs: an .npz file list, time index, level index, and option to plot all levels
    Outputs: a contour fill plot of the haze distribution, centered at the 
    substellar point for tidally locked planets"""
    
    mmr = data['mmr'] # Extract mmr cube from file
    levs = data['lev'] # Extract levels in units of sigma (ratio of pressure to surface pressure)
    heights = np.round(levs*999.99,0)  # Convert to pressure in mbar

    fig,ax = plt.subplots(figsize = (10,5)) # Create figure and axi
    im = ax.contourf(mmr[time_slice,level,:,:], cmap=plasma,extend="min") # Use a contour fill plot with plasma color map
    plt.title('Haze mass mixing ratio at %s mbar' %heights[level]) # Title and automatically label with the height
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.xticks((0,8,16,24,32,40,48,56,64),('180W','135W','90W','45W','0','45E','90E','135E','180E'))
    plt.yticks((0,8,16,24,32),('90S','45S','0','45N','90N')) 
    # plt.xticks((-32,-24,-16,-8,0,8,16,24,32),('180W','135W','90W','45W','0','45E','90E','135E','180E'))
    # plt.yticks((-16,-8,0,8,16),('90S','45S','0','45N','90N')) 
    cbar = plt.colorbar(im)
    cbar.set_label('kg/kg')
    plt.show() 
    
    if scroll==True: # Scroll option will make plots for every atmospheric layer, from top to bottom
        
        for layer in range(0,10):
            
            fig,ax = plt.subplots(figsize = (10,5))
            im = ax.contourf(mmr[time_slice,layer,:,:], cmap=plasma,extend="min")
            plt.title('Haze mass mixing ratio at %s mbar' %heights[layer])
            plt.xlabel('Longitude')
            plt.ylabel('Latitude')
            plt.xticks((0,8,16,24,32,40,48,56,64),('180W','135W','90W','45W','0','45E','90E','135E','180E'))
            plt.yticks((0,8,16,24,32),('90S','45S','0','45N','90N')) 
            # plt.xticks((-32,-24,-16,-8,0,8,16,24,32),('180W','135W','90W','45W','0','45E','90E','135E','180E'))
            # plt.yticks((-16,-8,0,8,16),('90S','45S','0','45N','90N')) 
            cbar = plt.colorbar(im)
            cbar.set_label('kg/kg')
            plt.show() 
        