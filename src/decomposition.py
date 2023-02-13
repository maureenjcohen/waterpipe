#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 13:36:06 2021

@author: Maureen Cohen
"""
import iris, windspharm
import numpy as np
import matplotlib.pyplot as plt



def decomposition(cubes, n=3, start=0, end=120, level=47):  
    
    """ Uses the windspharm package to perform a Helmholtz decomposition on an Iris cube
        Helmholtz composition splits the vector field into its divergent and rotational components
        Also plots the eddy rotational component (deviation from the zonal mean of the rotational component)
        
        Arguments: Iris CubeList, n relates to the density of the arrows (3 means plots every 3rd arrow), 
        time_slice is the time (defaults to last time output in the list), 
        level is the model level number (17 is around 0.4 bar)
        
        Outputs: 4 quiverplots of the wind vectors, divergent component, rotational, and eddy rotational components"""    
            
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'y_wind':
            y_wind = cube[start:end,:,:,:].copy()
        if cube.standard_name =='air_pressure':
            pressure = cube[start:end,:,:,:].copy()
            
    # Select the cubes we want from the cube list
 
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
    
    x_wind = x_wind.collapsed('time', iris.analysis.MEAN)
    y_wind = y_wind.collapsed('time', iris.analysis.MEAN)
    pressure = pressure.collapsed('time', iris.analysis.MEAN)
    
    height = [('level_height', x_wind.coord('level_height').points)]
    pressure = pressure.interpolate(height, iris.analysis.Linear())
    # Regrid so that all three cubes are on the same x, y, z grid. Uses the x_wind as reference for the others
    
    p_heights = np.round(pressure.data*1e-05,2)
    km_heights = np.round(pressure.coord('level_height').points*1e-03,2)
    # Extract pressure and kilometer values for labeling the plots
    
    winds = windspharm.iris.VectorWind(x_wind,y_wind)
    # Create a VectorWind data object from the x and y wind cubes
    uchi,vchi,upsi,vpsi = winds.helmholtz(truncation=21)
    # Calculate the Helmholtz decomposition. Truncation is set to 21 because 
    # this is what Hammond and Lewis 2021 used.
    
    zonal_upsi = upsi.collapsed('longitude', iris.analysis.MEAN)
    zonal_vpsi = vpsi.collapsed('longitude', iris.analysis.MEAN)
    # Calculate zonal means of the x and y components of the rotational component
    eddy_upsi = upsi - zonal_upsi
    eddy_vpsi = vpsi - zonal_vpsi 
    # Subtract zonal means from the original cubes to get the x and y eddy rotational components

    X,Y = np.meshgrid(np.arange(-72,72), np.arange(-45,45))
    # Create meshgrid with the spatial dimensions of the Iris cube
    
    fig1, ax1 = plt.subplots(figsize = (10,5)) 
    q1 = ax1.quiver(X[::n,::n], Y[::n,::n], np.roll(x_wind[level,::n,::n].data, 72, axis=1), np.roll(y_wind[level,::n,::n].data, 72, axis=1), angles='xy', scale_units='xy', scale=5)
    # Create a quiver plot. The np.roll function moves the cube data so the substellar point is centred.
    ax1.quiverkey(q1, X=0.9, Y=1.05, U=5, label='5 m/s', labelpos='E', coordinates='axes')
    # This creates the key with the arrow size. The value of U, the label string, and the value of scale in the previous line should all match.
    plt.title('Wind Vectors [m/s], h = %s km' %(km_heights[level]))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.yticks((-45,-30,-15,0,15,30,45),('90S','60S','30S','0','30N','60N','90N')) 
    # Labelled the longs and lats manually 
    # plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/laso/epsfigs/circulation.eps', format='eps')   
    plt.show()
    
    fig2, ax2 = plt.subplots(figsize = (10,5)) 
    q2 = ax2.quiver(X[::n,::n], Y[::n,::n], np.roll(uchi[level,::n,::n].data, 72, axis=1), np.roll(-vchi[level,::n,::n].data, 72, axis=1), angles='xy', scale_units='xy', scale=3)
    # Note -vchi is plotted. Same with the other y components below.
    ax2.quiverkey(q2, X=0.9, Y=1.05, U=3, label='3 m/s', labelpos='E', coordinates='axes')
    # Made the arrow size smaller since the divergent winds are weaker
    plt.title('Divergent component of wind [m s-1], h=%s bar, %s km' %(p_heights[level,0,0], km_heights[level]))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.yticks((-45,-30,-15,0,15,30,45),('90S','60S','30S','0','30N','60N','90N'))    
    plt.show() 
                
    fig3, ax3 = plt.subplots(figsize = (10,5)) 
    q3 = ax3.quiver(X[::n,::n], Y[::n,::n], np.roll(upsi[level,::n,::n].data, 72, axis=1), np.roll(-vpsi[level,::n,::n].data, 72, axis=1), angles='xy', scale_units='xy', scale=5)
    ax3.quiverkey(q3, X=0.9, Y=1.05, U=5, label='5 m/s', labelpos='E', coordinates='axes')
    plt.title('Rotational component of wind [m s-1], h=%s bar, %s km' %(p_heights[level,0,0], km_heights[level]))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.yticks((-45,-30,-15,0,15,30,45),('90S','60S','30S','0','30N','60N','90N'))    
    plt.show()
    
    fig4, ax4 = plt.subplots(figsize = (10,5)) 
    q4 = ax4.quiver(X[::n,::n], Y[::n,::n], np.roll(eddy_upsi[level,::n,::n].data, 72, axis=1), np.roll(-eddy_vpsi[level,::n,::n].data, 72, axis=1), angles='xy', scale_units='xy', scale=4)
    ax4.quiverkey(q4, X=0.9, Y=1.05, U=4, label='4 m/s', labelpos='E', coordinates='axes')
    plt.title('Eddy Rotational Component of Wind [m/s], h = %s km' %(km_heights[level]))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.yticks((-45,-30,-15,0,15,30,45),('90S','60S','30S','0','30N','60N','90N')) 
    # plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/laso/epsfigs/eddy_rot.eps', format='eps')   
    plt.show()
    
    