#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 12 15:22:45 2022

@author: Mo Cohen
"""
import iris, windspharm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import scipy as sp
from iris.coord_systems import GeogCS


plasma = mpl_cm.get_cmap('plasma')
redblu = mpl_cm.get_cmap('RdBu')


def wave_spectrum(cubes,start=500,end=520,nlat=90,nlon=144,level=8,omega=1.19e-05,g=9.12,radius=5797818,lat=73,n=3):  
    
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
        if cube.standard_name =='air_potential_temperature':
            theta = cube[start:end,:,:,:].copy()
            
    # Select the cubes we want from the cube list
 
    x_wind.coord('latitude').coord_system = GeogCS(radius)
    x_wind.coord('longitude').coord_system = GeogCS(radius)
    y_wind.coord('latitude').coord_system = GeogCS(radius)
    y_wind.coord('longitude').coord_system = GeogCS(radius)

    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
    
    height = [('level_height', x_wind.coord('level_height').points)]
    pressure = pressure.interpolate(height, iris.analysis.Linear())
    # Regrid so that all three cubes are on the same x, y, z grid. Uses the x_wind as reference for the others
    
    p_heights = np.round(pressure.data*1e-05,2)
    km_heights = np.round(pressure.coord('level_height').points*1e-03,2)
    # Extract pressure and kilometer values for labeling the plots
    
    winds = windspharm.iris.VectorWind(x_wind[:,level,:,:],y_wind[:,level,:,:])
    uchi,vchi,upsi,vpsi = winds.helmholtz(truncation=21)
   
    zonal_upsi = upsi.collapsed('longitude', iris.analysis.MEAN)
    zonal_vpsi = vpsi.collapsed('longitude', iris.analysis.MEAN)
    eddy_upsi = upsi - zonal_upsi
    eddy_vpsi = vpsi - zonal_vpsi
    magnitude = np.sqrt(eddy_upsi.data**2 + eddy_vpsi.data**2)
#    magnitude = np.sqrt(x_wind[:,level,:,:].data**2 + y_wind[:,level,:,:].data**2)
    eddies = windspharm.iris.VectorWind(eddy_upsi, eddy_vpsi)
    rel_vort = eddies.vorticity()
    
    fft2 = sp.fft.fftshift(sp.fftpack.fft2(sp.fft.ifftshift(rel_vort.data)))
    yfreqs = sp.fft.fftshift(sp.fft.fftfreq(fft2.shape[1],d=1./nlat))
    xfreqs = sp.fft.fftshift(sp.fft.fftfreq(fft2.shape[2],d=1./nlon))
    psd = np.abs(fft2)**2
    
    quadrant = psd[:,46:51,73:78] # Select positive zonal and meridional wavenumbers <= 5 
    xfreqs = xfreqs[73:78] # Select corresponding x-frequencies
    yfreqs = yfreqs[46:51] # Select corresponding y-frequencies
    print(yfreqs)
    
    for time in range(0,quadrant.shape[0]):
        
        X,Y = np.meshgrid(np.arange(0,nlon), np.arange(0,nlat)) 
        fig, ax = plt.subplots(figsize=(8.5,5))
        q1 = ax.quiver(X[::n,::n],Y[::n,::n], np.roll(x_wind[time,level,::n,::n].data,int(nlon/(2*n)),axis=1), 
                   np.roll(y_wind[time,level,::n,::n].data,int(nlon/(2*n)),axis=1),scale_units='xy',scale=5)
        ax.quiverkey(q1, X=0.9, Y=1.05, U=20, label='20 m/s', labelpos='E', coordinates='axes')
        plt.title('General circulation, h=%s km, day=%s' %(km_heights[level], start+time))
        plt.xticks((0,12,24,36,48,60,72,84,96,108,120,132,144),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
        plt.yticks((0,15,30,45,60,75,90),('90S','60S','30S','0','30N','60N','90N'))    
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        plt.show()
    
        X1,Y1 = np.meshgrid(np.arange(-72,72), np.arange(-45,45))
        # Create meshgrid with the spatial dimensions of the Iris cube
        
        fig1, ax1 = plt.subplots(figsize = (8.5,5)) 
        q = ax1.quiver(X1[::n,::n], Y1[::n,::n], np.roll(eddy_upsi[time,::n,::n].data, 72, axis=1), np.roll(-eddy_vpsi[time,::n,::n].data, 72, axis=1), angles='xy', scale_units='xy', scale=4)
        ax1.quiverkey(q, X=0.9, Y=1.05, U=4, label='4 m/s', labelpos='E', coordinates='axes')
        plt.title('Eddy Rotational Component of Wind [m/s], h=%s km, day=%s' %(km_heights[level], start+time))
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
        plt.yticks((-45,-30,-15,0,15,30,45),('90S','60S','30S','0','30N','60N','90N')) 
        # plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/laso/epsfigs/eddy_rot.eps', format='eps')   
        plt.show()
        
        # fig3,ax3 = plt.subplots(figsize = (8.5,5))
        # im2 = ax3.contourf(rel_vort[time,:,:].data, cmap=plasma)
        # # im = ax2.contourf(xfreqs,yfreqs,psd, cmap=plasma)
        # plt.title('Eddy rotational wind PV, h=%s km, day=%s' %(km_heights[level], start+time))
        # plt.xlabel('Longitude')
        # plt.ylabel('Latitude')
        # fig3.colorbar(im2)
        # plt.show()  
        
        # fig2,ax2 = plt.subplots(figsize = (8.5,5))
        # im = ax2.contourf(xfreqs,yfreqs,quadrant[time,:,:], cmap=plasma)
        # # im = ax2.contourf(xfreqs,yfreqs,psd, cmap=plasma)
        # plt.title('Power spectrum of eddy rotational wind PV, h=%s km, day=%s' %(km_heights[level], start+time))
        # plt.xlabel('Zonal wavenumber')
        # plt.ylabel('Meridional wavenumber')
        # fig2.colorbar(im)
        # plt.show()
        
    