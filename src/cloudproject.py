#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 09:42:00 2022

@author: Mo Cohen
"""

import os, iris
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import iris.plot as iplt
import iris.coords
from iris.coords import DimCoord
from matplotlib.colors import TwoSlopeNorm
import numpy as np
from numpy import unravel_index
import scipy as sp
#import windspharm
from iris.analysis import calculus
import pandas as pd

redblu = mpl_cm.get_cmap('coolwarm')
plasma = mpl_cm.get_cmap('plasma')
hot = mpl_cm.get_cmap('hot')
heat = mpl_cm.get_cmap('gist_heat')
brewer_bg = mpl_cm.get_cmap('brewer_PuBu_09')

def climgrid(datalist,nlat=90,nlon=144,nlev=38,level=15,start=0,end=300,ndata=5,n=4,save='no'):

    """ Create a full-page grid of plots displaying the vertical temperature profile, zonal
    mean zonal wind, vertical water vapour profile, and surface temp with wind flow vectors for
    each dataset in the input datalist.
    
    Inputs: Datalist - a list of Iris data lists or data sets
    Output: ndata x 4 plots in a grid layout in a single image """
    
    fig, ax = plt.subplots(figsize=(16,22), nrows=5, ncols=4)
    names = ['Control ProxB','Warm ProxB','Control TRAP1-e','Warm TRAP1-e','Dry TRAP1-e']
    
    for i in range(ndata):
        data = datalist[i]
        for cube in data:
            if cube.standard_name == 'air_potential_temperature':
                potential_temp = cube[start:end,:,:,:].copy()
            if cube.standard_name == 'air_pressure':
                air_pressure = cube[start:end,:,:,:].copy()
            if cube.standard_name == 'x_wind':
                x_wind = cube[start:end,:,:,:].copy()
            if cube.standard_name == 'y_wind':
                y_wind = cube[start:end,:,:,:].copy()
            if cube.standard_name == 'surface_temperature':
                surface_temp = cube[start:end,:,:].copy()
            if cube.standard_name == 'specific_humidity' and i < ndata-1:
                spec_humidity = cube[start:end,:,:,:].copy()
            elif cube.standard_name == 'specific_humidity' and i == ndata-1: 
                spec_humidity = cube.copy()*0.0
           

        y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
        
        vertical = [('level_height', x_wind.coord('level_height').points)]
        potential_temp = potential_temp.regrid(x_wind, iris.analysis.Linear())
        potential_temp = potential_temp.interpolate(vertical, iris.analysis.Linear())
        air_pressure = air_pressure.regrid(x_wind, iris.analysis.Linear())
        air_pressure = air_pressure.interpolate(vertical, iris.analysis.Linear())

        heights = np.round(x_wind.coord('level_height').points*1e-03,0)
        lats = np.round(x_wind.coord('latitude').points,0)
        lons = np.round(x_wind.coord('longitude').points,0)
        
        p0 = iris.coords.AuxCoord(100000.0, long_name='reference_pressure', units='Pa')
        p0.convert_units(air_pressure.units)
        absolute_temp = potential_temp*((air_pressure/p0)**(287.05/1005))
        absolute_temp = np.mean(absolute_temp.data,axis=0)
        
        zmzw = x_wind.collapsed('longitude',iris.analysis.MEAN)
        zmzw = np.mean(zmzw.data,axis=0)
        
        spec_humidity = spec_humidity.regrid(x_wind, iris.analysis.Linear())
        spec_humidity = spec_humidity.interpolate(vertical, iris.analysis.Linear())
        spec_humidity = np.mean(spec_humidity.data,axis=0)
        
        surface_temp = np.mean(surface_temp.data,axis=0)
        
        meaned_x = np.mean(x_wind.data,axis=0)
        meaned_y = np.mean(y_wind.data,axis=0)
        
#        X,Y = np.meshgrid(np.arange(-x_wind.shape[3]/2,x_wind.shape[3]/2), np.arange(-x_wind.shape[2]/2,x_wind.shape[2]/2))

        ax[i,0].plot(absolute_temp[:,45,0], heights, color='r', label='Substellar')
        ax[i,0].plot(absolute_temp[:,45,72], heights, color='b', label='Antistellar')
        ax[i,0].legend()
        ax[0,0].set_title('Temperature profile [K]',fontsize=12)
        ax[ndata-1,0].set_xlabel('Temperature [K]')
        ax[i,0].set_ylabel('%s \n Height [km]' %names[i],fontsize=12)
        
        cont = ax[i,1].contourf(lats,heights,zmzw,levels=np.arange(-60,140,20),cmap=redblu,norm=TwoSlopeNorm(0))
        ax[0,1].set_title('Zonal mean zonal wind [m/s]',fontsize=12)
        ax[ndata-1,1].set_xlabel('Latitude [degrees]')
        fig.colorbar(cont,ax=ax[i,1],orientation='vertical')
        
        
        ax[i,2].plot(spec_humidity[:,45,0], heights, color='r', label='Substellar')
        ax[i,2].plot(spec_humidity[:,45,72], heights, color='b', label='Antistellar')
        ax[i,2].legend()
        ax[0,2].set_title('Vapor profile [kg/kg]',fontsize=12)
        ax[ndata-1,2].set_xlabel('Specific humidity [kg/kg]')
        
        surf = ax[i,3].contourf(np.arange(-72,72)*2.5,lats,np.roll(surface_temp,72,axis=1),levels=np.arange(100,400,20),cmap=hot)
#        ax[i,3].streamplot(X,Y,np.roll(meaned_x[level,:,:],int(nlon/2),axis=1), 
#                       np.roll(meaned_y[level,:,:],int(nlon/(2)),axis=1),density=0.5,color='k')
        ax[0,3].set_title('Surface temperature [K]',fontsize=12)
        ax[ndata-1,3].set_xlabel('Longitude [degrees]')
        fig.colorbar(surf,ax=ax[i,3],orientation='vertical')



    if save == 'yes':
        plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/cloudproject/epsfigs/climgrid_tight.eps', format='eps',bbox_inches='tight')
    else:
        pass
    
    plt.show()


def composite(cubes, time_slice=500, nlat=90, nlon=144, nlev=38, level=8, n=4, cloudtype='both', fractype='mass', save='no'):

    """ Plot composites of the cloud cover and horizontal wind vectors
    
    Inputs: Iris CubeList, start and end times of data, number of latitudes, number of longitudes,
    number of levels, level to be plotted, meaning period (1 = no meaning), sparsity of quiver arrows,
    cloud type (options: ice, liq, both, none), fraction type (options: mass or volume), whether to save

    Outputs: Plot of horizontal wind vectors, with or without the cloud cover superimposed"""

    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube[time_slice,:,:,:].copy()
        if cube.standard_name == 'y_wind':
            y_wind = cube[time_slice,:,:,:].copy() 
        if cube.long_name == 'ice_cloud_volume_fraction_in_atmosphere_layer' and fractype=='volume' and cloudtype != 'none':
            ice = cube[time_slice,:,:,:].copy()
        if cube.long_name == 'liquid_cloud_volume_fraction_in_atmosphere_layer' and fractype=='volume' and cloudtype != 'none':
            liq = cube[time_slice,:,:,:].copy()
        if cube.standard_name == 'mass_fraction_of_cloud_ice_in_air' and fractype=='mass' and cloudtype != 'none':
            ice = cube[time_slice,:,:,:].copy()
        if cube.standard_name == 'mass_fraction_of_cloud_liquid_water_in_air' and fractype=='mass' and cloudtype != 'none':
            liq = cube[time_slice,:,:,:].copy()
            
        
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
    heights = np.round(x_wind.coord('level_height').points*1e-03,2)
            
    
    if cloudtype=='ice':
        cloud = ice
        titleterm = 'Ice cloud'
    elif cloudtype=='liq':
        cloud = liq
        titleterm = 'Liquid cloud'
    elif cloudtype=='both':
        cloud = ice + liq
        titleterm = 'Total cloud'
    elif cloudtype=='none':
        titleterm = 'Horizontal wind'
    else:
        print('Argument cloudtype must be ice, liq, both, or none. Default is both.')

    X,Y = np.meshgrid(np.arange(0,nlon), np.arange(0,nlat))   

    if cloudtype=='none':
        fig, ax = plt.subplots(figsize=(8.5,5))
        q1 = ax.quiver(X[::n,::n],Y[::n,::n], np.roll(x_wind[level,::n,::n].data,int(nlon/(2*n)),axis=1), 
                   np.roll(y_wind[level,::n,::n].data,int(nlon/(2*n)),axis=1),scale_units='xy',scale=5)
        ax.quiverkey(q1, X=0.9, Y=1.05, U=20, label='20 m/s', labelpos='E', coordinates='axes')
        plt.title('%s, day %s, h=%s km' %(titleterm, time_slice, heights[level]))
        plt.xticks((0,12,24,36,48,60,72,84,96,108,120,132,144),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
        plt.yticks((0,15,30,45,60,75,90),('90S','60S','30S','0','30N','60N','90N'))    
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
  
        if save == 'yes':
            plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/cloudproject/epsfigs/quiver_nocloud_%s_trap.eps' %(time_slice), format='eps')
        else:
            pass
        plt.show()

    else:
        fig, ax = plt.subplots(figsize=(10,5))
        plt.imshow(np.roll(cloud[level,:,:].data*1e4,int(nlon/2),axis=1), vmin=0, vmax=2,cmap=brewer_bg)
        cbar = plt.colorbar()

        if fractype=='mass':
            cbar.set_label('$10^{-4}$ kg/kg', loc='center')    

        q1 = ax.quiver(X[::n,::n],Y[::n,::n], np.roll(x_wind[level,::n,::n].data,int(nlon/(2*n)),axis=1), 
                np.roll(-y_wind[level,::n,::n].data,int(nlon/(2*n)),axis=1),scale_units='xy',scale=5)
        ax.quiverkey(q1, X=0.95, Y=1.05, U=20, label='20 m/s', labelpos='E', coordinates='axes')

        plt.title('%s and horizontal wind, day %s, h=%s km' %(titleterm, time_slice, heights[level]))
        plt.xticks((0,12,24,36,48,60,72,84,96,108,120,132,144),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
        plt.yticks((90,75,60,45,30,15,0),('90S','60S','30S','0','30N','60N','90N'))    
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        if save == 'yes':
            plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/cloudproject/epsfigs/quiver_withcloud_%s_prox.eps' %(time_slice), format='eps')
        else:
            pass            
        plt.show()

def rwave_velocity(datalist,start=500,end=600,nlat=90,nlon=144,level=8,omega=1.19e-05,g=9.12,radius=5797818,lat=73,save='no'):

    """ This function calculates the Rossby wave phase speed (including zonal wind) over time.

    The code performs a Helmholtz decomposition of the horizontal wind field at the input level.
    The wind vectors are converted to a wind speed (magnitude) and this 2-D data is input into a 
    2-D Fourier transform. The Fourier transform finds the spatial frequencies in the x- and y-directions.
    The function finds the maximum intensity wave and extracts the zonal and meridional wavenumbers.
    It then constructs a time series of the Rossby wave phase speed using the varying wavenumbers and 
    varying background zonal wind over time as inputs into the phase speed formula.

    At certain latitudes, the phase velocity alternates between positive and negative. This is the region
    where the cyclonic structure appears to travel back and forth periodically. """
    
    cphase_list = []
    cphase_hlist = []
    cgroup_list = []
    names = ['Control TRAP-1e', 'Dry TRAP-1e']
    
    for cubes in datalist:
        for cube in cubes:
            if cube.standard_name == 'x_wind':
                x_wind = cube[start:end,:,:,:].copy()
            if cube.standard_name == 'y_wind':
                y_wind = cube[start:end,:,:,:].copy()
            if cube.standard_name =='air_potential_temperature':
                theta = cube[start:end,:,:,:].copy()
            if cube.standard_name =='air_pressure':
                pressure = cube[start:end,:,:,:].copy()
            
        y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
        km_heights = np.round(x_wind.coord('level_height').points*1e-03,2)
        latitudes = x_wind.coord('latitude').points
        
        winds = windspharm.iris.VectorWind(x_wind[:,level,:,:],y_wind[:,level,:,:])
        uchi,vchi,upsi,vpsi = winds.helmholtz(truncation=21)
    
        
        zonal_upsi = upsi.collapsed('longitude', iris.analysis.MEAN)
        zonal_vpsi = vpsi.collapsed('longitude', iris.analysis.MEAN)
        eddy_upsi = upsi - zonal_upsi
        eddy_vpsi = vpsi - zonal_vpsi
        eddies = windspharm.iris.VectorWind(eddy_upsi, eddy_vpsi)
        rel_vort = eddies.vorticity()
#        magnitude = np.sqrt(eddy_upsi.data**2 + eddy_vpsi.data**2)
#        rel_vort = winds.vorticity()
        
        fft2 = sp.fft.fftshift(sp.fftpack.fft2(sp.fft.ifftshift(rel_vort.data)))
        yfreqs = sp.fft.fftshift(sp.fft.fftfreq(fft2.shape[1],d=1./nlat))
        xfreqs = sp.fft.fftshift(sp.fft.fftfreq(fft2.shape[2],d=1./nlon))
        psd = np.abs(fft2)**2
        x_zero = int(nlon/2) # Identify the midpoint of the x-axis
        y_zero = int(nlat/2) # Identify the midpoint of the y-axis
        
        quadrant = psd[:,y_zero:y_zero+5,x_zero:x_zero+5] # Select positive zonal and meridional wavenumbers <= 5 
        xfreqs = xfreqs[x_zero:x_zero+5] # Select corresponding x-frequencies
        yfreqs = yfreqs[y_zero:y_zero+5] # Select corresponding y-frequencies
        
        x_wavenumbers = []
        y_wavenumbers = []
        for time in range(0,psd.shape[0]):
            peak = unravel_index(np.argmax(quadrant[time,:,:], axis=None), quadrant[time,:,:].shape)
            print(peak)
            x_wavenumbers.append(1) # We add one here because the array index starts at 0, but index 0 is wavenumber 1
            y_wavenumbers.append(0)
            
        x_wavenumbers, y_wavenumbers = np.array(x_wavenumbers), np.array(y_wavenumbers)
        print(x_wavenumbers, y_wavenumbers)
        
        lat_deg = int(latitudes[lat]) # Convert input row number to latitude in degrees north
        print(lat_deg)
        zmzw = x_wind[:,level,lat].collapsed('longitude',iris.analysis.MEAN)
        zmzw = zmzw.data
        
        lat_rad = lat_deg*(np.pi/180) # Convert input latitude to radians
        beta = 2*omega*np.cos(lat_rad)/radius # Beta factor    
        circum = 2*np.pi*radius*np.cos(lat_rad) # Circumference in meters at input latitude    
        lambda_x = circum/x_wavenumbers # Wavelength in x-direction at input latitude
        x_num = 2*np.pi/lambda_x
        lambda_y = 2*np.pi*radius/y_wavenumbers # Wavelength in y-direction at input latitude
        y_num = 2*np.pi/lambda_y
        
        d_theta = iris.analysis.calculus.differentiate(theta, 'level_height')
        bv_freq = np.mean(np.sqrt(np.abs((g/theta[:,:-1,:,:].data)*d_theta.data)), axis=-1)
        Ld = (bv_freq[:,level,lat]*6800)/(2*omega*np.sin(lat_rad))
        
        c_phase = (zmzw - ((beta + (zmzw/Ld**2))/(x_num**2+y_num**2 + (1/Ld)**2)))                               
        cphase_list.append(c_phase)
        c_phase_h = (zmzw - (beta/(x_num**2+y_num**2))) 
        cphase_hlist.append(c_phase_h)
        
    plt.plot(cphase_list[0],color='r', label=names[0] + ', with kd')
    plt.plot(cphase_list[1], color='b', label=names[1] + ', with kd')
    plt.plot(cphase_hlist[0],color='r', linestyle='dashed', label=names[0] + ', no kd')
    plt.plot(cphase_hlist[1], color='b', linestyle='dashed', label=names[1] + ', no kd')

    plt.title('Rossby wave phase velocity at %s N' %(lat_deg))
    plt.xlabel('Time [days]')
    plt.ylabel('Velocity [m/s]')
#    plt.ylim(-35,25)
    plt.legend(fontsize='x-small')
    if save == 'yes':
        plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/cloudproject/epsfigs/rwave_vel_fixedaxis.eps', format='eps', bbox_inches='tight')
    else:
        pass  
    plt.show()

        

def cphase_vals(datalist, start=0, end=300, ndata=5, nlat=90, nlon=144, level=8, lat=73,save='no'):
    
    names = ['Control ProxB',
                'Warm ProxB','Control TRAP1-e','Warm TRAP1-e','Dry TRAP1-e']
    
    omegas = [6.46e-06, 7.93e-06,1.19e-05, 1.75e-05, 1.19e-05]
    radii = [7160000, 7160000, 5797818, 5797818, 5797818]
    
    output_list = []
    
    for i in range(ndata):
        
        data = datalist[i]
        name = names[i]
        omega = omegas[i]
        radius = radii[i]
        
        for cube in data:
            if cube.standard_name == 'x_wind':
                x_wind = cube[start:end,:,:,:].copy()
            if cube.standard_name == 'y_wind':
                y_wind = cube[start:end,:,:,:].copy()
            
        y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
        km_heights = np.round(x_wind.coord('level_height').points*1e-03,2)
        
        winds = windspharm.iris.VectorWind(x_wind[:,level,:,:],y_wind[:,level,:,:])
        uchi,vchi,upsi,vpsi = winds.helmholtz(truncation=21)
    
        
        zonal_upsi = upsi.collapsed('longitude', iris.analysis.MEAN)
        zonal_vpsi = vpsi.collapsed('longitude', iris.analysis.MEAN)
        eddy_upsi = upsi - zonal_upsi
        eddy_vpsi = vpsi - zonal_vpsi
        magnitude = np.sqrt(eddy_upsi.data**2 + eddy_vpsi.data**2)
        
        fft2 = sp.fft.fftshift(sp.fftpack.fft2(sp.fft.ifftshift(magnitude)))
        yfreqs = sp.fft.fftshift(sp.fft.fftfreq(fft2.shape[1],d=1./nlat))
        xfreqs = sp.fft.fftshift(sp.fft.fftfreq(fft2.shape[2],d=1./nlon))
        psd = np.abs(fft2)**2
        
        quadrant = psd[:,45:51,72:78] # Select positive zonal and meridional wavenumbers <= 5 
        xfreqs = xfreqs[72:78] # Select corresponding x- and y-frequencies
        yfreqs = yfreqs[45:51]
    
        x_wavenumbers = []
        y_wavenumbers = []
        for time in range(0,psd.shape[0]):
            peak = unravel_index(quadrant[time,:,:].argmax(), quadrant[time,:,:].shape)
            x_wavenumbers.append(peak[1])
            y_wavenumbers.append(peak[0])
            
        x_wavenumbers, y_wavenumbers = np.array(x_wavenumbers)+1, np.array(y_wavenumbers)+1
        
        lat_deg = np.abs(lat-45)*2 # Convert input row number to latitude in degrees north
        zmzw = x_wind[:,level,lat].collapsed('longitude',iris.analysis.MEAN)
        zmzw = zmzw.data
        mean_zmzw = np.mean(zmzw,axis=0)
        
        lat_rad = lat_deg*(np.pi/180) # Convert input latitude to radians
        beta = 2*omega*np.cos(lat_rad)/radius # Beta factor    
        circum = 2*np.pi*radius*np.cos(lat_rad) # Circumference in meters at input latitude    
        lambda_x = circum/x_wavenumbers # Wavelength in x-direction at input latitude
        x_num = 2*np.pi/lambda_x
        lambda_y = 2*np.pi*radius/y_wavenumbers # Wavelength in y-direction at input latitude
        y_num = 2*np.pi/lambda_y
        
        c_phase = (zmzw - (beta/(x_num**2+y_num**2)))
        mean_cphase = np.mean(c_phase, axis=0)
        print(mean_cphase)

        
        c_phase_df = pd.DataFrame(index=range(start,end), columns=['Simulation','Beta', 'X', 'Y', 'Factor','Mean cphase','ZMZW','Mean ZMZW','Lat'])
        c_phase_df['Simulation'] = name
        c_phase_df['Beta'] = (beta*10**12)*np.ones_like(x_num)
        c_phase_df['X'] = x_num
        c_phase_df['Y'] = y_num
        c_phase_df['Factor'] = (beta/(x_num**2+y_num**2))
        c_phase_df['Mean cphase'] = mean_cphase*np.ones_like(x_num)
        c_phase_df['ZMZW'] = zmzw
        c_phase_df['Mean ZMZW'] = mean_zmzw*np.ones_like(x_num)
        c_phase_df['Lat'] = lat_deg*np.ones_like(x_num)
        
        output_list.append(c_phase_df)
        
    output = pd.concat(output_list)
    print(output)
    if save=='yes':
        output.to_csv('/exports/csce/datastore/geos/users/s1144983/papers/cloudproject/c_phase_vals.csv')
        

def find_zero(datalist, ndata=5, start=0, end=300, nlat=90, nlon=144, level=8, save='no'):
    
    names = ['Control ProxB',
                'Warm ProxB','Control TRAP1-e','Warm TRAP1-e','Dry TRAP1-e']
    
    omegas = [6.46e-06, 7.93e-06,1.19e-05, 1.75e-05, 1.19e-05]
    radii = [7160000, 7160000, 5797818, 5797818, 5797818]
    
    output_list = []
    
    for i in range(ndata):
        
        data = datalist[i]
        name = names[i]
        omega = omegas[i]
        radius = radii[i]
        
        for cube in data:
            if cube.standard_name == 'x_wind':
                x_wind = cube[start:end,:,:,:].copy()
            if cube.standard_name == 'y_wind':
                y_wind = cube[start:end,:,:,:].copy()
            
        y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
        km_heights = np.round(x_wind.coord('level_height').points*1e-03,2)
        
        winds = windspharm.iris.VectorWind(x_wind[:,level,:,:],y_wind[:,level,:,:])
        uchi,vchi,upsi,vpsi = winds.helmholtz(truncation=21)
          
        zonal_upsi = upsi.collapsed('longitude', iris.analysis.MEAN)
        zonal_vpsi = vpsi.collapsed('longitude', iris.analysis.MEAN)
        eddy_upsi = upsi - zonal_upsi
        eddy_vpsi = vpsi - zonal_vpsi
        magnitude = np.sqrt(eddy_upsi.data**2 + eddy_vpsi.data**2)
        
        fft2 = sp.fft.fftshift(sp.fftpack.fft2(sp.fft.ifftshift(magnitude)))
        yfreqs = sp.fft.fftshift(sp.fft.fftfreq(fft2.shape[1],d=1./nlat))
        xfreqs = sp.fft.fftshift(sp.fft.fftfreq(fft2.shape[2],d=1./nlon))
        psd = np.abs(fft2)**2
        
        quadrant = psd[:,45:51,72:78] # Select positive zonal and meridional wavenumbers <= 5 
        xfreqs = xfreqs[72:78] # Select corresponding x- and y-frequencies
        yfreqs = yfreqs[45:51]

        x_wavenumbers = []
        y_wavenumbers = []
        for time in range(0,psd.shape[0]):
            peak = unravel_index(quadrant[time,:,:].argmax(), quadrant[time,:,:].shape)
            x_wavenumbers.append(peak[1])
            y_wavenumbers.append(peak[0])
            
        x_wavenumbers, y_wavenumbers = np.array(x_wavenumbers)+1, np.array(y_wavenumbers)+1
        
        lat_deg = x_wind.coord('latitude').points
        zmzw = x_wind[:,level,:,:].collapsed('longitude',iris.analysis.MEAN)
        zmzw = zmzw.data
        
        lat_rad = lat_deg*(np.pi/180) # Convert input latitude to radians
        beta = 2*omega*np.cos(lat_rad)/radius # Beta factor    
        circum = 2*np.pi*radius*np.cos(lat_rad) # Circumference in meters at input latitude

        lambda_x_list = []
        for i in range(len(x_wavenumbers)):          
            lambda_x_item = circum/x_wavenumbers[i] # Wavelength in x-direction at input latitude
            lambda_x_list.append(lambda_x_item)
        lambda_x = np.array(lambda_x_list)
        x_num = 2*np.pi/lambda_x
        lambda_y = 2*np.pi*radius/y_wavenumbers # Wavelength in y-direction at input latitude
        y_num = 2*np.pi/lambda_y
        
        c_phase_list = []
        for j in range(len(zmzw)):
            c_phase_item = (zmzw[j] - (beta/(x_num[j]**2+y_num[j]**2)))
            c_phase_list.append(c_phase_item)
        c_phase = np.array(c_phase_list)

        # mean = np.mean(c_phase[:,80:90])
        # textstr= r'Mean$=%.2f$' % mean
        # props = dict(boxstyle='square', facecolor='white', alpha=1.0)
        
        fig, ax = plt.subplots(figsize=(10,6))
        plt.contourf(c_phase.T, cmap=redblu, norm=TwoSlopeNorm(0))
        plt.title('Rossby wave phase velocity')
        plt.xlabel('Time [days]')
        plt.ylabel('Latitude [degrees]')
        plt.yticks((0,22, 45, 67, 90), ('90S', '45S', '0', '45N', '90N'))
        mbar = plt.colorbar(pad=0.1)
        # mbar.set_ticks(np.arange(-100, 100, 20))
        mbar.ax.set_title('m/s')
        # ax.text(0.05, 0.2, textstr, transform=ax.transAxes, fontsize=10, verticalalignment='top', bbox=props)
        
        if save == 'yes':
            plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/cloudproject/epsfigs/cphase_hov_%s.eps' %name, format='eps', bbox_inches='tight')
        else:
            pass  
        plt.show()
    
    
    

def zmzwplots(datalist, ndata=5, start=0, end=300, level=8, lats=(55,85), save='no'):
    
    names = ['Control ProxB',
                'Warm ProxB','Control TRAP1-e','Warm TRAP1-e','Dry TRAP1-e']
    
    output_list = []
    mean_list = []
    
    for i in range(ndata):
       
        data = datalist[i]
        name = names[i]
       
        for cube in data:
            if cube.standard_name == 'x_wind':
                x_wind = cube[start:end,level,:,:].copy()
               
        km_heights = np.round(x_wind.coord('level_height').points*1e-03,2)
#        lat_deg = np.abs(lat-45)*2 # Convert input row number to latitude in degrees north
        if x_wind.coord('latitude').bounds == None:
            x_wind.coord('latitude').guess_bounds()
        if x_wind.coord('longitude').bounds == None:
            x_wind.coord('longitude').guess_bounds()
            
        lat_band = x_wind.intersection(latitude=(lats[0],lats[1]))
        latpoints = np.round(lat_band.coord('latitude').points)
        lat_grid = iris.analysis.cartography.area_weights(lat_band)
    
        zmzw = lat_band.collapsed(['latitude','longitude'],iris.analysis.MEAN, weights=lat_grid)  
    
 #       zmzw = x_wind[:,level,lat].collapsed('longitude',iris.analysis.MEAN)
        zmzw = zmzw.data
        mean_zmzw = np.round(np.mean(zmzw,axis=0))
        mean_list.append(mean_zmzw)
        print(mean_zmzw)
        output_list.append(zmzw)
    
    plt.subplots(figsize=(12,8))
    plt.plot(output_list[0], color='b', label=names[0] + ', Mean: ' + str(mean_list[0]))
    plt.plot(output_list[1], color='r', label=names[1] + ', Mean: ' + str(mean_list[1]))
    plt.plot(output_list[2], color='g', label=names[2] + ', Mean: ' + str(mean_list[2]))
    plt.plot(output_list[3], color='k', label=names[3] + ', Mean: ' + str(mean_list[3]))
    plt.plot(output_list[4], color='y', label=names[4] + ', Mean: ' + str(mean_list[4]))
    plt.title('Zonal mean zonal wind from %sN to %sN' %(int(latpoints[0]), int(latpoints[-1])), fontsize=16)
    plt.ylabel('Wind velocity [m/s]')
    plt.xlabel('Time [days]')
    plt.legend()
    
    if save == 'yes':
        plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/cloudproject/epsfigs/zmzw_all.eps', format='eps', bbox_inches='tight')
    else:
        pass  
    plt.show()
    
def vertical_profile(cubes,start=500,end=600,level=8,top_level=30,select='absolute',save='no'):

    """ This function plots the vertical profile of the chosen data input over time.
     Possible inputs: air temperature ('absolute'), potential temperature ('potential'),
     air pressure ('pressure'), vertical wind ('z_wind'), cloud ('cloud') """

    for cube in cubes:
        if cube.standard_name == 'air_potential_temperature':
            theta = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'air_pressure':
            pressure = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'upward_air_velocity':
            z_wind = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'mass_fraction_of_cloud_ice_in_air':
            ice = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'mass_fraction_of_cloud_liquid_water_in_air':
            liq = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'x_wind':
            x_wind = cube[start:end,:,:,:].copy()

    heat = mpl_cm.get_cmap('gist_heat')
    blues = mpl_cm.get_cmap('Blues')
    redblu = mpl_cm.get_cmap('RdBu')

    p0 = iris.coords.AuxCoord(
        100000.0, long_name='reference_pressure', units='Pa')
    p0.convert_units(pressure.units)
    temperature = theta*((pressure/p0)**(287.05/1005)) # R and cp in J/kgK for 300K

    if select == 'absolute':
        datacube = temperature.copy()
        titleterm = 'absolute temperature'
        y_axis = 'Temperature [K]'
        colors = heat
        unit = 'K'
        norm = None
    elif select == 'potential':
        datacube = theta.copy()
        titleterm = 'potential temperature'
        y_axis = 'Temperature [K]'
        colors = heat
        unit = 'K'
        norm = None
    elif select == 'pressure':
        datacube = pressure.copy()
        titleterm = 'air pressure'
        y_axis = 'Pressure [Pa]'
        colors = blues
        unit = 'Pa'
        norm = None
    elif select == 'z_wind':
        datacube = (z_wind.copy())*1e4
        titleterm = 'vertical wind'
        y_axis = 'Wind speed [m/s]'
        colors = redblu
        unit = '$10^{-4}$ m/s'
        norm = TwoSlopeNorm(0)
    elif select == 'cloud':
        datacube = (ice.copy() + liq.copy())*1e6
        titleterm = 'total cloud (ice and liquid)'
        y_axis = 'Cloud mass [kg/kg]'
        colors = blues
        unit = '$10^{-6}$ kg/kg'
        norm = None
    elif select == 'x_wind':
        datacube = x_wind.copy()
        titleterm = 'zonal wind'
        y_axis = 'Wind speed [m/s]'
        colors = redblu
        unit = 'm/s'
        norm = TwoSlopeNorm(0)

    heights = np.round(datacube.coord('level_height').points*1e-03,0)
    time_axis = np.arange(0,datacube.shape[0])
    lats = datacube.coord('latitude')
    lons = datacube.coord('longitude')

    if lats.bounds == None:
        datacube.coord('latitude').guess_bounds()
    if lons.bounds == None:
        datacube.coord('longitude').guess_bounds()

    datacube = datacube.extract(iris.Constraint(longitude=lambda v: 270 < v <= 359 or 0 <= v <= 90, latitude=lambda v: -90 <= v <= 90))
    grid_areas = iris.analysis.cartography.area_weights(datacube)
    dayside_mean = datacube.collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=grid_areas)
    z_axis = heights[:top_level]
    dayside_time = dayside_mean[:,:top_level].data

    fig1, ax1 = plt.subplots(figsize=(8,5))
    ax1.set_xlabel('Time [days]')
    ax1.set_ylabel('%s' %y_axis)
    ax1.set_title('Dayside mean %s at h=%s km' %(titleterm, heights[level]))
    plt.plot(dayside_mean[:,level].data)
    plt.show()

    fig2, ax2 = plt.subplots(figsize=(8,5))
    ax2.set_xlabel('Time [days]')
    ax2.set_ylabel('Height [km]')
    ax2.set_title('Dayside mean %s' %titleterm)
    plt.contourf(time_axis, z_axis, dayside_time.T,cmap=colors,norm=norm)
    cbar = plt.colorbar(pad=0.1)
    cbar.ax.set_title('%s' %unit)

    if save == 'yes':
        plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/cloudproject/epsfigs/%s_%s_%s_trap.eps' %(select,start,end), format='eps')
    else:
        pass

    plt.show() 

def swresonance(cubes,start=0,end=100,level=0,save='no'):

    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'y_wind':
            y_wind = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'surface_net_downward_shortwave_flux':
            heat = cube[start:end,:,:,:].copy()

    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
    km_heights = np.round(x_wind.coord('level_height').points*1e-03,2)
    
    winds = windspharm.iris.VectorWind(x_wind[:,level,:,:],y_wind[:,level,:,:])
    # Create a VectorWind data object from the x and y wind cubes
    uchi,vchi,upsi,vpsi = winds.helmholtz(truncation=21)
    # Calculate the Helmholtz decomposition. Truncation is set to 21 because 
    # this is what Hammond and Lewis 2021 used.
    
    zonal_upsi = upsi.collapsed('longitude', iris.analysis.MEAN)
    zonal_vpsi = vpsi.collapsed('longitude', iris.analysis.MEAN)
    # Calculate zonal means of the x and y components of the rotational component
    eddy_upsi = upsi - zonal_upsi
    eddy_vpsi = vpsi - zonal_vpsi
    magnitude = np.sqrt(eddy_upsi.data**2 + eddy_vpsi.data**2)
    
    fft2 = sp.fft.fftshift(sp.fftpack.fft2(sp.fft.ifftshift(magnitude)))
    yfreqs = sp.fft.fftshift(sp.fft.fftfreq(fft2.shape[1],d=1./90))
    xfreqs = sp.fft.fftshift(sp.fft.fftfreq(fft2.shape[2],d=1./144))
    psd = np.abs(fft2)**2

    one_zero = psd[:,45,73]
    print(xfreqs[73],yfreqs[45])
    two_one = psd[:,46,74]
    two_two = psd[:,47,74]
    three_two = psd[:,47,75]
    one_one = psd[:,46,73]
    wave_sum = two_one + two_two + three_two + one_one

    lats = heat.coord('latitude')
    lons = heat.coord('longitude')

    if lats.bounds == None:
        heat.coord('latitude').guess_bounds()
    if lons.bounds == None:
        heat.coord('longitude').guess_bounds()

    grid_areas = iris.analysis.cartography.area_weights(heat)
    global_mean = heat.collapsed(['latitude','longitude'],iris.analysis.MEAN, weights=grid_areas)

    fig, ax1 = plt.subplots()
    ax1.set_xlabel('Time [days]')
    ax1.set_ylabel('Surface heating [W/m2]')
    ax1.plot(global_mean.data,color='k',label='Heating')

    ax2 = ax1.twinx()
    ax2.set_ylabel('PSD')
    ax2.plot(one_zero,color='r',label='1-0 wave')
    ax2.plot(wave_sum,color='b',label='Wave sum')
    ax2.ticklabel_format(axis='x',style='sci')
    plt.legend()

    plt.title('Mean SW heating and Rossby waves at h=%s km' %km_heights[level],y=1.05)


    if save == 'yes':
        plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/cloudproject/epsfigs/swresonance_dry.eps', format='eps')
    else:
        pass
    plt.show()

def cloud_bubble(cubes,lat=45,lon=72,start=0,end=120,final_height=35, save='no'):
    
    """ 25 km = 35 for ProxB, 25 for TRAP1-e"""
    
    for cube in cubes:
        if cube.standard_name == 'mass_fraction_of_cloud_ice_in_air':
            ice = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'mass_fraction_of_cloud_liquid_water_in_air':
            liq = cube[start:end,:,:,:].copy()
    
    ice = ice.collapsed('time',iris.analysis.MEAN)
    liq = liq.collapsed('time',iris.analysis.MEAN)
    
    heights = np.round(ice.coord('level_height').points*1e-03,0)
    lons = ice.coord('longitude').points
    lats = ice.coord('latitude').points
    total_cloud = ice + liq
    
    # if len(heights) > 39:
    #     final_height = 43
    # else:
    #     final_height = -1
    
    fig1, ax1 = plt.subplots(figsize=(5,5))
    plota = ax1.contourf(np.roll(lons,72), heights[:final_height], total_cloud[:final_height,lat,:].data*10**4, np.arange(0,3.5,0.1), cmap='Blues')
    ax1.set_title('Total cloud at latitude %s' %lats[lat])
    ax1.set_xlabel('Longitude [degrees]')
    ax1.set_xticks([0,90,180,270,360], ['180W','90W','0','90E','180E'])
    ax1.set_ylabel('Height [km]')
    cba = plt.colorbar(plota)
    cba.ax.set_title('$10^{-4}$ kg/kg', size=10)
    if save == 'yes':
        plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/cloudproject/epsfigs/cloud_bubble_lat_%s_prox.eps' %lats[lat], format='eps', bbox_inches='tight')
    else:
        pass
    plt.show()
    
    fig2, ax2 = plt.subplots(figsize=(5,5))
    plotb = ax2.contourf(lats, heights[:final_height], total_cloud[:final_height,:,lon+71].data*10**4, np.arange(0,3.5,0.1), cmap='Blues')
    ax2.set_title('Total cloud at longitude %s' %np.round(lons[lon+71],0))
    ax2.set_xlabel('Latitude [degrees]')
    ax2.set_xticks([-90,-60,-30,0,30,60,90], ['90S','60S','30S','0','30N','60N','90N'])
    ax2.set_ylabel('Height [km]')
    cbb = plt.colorbar(plotb)
    cbb.ax.set_title('$10^{-4}$ kg/kg', size=10)

    if save == 'yes':
        plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/cloudproject/epsfigs/cloud_bubble_lon_%s_prox.eps' %lons[lon], format='eps', bbox_inches='tight')
    else:
        pass
    plt.show()
    
    
def cloud_type(cubes, start=0,end=600, long1=36, long2=108, filtering=True, 
               ice_threshold = 0.1, liq_threshold = 0.1, save='no'):
    
    for cube in cubes:
        if cube.standard_name == 'mass_fraction_of_cloud_ice_in_air':
            ice_condensate_raw = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'mass_fraction_of_cloud_liquid_water_in_air':
            liquid_condensate_raw = cube[start:end,:,:,:].copy()
    
    if long1==36 and long2==108:
        titleloc = 'Terminators'
    elif long1==143 and long2==0:
        titleloc = 'Substellar Point'
    else:
        titleloc = '20W'
    
    ice_condensate = ice_condensate_raw.data
    liquid_condensate = liquid_condensate_raw.data

    y_axis = np.round(liquid_condensate_raw.coord('level_height').points*1e-03,0)
    x_axis = liquid_condensate_raw.coord('latitude').points
    time_axis = np.arange(0,liquid_condensate.shape[0])
    heights = np.array(ice_condensate_raw.coord('level_height').points)
        
    ice_east = np.mean(ice_condensate[:,:,:,long1], axis=2)
    ice_east = np.average(ice_east,axis=1,weights=heights)
    ice_west = np.mean(ice_condensate[:,:,:,long2], axis=2)
    ice_west = np.average(ice_west,axis=1,weights=heights)
    
    liq_east = np.mean(liquid_condensate[:,:,:,long1], axis=2)
    liq_east = np.average(liq_east,axis=1,weights=heights)
    liq_west = np.mean(liquid_condensate[:,:,:,long2], axis=2)
    liq_west = np.average(liq_west,axis=1,weights=heights)
    
    # ice_east = np.mean(ice_condensate[:,:,:,long1], axis=(1,2))
    # ice_west = np.mean(ice_condensate[:,:,:,long2], axis=(1,2))
    total_ice = (ice_east + ice_west)/2    
    
    # liq_east = np.mean(liquid_condensate[:,:,:,long1], axis=(1,2))
    # liq_west = np.mean(liquid_condensate[:,:,:,long2], axis=(1,2))
    total_liq = (liq_east + liq_west)/2
    
    if filtering == True:
        run_length = total_ice.shape[0]
        fft_ice = sp.fftpack.fft(total_ice)
        psd_ice = np.abs(fft_ice)**2
        # freqs_ice = sp.fftpack.fftfreq(len(psd_ice), 1./run_length)
        freqs_ice = sp.fftpack.fftfreq(len(psd_ice), 1./1)
        i = freqs_ice > 0
        lowpass_ice = fft_ice.copy()
        lowpass_ice[np.abs(freqs_ice) > ice_threshold] = 0
        total_ice = np.real(sp.fftpack.ifft(lowpass_ice))
        
        fft_liq = sp.fftpack.fft(total_liq)
        psd_liq = np.abs(fft_liq)**2
        # freqs_liq = sp.fftpack.fftfreq(len(psd_liq), 1./run_length)
        freqs_liq = sp.fftpack.fftfreq(len(psd_liq), 1./1)
        j = freqs_liq > 0
        lowpass_liq = fft_liq.copy()
        lowpass_liq[np.abs(freqs_liq) > liq_threshold] = 0
        total_liq = np.real(sp.fftpack.ifft(lowpass_liq))
        
        # fig, ax = plt.subplots(1,1,figsize=(8,4))
        # ax.plot(freqs_ice[i], psd_ice[i],color='r')
        # ax.plot(freqs_liq[j], psd_liq[j],color='b')
        # ax.set_xlim(0,150)
        # ax.set_xlabel('Frequency')
        # ax.set_ylabel('PSD')
        # ax.set_title('Periodicity of cloud cover')
        # print(np.argmax(psd_liq[i]))
        # print(np.argmax(psd_ice[i]))
    
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('Time [days]')
    ax1.set_ylabel('Ice condensate [kg/kg]')
    ax1.plot(time_axis, total_ice, color='b', label='Ice')
    ax1.tick_params(axis='y', labelcolor='b')
    
    ax2 = ax1.twinx()
    ax2.set_ylabel('Liquid condensate [kg/kg]')
    ax2.plot(time_axis, total_liq, color='r', label='Liquid')
    ax2.tick_params(axis='y', labelcolor='r')
    
    plt.title('Mean Ice and Liquid Condensate at %s' %titleloc)
    fig.tight_layout()
    plt.show()
    
    plt.plot(time_axis, total_ice, color='b', label='Ice')
    plt.plot(time_axis, total_liq, color='r', label='Liquid')
    plt.title('Mean Cloud Condensate at %s' %titleloc)
    plt.xlabel('Time [days]')
    plt.ylabel('Cloud condensate [kg/kg]')
    plt.legend()
    
    if save == 'yes':
        plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/cloudproject/epsfigs/cloudvari_lon_%s_%s_warmtrap.eps' %(start,end), format='eps')
    else:
        pass
    plt.show()
    
def hovmoeller_rwaves(cubes, start=0,end=100,level=8,lats=(55,85),title='trap',save='no'):

    for cube in cubes:
        if cube.standard_name == 'y_wind':
            y_wind = cube[start:end,level,:,:].copy()
    
    time_axis = np.arange(0,y_wind.shape[0])
    lons = y_wind.shape[2]/2
    if y_wind.coord('latitude').bounds == None:
        y_wind.coord('latitude').guess_bounds()
    if y_wind.coord('longitude').bounds == None:
        y_wind.coord('longitude').guess_bounds()
        
    lat_band = y_wind.intersection(latitude=(lats[0],lats[1]))
#    lat_grid = iris.analysis.cartography.area_weights(lat_band)
    
    band_mean = lat_band.collapsed('latitude',iris.analysis.MEAN)
    band_mean = band_mean.data
    
    plt.subplots(figsize=(8,6))
    plt.contourf(np.arange(-lons, lons), time_axis, np.roll(band_mean, 72, axis=1), levels=np.arange(-60, 61, 10), cmap=redblu, norm=TwoSlopeNorm(0))
    plt.title('Mean meridional wind from %s to %s N' %(lats[0],lats[1]))
    plt.xlabel('Longitude [degrees]')
    plt.ylabel('Time [days]')
    plt.xticks((-72, -48, -24, 0, 24, 48, 72), ('180W',
                 '120W', '60W', '0', '60E', '120E', '180E'))
    cbar = plt.colorbar(pad=0.1)
    cbar.set_ticks(np.arange(-60, 61, 10))
    cbar.ax.set_title('m/s')
    
    if save == 'yes':
        plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/cloudproject/epsfigs/hov_rwaves%sto%s_%s.eps' %(start,end, title), format='eps', bbox_inches='tight')
    else:
        pass
        
    plt.show()
    
    return


def rossby_no(datalist, ndata=4, lat=60):
    
    names = ['Control ProxB',
                'Warm ProxB','Control TRAP1-e','Warm TRAP1-e']
    
    omegas = [6.46e-06, 7.93e-06,1.19e-05, 1.75e-05]
    radii = [7160000, 7160000, 5797818, 5797818]
    
    rno_list = []
    
    for i in range(ndata):
       
        data = datalist[i]
        name = names[i]
        omega = omegas[i]
        radius = radii[i]
       
        for cube in data:
            if cube.standard_name == 'x_wind':
                x_wind = cube.copy()
                
        lats = x_wind.coord('latitude')
        lons = x_wind.coord('longitude')
        
        if lats.bounds == None:
            x_wind.coord('latitude').guess_bounds()
        if lons.bounds == None:
            x_wind.coord('longitude').guess_bounds()
        
        grid = iris.analysis.cartography.area_weights(x_wind)
        mean_wind = x_wind.collapsed(['longitude','latitude','level_height','time'], iris.analysis.MEAN, weights=grid)
        U = mean_wind.data
        L = np.pi*radius
        f = 2*omega*np.sin(lat*np.pi/180)
        
        rno = U/(L*f)
        print(name + f' has a Rossby number of {rno}')
        
        rno_list.append(rno)
    
    rno_df = pd.DataFrame(index=range(0,len(names)),columns = ['Name', 'RossbyNo'])
    rno_df['Name'] = names
    rno_df['RossbyNo'] = rno_list
    print(rno_df)
    rno_df.to_csv('/exports/csce/datastore/geos/users/s1144983/papers/cloudproject/rnos.csv')

    
def temp_gradient(datalist,start=0, end=600, ndata=4,lats=(0,90),level=(8,9)):
    
    names = ['Control ProxB',
                'Warm ProxB','Control TRAP1-e','Warm TRAP1-e']
    
    omegas = [6.46e-06, 7.93e-06,1.19e-05, 1.75e-05]
    radii = [7160000, 7160000, 5797818, 5797818]
    
    therms = []
    grads = []
    
    for i in range(ndata):
       
        data = datalist[i]
        name = names[i]
        omega = omegas[i]
        radius = radii[i]
       
        for cube in data:
            if cube.standard_name == 'air_potential_temperature':
                theta = cube[start:end,:,:,:].copy()
            if cube.standard_name == 'air_pressure':
                pressure = cube[start:end,:,:,:].copy()
            if cube.standard_name == 'x_wind':
                x_wind = cube[start:end,:,:,:].copy()
            if cube.standard_name == 'y_wind':
                y_wind = cube[start:end,:,:,:].copy()
                
        p0 = iris.coords.AuxCoord(
            100000.0, long_name='reference_pressure', units='Pa')
        p0.convert_units(pressure.units)
        temperature = theta*((pressure/p0)**(287.05/1005)) # R and cp in J/kgK for 300K
        temp_gradient = iris.analysis.calculus.differentiate(temperature, 'latitude')
        lat_band = temp_gradient.intersection(latitude=(lats[0],lats[1]))

        latitudes = lat_band.coord('latitude')
        longitudes = lat_band.coord('longitude')
        
        if latitudes.bounds == None:
            lat_band.coord('latitude').guess_bounds()
        if longitudes.bounds == None:
            lat_band.coord('longitude').guess_bounds()
        
        grid = iris.analysis.cartography.area_weights(lat_band)

        mean_grad = lat_band[:,level[0]:level[1],:,:].collapsed(['longitude','latitude','level_height','time'], iris.analysis.MEAN, weights=grid[:,level[0]:level[1],:,:])       
        mean_grad = mean_grad/(0.5*np.pi*radius/(lats[1]-lats[0]))
        flist = []
        for lat_item in range(lats[0], lats[1]):
            
            f = 2*omega*np.sin(lat_item*np.pi/180)
            flist.append(f)
        
        mean_f = np.mean(np.array(flist))
        
        platitudes = pressure.coord('latitude')
        plongitudes = pressure.coord('longitude')
        
        if platitudes.bounds == None:
            pressure.coord('latitude').guess_bounds()
        if plongitudes.bounds == None:
            pressure.coord('longitude').guess_bounds()
        
        pgrid = iris.analysis.cartography.area_weights(pressure)
        p1 = pressure[:,level[0]:level[1],:,:].collapsed(['longitude','latitude','level_height','time'], iris.analysis.MEAN, weights=pgrid[:,level[0]:level[1],:,:])
        p2 = pressure[:,level[0]+1:level[1]+1,:,:].collapsed(['longitude','latitude','level_height','time'], iris.analysis.MEAN, weights=pgrid[:,level[0]+1:level[1]+1,:,:])
        print(p1.data, p2.data, mean_f)
        thermal = -(287.05/mean_f)*(np.log(p1.data/p2.data))*mean_grad.data
        
        grads.append(mean_grad.data)
        therms.append(thermal)
        
        print(name + f': {mean_grad.data}')
        print(name + f': {thermal}')
        
    thermal_df = pd.DataFrame(index=range(0,len(names)), columns=['Name','Gradient','Thermal wind'])
    thermal_df['Name'] = names
    thermal_df['Gradient'] = grads
    thermal_df['Thermal wind'] = therms

    print(thermal_df)
    thermal_df.to_csv('/exports/csce/datastore/geos/users/s1144983/papers/cloudproject/thermal_wind.csv')
    return
    


def clim_stats(datalist,ndata=5,start=0,end=300):
    
    col_names = ['Name','Mean t anti', 'Max t anti','Min t anti', 
                 'Mean t sub', 'Max t sub', 'Min t sub',
                 'Mean ZMZW', 'Max ZMZW', 'Min ZMZW',
                 'Mean h anti', 'Max h anti', 'Min h anti',
                 'Mean h sub', 'Max h sub', 'Min h sub', 
                 'Mean surf t', 'Max surf t', 'Min surf t']

    stats = pd.DataFrame(columns=col_names)
    
    names = ['Control ProxB',
                'Warm ProxB','Control TRAP1-e','Warm TRAP1-e','Dry TRAP1-e']
    
    
    for i in range(ndata):
        
        data = datalist[i]
        name = names[i]
        
        for cube in data:
            if cube.standard_name == 'air_potential_temperature':
                potential_temp = cube[start:end,:,:,:].copy()
            if cube.standard_name == 'air_pressure':
                air_pressure = cube[start:end,:,:,:].copy()
            if cube.standard_name == 'x_wind':
                x_wind = cube[start:end,:,:,:].copy()
            if cube.standard_name == 'y_wind':
                y_wind = cube[start:end,:,:,:].copy()
            if cube.standard_name == 'surface_temperature':
                surface_temp = cube[start:end,:,:].copy()
            if cube.standard_name == 'specific_humidity' and i < ndata-1:
                spec_humidity = cube[start:end,:,:,:].copy()
            elif cube.standard_name == 'specific_humidity' and i == ndata-1: 
                spec_humidity = cube.copy()*0.0
           

        y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
        
        vertical = [('level_height', x_wind.coord('level_height').points)]
        potential_temp = potential_temp.regrid(x_wind, iris.analysis.Linear())
        potential_temp = potential_temp.interpolate(vertical, iris.analysis.Linear())
        air_pressure = air_pressure.regrid(x_wind, iris.analysis.Linear())
        air_pressure = air_pressure.interpolate(vertical, iris.analysis.Linear())

        heights = np.round(x_wind.coord('level_height').points*1e-03,0)
        lats = np.round(x_wind.coord('latitude').points,0)
        lons = np.round(x_wind.coord('longitude').points,0)
        
        p0 = iris.coords.AuxCoord(100000.0, long_name='reference_pressure', units='Pa')
        p0.convert_units(air_pressure.units)
        absolute_temp = potential_temp*((air_pressure/p0)**(287.05/1005))
        absolute_temp = np.mean(absolute_temp.data,axis=0)
        
        zmzw = x_wind.collapsed('longitude',iris.analysis.MEAN)
        zmzw = np.mean(zmzw.data,axis=0)
        
        spec_humidity = spec_humidity.regrid(x_wind, iris.analysis.Linear())
        spec_humidity = spec_humidity.interpolate(vertical, iris.analysis.Linear())
        spec_humidity = np.mean(spec_humidity.data,axis=0)
        
        surface_temp = np.mean(surface_temp.data,axis=0)
        
        meaned_x = np.mean(x_wind.data,axis=0)
        meaned_y = np.mean(y_wind.data,axis=0)
        
        sub_temp = absolute_temp[:,45,0]
        anti_temp = absolute_temp[:,45,72]
        sub_hum = spec_humidity[:,45,0]
        anti_hum = spec_humidity[:,45,72]
        
        zmzw_mean, zmzw_max, zmzw_min = np.mean(zmzw), np.max(zmzw), np.min(zmzw)
        st_mean, st_max, st_min = np.mean(surface_temp), np.max(surface_temp), np.min(surface_temp)
        
        values = [name, np.mean(anti_temp), np.max(anti_temp), np.min(anti_temp),
                  np.mean(sub_temp), np.max(sub_temp), np.min(sub_temp),
                  zmzw_mean, zmzw_max, zmzw_min,
                  np.mean(anti_hum), np.max(anti_hum), np.min(anti_hum),
                  np.mean(sub_hum), np.max(sub_hum), np.min(sub_hum),
                  st_mean, st_max, st_min]
        
        row = pd.DataFrame(data=values)
        print(row)
        
        stats.loc[len(stats.index)] = values
        

    print(stats)
    stats.to_csv('/exports/csce/datastore/geos/users/s1144983/papers/cloudproject/clim_stats.csv')
        
        
    return