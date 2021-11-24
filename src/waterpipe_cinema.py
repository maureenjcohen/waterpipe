#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 13:10:07 2021

@author: Mo Cohen

Pipeline for post-processing UM output data.
- Makes gifs of dynamical quantities

Model: University of Exeter Stand Alone model of Proxima Centauri b, UM vn11.8

"""

import iris
import iris.plot as iplt
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import gif
import windspharm
import numpy as np

brewer_red = mpl_cm.get_cmap('brewer_YlOrRd_09')
brewer_bg = mpl_cm.get_cmap('brewer_PuBu_09')
brewer_reds = mpl_cm.get_cmap('brewer_Reds_09')
brewer_redblu = mpl_cm.get_cmap('brewer_RdBu_11')

savepath = 'placeholder'

""" CODE BLOCK FOR ZONAL MEAN ZONAL WINDS"""

def calculate_zonal_wind(cubes):
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube.copy()
            
    dayside = x_wind.extract(iris.Constraint(longitude=lambda v: 270 < v <= 359 or 0 <= v <= 90, latitude=lambda v: -90 <= v <= 90))
    nightside = x_wind.extract(iris.Constraint(longitude=lambda v: 90 < v <= 270, latitude=lambda v: -90 <= v <= 90))
    
    dayside_zonal_mean = dayside.collapsed('longitude', iris.analysis.MEAN)
    nightside_zonal_mean = nightside.collapsed('longitude', iris.analysis.MEAN)
    
    return dayside_zonal_mean, nightside_zonal_mean

def calculate_120day_mean(cubes):
    
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube.copy()
            
    dayside = x_wind.extract(iris.Constraint(longitude=lambda v: 270 < v <= 359 or 0 <= v <= 90, latitude=lambda v: -90 <= v <= 90))
    nightside = x_wind.extract(iris.Constraint(longitude=lambda v: 90 < v <= 270, latitude=lambda v: -90 <= v <= 90))
    
    dayside = dayside.collapsed('longitude', iris.analysis.MEAN)
    nightside = nightside.collapsed('longitude', iris.analysis.MEAN)
    
    dayside_meaned = np.mean(dayside.data.reshape(-1,4,60,90),axis=1)
    nightside_meaned = np.mean(nightside.data.reshape(-1,4,60,90),axis=1)
    dayside[-15:,:,:].data = dayside_meaned
    nightside[-15:,:,:].data = nightside_meaned
    
    dayside_zonal_mean = dayside[-15:,:,:]
    nightside_zonal_mean = nightside[-15:,:,:]
    
    return dayside_zonal_mean, nightside_zonal_mean
    
@gif.frame
def gif_dayside_wind(i, dayside_zonal_mean):
        
    CS_day = iplt.contourf(dayside_zonal_mean[i,:,:], levels=np.linspace(-100,140,20), cmap=brewer_redblu)
#    iplt.contourf(dayside_zonal_mean[i,:,:], brewer_redblu.N, cmap=brewer_redblu)

    plt.title('Dayside Zonal Mean Zonal Wind [m s-1], 4month %s' %(i+1), y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.clabel(CS_day, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    
dayside_frames = []
for i in range(15):
    dayside_frame = gif_dayside_wind(i, dayside_zonal_mean)
    dayside_frames.append(dayside_frame)

gif.save(dayside_frames, str(savepath) + '/dayside_zonal_120daymeanwind_fixed.gif', duration = 30, unit = 's', between='startend')

@gif.frame
def gif_nightside_wind(i, nightside_zonal_mean):  
    
    CS_night = iplt.contourf(nightside_zonal_mean[i,:,:], levels=np.linspace(-100,140,20), cmap=brewer_redblu)
    plt.title('Nightside Zonal Mean Zonal Wind [m s-1], 4month %s' %(i+1), y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.clabel(CS_night, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)

nightside_frames = []
for i in range(15):
    nightside_frame = gif_nightside_wind(i, nightside_zonal_mean)
    nightside_frames.append(nightside_frame)
    
gif.save(nightside_frames, str(savepath) + '/nightside_zonal_120daymeanwind_fixed.gif', duration = 30, unit = 's', between='startend')

""" CODE BLOCK FOR ZONAL MEAN SPECIFIC HUMIDITY """

def calculate_zonal_humidity(cubes):
    
    for cube in cubes:
        if cube.standard_name == 'specific_humidity':
            spec_humidity = cube.copy() 
            
    dayside_humid = spec_humidity.extract(iris.Constraint(longitude=lambda v: 270 < v <= 359 or 0 <= v <= 90, latitude=lambda v: -90 <= v <= 90))
    nightside_humid = spec_humidity.extract(iris.Constraint(longitude=lambda v: 90 < v <= 270, latitude=lambda v: -90 <= v <= 90))
    
    dayside_zonal_humid = dayside_humid.collapsed('longitude', iris.analysis.MEAN)
    nightside_zonal_humid = nightside_humid.collapsed('longitude', iris.analysis.MEAN)
    
    return dayside_zonal_humid, nightside_zonal_humid

@gif.frame
def gif_dayside_humidity(i, dayside_zonal_humid):
        
    CS_day = iplt.contourf(dayside_zonal_humid[i,:,:], brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Dayside Zonal Mean Humidity [kg kg-1], day %s' %(i+1), y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.clabel(CS_day, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    
dayside_frames = []
for i in range(1710):
    dayside_frame = gif_dayside_humidity(i, dayside_zonal_humid)
    dayside_frames.append(dayside_frame)

gif.save(dayside_frames, str(savepath) + '/dayside_zonal_humidity.gif', duration = 900, unit = 's', between='startend')


@gif.frame
def gif_nightside_humidity(i, nightside_zonal_humid):
        
    CS_day = iplt.contourf(nightside_zonal_humid[i,:,:], brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Nightside Zonal Mean Humidity [kg kg-1], day %s' %(i+1), y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.clabel(CS_day, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    
nightside_frames = []
for i in range(1800):
    nightside_frame = gif_nightside_humidity(i, nightside_zonal_humid)
    nightside_frames.append(nightside_frame)

gif.save(nightside_frames, str(savepath) + '/nightside_zonal_humidity.gif', duration = 900, unit = 's', between='startend')

""" CODE BLOCK FOR ZONAL MEAN ABSOLUTE AIR TEMPERATURE """

def calculate_air_temp(cubes):
    
    for cube in cubes:
        if cube.standard_name == 'air_potential_temperature':
            potential_temp = cube.copy()
        if cube.standard_name == 'air_pressure':
            air_pressure = cube.copy()

    p0 = iris.coords.AuxCoord(100000.0, long_name='reference_pressure', units='Pa')
    p0.convert_units(air_pressure.units)
    absolute_temp = potential_temp*((air_pressure/p0)**(287.05/1005)) # R and cp in J/kgK for 300K
    
    dayside = absolute_temp.extract(iris.Constraint(longitude=lambda v: 270 < v <= 359 or 0 <= v <= 90, latitude=lambda v: -90 <= v <= 90))
    nightside = absolute_temp.extract(iris.Constraint(longitude=lambda v: 90 < v <= 270, latitude=lambda v: -90 <= v <= 90))
    
    dayside_absolute_temp = dayside.collapsed('longitude', iris.analysis.MEAN)
    nightside_absolute_temp = nightside.collapsed('longitude', iris.analysis.MEAN)
    
    return dayside_absolute_temp, nightside_absolute_temp

@gif.frame
def gif_dayside_temp(i, dayside_absolute_temp):
        
    CS_day = iplt.contourf(dayside_absolute_temp[i,:,:], brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Dayside Zonal Mean Air Temp [K], month %s' %(i+1), y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.clabel(CS_day, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    
dayside_frames = []
for i in range(60):
    dayside_frame = gif_dayside_temp(i, dayside_absolute_temp)
    dayside_frames.append(dayside_frame)

gif.save(dayside_frames, str(savepath) + '/dayside_zonal_airtemp.gif', duration = 30, unit = 's', between='startend')

@gif.frame
def gif_nightside_temp(i, nightside_absolute_temp):  
    
    CS_night = iplt.contourf(nightside_absolute_temp[i,:,:], brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Nightside Zonal Mean Air Temp [K], month %s' %(i+1), y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.clabel(CS_night, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)

nightside_frames = []
for i in range(60):
    nightside_frame = gif_nightside_temp(i, nightside_absolute_temp)
    nightside_frames.append(nightside_frame)
    
gif.save(nightside_frames, str(savepath) + '/nightside_zonal_airtemp.gif', duration = 30, unit = 's', between='startend')

""" CODE BLOCK FOR STREAMFUNCTION """

for cube in slow:
    if cube.standard_name == 'x_wind':
        x_wind = cube.copy()
    if cube.standard_name == 'y_wind':
        y_wind = cube.copy()
   
y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())

@gif.frame
def calculate_streamfunction(x_wind, y_wind, level, time_slice=-1):

    wind = windspharm.iris.VectorWind(x_wind, y_wind)
    streamfunction, velpotential = wind.sfvp()
    clevs = [-200, -180, -160, -120, -100, -80, -60, -40, -20, 0, 40, 80, 120, 160, 200]
    fig = iplt.contourf(streamfunction[time_slice,level,:,:]*1e-06, clevs, cmap=brewer_redblu, extend='both')
    plt.title('Streamfunction [$10^6$ m2 s-1], h = %s km' %(level+1), y=1.20)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.colorbar(orientation='horizontal')

frames = []
for level in range(38):
    frame = calculate_streamfunction(x_wind, y_wind, level)
    frames.append(frame)

gif.save(frames, str(savepath) + '/streamfunction.gif', duration = 19, unit = 's', between='startend')

""" CODE BLOCK FOR STREAMLINES """

for cube in control:
    if cube.standard_name == 'x_wind':
        x_wind = cube.copy()
    if cube.standard_name == 'y_wind':
        y_wind = cube.copy()
   
y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
speed = iris.analysis.maths.apply_ufunc(np.sqrt, (x_wind**2 + y_wind**2))
heights = np.round(x_wind.coord('level_height').points*1e-03,0)

x_wind = x_wind.data
y_wind = y_wind.data
speed = speed.data

@gif.frame
def gif_streamlines(x_wind, y_wind, heights, time_slice, level=15):   
    
    X,Y = np.meshgrid(np.arange(-72,72), np.arange(-45,45))
    fig = plt.figure(figsize = (10, 4)) 
    strm = plt.streamplot(X, Y, np.roll(x_wind[time_slice,level,:,:].data, 72, axis=1), np.roll(y_wind[time_slice,level,:,:].data, 72, axis=1), density = 0.5, color=np.roll(speed[time_slice,level,:,:].data, 72, axis=1), cmap=brewer_reds)
    # Since .data method extracts the numpy array and strips the metadata, the longitude/latitude information is lost.
    # To align plot so that (0,0) is at the center as in the Iris plots, use numpy.roll to shift columns (longitude) 180 degrees (72 places)
    cbar = plt.colorbar(strm.lines)
    cbar.set_ticks(np.arange(0,100,10))
    cbar.ax.set_title('m/s')
    plt.title('Wind speed and direction, h=%s km, t=%s months' %(heights[level], time_slice+1))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.xticks((-72,-52,-32,-12,0,12,32,52,72),('180W','140W','100W','60W','0','60E','100E','140E','180E'))
    plt.yticks((-45,-30,-15,0,15,30,45),('90S','60S','30S','0','30N','60N','90N'))    

frames = []
for time_slice in range(60):
    frame = gif_streamlines(x_wind, y_wind, heights, time_slice)
    frames.append(frame)

gif.save(frames, str(savepath) + '/streamlines_15_slow_fixed.gif', duration = 90, unit = 's', between='startend')


""" CODE BLOCK FOR OUTGOING LW RADIATION """

for cube in humid:
    if cube.standard_name == 'toa_outgoing_longwave_flux':
        outgoing_lw = cube.copy()
    if cube.standard_name == 'toa_outgoing_shortwave_flux':
        outgoing_sw = cube.copy()
        

@gif.frame
def gif_longwave(i, outgoing_lw):
    
    iplt.contourf(outgoing_lw[i,:,:], brewer_reds.N, cmap=brewer_reds)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Monthly Mean Outgoing LW Radiation, month %s' %(i+1), y=1.10)
    plt.ylabel('Latitude [degrees]')
    plt.xlabel('Longitude [degrees]')
    plt.colorbar(pad=0.1)
    
run_length = outgoing_lw.shape[0]    
frames = []
for i in range(run_length):
    frame = gif_longwave(i, outgoing_lw)
    frames.append(frame)
    
gif.save(frames, str(savepath) + '/longwave.gif', duration = 30, unit = 's', between='startend')

""" CODE BLOCK FOR COLUMN AIR TEMPERATURE """

for cube in humid:    
    if cube.standard_name == 'air_potential_temperature':
        potential_temp = cube.copy()
    if cube.standard_name == 'air_pressure':
        air_pressure = cube.copy()

p0 = iris.coords.AuxCoord(100000.0, long_name='reference_pressure', units='Pa')
p0.convert_units(air_pressure.units)
absolute_temp = potential_temp*((air_pressure/p0)**(287.05/1005)) # R and cp in J/kgK for 300K

mean_absolute_temp = absolute_temp.collapsed('model_level_number', iris.analysis.MEAN)
    

@gif.frame
def gif_airtemp(i, mean_absolute_temp):
    
    iplt.contourf(mean_absolute_temp[i,:,:], brewer_reds.N, cmap=brewer_reds)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Monthly Mean Air Temperature, month %s' %(i+1), y=1.10)
    plt.ylabel('Latitude [degrees]')
    plt.xlabel('Longitude [degrees]')
    plt.colorbar(pad=0.1)
    
run_length = mean_absolute_temp.shape[0]    
frames = []
for i in range(run_length):
    frame = gif_airtemp(i, mean_absolute_temp)
    frames.append(frame)
    
gif.save(frames, str(savepath) + '/airtemp.gif', duration = 30, unit = 's', between='startend')


""" CODE BLOCK FOR CLOUD COVER """

for cube in fast_only:    
    if cube.long_name == 'cloud_area_fraction_assuming_maximum_random_overlap':
        cloud_cover = cube.copy()
    
@gif.frame
def gif_cloudcover(i, cloud_cover):
    
    iplt.contourf(cloud_cover[i,:,:], brewer_bg.N, cmap=brewer_bg)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Cloud cover, day %s' %(i+1), y=1.10)
    plt.ylabel('Latitude [degrees]')
    plt.xlabel('Longitude [degrees]')
    plt.colorbar(pad=0.1)
    
#run_length = mean_absolute_temp.shape[0]    
frames = []
for i in range(1740,1800):
    frame = gif_cloudcover(i, cloud_cover)
    frames.append(frame)
    
gif.save(frames, str(savepath) + '/clouds_last60days.gif', duration = 30, unit = 's', between='startend')

""" CODE BLOCK FOR QBO """
savepath = '/exports/csce/datastore/geos/users/s1144983/um_data/ch1_control/qbo_gifs/'

for cube in control:
    if cube.standard_name == 'x_wind':
        x_wind = cube.copy()
        
longitudes = x_wind.shape[3]/2     
heights = np.round(x_wind.coord('level_height').points*1e-03,0)
equator = x_wind[:,:,45,:].data

@gif.frame            
def gif_equatorial_wind(i, equator, level):
    
    westerly = np.roll(equator[i,level,:], 72)
    
    plt.figure(figsize=(10,5))
    plt.plot(np.arange(-longitudes, longitudes), westerly)
    plt.title('Zonal Wind at Equator, t=%s months, h=%s km' %(i, heights[level]))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.axvline(x=-36, color='r', linestyle='--')
    plt.axvline(x=36, color='r', linestyle='--')
    plt.ylabel('Wind speed [m s-1]')
    
frames = []
for i in range(0,29):
    frame = gif_equatorial_wind(i, equator, level=47)
    frames.append(frame)
    
gif.save(frames, str(savepath) + '/equatorial_winds.gif', duration = 56, unit = 's', between='startend')
