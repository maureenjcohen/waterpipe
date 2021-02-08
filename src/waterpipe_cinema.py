#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 13:10:07 2021

@author: Mo Cohen

Pipeline for post-processing UM output data.
- Analyses model dynamics

Model: University of Exeter Stand Alone model of Proxima Centauri b, UM vn11.8

"""

import iris
import iris.plot as iplt
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import gif
import windspharm

brewer_red = mpl_cm.get_cmap('brewer_YlOrRd_09')
brewer_bg = mpl_cm.get_cmap('brewer_PuBu_09')
brewer_reds = mpl_cm.get_cmap('brewer_Reds_09')
brewer_redblu = mpl_cm.get_cmap('brewer_RdBu_11')



# def gif_zonal_wind(cubes):
    
#     for cube in cubes:
#         if cube.standard_name == 'x_wind':
#             x_wind = cube.copy()
            
#     dayside = x_wind.extract(iris.Constraint(longitude=lambda v: 270 < v <= 359 or 0 <= v <= 90, latitude=lambda v: -90 <= v <= 90))
#     nightside = x_wind.extract(iris.Constraint(longitude=lambda v: 90 < v <= 270, latitude=lambda v: -90 <= v <= 90))
    
#     dayside_zonal_mean = dayside.collapsed('longitude', iris.analysis.MEAN)
#     nightside_zonal_mean = nightside.collapsed('longitude', iris.analysis.MEAN)
    
#     dayside_frames = []
#     for i in range(dayside_zonal_mean.shape[0]):
#         CS_day = iplt.contourf(dayside_zonal_mean[i,:,:], brewer_redblu.N, cmap=brewer_redblu)
#         plt.title('Dayside Zonal Mean Zonal Wind [m s-1], month %s' %(i+1), y=1.05)
#         plt.ylabel('Height [m]')
#         plt.xlabel('Latitude [degrees]')
#         plt.clabel(CS_day, inline=False, colors='k', fmt='%1.1f')
#         plt.colorbar(pad=0.1)
        
#         dayside_frames.append(CS_day)
    
#     nightside_frames = []
#     for i in range(nightside_zonal_mean.shape[0]):
#         CS_night = iplt.contourf(nightside_zonal_mean[i,:,:], brewer_redblu.N, cmap=brewer_redblu)
#         plt.title('Nightside Zonal Mean Zonal Wind [m s-1], month %s' %(i+1), y=1.05)
#         plt.ylabel('Height [m]')
#         plt.xlabel('Latitude [degrees]')
#         plt.clabel(CS_day, inline=False, colors='k', fmt='%1.1f')
#         plt.colorbar(pad=0.1)    

#         nightside_frames.append(CS_night)
        
#     gif.save(dayside_frames, '/exports/csce/datastore/geos/users/s1144983/um_data/ch1_control/gifs/dayside_zonal_meanwind.gif', duration = 30, unit = 's', between='startend')
#     gif.save(nightside_frames, '/exports/csce/datastore/geos/users/s1144983/um_data/ch1_control/gifs/nightside_zonal_meanwind.gif', duration = 30, unit = 's', between='startend')


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
    
@gif.frame
def gif_dayside_wind(i, dayside_zonal_mean):
        
    CS_day = iplt.contourf(dayside_zonal_mean[i,:,:], brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Dayside Zonal Mean Zonal Wind [m s-1], month %s' %(i+1), y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.clabel(CS_day, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    
dayside_frames = []
for i in range(57):
    dayside_frame = gif_dayside_wind(i, dayside_zonal_mean)
    dayside_frames.append(dayside_frame)

gif.save(dayside_frames, '/exports/csce/datastore/geos/users/s1144983/um_data/ch1_5x/gifs/dayside_zonal_meanwind.gif', duration = 30, unit = 's', between='startend')

@gif.frame
def gif_nightside_wind(i, nightside_zonal_mean):  
    
    CS_night = iplt.contourf(nightside_zonal_mean[i,:,:], brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Nightside Zonal Mean Zonal Wind [m s-1], month %s' %(i+1), y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.clabel(CS_night, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)

nightside_frames = []
for i in range(57):
    nightside_frame = gif_nightside_wind(i, nightside_zonal_mean)
    nightside_frames.append(nightside_frame)
    
gif.save(nightside_frames, '/exports/csce/datastore/geos/users/s1144983/um_data/ch1_5x/gifs/nightside_zonal_meanwind.gif', duration = 30, unit = 's', between='startend')

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

gif.save(dayside_frames, '/exports/csce/datastore/geos/users/s1144983/um_data/ch1_5x/gifs/dayside_zonal_humidity.gif', duration = 900, unit = 's', between='startend')


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

gif.save(nightside_frames, '/exports/csce/datastore/geos/users/s1144983/um_data/ch1_control/gifs/nightside_zonal_humidity.gif', duration = 900, unit = 's', between='startend')

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
    plt.title('Dayside Zonal Mean Air Temp [m s-1], month %s' %(i+1), y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.clabel(CS_day, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    
dayside_frames = []
for i in range(57):
    dayside_frame = gif_dayside_temp(i, dayside_absolute_temp)
    dayside_frames.append(dayside_frame)

gif.save(dayside_frames, '/exports/csce/datastore/geos/users/s1144983/um_data/ch1_5x/gifs/dayside_zonal_airtemp.gif', duration = 30, unit = 's', between='startend')

@gif.frame
def gif_nightside_temp(i, nightside_absolute_temp):  
    
    CS_night = iplt.contourf(nightside_absolute_temp[i,:,:], brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Nightside Zonal Mean Air Temp [m s-1], month %s' %(i+1), y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.clabel(CS_night, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)

nightside_frames = []
for i in range(57):
    nightside_frame = gif_nightside_temp(i, nightside_absolute_temp)
    nightside_frames.append(nightside_frame)
    
gif.save(nightside_frames, '/exports/csce/datastore/geos/users/s1144983/um_data/ch1_5x/gifs/nightside_airtemp.gif', duration = 30, unit = 's', between='startend')

""" CODE BLOCK FOR STREAMFUNCTION """

for cube in control:
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
    iplt.contourf(streamfunction[time_slice,level,:,:]*1e-06, clevs, cmap=brewer_redblu, extend='both')
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.colorbar(orientation='horizontal')
    plt.title('Streamfunction [$10^6$ m2 s-1], h = %s km' %(level+1), y=1.20)

frames = []
for level in range(38):
    frame = calculate_streamfunction(x_wind, y_wind, level)
    frames.append(frame)

gif.save(frames, '/exports/csce/datastore/geos/users/s1144983/um_data/ch1_control/gifs/streamfunction.gif', duration = 19, unit = 's', between='startend')

    