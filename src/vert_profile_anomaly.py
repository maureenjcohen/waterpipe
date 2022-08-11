#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 19 16:30:25 2022

@author: Maureen Cohen
"""
import iris
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
from matplotlib.colors import TwoSlopeNorm
from iris.coord_systems import GeogCS
import numpy as np
from iris.analysis import calculus

# Import packages

redblu = mpl_cm.get_cmap('RdBu')
heat = mpl_cm.get_cmap('gist_heat')
blues = mpl_cm.get_cmap('Blues')

def vert_profile(cubes,start=500,end=600,level=8):

    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube[start:end,:,:,:].copy()

    heights = np.round(x_wind.coord('level_height').points*1e-03,0)
    time_axis = np.arange(0,x_wind.shape[0])
    lats = x_wind.coord('latitude')
    lons = x_wind.coord('longitude')

    if lats.bounds == None:
        x_wind.coord('latitude').guess_bounds()
    if lons.bounds == None:
        x_wind.coord('longitude').guess_bounds()

    grid_areas = iris.analysis.cartography.area_weights(x_wind)
    global_mean = x_wind.collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=grid_areas)

    fig1, ax1 = plt.subplots()
    ax1.set_xlabel('Time [days]')
    ax1.set_ylabel('Wind [m/s]')
    ax1.set_title('Global mean zonal wind at h=%s km' %heights[level])
    plt.plot(global_mean[:,level].data)
    plt.show()

    fig2, ax2 = plt.subplots()
    ax2.set_xlabel('Time [days]')
    ax2.set_ylabel('Height [km]')
    ax2.set_title('Global mean zonal wind')
    plt.contourf(time_axis, heights[:30], global_mean[:,:30].data.T, redblu.N, cmap=redblu, norm=TwoSlopeNorm(0))
    plt.colorbar(pad=0.1)
    plt.show()

def temp_profile(cubes,start=0,end=100,level=8, select='absolute'):

    for cube in cubes:
        if cube.standard_name == 'air_potential_temperature':
            theta = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'air_pressure':
            pressure = cube[start:end,:,:,:].copy()

    p0 = iris.coords.AuxCoord(
        100000.0, long_name='reference_pressure', units='Pa')
    p0.convert_units(pressure.units)
    # R and cp in J/kgK for 300K
    temperature = theta*((pressure/p0)**(287.05/1005))

    if select == 'absolute':
        temp = temperature.copy()
        titleterm = 'absolute temperature'
    elif select == 'potential':
        temp = theta.copy()
        titleterm = 'potential temperature'
    elif select == 'pressure':
        temp = pressure.copy()
        titleterm = 'air pressure'

    heights = np.round(temp.coord('level_height').points*1e-03,0)
    time_axis = np.arange(0,temp.shape[0])
    lats = temp.coord('latitude')
    lons = temp.coord('longitude')

    if lats.bounds == None:
        temp.coord('latitude').guess_bounds()
    if lons.bounds == None:
        temp.coord('longitude').guess_bounds()

    temp = temp.extract(iris.Constraint(longitude=lambda v: 270 < v <= 359 or 0 <= v <= 90, latitude=lambda v: -90 <= v <= 90))
    grid_areas = iris.analysis.cartography.area_weights(temp)
    global_mean = temp.collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=grid_areas)

    fig1, ax1 = plt.subplots()
    ax1.set_xlabel('Time [days]')
    ax1.set_ylabel('Temperature [K]')
    ax1.set_title('Dayside mean %s at h=%s km' %(titleterm, heights[level]))
    plt.plot(global_mean[:,level].data)
    plt.show()

    fig2, ax2 = plt.subplots()
    ax2.set_xlabel('Time [days]')
    ax2.set_ylabel('Height [km]')
    ax2.set_title('Dayside mean %s' %titleterm)
    plt.contourf(time_axis, heights[:30], global_mean[:,:30].data.T, heat.N, cmap=heat)
    plt.colorbar(pad=0.1)
    plt.show()

def w_profile(cubes,start=500,end=600,level=8):

    for cube in cubes:
        if cube.standard_name == 'upward_air_velocity':
            z_wind = cube[start:end,:,:,:].copy()

    heights = np.round(z_wind.coord('level_height').points*1e-03,0)
    time_axis = np.arange(0,z_wind.shape[0])
    lats = z_wind.coord('latitude')
    lons = z_wind.coord('longitude')

    if lats.bounds == None:
        z_wind.coord('latitude').guess_bounds()
    if lons.bounds == None:
        z_wind.coord('longitude').guess_bounds()

    z_wind = z_wind.extract(iris.Constraint(longitude=lambda v: 270 < v <= 359 or 0 <= v <= 90, latitude=lambda v: -90 <= v <= 90))
    grid_areas = iris.analysis.cartography.area_weights(z_wind)
    global_mean = z_wind.collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=grid_areas)

    fig1, ax1 = plt.subplots()
    ax1.set_xlabel('Time [days]')
    ax1.set_ylabel('Wind [m/s]')
    ax1.set_title('Dayside mean vertical wind at h=%s km' %heights[level])
    plt.plot(global_mean[:,level].data)
    plt.show()

    fig2, ax2 = plt.subplots()
    ax2.set_xlabel('Time [days]')
    ax2.set_ylabel('Height [km]')
    ax2.set_title('Dayside mean vertical wind')
    plt.contourf(time_axis, heights[:30], global_mean[:,:30].data.T, redblu.N, cmap=redblu, norm=TwoSlopeNorm(0))
    plt.colorbar(pad=0.1)
    plt.show()

def wind_shear(cubes,start=500,end=600,time_slice=-1,level=8):

    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube[start:end,:,:,:].copy()

    heights = np.round(x_wind.coord('level_height').points*1e-03,0)
    shear = iris.analysis.calculus.differentiate(x_wind, 'level_height')
    latitudes = x_wind.coord('latitude').points
    longitudes = x_wind.coord('longitude').points

    time_axis = np.arange(0,x_wind.shape[0])
    lats = shear.coord('latitude')
    lons = shear.coord('longitude')

    if lats.bounds == None:
        shear.coord('latitude').guess_bounds()
    if lons.bounds == None:
        shear.coord('longitude').guess_bounds()

    grid_areas = iris.analysis.cartography.area_weights(shear)
    global_mean = shear.collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=grid_areas)

    plt.contourf(np.roll(longitudes,72), latitudes, np.roll(shear[time_slice,level,:,:].data,72,axis=1), redblu.N, cmap=redblu)
    plt.xlabel('Longitude [degrees]')
    plt.ylabel('Latitude [degrees]')
    plt.title('Vertical shear of zonal wind at %s km' %heights[level])
    cbar = plt.colorbar(pad=0.1)
    plt.show()

    fig, ax = plt.subplots()
    ax.set_xlabel('Time [days]')
    ax.set_ylabel('Height [km]')
    ax.set_title('Global mean vertical wind shear')
    plt.contourf(time_axis, heights[:30], global_mean[:,:30].data.T, redblu.N, cmap=redblu, norm=TwoSlopeNorm(0))
    plt.colorbar(pad=0.1)
    plt.show()

def cloud_curve(cubes,start=500,end=600,level=8):

    for cube in cubes:
        if cube.standard_name == 'mass_fraction_of_cloud_ice_in_air':
            ice = cube[start:end,:,:,:].copy()
        if cube.standard_name == 'mass_fraction_of_cloud_liquid_water_in_air':
            liq = cube[start:end,:,:,:].copy()

    total_cloud = ice + liq
    heights = np.round(total_cloud.coord('level_height').points*1e-03,0)
    time_axis = np.arange(0,total_cloud.shape[0])
    lats = total_cloud.coord('latitude')
    lons = total_cloud.coord('longitude')

    if lats.bounds == None:
        total_cloud.coord('latitude').guess_bounds()
    if lons.bounds == None:
        total_cloud.coord('longitude').guess_bounds()

    total_cloud = total_cloud.extract(iris.Constraint(longitude=lambda v: 270 < v <= 359 or 0 <= v <= 90, latitude=lambda v: -90 <= v <= 90))
    grid_areas = iris.analysis.cartography.area_weights(total_cloud)
    global_mean = total_cloud.collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=grid_areas)

    fig1, ax1 = plt.subplots()
    ax1.set_xlabel('Time [days]')
    ax1.set_ylabel('Wind [m/s]')
    plt.title('Dayside mean total cloud at h=%s km' %heights[level])
    plt.plot(global_mean[:,level].data)
    plt.show()

    fig2, ax2 = plt.subplots()
    ax2.set_xlabel('Time [days]')
    ax2.set_ylabel('Height [km]')
    ax2.set_title('Dayside mean total cloud')
    plt.contourf(time_axis, heights[:15], global_mean[:,:15].data.T, blues.N, cmap=blues)
    plt.colorbar(pad=0.1)
    plt.show()
