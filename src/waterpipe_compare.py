# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 13:31:24 2020

@author: Mo Cohen

Pipeline for post-processing UM output data. 
- Compares model data to a control run

Model: University of Exeter Stand Alone model of Proxima Centauri b, UM vn11.8

"""
import os, iris
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import iris.plot as iplt
import iris.coords
from iris.coords import DimCoord
import numpy as np
import scipy as sp
import windspharm
# Import packages


brewer_red = mpl_cm.get_cmap('brewer_YlOrRd_09')
brewer_bg = mpl_cm.get_cmap('brewer_PuBu_09')
brewer_reds = mpl_cm.get_cmap('brewer_Reds_09')
brewer_redblu = mpl_cm.get_cmap('brewer_RdBu_11')
# Load colormaps to use in plots

def compare_standard(cubes1, cubes2, time_slice=-1):
    
    """ Find differences in matching cubes from two CubeLists 
        Compare standard 2D outputs of 30-day means from the UM:
        Mean surface temperature
        Mean precipitation flux
        Mean downward SW flux at surface                    """
        
    for cube in cubes1:
        if cube.standard_name == 'surface_temperature':
            if len(cube.cell_methods) != 0 and cube.cell_methods[0].method == 'mean':
                surface_temp_1 = cube
        if cube.standard_name == 'precipitation_flux':
            precip_flux_1 = cube
        if cube.standard_name == 'surface_net_downward_shortwave_flux':
            downward_sw_flux_1 = cube

    for cube in cubes2:
        if cube.standard_name == 'surface_temperature':
            if len(cube.cell_methods) != 0 and cube.cell_methods[0].method == 'mean':
                surface_temp_2 = cube
        if cube.standard_name == 'precipitation_flux':
            precip_flux_2 = cube
        if cube.standard_name == 'surface_net_downward_shortwave_flux':
            downward_sw_flux_2 = cube
           
    surface_temp_diff = surface_temp_1.data - surface_temp_2.data
    surface_temp = surface_temp_1
    surface_temp.data = surface_temp_diff
    precip_diff = precip_flux_1.data - precip_flux_2.data
    precip_flux = precip_flux_1
    precip_flux.data = precip_diff
    downward_sw_flux_diff = downward_sw_flux_1.data - downward_sw_flux_2.data
    downward_sw_flux = downward_sw_flux_1
    downward_sw_flux.data = downward_sw_flux_diff
    
    
    iplt.contourf(surface_temp[time_slice,:,:], brewer_red.N, cmap=brewer_red)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Difference in Mean Surface Temp [K]', y=1.20)
    plt.colorbar(pad=0.1)
    plt.show()
    
    iplt.contourf(precip_flux[time_slice,:,:]*86400, brewer_bg.N, cmap=brewer_bg)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Difference in Precipitation Flux [mm day-1]', y=1.20)
    plt.colorbar(pad=0.1)
    plt.show()
            
    iplt.contourf(downward_sw_flux[time_slice,:,:], brewer_bg.N, cmap=brewer_bg)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Difference in Downward SW Flux [W m-2]', y=1.20)
    plt.colorbar(pad=0.1)
    plt.show()
    

def compare_daily(cubes1, cubes2, time_slice=-1):
     
    """ Find differences in matching cubes from two CubeLists 
        Compare standard 2D daily outputs from the UM:
        Cloud area fraction                          """
    
    for cube in cubes1:
        if cube.long_name == 'cloud_area_fraction_assuming_maximum_random_overlap':
            cloud_frac_1 = cube
            
    for cube in cubes2:
        if cube.long_name == 'cloud_area_fraction_assuming_maximum_random_overlap':
            cloud_frac_2 = cube


    cloud_frac_diff = cloud_frac_1.data - cloud_frac_2.data
    cloud_frac = cloud_frac_1
    cloud_frac.data = cloud_frac_diff
    
    iplt.contourf(cloud_frac[time_slice,:,:], brewer_bg.N, cmap=brewer_bg)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Difference in Cloud Area Fraction', y=1.20)
    plt.colorbar(pad=0.1)
    plt.show()
    
    
def compare_profile(cubes1, cubes2, time_slice=-1, experiment='Experiment'):
    
    """ Compare temperature profile of control and experiment
        Compare specific humidity profile of control and experiment """
    
    for cube in cubes1:
        if cube.standard_name == 'air_potential_temperature':
            potential_temp_1 = cube
        if cube.standard_name == 'air_pressure':
            air_pressure_1 = cube
        if cube.standard_name == 'specific_humidity':
            spec_humidity_1 = cube 
            
    for cube in cubes2:
        if cube.standard_name == 'air_potential_temperature':
            potential_temp_2 = cube
        if cube.standard_name == 'air_pressure':
            air_pressure_2 = cube
        if cube.standard_name == 'specific_humidity':
            spec_humidity_2 = cube

            
    p0 = iris.coords.AuxCoord(100000.0, long_name='reference_pressure', units='Pa')
    p0.convert_units(air_pressure_1.units)
    absolute_temp_1 = potential_temp_1*((air_pressure_1/p0)**(287.05/1005))
    absolute_temp_2 = potential_temp_2*((air_pressure_2/p0)**(287.05/1005))

    
    plt.plot(absolute_temp_1[time_slice,:,0,0].data, np.arange(0,39), linestyle='-', color='b', label='Control, sub')
    plt.plot(absolute_temp_2[time_slice,:,0,0].data, np.arange(0,39), linestyle='-', color='r', label='%s, sub' %(experiment))
    plt.plot(absolute_temp_1[time_slice,:,0,72].data, np.arange(0,39), linestyle='--', color='b', label='Control, anti')
    plt.plot(absolute_temp_2[time_slice,:,0,72].data, np.arange(0,39), linestyle='--', color='r', label='%s, anti' %(experiment))
    plt.title('Temperature Profiles at Substellar and Antistellar Point')
    plt.xlabel('Temperature [K]')
    plt.ylabel('Height [km]')
    plt.legend()
    plt.show()
    
    plt.plot(spec_humidity_1[time_slice,:,0,0].data*1000, np.arange(0,39), linestyle='-', color='b', label='Control, sub')
    plt.plot(spec_humidity_2[time_slice,:,0,0].data*1000, np.arange(0,39), linestyle='-', color='r', label='%s, sub' %(experiment))
    plt.plot(spec_humidity_1[time_slice,:,0,72].data*1000, np.arange(0,39), linestyle='--', color='b', label='Control, anti')
    plt.plot(spec_humidity_2[time_slice,:,0,72].data*1000, np.arange(0,39), linestyle='--', color='r', label='%s, anti' %(experiment))
    plt.title('Humidity Profile at Substellar and Antistellar Point')
    plt.xlabel('Specific Humidity [g kg-1]')
    plt.ylabel('Height [km]')
    plt.legend()
    plt.show()
    
    

    
    
