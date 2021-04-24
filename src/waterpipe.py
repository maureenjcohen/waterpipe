# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 14:54:07 2020

@author: Mo Cohen

Pipeline for post-processing UM output data.
- Analyses model data and outputs plots and calculated data

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
# Import packages


brewer_red = mpl_cm.get_cmap('brewer_YlOrRd_09')
brewer_bg = mpl_cm.get_cmap('brewer_PuBu_09')
brewer_reds = mpl_cm.get_cmap('brewer_Reds_09')
brewer_redblu = mpl_cm.get_cmap('brewer_RdBu_11')
# Load colormaps to use in plots

# Block 1: Data management


def load_files(directory):
    
    """ Read PP files from a directory and load them into an Iris CubeList"""
    
    filenames = os.listdir(directory)
    # Create list of names of files in source directory

    allfiles = []
    for entry in filenames:
        fullpath = os.path.join(directory, entry)
        allfiles.append(fullpath)
    # Add full address to file name and append to a new list
    cubes= iris.load(allfiles)
    # Load files into an Iris CubeList
    
    return cubes


def load_loop(directory):
    
    """ Load model data dumps into an Iris cube individually
        Concatenate individual cubes into a CubeList           """
        
    filenames = os.listdir(directory)
    print(filenames)
    
    allfiles = []
    for entry in filenames:
        fullpath = os.path.join(directory, entry)
        allfiles.append(fullpath)
    print(allfiles)
    
    dictionary = {}
    for entry in allfiles:
        key = entry[37:53]
#        key = entry[69:83]
        print(key)
        group = dictionary.get(key, [])
        group.append(entry)
        dictionary[key] = group
     
    list_of_cubes = []    
    for key in dictionary:
        model_dump = []
        for entry in allfiles:
            if entry[37:53] == key:
                model_dump.append(entry)
        cube = iris.load(model_dump)
        list_of_cubes.append(cube)
        print(str('cubes: ') + str(len(list_of_cubes)))
        
    cubes = list_of_cubes
    
    return cubes
                


def convert_lazy_data(cubes):
    
    """ Convert cube data from lazy to real"""
    
    for cube in cubes:
        
        if cube.has_lazy_data() == True:
            cube.data
        else:
            pass
     # Convert lazy data to real data
     
    return cubes


def convert_names(cubes):
    
    """ Iris doesn't seem to recognise some STASH indices. This function adds
    an intelligible long_name attribute to cubes identified only by index."""
    
    for cube in cubes:
        if cube.name() == 'm01s03i026':
            cube.long_name = 'roughness_length_after_boundary_layer'
        if cube.name() == 'm01s03i027':
            cube.long_name = 'effective_roughness_length_for_scalars'
        if cube.name() == 'm01s03i028':
            cube.long_name = 'roughness_length_for_momentum_without_orography'
        if cube.name() == 'm01s03i464':
            cube.long_name = 'obukhov_length'
        if cube.name() == 'm01s03i465':
            cube.long_name = 'explicit_friction_velocity'
        if cube.name() == 'm01s09i202':
            cube.long_name = 'very_low_cloud_amount'
        if cube.name() == 'm01s30i407':
            cube.long_name = 'water_vapour_flux_x'
        if cube.name() == 'm01s30i408':
            cube.long_name = 'water_vapour_flux_y'
        if cube.name() == 'm01s30i409':
            cube.long_name = 'water_vapour_flux_z'
        if cube.name() == 'm01s00i253':
            cube.long_name = 'density_r_r'
        if cube.name() == 'm01s09i181':
            cube.long_name = 'change_over_time_in_air_temperature_due_to_boundary_layer_ls_cloud'
            
    return cubes
            

def load_and_save(directory, filename):
    
    cubes = load_files(directory)
    print(cubes)
    cubes = convert_names(cubes)
    print(cubes)
    cubes = convert_lazy_data(cubes)
    
    iris.save(cubes, directory+filename)
    print('Completed save of ' + filename)
    
def convert_and_save(cubes,directory,filename):
    
    cubes = convert_lazy_data(cubes)
    
    iris.save(cubes, directory+filename)
    print('Completed save of ' + filename)
           

def check_difference(cubes1, cubes2):
    
 """ Check whether two Iris CubeLists are identical """
 
 for i in range(0,len(cubes1)):
    if np.sum(cubes1[i].data - cubes2[i].data) > 0.0:
        print(cubes1[i].name() + 'is different')
    else:
        print('No difference in ' + cubes1[i].name())
 # Check if there is a difference between two cubelists

        

# Block 2: Plotting functions

def plot_temperature(cubes, time_slice=-1):
    
    """ Plot mean surface temperature from final output data dump
        Plot development of average daily surface temperature over time:
            dayside, nightside, and global                             """ 
        
    for cube in cubes:
        if cube.standard_name == 'surface_temperature':
            surface_temp = cube.copy()
            
            if len(surface_temp.cell_methods) != 0 and surface_temp.cell_methods[0].method == 'mean':
            # Extract mean surface temperature from CubeList
                
                iplt.contourf(surface_temp[time_slice,:,:], brewer_red.N, cmap=brewer_red)
                ax = plt.gca()
                ax.gridlines(draw_labels=True)
                plt.title('Mean Surface Temp [K]', y=1.20)
                plt.colorbar(pad=0.1)
                plt.show()
            # Plot mean surface temperature of final data dump
            
            elif len(surface_temp.cell_methods) == 0:
            # Extract daily surface temperature from CubeList
            
                lats = surface_temp.coord('latitude')
                longs = surface_temp.coord('longitude')

                if lats.bounds == None:
                    surface_temp.coord('latitude').guess_bounds()
                if longs.bounds == None:
                    surface_temp.coord('longitude').guess_bounds()
                
                dayside = surface_temp.extract(iris.Constraint(longitude=lambda v: 270 < v <= 359 or 0 <= v <= 90, latitude=lambda v: -90 <= v <= 90))
                day_grid = iris.analysis.cartography.area_weights(dayside)
                nightside = surface_temp.extract(iris.Constraint(longitude=lambda v: 90 < v <= 270, latitude=lambda v: -90 <= v <= 90))
                night_grid = iris.analysis.cartography.area_weights(nightside)
                global_grid = iris.analysis.cartography.area_weights(surface_temp)
                
                dayside_temp = dayside.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=day_grid)
                
                nightside_temp = nightside.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=night_grid)
                
                global_temp = surface_temp.collapsed(['latitude','longitude'],iris.analysis.MEAN, weights=global_grid)
                
                run_length = surface_temp.shape[0]
                days = DimCoord((np.arange(0,run_length)),standard_name='time', units='days')
                
                iplt.plot(days,dayside_temp)
                plt.title('Average Dayside Temperature')
                plt.ylabel('Temperature [K]')
                plt.xlabel('Time [days]')
                plt.show()
                
                iplt.plot(days,nightside_temp)
                plt.title('Average Nightside Temperature')
                plt.ylabel('Temperature [K]')
                plt.xlabel('Time [days]')
                plt.show()
                
                iplt.plot(days,global_temp)
                plt.title('Average Global Temperature')
                plt.ylabel('Temperature [K]')
                plt.xlabel('Time [days]')
                plt.show()
                    

def plot_temp_profile(cubes, time_slice=-1):
    
    """ Plot temperature profile with height at substellar point
        Plot temperature profile with height at antistellar point
        Plot air temperature at 1.5 m for all longs/lats

        Returns: absolute air temperature as Iris cube       """
        
    for cube in cubes:
        if cube.standard_name == 'air_potential_temperature':
            potential_temp = cube.copy()
        if cube.standard_name == 'air_pressure':
            air_pressure = cube.copy()
        if cube.standard_name == 'air_temperature':
            if len(cube.cell_methods) != 0 and cube.cell_methods[0].method == 'mean':
                air_temp_bl = cube.copy()
            
    p0 = iris.coords.AuxCoord(100000.0, long_name='reference_pressure', units='Pa')
    p0.convert_units(air_pressure.units)
    absolute_temp = potential_temp*((air_pressure/p0)**(287.05/1005)) # R and cp in J/kgK for 300K
    
    plt.plot(absolute_temp[time_slice,:,45,0].data, np.arange(0,absolute_temp.shape[1]))
    plt.title('Temperature Profile at Substellar Point')
    plt.xlabel('Temperature [K]')
    plt.ylabel('Height [km]')
    plt.show()
    
    plt.plot(absolute_temp[time_slice,:,45,72].data, np.arange(0,absolute_temp.shape[1]))
    plt.title('Temperature Profile at Antistellar Point')
    plt.xlabel('Temperature [K]')
    plt.ylabel('Height [km]')
    plt.show()
    
    plt.plot(absolute_temp[time_slice,:,45,0].data, np.arange(0,absolute_temp.shape[1]), color='r', label='Substellar')
    plt.plot(absolute_temp[time_slice,:,45,72].data, np.arange(0,absolute_temp.shape[1]), color='b', label='Antistellar')
    plt.title('Temperature Profiles at Substellar and Antistellar Point')
    plt.xlabel('Temperature [K]')
    plt.ylabel('Height [km]')
    plt.legend()
    plt.show()
    
    iplt.contourf(air_temp_bl[time_slice,:,:], brewer_red.N, cmap=brewer_red)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Air Temp at 1.5 m [K]', y=1.20)
    plt.xlabel('Latitude [degrees]')
    plt.ylabel('Longitude [degrees]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    return absolute_temp


def plot_xtemp(cubes, time_slice=-1):    
    
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
    
    dayside_zonal_mean = dayside.collapsed('longitude', iris.analysis.MEAN)
    nightside_zonal_mean = nightside.collapsed('longitude', iris.analysis.MEAN)
    
    CS_day = iplt.contourf(dayside_zonal_mean[time_slice,:,:], brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Dayside Zonal Mean Air Temp [K]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.clabel(CS_day, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    plt.show()
    
    CS_night = iplt.contourf(nightside_zonal_mean[time_slice,:,:], brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Nightside Zonal Mean Air Temp [K]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.clabel(CS_night, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    plt.show()
    
    dayside_meridional_mean = dayside.collapsed('latitude', iris.analysis.MEAN)
    nightside_meridional_mean = nightside.collapsed('latitude', iris.analysis.MEAN)
    
    CS_day = iplt.contourf(dayside_meridional_mean[time_slice,:,:], brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Dayside Meridional Mean Air Temp [K]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Longitude [degrees]')
    plt.clabel(CS_day, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    plt.show()
    
    CS_night = iplt.contourf(nightside_meridional_mean[time_slice,:,:], brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Nightside Meridional Mean Air Temp [K]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Longitude [degrees]')
    plt.clabel(CS_night, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    plt.show()
    
    return absolute_temp


def plot_stratosphere(cube, level=50):
    
    """ Plot average stratospheric temperature over time
        Input: absolute temperature cube, level of tropopause       """

    if cube.shape[1] < level:
        raise Exception('Cube height is less than stratosphere height. Get a taller atmosphere.')
    
    absolute_temp = cube.copy()
    stratosphere = absolute_temp[:,level:-1,:,:]
    
    lats = stratosphere.coord('latitude')
    longs = stratosphere.coord('longitude')

    if lats.bounds == None:
        stratosphere.coord('latitude').guess_bounds()
    if longs.bounds == None:
        stratosphere.coord('longitude').guess_bounds()

    grid_areas = iris.analysis.cartography.area_weights(stratosphere)  
    stratospheric_temp = stratosphere.collapsed(['latitude', 'longitude', 'model_level_number'],iris.analysis.MEAN, weights=grid_areas)
    
    run_length = stratosphere.shape[0]
    months = DimCoord((np.arange(0,run_length)), standard_name='time', units='months')    
    
    iplt.plot(months, stratospheric_temp)
    plt.title('Average Stratospheric Temperature [K]')
    plt.xlabel('Time [months]')
    plt.ylabel('Temperature [K]')
    plt.show()
                
        
def plot_radiation(cubes, time_slice=-1):
    
    """ Plot radiation balance at top of atmosphere over time
        Plot surface net shortwave flux from final output data dump """
        
    for cube in cubes:
        if cube.standard_name == 'toa_incoming_shortwave_flux':
            incoming_rad = cube.copy()
        if cube.standard_name == 'toa_outgoing_longwave_flux_assuming_clear_sky':
            outgoing_lw_clear = cube.copy()
        if cube.standard_name == 'toa_outgoing_shortwave_flux_assuming_clear_sky':
            outgoing_sw_clear = cube.copy()
        if cube.standard_name == 'toa_outgoing_longwave_flux':
            outgoing_lw = cube.copy()
        if cube.standard_name == 'toa_outgoing_shortwave_flux':
            outgoing_sw = cube.copy()
        if cube.standard_name == 'surface_net_downward_shortwave_flux':
            surface_sw_rad = cube.copy()
    # Extract desired cubes from CubeList
        
    run_length = incoming_rad.shape[0]
    
    lats = incoming_rad.coord('latitude')
    longs = incoming_rad.coord('longitude')

    if lats.bounds == None:
        incoming_rad.coord('latitude').guess_bounds()
    if longs.bounds == None:
        incoming_rad.coord('longitude').guess_bounds()
    
    grid_areas = iris.analysis.cartography.area_weights(incoming_rad)
    radiation_difference = incoming_rad.data - outgoing_lw.data - outgoing_sw.data
    radiation_difference_watts = radiation_difference*grid_areas
    
    summed = []
    # incoming_list = []
    # outgoing_lw_list = []
    # outgoing_sw_list = []
#   outgoing_lw_clear_list = []
#   outgoing_sw_clear_list = []

    for i in range(0,run_length):
        total = np.sum(radiation_difference_watts[i,:,:])
        summed.append(total)
        # total_incoming = np.sum(incoming_rad[i,:,:].data)
        # incoming_list.append(total_incoming)
        # total_outgoing_lw = np.sum(outgoing_lw[i,:,:].data)
        # outgoing_lw_list.append(total_outgoing_lw)
        # total_outgoing_sw = np.sum(outgoing_sw[i,:,:].data)
        # outgoing_sw_list.append(total_outgoing_sw)
        
    incoming_mean = incoming_rad.collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=grid_areas)
    outgoing_lwmean = outgoing_lw.collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=grid_areas)
    outgoing_swmean = outgoing_sw.collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=grid_areas)
    radiation_difference_means = incoming_mean.data - outgoing_lwmean.data - outgoing_swmean.data 
            
    plt.plot(np.arange(0,run_length), np.array(summed))
    plt.xlabel('Time [months]')
    plt.ylabel('TOA radiation balance [W]')
    plt.title('Radiation balance (summed Watts)')
    plt.show()
    
    plt.plot(np.arange(0,run_length), radiation_difference_means)
    plt.xlabel('Time [months]')
    plt.ylabel('TOA radiation balance [W m-2]')
    plt.title('Radiation balance (meaned Watts)')
    plt.show()    
    
    # plt.plot(np.arange(0,run_length), incoming_list, linestyle='-', color='k', label='Incoming SW')
    # plt.plot(np.arange(0,run_length), outgoing_lw_list, linestyle='--', color='b', label='Outgoing LW')
    # plt.plot(np.arange(0,run_length), outgoing_sw_list, linestyle='--', color='r', label='Outgoing SW')
    # # plt.plot(np.arange(0,run_length), outgoing_lw_clear_list, linestyle=':', color='b', label='Outgoing LW, clear')
    # # plt.plot(np.arange(0,run_length), outgoing_sw_clear_list, linestyle=':', color='r', label='Outgoing SW, clear')
    # plt.title('Radiation Balance Elements')
    # plt.xlabel('Time [months]')
    # plt.ylabel('Radiation [W]')
    # plt.legend()
    # plt.show()    
      
    iplt.contourf(outgoing_lw[time_slice,:,:], brewer_red.N, cmap=brewer_red)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('TOA Outgoing Longwave Flux [W m-2]', y=1.20)
    plt.ylabel('Longitude [degrees]')
    plt.xlabel('Latitude [degrees]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    iplt.contourf(outgoing_sw[time_slice,:,:], brewer_red.N, cmap=brewer_red)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('TOA Outgoing Shortwave Flux [W m-2]', y=1.20)
    plt.ylabel('Longitude [degrees]')
    plt.xlabel('Latitude [degrees]')
    plt.colorbar(pad=0.1)
    plt.show()
  
    iplt.contourf(surface_sw_rad[time_slice,:,:], brewer_red.N, cmap=brewer_red)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Surface Net Downward Shortwave Flux [W m-2]', y=1.20)
    plt.ylabel('Longitude [degrees]')
    plt.xlabel('Latitude [degrees]')
    plt.colorbar(pad=0.1)
    plt.show()
    

def plot_outgoing(cubes):
    
    for cube in cubes:
        if cube.standard_name == 'toa_outgoing_longwave_flux':
            outgoing_lw = cube.copy()

    
    global_lw = outgoing_lw[24:,:,:].collapsed('time',iris.analysis.MEAN)
    
    iplt.contourf(global_lw, brewer_red.N, cmap=brewer_red)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Mean TOA Outgoing Longwave Flux [W m-2]', y=1.20)
    plt.ylabel('Longitude [degrees]')
    plt.xlabel('Latitude [degrees]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    data = outgoing_lw[24:,:,:].data
    spatial_stdev = np.std(data, axis=0)
    
    plt.figure(figsize=(10,4))
    plt.contourf(np.roll(spatial_stdev, 72, axis=1), cmap=brewer_red)
    plt.title('Standard Deviation in TOA Outgoing LW [W m-2]')
    plt.ylabel('Longitude [degrees]')
    plt.xlabel('Latitude [degrees]')
    plt.yticks((0,45,89),('90S', '0', '90N'))
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.colorbar(pad=0.1)
    plt.show()
    
    plt.figure(figsize=(10,4))
    plt.contourf(spatial_stdev, cmap=brewer_red)
    plt.title('Standard Deviation in TOA Outgoing LW [W m-2]')
    plt.ylabel('Longitude [degrees]')
    plt.xlabel('Latitude [degrees]')
    plt.yticks((0,45,89),('90S', '0', '90N'))
    plt.xticks((0,36, 72, 104, 143),('0', '90E', '180E/W', '90W','0'))
    plt.colorbar(pad=0.1)
    plt.show()


def plot_humidity(cubes, level=14, time_slice=-1):
    
    """ Plot dayside atmospheric humidity over time
        Plot nightside atmospheric humidity over time
        Plot specific humidity profile with height at substellar point
        Plot specific humidity distribution at 15 km (~ tropopause) 
        Plot averaged specific humidity at terminators
        Plot specific humidity at longitude 0                        """
    
    for cube in cubes:
        if cube.standard_name == 'specific_humidity':
            spec_humidity = cube.copy()         
    
    lats = spec_humidity.coord('latitude')
    longs = spec_humidity.coord('longitude')
    heights = np.round(spec_humidity.coord('level_height').points,0)


    if lats.bounds == None:
        spec_humidity.coord('latitude').guess_bounds()
    if longs.bounds == None:
        spec_humidity.coord('longitude').guess_bounds()
        
    dayside = spec_humidity.extract(iris.Constraint(longitude=lambda v: 270 < v <= 359 or 0 <= v <= 90, latitude=lambda v: -90 <= v <= 90))
    day_grid = iris.analysis.cartography.area_weights(dayside)
    dayside_humidity = dayside.collapsed(['latitude', 'longitude', 'model_level_number'],iris.analysis.MEAN, weights=day_grid)

    nightside = spec_humidity.extract(iris.Constraint(longitude=lambda v: 90 < v <= 270, latitude=lambda v: -90 <= v <= 90))
    night_grid = iris.analysis.cartography.area_weights(nightside)
    nightside_humidity = nightside.collapsed(['latitude', 'longitude', 'model_level_number'],iris.analysis.MEAN, weights=night_grid)
    
    global_grid = iris.analysis.cartography.area_weights(spec_humidity)
    global_humidity = spec_humidity.collapsed(['latitude','longitude', 'model_level_number'],iris.analysis.MEAN, weights=global_grid)    
    
    run_length = spec_humidity.shape[0]
    days = DimCoord((np.arange(0,run_length)),standard_name='time', units='days')
        
    iplt.plot(days, dayside_humidity)
    plt.xlabel('Time [days]')
    plt.ylabel('Specific Humidity [kg kg-1]')
    plt.title('Average Dayside Specific Humidity')
    plt.show()

    iplt.plot(days, nightside_humidity)
    plt.xlabel('Time [days]')
    plt.ylabel('Specific Humidity [kg kg-1]')
    plt.title('Average Nightside Specific Humidity')
    plt.show() 
        
    iplt.plot(days, global_humidity)
    plt.xlabel('Time [days]')
    plt.ylabel('Specific Humidity [kg kg-1]')
    plt.title('Average Global Specific Humidity')
    plt.show() 
    
    plt.plot(spec_humidity[time_slice,:,45,0].data, np.arange(0,39))
    plt.title('Humidity Profile at Substellar Point')
    plt.xlabel('Specific Humidity [kg kg-1]')
    plt.ylabel('Height [km]')
    plt.show()
    
    iplt.contourf(spec_humidity[time_slice,level,:,:], brewer_bg.N, cmap=brewer_bg)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Specific Humidity at h=%s km [kg kg-1]'%(heights[level]), y=1.20)
    plt.colorbar(pad=0.1)
    plt.show()
    
    iplt.contourf((spec_humidity[time_slice,:,:,36]+spec_humidity[time_slice,:,:,108])/2, brewer_bg.N, cmap=brewer_bg)
    plt.title('Specific Humidity at Terminators, final day [kg kg-1]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude')
    plt.colorbar(pad=0.1)
    plt.show()
    
    iplt.contourf(spec_humidity[time_slice,:,:,0], brewer_bg.N, cmap=brewer_bg)
    plt.title('Specific Humidity at Longitude 0 [kg kg-1]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude')
    plt.colorbar(pad=0.1)
    plt.show()


def plot_xhumidity(cubes, time_slice=-1):
    
    for cube in cubes:
        if cube.standard_name == 'specific_humidity':
            spec_humidity = cube.copy()         

            
    dayside = spec_humidity.extract(iris.Constraint(longitude=lambda v: 270 < v <= 359 or 0 <= v <= 90, latitude=lambda v: -90 <= v <= 90))
    nightside = spec_humidity.extract(iris.Constraint(longitude=lambda v: 90 < v <= 270, latitude=lambda v: -90 <= v <= 90))
    
    dayside_zonal_mean = dayside.collapsed('longitude', iris.analysis.MEAN)
    nightside_zonal_mean = nightside.collapsed('longitude', iris.analysis.MEAN)
    
    CS_day = iplt.contourf(dayside_zonal_mean[time_slice,:,:], brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Dayside Zonal Mean Specific Humidity [kg kg-1]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
#    plt.clabel(CS_day, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    plt.show()
    
    CS_night = iplt.contourf(nightside_zonal_mean[time_slice,:,:], brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Nightside Zonal Mean Specific Humidity [kg kg-1]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
#    plt.clabel(CS_night, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    plt.show()

    dayside_meridional_mean = dayside.collapsed('latitude', iris.analysis.MEAN)
    nightside_meridional_mean = nightside.collapsed('latitude', iris.analysis.MEAN)
    
    CS_day = iplt.contourf(dayside_meridional_mean[time_slice,:,:], brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Dayside Meridional Mean Specific Humidity [kg kg-1]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Longitude [degrees]')
#    plt.clabel(CS_day, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    plt.show()
    
    CS_night = iplt.contourf(nightside_meridional_mean[time_slice,:,:], brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Nightside Meridional Mean Specific Humidity [kg kg-1]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Longitude [degrees]')
#    plt.clabel(CS_night, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    plt.show()


def plot_precipitation(cubes, time_slice=-1):
    
    """ Plot dayside precipitation rate over time
        Plot nightside precipitation rate over time
        Plot mean precipitation of final output data dump  """
        
    for cube in cubes:
        if cube.standard_name == 'precipitation_flux':
            precipitation = cube.copy()
        
    lats = precipitation.coord('latitude')
    longs = precipitation.coord('longitude')

    if lats.bounds == None:
        precipitation.coord('latitude').guess_bounds()
    if longs.bounds == None:
        precipitation.coord('longitude').guess_bounds()
    
    dayside = precipitation.extract(iris.Constraint(longitude=lambda v: 270 < v <= 359 or 0 <= v <= 90, latitude=lambda v: -90 <= v <= 90))
    day_grid = iris.analysis.cartography.area_weights(dayside)
    dayside_precip = dayside.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=day_grid)

    nightside = precipitation.extract(iris.Constraint(longitude=lambda v: 90 < v <= 270, latitude=lambda v: -90 <= v <= 90))
    night_grid = iris.analysis.cartography.area_weights(nightside)
    nightside_precip = nightside.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=night_grid)
    
    global_grid = iris.analysis.cartography.area_weights(precipitation)
    global_precip = precipitation.collapsed(['latitude','longitude'],iris.analysis.MEAN, weights=global_grid)
    
    substellar = precipitation.extract(iris.Constraint(longitude=lambda v: 345 < v <= 359 or 0 <= v <= 15, latitude=lambda v: -15 <= v <= 15))
    substellar_grid = iris.analysis.cartography.area_weights(substellar)
    substellar_precip = substellar.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=substellar_grid)
    
    run_length = precipitation.shape[0]
    months = DimCoord((np.arange(0,run_length)), standard_name='time', units='months')
    
    iplt.plot(months, dayside_precip*86400)
    plt.xlabel('Time [months]')
    plt.ylabel('Precipitation Flux [mm day-1]')
    plt.title('Average Dayside Precipitation Flux')
    plt.show()

    iplt.plot(months, nightside_precip*86400)
    plt.xlabel('Time [months]')
    plt.ylabel('Precipitation Flux [mm day-1]')
    plt.title('Average Nightside Precipitation Flux')
    plt.show()
    
    iplt.plot(months, global_precip*86400)
    plt.xlabel('Time [months]')
    plt.ylabel('Precipitation Flux [mm day-1]')
    plt.title('Average Global Precipitation Flux')
    plt.show()
    
    iplt.plot(months, substellar_precip*86400)
    plt.xlabel('Time [months]')
    plt.ylabel('Precipitation Flux [mm day-1]')
    plt.title('Average Substellar Region Precipitation Flux')
    plt.show()
    
    iplt.contourf(precipitation[time_slice,:,:]*86400, brewer_bg.N, cmap=brewer_bg)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Mean Precipitation Flux [mm day-1]', y=1.20)
    plt.xlabel('Latitude [degrees]')
    plt.ylabel('Longitude [degrees]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    precipitation = precipitation*86400
    precipitation.units = 'mm.day-1'
    
    return precipitation
    

def plot_evaporation(cubes, time_slice=-1):
    
    """ Plot surface evaporation 
        Calculated from surface latent heat flux 
        
        Returns: evaporation as Iris cube       """
    
    for cube in cubes:
        if cube.standard_name == 'surface_upward_latent_heat_flux':
            latent_heat = cube.copy()
    
    density = 999.7 #kg/m^3, density of water at 283.15 K
    l_v = 2477300 #J/kg, latent heat of vaporisation of water at 283.15 K
    evaporation = (latent_heat/(density*l_v))*86400*1000 # Calculate evap and convert to mm/day
    
    iplt.contourf(evaporation[time_slice,:,:], brewer_red.N, cmap=brewer_red)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Evaporation [mm day-1]', y=1.20)
    plt.ylabel('Longitude [degrees]')
    plt.xlabel('Latitude [degrees]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    evaporation.units = 'mm.day-1'
    
    return evaporation

def plot_water_balance(evap, precip, time_slice=-1):
    
    """ Plot difference between evaporation and precipitation
        Inputs are Iris cubes calculated by plot_evaporation and plot_precipitation """
        
    if evap.shape != precip.shape:
        raise Exception('Cubes must be same shape')
    if evap.units != precip.units:
        raise Exception('Cubes must have same units')

    evaporation = evap.copy()
    precipitation = precip.copy()
    
    evaporation = evaporation.regrid(precipitation, iris.analysis.Linear())    
    
    difference = evaporation - precipitation
    
    lats = precipitation.coord('latitude')
    longs = precipitation.coord('longitude')

    if lats.bounds.all() == None:
        precipitation.coord('latitude').guess_bounds()
    if longs.bounds.all() == None:
        precipitation.coord('longitude').guess_bounds()
    
    global_grid = iris.analysis.cartography.area_weights(precipitation)

    iplt.contourf(difference[time_slice,:,:], brewer_bg.N, cmap=brewer_bg)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Evaporation Minus Precipitation [mm day-1]', y=1.20)
    plt.ylabel('Longitude [degrees]')
    plt.xlabel('Latitude [degrees]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    run_length = precipitation.shape[0]
    water_balance = []
    for i in range(0,run_length):
        water_difference = evaporation[i,:,:] - precipitation[i,:,:]
        total = np.sum(water_difference.data*global_grid)
        water_balance.append(total)
        
    plt.plot(np.arange(0,run_length), water_balance)
    plt.xlabel('Time [months]')
    plt.ylabel('Evaporation Minus Precipitation [mm day-1]')
    plt.title('Evaporation Minus Precipitation')
    plt.show()  
                
    evap_mean = evap.collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=global_grid)
    precip_mean = precip.collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=global_grid)
    mean_difference = evap_mean.data - precip_mean.data
    
    plt.plot(np.arange(0,run_length), mean_difference)
    plt.xlabel('Time [months]')
    plt.ylabel('Evaporation Minus Precipitation [mm day-1 m-2]')
    plt.title('Evaporation Minus Precipitation')
    plt.show() 


def plot_clouds(cubes, time_slice=-1, periodicity=False):
    
    """ Plot cloud area fraction for all longs/lats on final output data dump
        Plot cloud volume fraction at longitude 0, final output
        Plot cloud volume fraction at terminators (averaged), final output 
        Plot percentage of latitudes with clear skies at terminators over time 
        If periodicity=True, plot periodicity of clear sky days      
        
        Returns: numpy array of percentage of latitudes with clear skies at terminators """

        
    for cube in cubes:
        if cube.long_name == 'cloud_area_fraction_assuming_maximum_random_overlap':
            cloud_cover = cube.copy()
        if cube.long_name == 'cloud_volume_fraction_in_atmosphere_layer':
            cloud_volume = cube.copy()
        if cube.long_name == 'liquid_cloud_volume_fraction_in_atmosphere_layer':
            liquid_cloud = cube.copy()
        if cube.long_name == 'ice_cloud_volume_fraction_in_atmosphere_layer':
            ice_cloud = cube.copy()
        if cube.standard_name == 'mass_fraction_of_cloud_ice_in_air':
            ice_condensate = cube.copy()
        if cube.standard_name == 'mass_fraction_of_cloud_liquid_water_in_air':
            liquid_condensate = cube.copy()
    
    iplt.contourf(cloud_cover[time_slice,:,:], brewer_bg.N, cmap=brewer_bg)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Cloud Area Fraction', y=1.20)
    plt.colorbar(pad=0.1)
    plt.show()

    iplt.contourf(cloud_volume[time_slice,:,:,0], brewer_bg.N, cmap=brewer_bg)
    plt.title('Cloud Volume Fraction at Longitude 0', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    iplt.contourf((cloud_volume[time_slice,:,:,36]+cloud_volume[time_slice,:,:,108])/2, brewer_bg.N, cmap=brewer_bg)
    plt.title('Cloud Volume Fraction at Terminators [kg kg-1]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    clear_list = []
    for day in range(0,cloud_cover.shape[0]):
        data = (cloud_cover[day,:,36].data + cloud_cover[day,:,108].data)/2
        clear = data[np.where( data <= 0.2)]
        count = len(clear.tolist())
        clear_list.append(count)
    
    average = round(np.mean(np.array(clear_list)), 3)
    std_dev = round(np.std(np.array(clear_list)), 3)
    
    time_axis = np.arange(0,cloud_cover.shape[0])
    plt.plot(time_axis,(np.array(clear_list)/90)*100)
    plt.title('Percent of latitudes with clear skies (cloud frac <= 0.2)')
    plt.xlabel('Days')
    plt.ylabel('Percent')
    plt.text(x=0,y=0, s="Mean = %s, std = %s"%(average, std_dev))
    plt.show()
    # Get an overview of what percentage of days might allow observation
    # of atmospheric spectrum in transmission spectroscopy
    
    if periodicity == True:
        
        daily_cloud = np.array(clear_list)
        run_length = cloud_cover.shape[0]
        cloud_fft = sp.fftpack.fft(daily_cloud)
        cloud_psd = np.abs(cloud_fft)**2
        cloudfreq = sp.fftpack.fftfreq(len(cloud_psd), 1./run_length)
        i = cloudfreq > 0
        
        fig, ax = plt.subplots(1,1,figsize=(8,4))
        ax.plot(cloudfreq[i], cloud_psd[i])
        ax.set_xlim(0,15)
        ax.set_xlabel('Frequency [1/run_length]')
        ax.set_ylabel('PSD')
        ax.set_title('Periodicity of clear sky')
        # Check if there is periodicity in cloud cover over time
        # Uses scipy discrete Fast Fourier Transform. Units of 1/run_length as
        # there is no meaningful time period in a tidally-locked planet without
        # eccentricity or obliquity. The plot tells you the amplitude of a cycle
        # and how many times it repeated during the model run time.
    
    return np.array(clear_list)

def quickview_clouds(cubes, time_slice=-1):
    
    for cube in cubes:
        if cube.long_name == 'cloud_area_fraction_assuming_maximum_random_overlap':
            cloud_cover = cube.copy()
        if cube.long_name == 'cloud_volume_fraction_in_atmosphere_layer':
            cloud_volume = cube.copy()
        if cube.long_name == 'liquid_cloud_volume_fraction_in_atmosphere_layer':
            liquid_cloud = cube.copy()
        if cube.long_name == 'ice_cloud_volume_fraction_in_atmosphere_layer':
            ice_cloud = cube.copy()
        if cube.standard_name == 'mass_fraction_of_cloud_ice_in_air':
            ice_condensate = cube.copy()
        if cube.standard_name == 'mass_fraction_of_cloud_liquid_water_in_air':
            liquid_condensate = cube.copy()
    
    iplt.contourf(cloud_cover[time_slice,:,:], brewer_bg.N, cmap=brewer_bg)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Cloud Area Fraction', y=1.20)
    plt.colorbar(pad=0.1)
    plt.show()

    iplt.contourf(cloud_volume[time_slice,:,:,0], brewer_bg.N, cmap=brewer_bg)
    plt.title('Cloud Volume Fraction at Longitude 0', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    iplt.contourf(liquid_cloud[time_slice,:,:,0], brewer_bg.N, cmap=brewer_bg)
    plt.title('Liquid Cloud Volume Fraction at Longitude 0', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    iplt.contourf(ice_cloud[time_slice,:,:,0], brewer_bg.N, cmap=brewer_bg)
    plt.title('Ice Cloud Volume Fraction at Longitude 0', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    iplt.contourf(ice_condensate[time_slice,:,:,0], brewer_bg.N, cmap=brewer_bg)
    plt.title('Ice Cloud Mass Fraction at Longitude 0 [kg kg-1]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    iplt.contourf(liquid_condensate[time_slice,:,:,0], brewer_bg.N, cmap=brewer_bg)
    plt.title('Liquid Cloud Mass Fraction at Longitude 0 [kg kg-1]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    iplt.contourf((cloud_volume[time_slice,:,:,36]+cloud_volume[time_slice,:,:,108])/2, brewer_bg.N, cmap=brewer_bg)
    plt.title('Cloud Volume Fraction at Terminators', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    iplt.contourf((liquid_cloud[time_slice,:,:,36]+liquid_cloud[time_slice,:,:,108])/2, brewer_bg.N, cmap=brewer_bg)
    plt.title('Liquid Cloud Volume Fraction at Terminators', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    iplt.contourf((ice_cloud[time_slice,:,:,36]+ice_cloud[time_slice,:,:,108])/2, brewer_bg.N, cmap=brewer_bg)
    plt.title('Ice Cloud Volume Fraction at Terminators', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    iplt.contourf((ice_condensate[time_slice,:,:,36]+ice_condensate[time_slice,:,:,108])/2, brewer_bg.N, cmap=brewer_bg)
    plt.title('Ice Cloud Mass Fraction at Terminators [kg kg-1]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    iplt.contourf((liquid_condensate[time_slice,:,:,36]+liquid_condensate[time_slice,:,:,108])/2, brewer_bg.N, cmap=brewer_bg)
    plt.title('Liquid Cloud Mass Fraction at Terminators [kg kg-1]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    
    
def check_periodicity(data):
    
    """ Check for periodicity in a 1-D numpy array 
        Plot Fast Fourier Transform of data        """
    
    run_length = data.shape[0]
    fft = sp.fftpack.fft(data)
    psd = np.abs(fft)**2
    fftfreq = sp.fftpack.fftfreq(len(psd), 1./run_length)
    i = fftfreq > 0
    
    fig, ax = plt.subplots(1,1,figsize=(8,4))
    ax.plot(fftfreq[i], psd[i])
    ax.set_xlim(0,15)
    ax.set_xlabel('Frequency [1/run_length]')
    ax.set_ylabel('PSD')
    ax.set_title('Periodicity in data')


def plot_vapour(cubes, time_slice=-1):
    
    """ Plot zonal vapour flux
        Plot meridional vapour flux
        Plot upward vapour flux
        Plot streamplot of 2-D vapour fluxes """
    
    for cube in cubes:
        if cube.long_name == 'water_vapour_flux_x':
            x_vapour = cube.copy()
        if cube.long_name == 'water_vapour_flux_y':
            y_vapour = cube.copy()
        if cube.long_name == 'water_vapour_flux_z':
            z_vapour = cube.copy()
    
    iplt.contourf(x_vapour[time_slice,:,:], brewer_bg.N, cmap=brewer_bg)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Zonal Vapour Flux [kg m-2 s-1]', y=1.20)
    plt.colorbar(pad=0.1)
    plt.show()
    
    iplt.contourf(y_vapour[time_slice,:,:], brewer_bg.N, cmap=brewer_bg)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Meridional Vapour Flux [kg m-2 s-1]', y=1.20)
    plt.colorbar(pad=0.1)
    plt.show()
    
    iplt.contourf(z_vapour[time_slice,:,:], brewer_bg.N, cmap=brewer_bg)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Upward Vapour Flux [kg m-2 s-1]', y=1.20)
    plt.colorbar(pad=0.1)
    plt.show()
    
    magnitude = iris.analysis.maths.apply_ufunc(np.sqrt, (x_vapour**2 + y_vapour**2))
    
    X,Y = np.meshgrid(np.arange(-72,72), np.arange(-45,45))
    fig = plt.figure(figsize = (12, 7)) 
    strm = plt.streamplot(X, Y, np.roll(x_vapour[time_slice,:,:].data, 72, axis=1), np.roll(y_vapour[time_slice,:,:].data, 72, axis=1), density = 0.5, color=np.roll(magnitude[time_slice,:,:].data, 72, axis=1), cmap=brewer_bg)
    # Since .data method extracts the numpy array and strips the metadata, the longitude/latitude information is lost.
    # To align plot so that (0,0) is at the center as in the Iris plots, use numpy.roll to shift columns (longitude) 180 degrees (72 places)
    fig.colorbar(strm.lines)
    plt.title('Vapour flux magnitude and direction [kg m-2 s-1]')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.yticks((-45,-30,-15,0,15,30,45),('90S','60S','30S','0','30N','60N','90N'))    
    plt.show() 
    
    magnitude3D = iris.analysis.maths.apply_ufunc(np.sqrt, (x_vapour**2 + y_vapour**2 + z_vapour**2))
    
    lats = x_vapour.coord('latitude')
    longs = x_vapour.coord('longitude')

    if lats.bounds == None:
        x_vapour.coord('latitude').guess_bounds()
    if longs.bounds == None:
        x_vapour.coord('longitude').guess_bounds()
    
    run_length = x_vapour.shape[0]
    global_grid = iris.analysis.cartography.area_weights(x_vapour)
    global_avg = magnitude3D.collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=global_grid)
    months = DimCoord((np.arange(0,run_length)), standard_name='time', units='months')

    
    # global_avgs = []
    # for day in range(0,x_vapour.shape[0]):
    #     global_avg = np.mean(magnitude3D[day,:,:].data*global_grid)
    #     global_avgs.append(global_avg)

    iplt.plot(months, global_avg)
    plt.xlabel('Time [months]')
    plt.ylabel('Average Total Vapour Flux [kg s-1]')
    plt.title('Vapour Flux Over Time')
    plt.show() 


def plot_efficiency(cubes):    
    
    """ Plot day-night heat redistribution efficiency, defined in Leconte et al. 2013
        as the ratio of mean nightside outgoing longwave radiation over mean dayside 
        outgoing longwave radiation                                                  """
    
    for cube in cubes:
        if cube.standard_name == 'toa_outgoing_longwave_flux':
            outgoing_lw = cube.copy()
            
    lats = outgoing_lw.coord('latitude')
    longs = outgoing_lw.coord('longitude')

    if lats.bounds == None:
        outgoing_lw.coord('latitude').guess_bounds()
    if longs.bounds == None:
        outgoing_lw.coord('longitude').guess_bounds()
    
    dayside = outgoing_lw.extract(iris.Constraint(longitude=lambda v: 270 < v <= 359 or 0 <= v <= 90, latitude=lambda v: -90 <= v <= 90))
    day_grid = iris.analysis.cartography.area_weights(dayside)
    dayside_avg = dayside.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=day_grid)

    nightside = outgoing_lw.extract(iris.Constraint(longitude=lambda v: 90 < v <= 270, latitude=lambda v: -90 <= v <= 90))
    night_grid = iris.analysis.cartography.area_weights(nightside)
    nightside_avg = nightside.collapsed(['latitude', 'longitude'],iris.analysis.MEAN, weights=night_grid)    

    run_length = outgoing_lw.shape[0]
    months = DimCoord((np.arange(0,run_length)),standard_name='time', units='months')
    
    iplt.plot(months,(nightside_avg/dayside_avg))
    plt.title('Day-Night Heat Redistribution Efficiency')
    plt.ylabel('Ratio')
    plt.xlabel('Time [months]')
    plt.show()
    

def plot_cloud_variability(cubes):
    
    """ Plot spatial standard deviation in cloud area fraction"""
    
    for cube in cubes:
        if cube.long_name == 'cloud_area_fraction_assuming_maximum_random_overlap':
            cloud_cover = cube.copy()
    
    data = cloud_cover[1000:,:,:].data
    spatial_stdev = np.std(data, axis=0)
    
    plt.figure(figsize=(10,4))
    plt.contourf(np.roll(spatial_stdev, 72, axis=1), cmap=brewer_bg)
    plt.title('Standard Deviation in Cloud Cover [Area Fraction]')
    plt.ylabel('Latitude [degrees]')
    plt.xlabel('Longitude [degrees]')
    plt.yticks((0,45,89),('90S', '0', '90N'))
    plt.xticks((0,36, 72, 104, 143),('180W', '90W', '0', '90E', '180E'))
    plt.colorbar(pad=0.1)
    plt.show()  


def plot_surftemp_variability(cubes):
    
    """ Plot spatial standard deviation in mean surface temperature"""
    
    for cube in cubes:
        if cube.standard_name == 'surface_temperature':
            surface_temp = cube.copy()
            
            if len(surface_temp.cell_methods) != 0 and surface_temp.cell_methods[0].method == 'mean':
    
                data = surface_temp[24:,:,:].data
                spatial_stdev = np.std(data, axis=0)
            
                plt.figure(figsize=(10,4))
                plt.contourf(np.roll(spatial_stdev, 72, axis=1), cmap=brewer_reds)
                plt.title('Standard Deviation in Mean Surface Temperature [K]')
                plt.ylabel('Latitude [degrees]')
                plt.xlabel('Longitude [degrees]')
                plt.yticks((0,45,89),('90S', '0', '90N'))
                plt.xticks((0,36, 72, 104, 143),('180W', '90W', '0', '90E', '180E'))
                plt.colorbar(pad=0.1)
                plt.show()
                

def plot_airtemp_variability(cubes, level=35):
    
    """ Plot spatial standard deviation in mean air temperature at levels"""
    
    for cube in cubes:
        if cube.standard_name == 'air_potential_temperature':
            potential_temp = cube.copy()
        if cube.standard_name == 'air_pressure':
            air_pressure = cube.copy()
            
    heights = np.round(air_pressure.coord('level_height').points,0)
    p0 = iris.coords.AuxCoord(100000.0, long_name='reference_pressure', units='Pa')
    p0.convert_units(air_pressure.units)
    absolute_temp = potential_temp*((air_pressure/p0)**(287.05/1005)) # R and cp in J/kgK for 300K
    
    data = absolute_temp[24:,:,:,:].data
    spatial_stdev = np.std(data, axis=0)
    
    plt.figure(figsize=(10,4))
    plt.contourf(np.roll(spatial_stdev[level[0],:,:], 72, axis=1), cmap=brewer_reds)
    plt.title('Standard Deviation in Air Temperature [K], height=%s km' %(heights[level]))
    plt.ylabel('Longitude [degrees]')
    plt.xlabel('Latitude [degrees]')
    plt.yticks((0,45,89),('90S', '0', '90N'))
    plt.xticks((0,36, 72, 104, 143),('180W', '90W', '0', '90E', '180E'))
    plt.colorbar(pad=0.1)
    plt.show()
    
    plt.figure(figsize=(10,4))
    plt.contourf(np.roll(spatial_stdev[level[1],:,:], 72, axis=1), cmap=brewer_reds)
    plt.title('Standard Deviation in Air Temperature [K], height=%s km' %(heights[level]))
    plt.ylabel('Longitude [degrees]')
    plt.xlabel('Latitude [degrees]')
    plt.yticks((0,45,89),('90S', '0', '90N'))
    plt.xticks((0,36, 72, 104, 143),('180W', '90W', '0', '90E', '180E'))
    plt.colorbar(pad=0.1)
    plt.show()
    
    plt.figure(figsize=(10,4))
    plt.contourf(np.roll(spatial_stdev[level[2],:,:], 72, axis=1), cmap=brewer_reds)
    plt.title('Standard Deviation in Air Temperature [K], height=%s km' %(heights[level]))
    plt.ylabel('Longitude [degrees]')
    plt.xlabel('Latitude [degrees]')
    plt.yticks((0,45,89),('90S', '0', '90N'))
    plt.xticks((0,36, 72, 104, 143),('180W', '90W', '0', '90E', '180E'))
    plt.colorbar(pad=0.1)
    plt.show()
    

def plot_humidity_variability(cubes, level=35):
    
    """ Plot spatial standard deviation in specific humidity at levels"""
    
    for cube in cubes:
        if cube.standard_name == 'specific_humidity':
            spec_humid = cube.copy()
    
    heights = np.round(spec_humid.coord('level_height').points,0)

    spec_humid = spec_humid[500:,:,:,:]
    flattened = spec_humid.collapsed('model_level_number', iris.analysis.MEAN)
    data = spec_humid.data
    flattened_data = flattened.data
    spatial_stdev = np.std(data, axis=0)
    flattened_stdev = np.std(flattened_data, axis=0)
    
    plt.figure(figsize=(10,4))
    plt.contourf(np.roll(flattened_stdev[:,:], 72, axis=1), cmap=brewer_reds)
    plt.title('Standard Deviation in Specific Humidity Column [kg kg-1]')
    plt.ylabel('Longitude [degrees]')
    plt.xlabel('Latitude [degrees]')
    plt.yticks((0,45,89),('90S', '0', '90N'))
    plt.xticks((0,36, 72, 104, 143),('180W', '90W', '0', '90E', '180E'))
    plt.colorbar(pad=0.1)
    plt.show()    
    
    plt.figure(figsize=(10,4))
    plt.contourf(np.roll(spatial_stdev[level[0],:,:], 72, axis=1), cmap=brewer_reds)
    plt.title('Standard Deviation in Specific Humidity [kg kg-1], height=%s km' %(heights[level]))
    plt.ylabel('Longitude [degrees]')
    plt.xlabel('Latitude [degrees]')
    plt.yticks((0,45,89),('90S', '0', '90N'))
    plt.xticks((0,36, 72, 104, 143),('180W', '90W', '0', '90E', '180E'))
    plt.colorbar(pad=0.1)
    plt.show()
    
    plt.figure(figsize=(10,4))
    plt.contourf(np.roll(spatial_stdev[level[1],:,:], 72, axis=1), cmap=brewer_reds)
    plt.title('Standard Deviation in Specific Humidity [kg kg-1], height=%s km' %(heights[level]))
    plt.ylabel('Longitude [degrees]')
    plt.xlabel('Latitude [degrees]')
    plt.yticks((0,45,89),('90S', '0', '90N'))
    plt.xticks((0,36, 72, 104, 143),('180W', '90W', '0', '90E', '180E'))
    plt.colorbar(pad=0.1)
    plt.show()
    
    plt.figure(figsize=(10,4))
    plt.contourf(np.roll(spatial_stdev[level[2],:,:], 72, axis=1), cmap=brewer_reds)
    plt.title('Standard Deviation in Specific Humidity [kg kg-1], height=%s km' %(heights[level]))
    plt.ylabel('Longitude [degrees]')
    plt.xlabel('Latitude [degrees]')
    plt.yticks((0,45,89),('90S', '0', '90N'))
    plt.xticks((0,36, 72, 104, 143),('180W', '90W', '0', '90E', '180E'))
    plt.colorbar(pad=0.1)
    plt.show()
    

def zonal_heat_transport(cubes, time_slice=-1):
    
    """ Total eastward heat transport at longitude theta """
            
    for cube in cubes:
        if cube.standard_name == 'toa_incoming_shortwave_flux':
            incoming_rad = cube.copy()
        if cube.standard_name == 'toa_outgoing_longwave_flux':
            outgoing_lw = cube.copy()
        if cube.standard_name == 'toa_outgoing_shortwave_flux':
            outgoing_sw = cube.copy()
        if cube.standard_name == 'surface_net_downward_shortwave_flux':
            surface_sw_rad = cube.copy()   
            
    
    lats = incoming_rad.coord('latitude')
    longs = incoming_rad.coord('longitude')

    if lats.bounds == None:
        incoming_rad.coord('latitude').guess_bounds()
    if longs.bounds == None:
        incoming_rad.coord('longitude').guess_bounds()
    
    grid_areas = iris.analysis.cartography.area_weights(incoming_rad)
    print(grid_areas.shape)
    radiation_difference = incoming_rad[time_slice,:,:].data - outgoing_lw[time_slice,:,:].data - outgoing_sw[time_slice,:,:].data
    radiation_difference_watts = radiation_difference*grid_areas[time_slice,:,:]
    print(radiation_difference_watts.shape)
            
    plt.figure(figsize=(10,4))
    plt.contourf(np.roll(radiation_difference_watts[:,:], 72, axis=1), cmap=brewer_redblu)
    plt.title('Radiation imbalance [W]')
    plt.ylabel('Longitude [degrees]')
    plt.xlabel('Latitude [degrees]')
    plt.yticks((0,45,89),('90S', '0', '90N'))
    plt.xticks((0,36, 72, 104, 143),('180W', '90W', '0', '90E', '180E'))
    plt.colorbar(pad=0.1)
    plt.show()
    

def plot_geopotential(cubes, g=10.9, R=7160000):
    
    for cube in cubes:
        if cube.standard_name =='air_pressure':
            pressure = cube.copy()
            
    heights = pressure.coord('level_height').points
    geopotential = g*R*heights/(R+heights)
    

# if __name__ == '__main__':
#     print('Source folder:')
#     directory = input()

