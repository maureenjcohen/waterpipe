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
import windspharm
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
#        key = entry[63:76]
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
            if len(cube.cell_methods) != 0 and cube.cell_methods[0].method == 'mean':
            # Extract mean surface temperature from CubeList
                
                iplt.contourf(cube[time_slice,:,:], brewer_red.N, cmap=brewer_red)
                ax = plt.gca()
                ax.gridlines(draw_labels=True)
                plt.title('Mean Surface Temp [K]', y=1.20)
                plt.colorbar(pad=0.1)
                plt.show()
            # Plot mean surface temperature of final data dump
            
            elif len(cube.cell_methods) == 0:
            # Extract daily surface temperature from CubeList
    
                dayside = cube.extract(iris.Constraint(longitude=lambda v: 270 < v < 360 or 0 < v < 90, latitude=lambda v: -90 < v < 90))
                nightside = cube.extract(iris.Constraint(longitude=lambda v: 90 < v < 270, latitude=lambda v: -90 < v < 90))
                
                dayside_temp = dayside.collapsed('latitude',iris.analysis.MEAN)
                dayside_temp = dayside_temp.collapsed('longitude',iris.analysis.MEAN)
                
                nightside_temp = nightside.collapsed('latitude',iris.analysis.MEAN)
                nightside_temp = nightside_temp.collapsed('longitude',iris.analysis.MEAN)
                
                global_temp = cube.collapsed('longitude',iris.analysis.MEAN)
                global_temp = global_temp.collapsed('latitude',iris.analysis.MEAN)
                
                run_length = cube.shape[0]
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
            potential_temp = cube
        if cube.standard_name == 'air_pressure':
            air_pressure = cube
        if cube.standard_name == 'air_temperature':
            if len(cube.cell_methods) != 0 and cube.cell_methods[0].method == 'mean':
                air_temp_bl = cube
            
    p0 = iris.coords.AuxCoord(100000.0, long_name='reference_pressure', units='Pa')
    p0.convert_units(air_pressure.units)
    absolute_temp = potential_temp*((air_pressure/p0)**(287.05/1005))
    
    plt.plot(absolute_temp[time_slice,:,45,0].data, np.arange(0,39))
    plt.title('Temperature Profile at Substellar Point')
    plt.xlabel('Temperature [K]')
    plt.ylabel('Height [km]')
    plt.show()
    
    plt.plot(absolute_temp[time_slice,:,45,72].data, np.arange(0,39))
    plt.title('Temperature Profile at Antistellar Point')
    plt.xlabel('Temperature [K]')
    plt.ylabel('Height [km]')
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
                
        
def plot_radiation(cubes, time_slice=-1):
    
    """ Plot radiation balance at top of atmosphere over time
        Plot surface net shortwave flux from final output data dump """
        
    for cube in cubes:
        if cube.standard_name == 'toa_incoming_shortwave_flux':
            incoming_rad = cube
        if cube.standard_name == 'toa_outgoing_longwave_flux':
            outgoing_lw = cube
        if cube.standard_name == 'toa_outgoing_shortwave_flux':
            outgoing_sw = cube
        if cube.standard_name == 'surface_net_downward_shortwave_flux':
            surface_sw_rad = cube
    # Extract desired cubes from CubeList
    
    run_length = incoming_rad.shape[0]
    gridbox_area = 300000*222390
    radiation_balance = []
    for i in range(0,run_length-1):
        radiation_difference = incoming_rad[i,:,:] - outgoing_lw[i,:,:] - outgoing_sw[i,:,:]
        total = np.sum(radiation_difference.data)*gridbox_area
        radiation_balance.append(total)
        
    plt.plot(np.arange(0,run_length-1), radiation_balance)
    plt.xlabel('Time [months]')
    plt.ylabel('TOA radiation balance [W]')
    plt.title('Radiation balance')
    plt.show()    
      
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


def plot_humidity(cubes, level=14, time_slice=-1):
    
    """ Plot dayside atmospheric humidity over time
        Plot nightside atmospheric humidity over time
        Plot specific humidity profile with height at substellar point
        Plot specific humidity distribution at 15 km (~ tropopause) 
        Plot averaged specific humidity at terminators
        Plot specific humidity at longitude 0                        """
    
    for cube in cubes:
        if cube.standard_name == 'specific_humidity':
            spec_humidity = cube          
    
    # dayside = spec_humidity.extract(iris.Constraint(longitude=lambda v: 270 < v < 360 or 0 < v < 90, latitude=lambda v: -90 < v < 90))
    # nightside = spec_humidity.extract(iris.Constraint(longitude=lambda v: 90 < v < 270, latitude=lambda v: -90 < v < 90))
    
    # dayside_avg = dayside.collapsed('model_level_number', iris.analysis.MEAN)
    # dayside_avg = dayside_avg.collapsed('longitude', iris.analysis.MEAN)
    # dayside_avg = dayside_avg.collapsed('latitude', iris.analysis.MEAN)
    
    # nightside_avg = nightside.collapsed('model_level_number', iris.analysis.MEAN)
    # nightside_avg = nightside_avg.collapsed('longitude', iris.analysis.MEAN)
    # nightside_avg = nightside_avg.collapsed('latitude', iris.analysis.MEAN)
    
    # run_length = spec_humidity.shape[0]
    # days = DimCoord((np.arange(0,run_length)),standard_name='time', units='days')
    
    # iplt.plot(days,dayside_avg)
    # plt.title('Average Dayside Specific Humidity [kg kg-1]')
    # plt.ylabel('Humidity [kg kg-1]')
    # plt.xlabel('Time [days]')
    # plt.show()
    
    # iplt.plot(days,nightside_avg)
    # plt.title('Average Nightside Specific Humidity [kg kg-1]')
    # plt.ylabel('Humidity [kg kg-1]')
    # plt.xlabel('Time [days]')
    # plt.show()
    
    
    # dayside_avgs = []
    # nightside_avgs = []
    # for day in range(0,dayside.shape[0]):
    #     dayside_avg = np.mean(dayside[day,:,:,:].data)
    #     dayside_avgs.append(dayside_avg)
    #     nightside_avg = np.mean(nightside[day,:,:,:].data)
    #     nightside_avgs.append(nightside_avg)
    
        
    # plt.plot(np.arange(0,spec_humidity.shape[0]),np.array(dayside_avgs))
    # plt.xlabel('Time [days]')
    # plt.ylabel('Specific Humidity [kg kg-1]')
    # plt.title('Average Dayside Specific Humidity')
    # plt.show()

    # plt.plot(np.arange(0,spec_humidity.shape[0]),np.array(nightside_avgs))
    # plt.xlabel('Time [days]')
    # plt.ylabel('Specific Humidity [kg kg-1]')
    # plt.title('Average Nightside Specific Humidity')
    # plt.show() 
    
    plt.plot(spec_humidity[time_slice,:,45,0].data, np.arange(0,39))
    plt.title('Humidity Profile at Substellar Point')
    plt.xlabel('Specific Humidity [kg kg-1]')
    plt.ylabel('Height [km]')
    plt.show()
    
    iplt.contourf(spec_humidity[time_slice,14,:,:], brewer_bg.N, cmap=brewer_bg)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Specific Humidity at h=%s km [kg kg-1]'%(level+1), y=1.20)
    plt.colorbar(pad=0.1)
    plt.show()
    
    iplt.contourf((spec_humidity[time_slice,:,:,36]+spec_humidity[-1,:,:,108])/2, brewer_bg.N, cmap=brewer_bg)
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


def plot_precipitation(cubes, time_slice=-1):
    
    """ Plot dayside precipitation rate over time
        Plot nightside precipitation rate over time
        Plot mean precipitation of final output data dump  """
        
    for cube in cubes:
        if cube.standard_name == 'precipitation_flux':
            precipitation = cube
        
    dayside_precip = precipitation.extract(iris.Constraint(longitude=lambda v: 270 < v < 360 or 0 < v < 90, latitude=lambda v: -90 < v < 90))
    nightside_precip = precipitation.extract(iris.Constraint(longitude=lambda v: 90 < v < 270, latitude=lambda v: -90 < v < 90))
    
    dayside_precips = []
    nightside_precips = []
    for day in range(0,dayside_precip.shape[0]):
        dayside_flux = np.mean(dayside_precip[day,:,:].data)
        dayside_precips.append(dayside_flux)
        nightside_flux = np.mean(nightside_precip[day,:,:].data)
        nightside_precips.append(nightside_flux)
        
    plt.plot(np.arange(0,precipitation.shape[0]),np.array(dayside_precips)*86400)
    plt.xlabel('Time [months]')
    plt.ylabel('Precipitation Flux [mm day-1]')
    plt.title('Average Dayside Precipitation Flux')
    plt.show()

    plt.plot(np.arange(0,precipitation.shape[0]),np.array(nightside_precips)*86400)
    plt.xlabel('Time [months]')
    plt.ylabel('Precipitation Flux [mm day-1]')
    plt.title('Average Nightside Precipitation Flux')
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
            latent_heat = cube
    
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

def plot_water_balance(evaporation, precipitation, time_slice=-1):
    
    """ Plot difference between evaporation and precipitation
        Inputs are Iris cubes calculated by plot_evaporation and plot_precipitation """
        
    if evaporation.shape != precipitation.shape:
        raise Exception('Cubes must be same shape')
    if evaporation.units != precipitation.units:
        raise Exception('Cubes must have same units')
        
    difference = evaporation - precipitation
    
    iplt.contourf(difference[time_slice,:,:], brewer_bg.N, cmap=brewer_bg)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Evaporation Minus Precipitation [mm day-1]', y=1.20)
    plt.ylabel('Longitude [degrees]')
    plt.xlabel('Latitude [degrees]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    run_length = evaporation.shape[0]
    water_balance = []
    for i in range(0,run_length-1):
        water_difference = evaporation[i,:,:] - precipitation[i,:,:]
        total = np.sum(water_difference.data)
        water_balance.append(total)
        
    plt.plot(np.arange(0,run_length-1), water_balance)
    plt.xlabel('Time [months]')
    plt.ylabel('Evaporation Minus Precipitation [mm day-1]')
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
            cloud_cover = cube
        if cube.long_name == 'cloud_volume_fraction_in_atmosphere_layer':
            cloud_volume = cube
        if cube.long_name == 'liquid_cloud_volume_fraction_in_atmosphere_layer':
            liquid_cloud = cube
        if cube.long_name == 'ice_cloud_volume_fraction_in_atmosphere_layer':
            ice_cloud = cube
        if cube.standard_name == 'mass_fraction_of_cloud_ice_in_air':
            ice_condensate = cube
        if cube.standard_name == 'mass_fraction_of_cloud_liquid_water_in_air':
            liquid_condensate = cube
    
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
    
    iplt.contourf((cloud_volume[time_slice,:,:,36]+cloud_volume[-1,:,:,108])/2, brewer_bg.N, cmap=brewer_bg)
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
    

def plot_winds(cubes, level=14, time_slice=-1, wpharm=False, omega=0.64617667):
    
    """ Inputs: Iris CubeList and model level you want to look at
        wpharm=True allows windspharm library to be used on Linux machines
        omega is the rotation rate of the planet in 10^-5 rad/s. Default is Proxima Centauri b.
        
        Plot zonal wind speed at level
        Plot streamplot at level
        Plot zonal mean zonal winds for planet
        
        If wpharm=True, plot stream function and relative vorticity    """
        
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube
        if cube.standard_name == 'y_wind':
            y_wind = cube
        if cube.standard_name == 'upward_air_velocity':
            z_wind = cube
 
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
   
    speed = iris.analysis.maths.apply_ufunc(np.sqrt, (x_wind**2 + y_wind**2))
    
    iplt.contourf(x_wind[time_slice,level,:,:], brewer_reds.N, cmap=brewer_reds)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Zonal Wind Speed, h=%s km [m s-1]' %(level+1), y=1.20)
    plt.colorbar(pad=0.1)
    plt.show()

    X,Y = np.meshgrid(np.arange(-72,72), np.arange(-45,45))
    fig = plt.figure(figsize = (12, 7)) 
    strm = plt.streamplot(X, Y, np.roll(x_wind[time_slice,level,:,:].data, 72, axis=1), np.roll(y_wind[time_slice,level,:,:].data, 72, axis=1), density = 0.5, color=np.roll(speed[time_slice,level,:,:].data, 72, axis=1), cmap=brewer_reds)
    # Since .data method extracts the numpy array and strips the metadata, the longitude/latitude information is lost.
    # To align plot so that (0,0) is at the center as in the Iris plots, use numpy.roll to shift columns (longitude) 180 degrees (72 places)
    fig.colorbar(strm.lines)
    plt.title('Wind speed and direction [m s-1], h=%s km' %(level+1))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.xticks((-72,-52,-32,-12,0,12,32,52,72),('180W','140W','100W','60W','0','60E','100E','140E','180E'))
    plt.yticks((-45,-30,-15,0,15,30,45),('90S','60S','30S','0','30N','60N','90N'))    
    plt.show() 
        
    zonal_mean = x_wind.collapsed('longitude', iris.analysis.MEAN)
    
    CS = iplt.contourf(zonal_mean[time_slice,:,:], brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Zonal Mean Zonal Wind [m s-1]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.clabel(CS, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    plt.show()
    
    if wpharm == True:

        wind = windspharm.iris.VectorWind(x_wind, y_wind)
        streamfunction, velpotential = wind.sfvp()
        clevs = [-200, -180, -160, -120, -100, -80, -60, -40, -20, 0, 40, 80, 120, 160, 200]
        iplt.contourf(streamfunction[time_slice,level,:,:]*1e-06, clevs, cmap=brewer_redblu, extend='both')
        ax = plt.gca()
        ax.gridlines(draw_labels=True)
        plt.colorbar(orientation='horizontal')
        plt.title('Streamfunction [$10^6$ m2 s-1], h = %s km' %(level+1), y=1.20)
        plt.show()

        planet_vort = wind.planetaryvorticity(omega=omega)
        relative_vort = wind.vorticity()
        absolute_vort = wind.absolutevorticity()
        
        iplt.contourf(relative_vort[time_slice,level,:,:], brewer_bg.N, cmap=brewer_bg)
        plt.title('Relative Vorticity, h = %s km' %(level+1), y=1.20)
        plt.ylabel('Latitude [degrees]')
        plt.xlabel('Longitude [degrees]')
        ax = plt.gca()
        ax.gridlines(draw_labels=True)
        plt.colorbar(orientation='horizontal')
        plt.show()


def plot_vapour(cubes, time_slice=-1):
    
    """ Plot zonal vapour flux
        Plot meridional vapour flux
        Plot upward vapour flux
        Plot streamplot of 2-D vapour fluxes """
    
    for cube in cubes:
        if cube.long_name == 'water_vapour_flux_x':
            x_vapour = cube
        if cube.long_name == 'water_vapour_flux_y':
            y_vapour = cube
        if cube.long_name == 'water_vapour_flux_z':
            z_vapour = cube
    
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
    plt.xticks((-72,-52,-32,-12,0,12,32,52,72),('180W','140W','100W','60W','0','60E','100E','140E','180E'))
    plt.yticks((-45,-30,-15,0,15,30,45),('90S','60S','30S','0','30N','60N','90N'))    
    plt.show() 
    
    magnitude3D = iris.analysis.maths.apply_ufunc(np.sqrt, (x_vapour**2 + y_vapour**2 + z_vapour**2))
    
    global_avgs = []
    for day in range(0,x_vapour.shape[0]):
        global_avg = np.mean(magnitude3D[day,:,:].data)
        global_avgs.append(global_avg)

    plt.plot(np.arange(0,x_vapour.shape[0]), np.array(global_avgs))
    plt.xlabel('Time [months]')
    plt.ylabel('Average Total Vapour Flux [kg s-1]')
    plt.title('Vapour Flux Over Time')
    plt.show() 


def plot_efficiency(cubes):
    
    
    """ Plot day-night heat redistribution efficiency, defined in Leconte et al. 2013
        as the ratio of mean night-side outgoing longwave radiation over mean dayside 
        outgoing longwave radiation                                                  """
    
    for cube in cubes:
        if cube.standard_name == 'toa_outgoing_longwave_flux':
            outgoing_lw = cube
            
    dayside = outgoing_lw.extract(iris.Constraint(longitude=lambda v: 270 < v < 360 or 0 < v < 90, latitude=lambda v: -90 < v < 90))
    nightside = outgoing_lw.extract(iris.Constraint(longitude=lambda v: 90 < v < 270, latitude=lambda v: -90 < v < 90))
                
    dayside_avg = dayside.collapsed('latitude',iris.analysis.MEAN)
    dayside_avg = dayside_avg.collapsed('longitude',iris.analysis.MEAN)
    
    nightside_avg = nightside.collapsed('latitude',iris.analysis.MEAN)
    nightside_avg = nightside_avg.collapsed('longitude',iris.analysis.MEAN)
    
    run_length = cube.shape[0]
    months = DimCoord((np.arange(0,run_length)),standard_name='time', units='months')
    
    iplt.plot(months,(nightside_avg/dayside_avg))
    plt.title('Day-Night Heat Redistribution Efficiency')
    plt.ylabel('Ratio')
    plt.xlabel('Time [months]')
    plt.show()

                

    
    

# if __name__ == '__main__':
#     print('Source folder:')
#     directory = input()

