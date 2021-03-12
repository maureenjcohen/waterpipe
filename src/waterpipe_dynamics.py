# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 14:34:14 2021

@author: Mo Cohen

Pipeline for post-processing UM output data.
- Analyses model dynamics

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

def plot_zonal_wind(cubes, time_slice=-1):
    
    """ Inputs: Iris CubeList and model level you want to look at
        Plots dayside, nightside, and global zonal mean winds     """
        
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube.copy()
    
    dayside = x_wind.extract(iris.Constraint(longitude=lambda v: 270 < v <= 359 or 0 <= v <= 90, latitude=lambda v: -90 <= v <= 90))
    nightside = x_wind.extract(iris.Constraint(longitude=lambda v: 90 < v <= 270, latitude=lambda v: -90 <= v <= 90))

    x_end = x_wind[-4:,:,:,:]
    x_mean = x_end.collapsed('time',iris.analysis.MEAN)
    
    dayside_zonal_mean = dayside.collapsed('longitude', iris.analysis.MEAN)
    nightside_zonal_mean = nightside.collapsed('longitude', iris.analysis.MEAN)    
    dayside_zonal_end = dayside_zonal_mean[-4:,:,:]
    nightside_zonal_end = nightside_zonal_mean[-4:,:,:]
    dayside_zonal_time = dayside_zonal_end.collapsed('time', iris.analysis.MEAN)
    nightside_zonal_time = nightside_zonal_end.collapsed('time', iris.analysis.MEAN)
    
    CS_day = iplt.contourf(dayside_zonal_mean[time_slice,:,:], levels=np.linspace(-170,130,20), cmap=brewer_redblu)
    plt.title('Dayside Zonal Mean Zonal Wind [m s-1]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.clabel(CS_day, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    plt.show()
    
    CS_night = iplt.contourf(nightside_zonal_mean[time_slice,:,:], levels=np.linspace(-170,130,20), cmap=brewer_redblu)
    plt.title('Nightside Zonal Mean Zonal Wind [m s-1]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.clabel(CS_night, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    plt.show()
    
    CS_day = iplt.contourf(x_mean[:,:,0], levels=np.linspace(-170,130,20), cmap=brewer_redblu)
    plt.title('Mean Zonal Wind, Long 0 [m s-1]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.clabel(CS_day, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    plt.show()
    
    CS_night = iplt.contourf(x_mean[:,:,72], levels=np.linspace(-170,130,20), cmap=brewer_redblu)
    plt.title('Mean Zonal Wind, Long 180 [m s-1]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
    plt.clabel(CS_night, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    plt.show()

    
def plot_zonal_line(cubes):
    
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube.copy()
    
    height = x_wind.shape[1]
    lats = x_wind.coord('latitude')
    longs = x_wind.coord('longitude')

    if lats.bounds == None:
        x_wind.coord('latitude').guess_bounds()
    if longs.bounds == None:
        x_wind.coord('longitude').guess_bounds()
        
    dayside = x_wind.extract(iris.Constraint(longitude=lambda v: 270 < v <= 359 or 0 <= v <= 90, latitude=lambda v: -90 <= v <= 90))
    nightside = x_wind.extract(iris.Constraint(longitude=lambda v: 90 < v <= 270, latitude=lambda v: -90 <= v <= 90))
    day_grid = iris.analysis.cartography.area_weights(dayside[-12:,:,:,:])
    night_grid = iris.analysis.cartography.area_weights(nightside[-12:,:,:,:])
    dayside = dayside[-12:,:,:,:]
    nightside = nightside[-12:,:,:,:]
    
    dayside_mean = dayside.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=day_grid)
    nightside_mean = nightside.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=night_grid)
    dayside_data = dayside_mean.data
    dayside_stdev = np.std(dayside_data, axis=0)
    nightside_data = nightside_mean.data
    nightside_stdev = np.std(nightside_data, axis=0)
    
    dayside_mean = dayside_mean.collapsed('time',iris.analysis.MEAN)
    nightside_mean = nightside_mean.collapsed('time',iris.analysis.MEAN)
                                 
    plt.plot(dayside_mean.data, np.arange(0,height))
    plt.title('Dayside Mean Zonal Wind [m s-1]')
    plt.ylabel('Level')
    plt.xlabel('Speed [m s-1]')  
    plt.show()
    
    plt.plot(nightside_mean.data, np.arange(0,height))
    plt.title('Nightside Mean Zonal Wind [m s-1]')
    plt.ylabel('Level')
    plt.xlabel('Speed [m s-1]')  
    plt.show()
    
    plt.plot(dayside_stdev, np.arange(0,height))
    plt.title('Dayside Mean Zonal Wind StDev [m s-1]')
    plt.ylabel('Level')
    plt.xlabel('Speed [m s-1]')  
    plt.show()
    
    plt.plot(nightside_stdev, np.arange(0,height))
    plt.title('Nightside Mean Zonal Wind StDev [m s-1]')
    plt.ylabel('Level')
    plt.xlabel('Speed [m s-1]')  
    plt.show()
    
    
def plot_zonal_stdev(cubes):
    
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube.copy()    
    
    data = x_wind[-4:,:,:,:].data
    spatial_stdev = np.std(data, axis=0)
    
    plt.figure(figsize=(10,4))
    plt.contourf(spatial_stdev[:,:,0], cmap=brewer_bg)
    plt.title('Standard Deviation in Zonal Wind, Long 0')
    plt.ylabel('Model level')
    plt.xlabel('Latitude [degrees]')
    plt.yticks((0,40,60),('0', '40', '85'))
    plt.xticks((-90,-45,0,45,90),('90S', '45S', '0', '45N', '90N'))
    plt.colorbar(pad=0.1)
    plt.show()  


def plot_streamfunction(cubes, level=14, time_slice=-1, omega=0.64617667):
    
    """ Uses windspharm package to plot the streamfunction       """
    
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube.copy()
        if cube.standard_name == 'y_wind':
            y_wind = cube.copy()
        if cube.standard_name == 'upward_air_velocity':
            z_wind = cube.copy()    
    
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())

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

        

def plot_streamlines(cubes, level=14, time_slice=-1):
    
    """ Inputs: Iris CubeList and model level you want to look at
        wpharm=True allows windspharm library to be used on Linux machines
        omega is the rotation rate of the planet in 10^-5 rad/s. Default is Proxima Centauri b.
        
        Plot zonal wind speed at level
        Plot streamplot at level
        Plot zonal mean zonal winds for planet
        
        If wpharm=True, plot stream function and relative vorticity    """
        
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube.copy()
        if cube.standard_name == 'y_wind':
            y_wind = cube.copy()
 
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


def plot_zwind(cubes, time_slice=-1):
    
    for cube in cubes:
        if cube.standard_name == 'upward_air_velocity':
            z_wind = cube.copy()
            
    dayside = z_wind.extract(iris.Constraint(longitude=lambda v: 270 < 359 or 0 <= v <= 90, latitude=lambda v: -90 <= v <= 90))
    nightside = z_wind.extract(iris.Constraint(longitude=lambda v: 90 < v <= 270, latitude=lambda v: -90 <= v <= 90))
    
    dayside_zonal_mean = dayside.collapsed('longitude', iris.analysis.MEAN)
    nightside_zonal_mean = nightside.collapsed('longitude', iris.analysis.MEAN)
    
    CS_day = iplt.contourf(dayside_zonal_mean[time_slice,:,:], brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Dayside Zonal Mean Z Wind [m s-1]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
#    plt.clabel(CS_day, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    plt.show()
    
    CS_night = iplt.contourf(nightside_zonal_mean[time_slice,:,:], brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Nightside Zonal Mean Z Wind [m s-1]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Latitude [degrees]')
#   plt.clabel(CS_night, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    plt.show()

    dayside_meridional_mean = dayside.collapsed('latitude', iris.analysis.MEAN)
    nightside_meridional_mean = nightside.collapsed('latitude', iris.analysis.MEAN)
    
    CS_day = iplt.contourf(dayside_meridional_mean[time_slice,:,:], brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Dayside Meridional Mean Z Wind [m s-1]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Longitude [degrees]')
#    plt.clabel(CS_day, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    plt.show()
    
    CS_night = iplt.contourf(nightside_meridional_mean[time_slice,:,:], brewer_redblu.N, cmap=brewer_redblu)
    plt.title('Nightside Meridional Mean Z Wind [m s-1]', y=1.05)
    plt.ylabel('Height [m]')
    plt.xlabel('Longitude [degrees]')
#    plt.clabel(CS_night, inline=False, colors='k', fmt='%1.1f')
    plt.colorbar(pad=0.1)
    plt.show()
    

def plot_hovmoellerx(cubes):
    
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube.copy()
    
    height = x_wind.shape[1]
    run_length = x_wind.shape[0]
    lats = x_wind.coord('latitude')
    longs = x_wind.coord('longitude')

    if lats.bounds == None:
        x_wind.coord('latitude').guess_bounds()
    if longs.bounds == None:
        x_wind.coord('longitude').guess_bounds()
        
    dayside = x_wind.extract(iris.Constraint(longitude=lambda v: 270 < v <= 359 or 0 <= v <= 90, latitude=lambda v: -10<= v <= 10))
    nightside = x_wind.extract(iris.Constraint(longitude=lambda v: 90 < v <= 270, latitude=lambda v: -10 <= v <= 10))
    day_grid = iris.analysis.cartography.area_weights(dayside)
    night_grid = iris.analysis.cartography.area_weights(nightside)
    
    dayside_mean = dayside.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=day_grid)
    nightside_mean = nightside.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=night_grid)
    
    # dayside_zonal_mean = dayside.collapsed('longitude', iris.analysis.MEAN)
    # nightside_zonal_mean = nightside.collapsed('longitude', iris.analysis.MEAN)  
    
    # months = DimCoord((np.arange(0,run_length)), standard_name='time', units='months')

    iplt.contourf(dayside_mean, levels=np.linspace(-80,130,20), cmap=brewer_redblu)
    plt.title('Dayside Mean Zonal Equatorial Wind [m s-1]')
    plt.xlabel('Time')
    plt.ylabel('Height [m]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    iplt.contourf(nightside_mean, levels=np.linspace(-80,130,20), cmap=brewer_redblu)
    plt.title('Nightside Mean Zonal Equatorial Wind [m s-1]')
    plt.xlabel('Time')
    plt.ylabel('Height [m]')
    plt.colorbar(pad=0.1)
    plt.show()    
    
    # iplt.contourf(dayside_zonal_mean[:,40], levels=np.linspace(-80,130,20), cmap=brewer_redblu)
    # plt.title('Dayside Mean Zonal Wind at 40 km [m s-1]')
    # plt.xlabel('Time')
    # plt.ylabel('Latitude')
    # plt.colorbar(pad=0.1)
    # plt.show()
    
    # iplt.contourf(nightside_zonal_mean[:,40], levels=np.linspace(-80,130,20), cmap=brewer_redblu)
    # plt.title('Nightside Mean Zonal Wind at 40 km [m s-1]')
    # plt.xlabel('Time')
    # plt.ylabel('Latitude')
    # plt.colorbar(pad=0.1)
    # plt.show()
    
def plot_hovmoellerz(cubes):
    
    for cube in cubes:
        if cube.standard_name == 'upward_air_velocity':
            z_wind = cube.copy()
    
    height = z_wind.shape[1]
    run_length = z_wind.shape[0]
    lats = z_wind.coord('latitude')
    longs = z_wind.coord('longitude')

    if lats.bounds == None:
        z_wind.coord('latitude').guess_bounds()
    if longs.bounds == None:
        z_wind.coord('longitude').guess_bounds()
        
    dayside = z_wind.extract(iris.Constraint(longitude=lambda v: 270 < v <= 359 or 0 <= v <= 90, latitude=lambda v: -2 <= v <= 2))
    nightside = z_wind.extract(iris.Constraint(longitude=lambda v: 90 < v <= 270, latitude=lambda v: -2 <= v <= 2))
    day_grid = iris.analysis.cartography.area_weights(dayside)
    night_grid = iris.analysis.cartography.area_weights(nightside)
    global_grid = iris.analysis.cartography.area_weights(z_wind)
    
    dayside_mean = dayside.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=day_grid)
    nightside_mean = nightside.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=night_grid)
    global_mean = z_wind.collapsed(['latitude','longitude'], iris.analysis.MEAN, weights=global_grid)
    months = DimCoord((np.arange(0,run_length)), standard_name='time', units='months')

    iplt.contourf(dayside_mean, brewer_reds.N, cmap=brewer_reds)
    plt.title('Dayside Mean Equatorial Z-Wind [m s-1]')
    plt.xlabel('Time')
    plt.ylabel('Height [m]')
    plt.colorbar(pad=0.1)
    plt.show()
    
    iplt.contourf(nightside_mean, brewer_reds.N, cmap=brewer_reds)
    plt.title('Nightside Mean Equatorial Z-Wind [m s-1]')
    plt.xlabel('Time')
    plt.ylabel('Height [m]')
    plt.colorbar(pad=0.1)
    plt.show()



def calculate_timescales(cubes, stellar_constant=837):
    
    """ Implements formula from Koll & Abbot 2016 for:
        Timescale for equatorial waves to transport energy across the planet
        Atmosphere's radiative cooling time
        Stellar constant for FAST experiment: 1100
        Stellar constant for SLOW experiment: 420                        """

# I made every conceivable kind of error while writing this function    

    for cube in cubes:
        if cube.long_name == 'surface_diffuse_albedo_assuming_no_snow': # long name not standard name
            albedo = cube.copy()
        # if cube.standard_name == 'specific_humidity':
        #     q = cube.copy()
            
    R = 297.0 # gas constant for air in J/kgK, from UM Planet Constants
    a = 7160000 # planet's radius in m, from UM Planet Constants
    g = 10.9 # mean gravity in m/s2, from UM Planet Constants
    sigma = 5.67e-8 # Stefan-Boltzmann constant, SI units NOT BOLTZMANN'S CONSTANT GODDAMN ONLY 20 ORDERS OF MAGNITUDE OFF
    pressure = 100000 # surface pressure in Pa, from UM Planet Constants DO NOT PUT A COMMA IN HERE!! NUMPY NOT LIKE!!
    cp = 1038 # specific heat capacity at constant pressure for N2, J/kgK
#    solar_constant = 0.0017*3.828e26 # luminosity of Proxima Centauri in W, wrong constant, WHY DID U USE THE SYMBOL FOR STELLAR LUMINOSITY THEN DANIEL AND DORIAN
    
    mean_albedo = albedo.collapsed(['time','pseudo_level'], iris.analysis.MEAN)
    iplt.contourf(mean_albedo, brewer_bg.N, cmap=brewer_bg)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Mean Surface Albedo [K]', y=1.20)
    plt.colorbar(pad=0.1)
    plt.show()
    
    lats = albedo.coord('latitude')
    longs = albedo.coord('longitude')

    if lats.bounds == None:
        albedo.coord('latitude').guess_bounds()
    if longs.bounds == None:
        albedo.coord('longitude').guess_bounds()
    
    global_grid = iris.analysis.cartography.area_weights(albedo)
    
    alpha = albedo.collapsed(['latitude', 'longitude', 'time', 'pseudo_level'], iris.analysis.MEAN, weights=global_grid)
    alpha = alpha.data   
    
    Teq = ((stellar_constant*(1-alpha))/(4*sigma))**(1/4) # brackets matter
    cwave = R*np.sqrt(Teq/cp)  
    twave = a/cwave
    trad = (cp*pressure)/(g*sigma*(Teq**3)) 
    ratio = twave/trad
    
    print('Ratio of wave to radiative timescale is ' + str(ratio))
    print('Wave timescale is ' + str(twave) + ' seconds')
    print('Radiative timescale is ' + str(trad) + ' seconds')
    
    return ratio

    
def qbo_period(cubes, periodicity=False):
    
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube.copy()
    
    height = x_wind.shape[1]
    run_length = x_wind.shape[0]
    lats = x_wind.coord('latitude')
    longs = x_wind.coord('longitude')

    if lats.bounds == None:
        x_wind.coord('latitude').guess_bounds()
    if longs.bounds == None:
        x_wind.coord('longitude').guess_bounds()
        
    strat = x_wind.extract(iris.Constraint(longitude=lambda v: 355 < v <= 359 or 0 <= v <= 4, latitude=lambda v: -4 <= v <= 4, model_level_number=40))
    trop = x_wind.extract(iris.Constraint(longitude=lambda v: 355 < v <= 359 or 0 <= v <= 4, latitude=lambda v: -4 <= v <= 4, model_level_number=24))
    high = x_wind.extract(iris.Constraint(longitude=lambda v: 355 < v <= 359 or 0 <= v <= 4, latitude=lambda v: -4 <= v <= 4, model_level_number=50))
    low = x_wind.extract(iris.Constraint(longitude=lambda v: 355 < v <= 359 or 0 <= v <= 4, latitude=lambda v: -4 <= v <= 4, model_level_number=29))
    strat_grid = iris.analysis.cartography.area_weights(strat)
    trop_grid = iris.analysis.cartography.area_weights(trop)
    high_grid = iris.analysis.cartography.area_weights(high)
    low_grid = iris.analysis.cartography.area_weights(low)
    
    strat_mean = strat.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=strat_grid)
    trop_mean = trop.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=trop_grid)
    
    high_mean = high.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=high_grid)
    low_mean = low.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=low_grid)
    
    plt.plot(np.arange(0,run_length), high_mean.data, linestyle='--', color='r', label='50 km')
    plt.plot(np.arange(0,run_length), strat_mean.data, linestyle='-', color='r', label='40 km')
    plt.plot(np.arange(0,run_length), low_mean.data, linestyle='--', color='b', label='30 km')
    plt.plot(np.arange(0,run_length), trop_mean.data, linestyle='-', color='b', label='25 km')
    plt.title('Mean Substellar Zonal Wind')
    plt.xlabel('Time [months]') 
    plt.ylabel('Velocity [m s-1]')
    plt.legend()
    plt.show()
    
    if periodicity == True:
        
        data = np.array(high_mean.data)
        fft = sp.fftpack.fft(data)
        psd = np.abs(fft)**2
        freq = sp.fftpack.fftfreq(len(psd), 1./run_length)
        i = freq > 0
        
        data2 = np.array(strat_mean.data)
        fft2 = sp.fftpack.fft(data2)
        psd2 = np.abs(fft2)**2
        freq2 = sp.fftpack.fftfreq(len(psd2), 1./run_length)
        i2 = freq2 > 0
        
        periods_50km = np.round(run_length/freq[i], 2)
        fig, ax = plt.subplots(1,1,figsize=(8,4))
        ax.plot(periods_50km, psd[i])
        ax.set_xlabel('Period [months]')
        ax.set_ylabel('PSD')
        ax.set_title('Periodicity of mean substellar zonal wind at 50 km')
        psd_sorted = -np.sort(-psd[i])
        top_three = psd_sorted[:3]
        for i,j in zip(periods_50km,psd[i]):
            if j in top_three:
                ax.annotate('%s' %i, xy=(i,j), xytext=(2, 5), textcoords='offset points')
        # plt.annotate(str(maxpsd), location)
        plt.show()
        
        periods_40km = np.round(run_length/freq2[i2], 2)
        fig, ax = plt.subplots(1,1,figsize=(8,4))
        ax.plot(periods_40km, psd2[i2])
        ax.set_xlabel('Period [months]')
        ax.set_ylabel('PSD')
        ax.set_title('Periodicity of mean substellar zonal wind at 40 km')
        psd_sorted2 = -np.sort(-psd2[i2])
        top_three2 = psd_sorted2[:3]
        for i,j in zip(periods_40km,psd2[i2]):
            if j in top_three2:
                ax.annotate('%s' %i, xy=(i,j), xytext=(2, 5), textcoords='offset points')
        # plt.annotate(str(maxpsd), location)
        plt.show()
        
    
    