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
    heights = np.round(x_wind.coord('level_height').points,0)

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
                                 
    plt.plot(dayside_mean.data, np.array(heights))
    plt.title('Dayside Mean Zonal Wind [m s-1]')
    plt.ylabel('Height [m]')
    plt.xlabel('Speed [m s-1]')  
    plt.show()
    
    plt.plot(nightside_mean.data, np.array(heights))
    plt.title('Nightside Mean Zonal Wind [m s-1]')
    plt.ylabel('Height [m]')
    plt.xlabel('Speed [m s-1]')  
    plt.show()
    
    plt.plot(dayside_stdev, np.array(heights))
    plt.title('Dayside Mean Zonal Wind StDev [m s-1]')
    plt.ylabel('Height [m]')
    plt.xlabel('Speed [m s-1]')  
    plt.show()
    
    plt.plot(nightside_stdev, np.array(heights))
    plt.title('Nightside Mean Zonal Wind StDev [m s-1]')
    plt.ylabel('Height [m]')
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
    heights = np.round(x_wind.coord('level_height').points,0)

    wind = windspharm.iris.VectorWind(x_wind, y_wind)
    streamfunction, velpotential = wind.sfvp()
    clevs = [-200, -180, -160, -120, -100, -80, -60, -40, -20, 0, 40, 80, 120, 160, 200]
    iplt.contourf(streamfunction[time_slice,level,:,:]*1e-06, clevs, cmap=brewer_redblu, extend='both')
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.colorbar(orientation='horizontal')
    plt.title('Streamfunction [$10^6$ m2 s-1], h = %s km' %(heights[level]), y=1.20)
    plt.show()

    planet_vort = wind.planetaryvorticity(omega=omega)
    relative_vort = wind.vorticity()
    absolute_vort = wind.absolutevorticity()
    
    iplt.contourf(relative_vort[time_slice,level,:,:], brewer_bg.N, cmap=brewer_bg)
    plt.title('Relative Vorticity, h = %s km' %(heights[level]), y=1.20)
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
    heights = np.round(x_wind.coord('level_height').points,0)
   
    speed = iris.analysis.maths.apply_ufunc(np.sqrt, (x_wind**2 + y_wind**2))
    
    iplt.contourf(x_wind[time_slice,level,:,:], brewer_reds.N, cmap=brewer_reds)
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.title('Zonal Wind Speed, h=%s km [m s-1]' %(heights[level]), y=1.20)
    plt.colorbar(pad=0.1)
    plt.show()

    X,Y = np.meshgrid(np.arange(-72,72), np.arange(-45,45))
    fig = plt.figure(figsize = (12, 7)) 
    strm = plt.streamplot(X, Y, np.roll(x_wind[time_slice,level,:,:].data, 72, axis=1), np.roll(y_wind[time_slice,level,:,:].data, 72, axis=1), density = 0.5, color=np.roll(speed[time_slice,level,:,:].data, 72, axis=1), cmap=brewer_reds)
    # Since .data method extracts the numpy array and strips the metadata, the longitude/latitude information is lost.
    # To align plot so that (0,0) is at the center as in the Iris plots, use numpy.roll to shift columns (longitude) 180 degrees (72 places)
    fig.colorbar(strm.lines)
    plt.title('Wind speed and direction [m s-1], h=%s km' %(heights[level]))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.yticks((-45,-30,-15,0,15,30,45),('90S','60S','30S','0','30N','60N','90N'))    
    plt.show() 


def plot_zwind(cubes, time_slice=-1):
    
    for cube in cubes:
        if cube.standard_name == 'upward_air_velocity':
            z_wind = cube.copy()
            
    # dayside = z_wind.extract(iris.Constraint(longitude=lambda v: 270 < 359 or 0 <= v <= 90, latitude=lambda v: -90 <= v <= 90))
    # nightside = z_wind.extract(iris.Constraint(longitude=lambda v: 90 < v <= 270, latitude=lambda v: -90 <= v <= 90))
    dayside = z_wind.intersection(longitude=(-90,90), latitude=(-90,90))    
    nightside = z_wind.intersection(longitude=(90,270), latitude=(-90,90))

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
    

def rossby_source(cubes, time_slice=-1, level=0):  
            
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube.copy()
        if cube.standard_name == 'y_wind':
            y_wind = cube.copy()
 
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
    
    winds = windspharm.iris.VectorWind(x_wind,y_wind)
    
    eta = winds.absolutevorticity()
    div = winds.divergence()
    uchi, vchi = winds.irrotationalcomponent()
    etax, etay = winds.gradient(eta)
    etax.units = 'm**-1 s**-1'
    etay.units = 'm**-1 s**-1'
    
    S = eta*-1.*div-(uchi*etax+vchi*etay)
    S.coord('longitude').attributes['circular'] = True
    print(S.shape)
    print(S)
    
    clevs = [-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30]
    iplt.contourf(S[time_slice,level,:,:]*1e11,clevs,cmap=brewer_redblu, extend='both')
    ax = plt.gca()
    ax.gridlines(draw_labels=True)
    plt.colorbar(orientation='horizontal')
    plt.title('Rossby Wave Source ($10^{-11}$s$^{-1}$)')
    plt.show()
  

def decomposition(cubes, time_slice=-1, level=0, vapour=False):  
            
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube.copy()
        if cube.standard_name == 'y_wind':
            y_wind = cube.copy()
        if cube.long_name == 'water_vapour_flux_x':
            x_vapour = cube.copy()
        if cube.long_name == 'water_vapour_flux_y':
            y_vapour = cube.copy()
 
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
    y_vapour = y_vapour.regrid(x_vapour, iris.analysis.Linear())
    heights = np.round(x_wind.coord('level_height').points,0)
    
    winds = windspharm.iris.VectorWind(x_wind,y_wind)  
    vapour = windspharm.iris.VectorWind(x_vapour,y_vapour)
    uchi,vchi,upsi,vpsi = winds.helmholtz()
    u2chi,v2chi,u2psi,v2psi = vapour.helmholtz()
    
    chi_magnitude = iris.analysis.maths.apply_ufunc(np.sqrt, (uchi**2 + vchi**2))
    psi_magnitude =  iris.analysis.maths.apply_ufunc(np.sqrt, (upsi**2 + vpsi**2))
    chi2_magnitude = iris.analysis.maths.apply_ufunc(np.sqrt, (u2chi**2 + v2chi**2))
    psi2_magnitude =  iris.analysis.maths.apply_ufunc(np.sqrt, (u2psi**2 + v2psi**2))

    X,Y = np.meshgrid(np.arange(-72,72), np.arange(-45,45))
    fig = plt.figure(figsize = (10,4)) 
    strm = plt.streamplot(X, Y, np.roll(uchi[time_slice,level,:,:].data, 72, axis=1), np.roll(vchi[time_slice,level,:,:].data, 72, axis=1), density = 0.5, color=np.roll(chi_magnitude[time_slice,level,:,:].data, 72, axis=1), cmap=brewer_reds)
    # Since .data method extracts the numpy array and strips the metadata, the longitude/latitude information is lost.
    # To align plot so that (0,0) is at the center as in the Iris plots, use numpy.roll to shift columns (longitude) 180 degrees (72 places)
    fig.colorbar(strm.lines)
    plt.title('Divergent component of wind [m s-1], h=%s m' %(heights[level]))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.yticks((-45,-30,-15,0,15,30,45),('90S','60S','30S','0','30N','60N','90N'))    
    plt.show() 
    
    X,Y = np.meshgrid(np.arange(-72,72), np.arange(-45,45))
    fig = plt.figure(figsize = (10,4)) 
    strm = plt.streamplot(X, Y, np.roll(upsi[time_slice,level,:,:].data, 72, axis=1), np.roll(vpsi[time_slice,level,:,:].data, 72, axis=1), density = 0.5, color=np.roll(psi_magnitude[time_slice,level,:,:].data, 72, axis=1), cmap=brewer_reds)
    # Since .data method extracts the numpy array and strips the metadata, the longitude/latitude information is lost.
    # To align plot so that (0,0) is at the center as in the Iris plots, use numpy.roll to shift columns (longitude) 180 degrees (72 places)
    fig.colorbar(strm.lines)
    plt.title('Rotational component of wind [m s-1], h=%s m' %(heights[level]))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.yticks((-45,-30,-15,0,15,30,45),('90S','60S','30S','0','30N','60N','90N'))    
    plt.show() 
    
    if vapour==True:
    
        X,Y = np.meshgrid(np.arange(-72,72), np.arange(-45,45))
        fig = plt.figure(figsize = (10,4)) 
        strm = plt.streamplot(X, Y, np.roll(u2chi[time_slice,:,:].data, 72, axis=1), np.roll(v2chi[time_slice,:,:].data, 72, axis=1), density = 0.5, color=np.roll(chi2_magnitude[time_slice,:,:].data, 72, axis=1), cmap=brewer_reds)
        # Since .data method extracts the numpy array and strips the metadata, the longitude/latitude information is lost.
        # To align plot so that (0,0) is at the center as in the Iris plots, use numpy.roll to shift columns (longitude) 180 degrees (72 places)
        fig.colorbar(strm.lines)
        plt.title('Divergent component of vapour flux')
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
        plt.yticks((-45,-30,-15,0,15,30,45),('90S','60S','30S','0','30N','60N','90N'))    
        plt.show() 
        
        X,Y = np.meshgrid(np.arange(-72,72), np.arange(-45,45))
        fig = plt.figure(figsize = (10,4)) 
        strm = plt.streamplot(X, Y, np.roll(u2psi[time_slice,:,:].data, 72, axis=1), np.roll(v2psi[time_slice,:,:].data, 72, axis=1), density = 0.5, color=np.roll(psi2_magnitude[time_slice,:,:].data, 72, axis=1), cmap=brewer_reds)
        # Since .data method extracts the numpy array and strips the metadata, the longitude/latitude information is lost.
        # To align plot so that (0,0) is at the center as in the Iris plots, use numpy.roll to shift columns (longitude) 180 degrees (72 places)
        fig.colorbar(strm.lines)
        plt.title('Rotational component of vapour flux')
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
        plt.yticks((-45,-30,-15,0,15,30,45),('90S','60S','30S','0','30N','60N','90N'))    
        plt.show() 


def plot_coriolis(cubes, omega=0.64617667e-05):
    
    for cube in cubes:
        if cube.standard_name == 'air_potential_temperature':
            theta = cube.copy()
    
    earth_omega = 7.2921159e-05
    lats = theta.coord('latitude').points
    lats_rad = (theta.coord('latitude').points)*(np.pi/180) # Latitudes in radians
    coriolis = 2*omega*np.sin(lats_rad) # Coriolis force
    coriolis_earth = 2*earth_omega*np.sin(lats_rad)
    
    plt.plot(lats, coriolis, color='b', label='Proxima Centauri b')
    plt.plot(lats, coriolis_earth, color='r', label='Earth')
    plt.title('Coriolis parameter')
    plt.xlabel('Latitude [degrees]')
    plt.ylabel('Coriolis parameter')
    plt.annotate(r"${0:0.2e}$".format(coriolis[-1]), xy=(65,-0.000005), xytext=(0, 0), textcoords='offset points')
    plt.annotate(r"${0:0.2e}$".format(coriolis_earth[-1]), xy=(65,0.00010), xytext=(0, 3), textcoords='offset points')
    plt.legend()
    plt.show()
    

def plot_pressure_force(cubes, time_slice=10, low=47, high=53):
    
    for cube in cubes:
        if cube.standard_name == 'air_pressure':
            pressure = cube.copy()
        if cube.standard_name == 'air_potential_temperature':
            theta = cube.copy()
        if cube.long_name == 'density_r_r':
            density_r = cube.copy()
            
    x_gradient = iris.analysis.calculus.differentiate(pressure, 'longitude')
    longitudes = pressure.shape[3]/2     
    heights = np.round(pressure.coord('level_height').points*1e-03,0)
    
    p0 = iris.coords.AuxCoord(100000.0, long_name='reference_pressure', units='Pa')
    p0.convert_units(pressure.units)
    temperature = theta*((pressure/p0)**(287.05/1005)) # R and cp in J/kgK for 300K
    
    R = 8.314 # ideal gas constant in J/K*mol
    M = 28.0134e-03 # molar mass of N2 in kg/mol
    
    density = (pressure*R*temperature)/M
#    density = density_r/(R*R)
    pressure_force = -(1/density.data)*x_gradient.data
    
    plt.figure(figsize=(10,5))
    plt.plot(np.arange(-longitudes, longitudes), np.roll(pressure_force[time_slice,38,45,:], 72))
    plt.title('Longitudinal Pressure Gradient Acceleration at Equator, t=%s months, h=%s km' %(time_slice+1, heights[38]))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.axvline(x=-36, color='r', linestyle='--')
    plt.axvline(x=36, color='r', linestyle='--')
    plt.ylabel('Pressure gradient acceleration [m s-2]')
    plt.show()
    
    plt.figure(figsize=(10,5))
    plt.plot(np.arange(-longitudes, longitudes), np.roll(pressure_force[time_slice,low,45,:], 72))
    plt.title('Longitudinal Pressure Gradient Acceleration at Equator, t=%s months, h=%s km' %(time_slice+1, heights[low]))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.axvline(x=-36, color='r', linestyle='--')
    plt.axvline(x=36, color='r', linestyle='--')
    plt.ylabel('Pressure gradient acceleration [m s-2]')
    plt.show()
    
    plt.figure(figsize=(10,5))
    plt.plot(np.arange(-longitudes, longitudes), np.roll(pressure_force[time_slice,high,45,:], 72))
    plt.title('Longitudinal Pressure Gradient Acceleration at Equator, t=%s months, h=%s km' %(time_slice+1, heights[high]))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72,-60,-48,-36,-24,-12,0,12,24,36,48,60,72),('180W','150W','120W','90W','60W','30W','0','30E','60E','90E','120E','150E','180E'))
    plt.axvline(x=-36, color='r', linestyle='--')
    plt.axvline(x=36, color='r', linestyle='--')
    plt.ylabel('Pressure gradient acceleration [m s-2]')
    plt.show()
    