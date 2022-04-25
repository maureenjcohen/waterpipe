#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 12:00:11 2021

@author: Mo Cohen
"""

import iris
import pywt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
from matplotlib import ticker
from matplotlib.colors import TwoSlopeNorm
import scipy as sp
from iris.coord_systems import GeogCS
from iris.analysis import calculus



brewer_redblu = mpl_cm.get_cmap('RdBu_r')
reds = mpl_cm.get_cmap('Reds')
redblu = mpl_cm.get_cmap('RdBu_r')
magma = mpl_cm.get_cmap('magma')


# Import packages used to postprocess data
# The iris package is used for I/O only

""" The project data is stored in two NETCDF4 files, qbo.nc and extras.nc.
The file qbo.nc contains six-hourly data from the 'default' Unified Model
vn 11.9 configured for Proxima Centauri b. 

qbo.nc contains the data cubes:
    Air pressure
    Potential temperature
    Specific humidity
    Eastward wind
    Northward wind
    Upward wind
    
The air pressure and potential temperature were used to calculate the absolute
temperature oscillation; the eastward wind for Hovmoeller plots, latitudinal
extent, u-prime, and asymmetric westward/eastward transitions; the eastward
and upward wind were used for the wave momentfum flux; and the eastward and
northward wind for the Lagrangian Rossby number.

Six hourly data was required to visualise the resolved gravity waves.

extras.nc contains the data cubes:
    Deep convection indicator (30 day average)
    Specific humidity (6-hourly, from a later part of the run than qbo.nc)
    X_wind for horizontal resolution increased to 180 columns
    X_wind for timestep reduced from 20 to 2 minutes
    X_wind, y_wind, upward_wind for precipitation-based parameterised gravity
    wave source
    
    """


def plot_hovmoellerx(cubes, radius=7160000, time='days'):
    
    """ Plot Hovmoeller diagrams of zonal wind 
        Arguments: CubeList, time (can be '6-hours' or 'days')
        Outputs: 1) Hovmoeller plot of mean equatorial wind between -10 and 10 
                lat, axes time vs. height
                 2) Hovmoeller plot of zonal mean wind at 40 km, axes time vs. 
                latitude (top view) """

    for cube in cubes:
        if cube.standard_name == 'eastward_wind' or cube.standard_name == 'x_wind':
            x_wind = cube.copy()

    x_wind.coord('latitude').coord_system = GeogCS(radius)
    x_wind.coord('longitude').coord_system = GeogCS(radius)

    run_length, lats, longs = x_wind.shape[0], x_wind.coord(
        'latitude'), x_wind.coord('longitude')

    for coord in x_wind.coords():
        if coord.long_name == 'Hybrid height':
            heights = np.round(x_wind.coord('Hybrid height').points*1e-03, 0)
        elif coord.long_name == 'level_height':
            heights = np.round(x_wind.coord('level_height').points*1e-03, 0)

    if lats.bounds == None:
        x_wind.coord('latitude').guess_bounds()
    if longs.bounds == None:
        x_wind.coord('longitude').guess_bounds()

    dayside, nightside = x_wind.intersection(
        longitude=(-90, 89), latitude=(-10, 10)), x_wind.intersection(longitude=(90, 270), latitude=(-10, 10))
    dayside_full, nightside_full = x_wind.intersection(
        longitude=(-90, 89)), x_wind.intersection(longitude=(90, 270))
    day_grid, night_grid = iris.analysis.cartography.area_weights(
        dayside), iris.analysis.cartography.area_weights(nightside)

    dayside_mean = dayside.collapsed(
        ['latitude', 'longitude'], iris.analysis.MEAN, weights=day_grid)
    nightside_mean = nightside.collapsed(
        ['latitude', 'longitude'], iris.analysis.MEAN, weights=night_grid)

    dayside_zonal_mean = dayside_full.collapsed(
        'longitude', iris.analysis.MEAN)
    nightside_zonal_mean = nightside_full.collapsed(
        'longitude', iris.analysis.MEAN)

    for cube in (dayside_mean, nightside_mean):
        if cube == dayside_mean:
            side = 'Dayside'
        else:
            side = 'Nightside'
        plt.contourf(np.arange(0, run_length)*0.25, np.array(heights), cube.data.T,
                     levels=np.arange(-120, 121, 20), cmap=brewer_redblu, norm=TwoSlopeNorm(0))
        plt.title('%s Mean Zonal Equatorial Wind' % (side))
        plt.xlabel('Time [%s]' % (time))
        plt.ylabel('Height [km]')
        cbar = plt.colorbar(pad=0.1)
        cbar.set_ticks(np.arange(-120, 121, 20))
        cbar.ax.set_title('m/s')
        # plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/laso/epsfigs/hovmoeller_%s_new.eps' %(side), format='eps')
        plt.show()

    for cube in (dayside_zonal_mean, nightside_zonal_mean):
        if cube == dayside_zonal_mean:
            side = 'Dayside'
        else:
            side = 'Nightside'
        plt.contourf(np.arange(0, run_length)*0.25, np.arange(-45, 45),
                     cube[:, 47, :].data.T, np.arange(-120, 121, 20), cmap=brewer_redblu, norm=TwoSlopeNorm(0))
        plt.title('%s Mean Zonal Wind at 41 km' % (side))
        plt.xlabel('Time [%s]' % (time))
        plt.ylabel('Latitude')
        plt.yticks((-45, 0, 45), ('90S', '0', '90N'))
        mbar = plt.colorbar(pad=0.1)
        mbar.set_ticks(np.arange(-120, 121, 20))
        mbar.ax.set_title('m/s')
        # plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/laso/epsfigs/topview_%s_new.eps' %(side), format='eps')
        plt.show()


def plot_temp_anomaly(cubes, period=(0, 220), lat=45, level=47):
    """ Plot temperature anomaly, temperature minus the local time-averaged temperature
        Default latitude is the equator
        Arguments: CubeList, latitude, and atmospheric level"""

    for cube in cubes:
        if cube.standard_name == 'air_potential_temperature':
            theta = cube[period[0]:period[1], :, lat, :].copy()
        if cube.standard_name == 'air_pressure':
            pressure = cube[period[0]:period[1], :, lat, :].copy()

    run_length, longitudes = theta.shape[0], theta.shape[2]/2

    for coord in theta.coords():
        if coord.long_name == 'Hybrid height':
            heights = np.round(theta.coord('Hybrid height').points*1e-03, 0)
        elif coord.long_name == 'level_height':
            heights = np.round(theta.coord('level_height').points*1e-03, 0)

    p0 = iris.coords.AuxCoord(
        100000.0, long_name='reference_pressure', units='Pa')
    p0.convert_units(pressure.units)
    # R and cp in J/kgK for 300K
    temperature = theta*((pressure/p0)**(287.05/1005))

    temp_time_mean = temperature.collapsed('time', iris.analysis.MEAN)
    anomaly = temperature - temp_time_mean

    plt.figure(figsize=(10, 5))
    plt.contourf(np.arange(-longitudes, longitudes), np.arange(period[0], period[1]), np.roll(
        anomaly[:, level, :].data, 72, axis=1), np.arange(-20, 21, 5), cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    plt.title('Temperature Anomaly at Equator, h=%s km' % (heights[level]))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-144, -120, -96, -72, -48, -24, 0, 24, 48, 72, 96, 120, 144), ('180W', '150W',
                '120W', '90W', '60W', '30W', '0', '30E', '60E', '90E', '120E', '150E', '180E'))
    # plt.xticks((-72, -60, -48, -36, -24, -12, 0, 12, 24, 36, 48, 60, 72), ('180W', '150W',
    #            '120W', '90W', '60W', '30W', '0', '30E', '60E', '90E', '120E', '150E', '180E'))
    plt.ylabel('Time [months]')
    cbar = plt.colorbar(pad=0.1)
    cbar.set_ticks(np.arange(-20, 21, 5))
    cbar.ax.set_title('K')
    # plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/laso/epsfigs/abs_temp_anomaly_%s_new.eps' %(heights[level]), format='eps')
    plt.show()


def wave_acceleration(cubes, density,hlevel=47, lat=45, long=0, start=2880, end=3240, plot=False):
    
    """ This function plots several things:
        1) The high-pass filtered u-prime. We need to filter u-prime because
    the standing Rossby waves create an asymmetrical background flow. Taking
    the zonal anomaly without filtering out long Rossby waves results in a 
    zonal anomaly that represents these waves instead of the transient gravity
    waves we're looking for. We filter out wavelengths of wave number < 5.1.
        2) The time anomaly of the upward wind, w'. We use the time anomaly
    instead of zonal anomaly because the asymmetrical dayside/nightside 
    rising/subsiding air pattern creates the same problem as the standing 
    Rossby waves do for the zonal wind.
        3) The wave momentum flux, calculated by multiplying 1) and 2) and
    differentiating in the vertical direction. We plot the latitudinal 
    cross-section of the zonal mean wave momentum flux together with the change
    in zonal mean zonal wind to see if the direction of acceleration matches the
    change in direction of the zonal wind. We also plot the longitudinal 
    cross-section along the equator (not the meridional mean, just a slice) and
    overlay the longitudinal cross-section of the zonal wind to show that the
    areas of high wave-induced acceleration are collocated with the jet exit region
    at the western terminator and the deep convection zone to the west of the 
    substellar point.
    
    Note: The wave-induced acceleration is calculated ONLY for the resolved,
    larger-than-gridbox gravity waves and does not include the acceleration
    from parameterised gravity waves.
    
    You must use the six-hourly data with this function.
    
    Paper plot time spans: 2880-3240 and 3340-3680"""

    for cube in cubes:
        if cube.standard_name == 'eastward_wind':
            x_wind = cube[start:end, :, :, :].copy()
        if cube.standard_name == 'upward_air_velocity':
            z_wind = cube[start:end, :, :, :].copy()

    z_wind = z_wind.regrid(x_wind, iris.analysis.Linear())

    vertical = [('Hybrid height', x_wind.coord('Hybrid height').points)]
    z_wind = z_wind.interpolate(vertical, iris.analysis.Linear())

    longitudes = x_wind.shape[3]/2
    latitudes = x_wind.shape[2]/2
    heights = x_wind.coord('Hybrid height').points
    heights_km = np.round(x_wind.coord('Hybrid height').points*1e-03, 0)

    """ Filter zonal wind, remove longer wavelengths, calculate u-prime"""

    x_data = x_wind.data

    u_prime_list = []

    for time in range(0, x_wind.shape[0]):
        time_list = []

        for level in range(0, x_wind.shape[1]):
            level_list = []

            for latx in range(0, x_wind.shape[2]):

                u_fft = sp.fftpack.fft(x_data[time, level, latx, :])
                u_psd = np.abs(u_fft)**2
                u_freq = sp.fftpack.fftfreq(len(u_psd), 1./144)

                highpass = u_fft.copy()
                highpass[np.abs(u_freq) < 5.1] = 0

                u_cleaned = np.real(sp.fftpack.ifft(highpass))
                u_bar = np.mean(u_cleaned)
                u_prime = u_cleaned - u_bar
                level_list.append(u_prime)

            time_list.append(level_list)
        u_prime_list.append(time_list)

    u_prime = np.array(u_prime_list)
    print(u_prime.shape)

    if plot == True:
        
        zonal_mean = x_wind.collapsed('longitude', iris.analysis.MEAN)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16,8), gridspec_kw={'width_ratios': [4,1]}, sharey=True)
        up = ax1.contourf(np.arange(-longitudes, longitudes), np.array(heights_km), np.roll(
            u_prime[-1, :, lat, :], 72, axis=1), np.linspace(-35, 35, 70), cmap=brewer_redblu, norm=TwoSlopeNorm(0))
        ax1.set_title('$U^{\prime}$ at Equator, t=%s days' % (end/4))
        ax1.set_xlabel('Longitude [degrees]')
        ax1.set_xticks((-72, -48, -24, 0, 24, 48, 72), ('180W',
                   '120W', '60W', '0', '60E', '120E', '180E'))
        ax1.set_ylabel('Height [km]')
        # ax1.annotate('', xy=(0.177, 0.55), xytext=(0.11, 0.55), xycoords='figure fraction', arrowprops=dict(facecolor='black',shrink=0.005),
        #              horizontalalignment='left', verticalalignment='top')
        cbar = plt.colorbar(up, pad=0.05, ax=ax1)
        cbar.locator = ticker.AutoLocator()
        cbar.update_ticks()
        cbar.ax.set_title('m/s')
        
        # plt.figure(figsize=(10, 5))
        # plt.contourf(np.arange(-longitudes, longitudes), np.array(heights_km), np.roll(
        #     u_prime[-1, :, lat, :], 72, axis=1), np.linspace(-35, 35, 70), cmap=brewer_redblu, norm=TwoSlopeNorm(0))
        # plt.title('$U^{\prime}$ at Equator, t=%s days' % (end/4))
        # plt.xlabel('Longitude [degrees]')
        # plt.xticks((-72, -60, -48, -36, -24, -12, 0, 12, 24, 36, 48, 60, 72), ('180W', '150W',
        #            '120W', '90W', '60W', '30W', '0', '30E', '60E', '90E', '120E', '150E', '180E'))
        # plt.ylabel('Height [km]')
        # cbar = plt.colorbar(pad=0.1)
        # cbar.locator = ticker.AutoLocator()
        # cbar.update_ticks()
        # cbar.ax.set_title('m/s')

        # # plt.savefig(
        #     # '/exports/csce/datastore/geos/users/s1144983/papers/laso/epsfigs/uprime_%s_ticks.eps' % (end), format='eps')
        # plt.show()
        
        
        # plt.figure(figsize=(2.5, 5))
        zmzw = ax2.contour(np.arange(-latitudes,latitudes)*2, np.array(heights_km), zonal_mean[-1,:,:].data, np.arange(-80, 81, 20), colors='black', linewidths=2.0)
        ax2.set_title('Zonal Mean Zonal Wind, t=%s days' % (end/4))
        ax2.set_xlabel('Latitude [degrees]')
        # ax2.set_ylabel('Height [km]')
        ax2.clabel(zmzw, inline=False, colors='k', fmt='%1.1f')
        fig.tight_layout()
        # fig.suptitle('$U^{\prime}$ at Equator and Zonal Mean Zonal Wind, t=%s days' %(end/4), fontsize=18)
    
        # plt.savefig(
            # '/exports/csce/datastore/geos/users/s1144983/papers/laso/epsfigs/uprime_subplots_%s_right.eps' % (end), format='eps')

        plt.show()
        

    """ Calculate w-prime"""
    
    z_mean = z_wind.collapsed('t', iris.analysis.MEAN)
    z_anomaly = z_wind - z_mean

    if plot == True:

        plt.figure(figsize=(10, 5))
        plt.contourf(np.arange(-longitudes, longitudes), np.array(heights_km), np.roll(
            z_anomaly[-1, :, lat, :].data, 72, axis=1), np.linspace(-0.09, 0.09, 20), cmap=brewer_redblu, norm=TwoSlopeNorm(0))
        plt.title(
            'Vertical Wind Anomaly at Equator [m s-1], t = %s days' % (end/4))
        plt.xlabel('Longitude [degrees]')
        plt.xticks((-72, -60, -48, -36, -24, -12, 0, 12, 24, 36, 48, 60, 72), ('180W', '150W',
                   '120W', '90W', '60W', '30W', '0', '30E', '60E', '90E', '120E', '150E', '180E'))
        plt.ylabel('Height [km]')
        plt.colorbar(pad=0.1)
        plt.show()

    """ Differentiate wrt height"""

    array = u_prime*z_anomaly.data

    h = np.array([heights])
    h = np.repeat(h[:, np.newaxis], (end-start), axis=0)
    h = np.reshape(h, ((end-start), 60))
    h = np.repeat(h[:, :, np.newaxis], 90, axis=2)
    h = np.repeat(h[:, :, :, np.newaxis], 144, axis=3)

    print(array.shape, h.shape)

    acceleration = array.copy()*density[start:end,:-1,:,:]
    acceleration[:, 0, ...] = (
        array[:, 1, ...]-array[:, 0, ...])/(h[:, 1, ...]-h[:, 0, ...])
    acceleration[:, -1, ...] = (array[:, -1, ...] -
                                array[:, -2, ...])/(h[:, -1, ...]-h[:, -2, ...])
    acceleration[:, 1:-1, ...] = (array[:, 2:, ...] -
                                  array[:, 0:-2, ...])/(h[:, 2:, ...]-h[:, 0:-2, ...])
    acceleration = -acceleration
    # Negative because Plumb 1977 paper formula has a negative, just so you don't forget this and freak out

    """ Plot """
    zonal_acc = np.mean(
        acceleration, axis=3)  
    net_acc = np.mean(zonal_acc, axis=0)

    x_mean = x_wind.collapsed('longitude', iris.analysis.MEAN)
    x_avg = x_wind.collapsed('t', iris.analysis.MEAN)
    mean_acc = np.mean(acceleration, axis=0)

    plt.figure(figsize=(10, 5))
    plt.contourf(np.arange(-longitudes, longitudes), np.array(heights_km),
                 np.roll(mean_acc[:, lat, :], 72, axis=1), np.linspace(-8e-05, 8e-05, 20), cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    mbar = plt.colorbar(pad=0.1)
    mbar.locator = ticker.AutoLocator()
    mbar.update_ticks()
    mbar.set_label('$m/s^2$')
    contours = plt.contour(np.arange(-longitudes, longitudes), np.array(heights_km), np.roll(
        x_avg[:, lat, :].data, 72, axis=1), np.arange(-80, 81, 20), colors='black', linewidths=0.3)
    plt.title('Mean Wave-Induced Acceleration at Equator, t=%s to %s days' %
              (start/4, end/4))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72, -60, -48, -36, -24, -12, 0, 12, 24, 36, 48, 60, 72), ('180W', '150W',
               '120W', '90W', '60W', '30W', '0', '30E', '60E', '90E', '120E', '150E', '180E'))
    plt.ylabel('Height [km]')
    plt.clabel(contours, inline=False, colors='k', fmt='%1.1f')
    # plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/laso/epsfigs/jetexit_%s_ticks_by20s.eps' %(end), format='eps')
    plt.show()


    x_axis = x_wind.coord('latitude').points
    y_axis = np.round(x_wind.coord('Hybrid height').points*1e-03, 0)
    plt.figure(figsize=(8, 10))
    plt.contourf(x_axis, y_axis, net_acc, brewer_redblu.N,
                 cmap=brewer_redblu, norm=TwoSlopeNorm(0))
    wbar = plt.colorbar(pad=0.1)
    wbar.locator = ticker.AutoLocator()
    wbar.update_ticks()
    wbar.set_label('$m/s^2$', rotation=90)
    CS = plt.contour(
        x_axis, y_axis, (x_mean[-1, :, :].data - x_mean[0, :, :].data), colors='black', linewidths=1.5)
    plt.title('Mean Zonal Mean Wave-Induced Acceleration, t=%s to %s days' %
              (start/4, end/4))
    plt.xlabel('Latitude [degrees]')
    plt.ylabel('Height [km]')
    plt.clabel(CS, inline=False, colors='k', fmt='%1.1f')
    # plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/laso/epsfigs/waveinducedacc_%s_ticks.eps' %(end), format='eps')
    plt.show()
    
    zonal_acc_day = np.mean(zonal_acc[-4:-1,:,:], axis=0)
    latmean = np.mean(zonal_acc_day[:,40:51], axis=1)
    latwind = x_mean[-4:-1,:,40:51].collapsed(['latitude', 't'], iris.analysis.MEAN)
    
    fig, ax1 = plt.subplots(figsize=(8, 10))
    ax1.set_xlabel('Acceleration [m/s$^2$]')
    ax1.set_ylabel('Height [km]')
    ax1.plot(latmean, y_axis, color='b', label='Acc')
    ax1.tick_params(axis='x', labelcolor='b')
    
    ax2 = ax1.twiny()
    ax2.set_xlabel('Zonal mean equatorial wind [m/s]')
    ax2.plot(latwind.data, y_axis, color='r', label='Wind')
    ax2.tick_params(axis='x', labelcolor='r')
    
    plt.title('Mean Zonal Mean Equatorial Wave-Induced Acceleration, t=%s to %s days' %
              (start/4, end/4))
    # fig.tight_layout()
    plt.show()



def deep_convection(cubes, time=39):
    
    """ Simple plot of the deep convection indicator.
    This cube is found in the extras.nc file. 
    The 30-day period plotted should match the 30-day period plotted
    in the Lagrangian Rossby wave function, i.e. month 39 is equivalent to
    days 1170-1200 of the 1800-day run."""

    for cube in cubes:
        if cube.long_name == 'deep convection indicator':
            data = cube.copy()

    run_length = np.arange(0, data.shape[0])
    longitudes = data.shape[2]
    latitudes = data.shape[1]

    plt.figure(figsize=(10, 5))
    plt.contourf(np.arange(-longitudes/2, longitudes/2), np.arange(-latitudes/2, latitudes/2),
                 np.roll(data[time, :, :].data, 144, axis=1), reds.N, cmap=reds)
    plt.title('Deep convection indicator')
    plt.xlabel('Longitude [degrees]')
    plt.ylabel('Latitude [degrees]')
    plt.xticks((-144, -120, -96, -72, -48, -24, 0, 24, 48, 72, 96, 120, 144), ('180W', '150W',
                '120W', '90W', '60W', '30W', '0', '30E', '60E', '90E', '120E', '150E', '180E'))
    plt.yticks((-90,-60, -30, 0, 30, 60, 90), ('90S', '60S', '30S', '0', '30N', '60N', '90N'))
    # plt.xticks((-72, -60, -48, -36, -24, -12, 0, 12, 24, 36, 48, 60, 72), ('180W', '150W',
    #            '120W', '90W', '60W', '30W', '0', '30E', '60E', '90E', '120E', '150E', '180E'))
    # plt.yticks((-45, -30, -15, 0, 15, 30, 45),
    #            ('90S', '60S', '30S', '0', '30N', '60N', '90N'))
    cbar = plt.colorbar(pad=0.1)
    # plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/laso/epsfigs/deepconvection.eps', format='eps')
    plt.show()


def plot_number(cubes, radius=7160000, omega=0.64617667e-05, g=10.9, start=1080, end=1200, time_slice=-1, lat_slice=45, level=47):
    
    """ Level 47 is 41 km, level 20 is 8 km
    
    The Lagrangian Rossby number diagnoses areas of instability in the flow
    which generate gravity waves. These areas are, in particular, the area of
    deep convection to the west of the substellar point and the jet exit
    region at the western terminator. 
    
    Modifications made to adapt the RoL to a slowly rotating planet are 
    described in the LASO paper in the Methods section."""

    for cube in cubes:
        if cube.standard_name == 'eastward_wind' or cube.standard_name == 'x_wind':
            x_wind = cube[start:end, :-5, 25:66, :].copy()
        if cube.standard_name == 'northward_wind' or cube.standard_name == 'y_wind':
            y_wind = cube[start:end, :-5, 25:66, :].copy()
        # if cube.standard_name =='air_potential_temperature':
        #     theta = cube[start:end,:,:,:].copy()

    vertical, latitudes, longitudes = x_wind.shape[1], x_wind.shape[2], x_wind.shape[3]
    x_axis, y_axis = x_wind.coord(
        'longitude').points, x_wind.coord('latitude').points

    for coord in x_wind.coords():
        if coord.long_name == 'Hybrid height':
            heights = np.round(x_wind.coord('Hybrid height').points*1e-03, 0)
            tcoord = 't'
            hcoord = 'Hybrid height'
        elif coord.long_name == 'level_height':
            heights = np.round(x_wind.coord('level_height').points*1e-03, 0)
            tcoord = 'time'
            hcoord = 'level_height'

    # d_theta = iris.analysis.calculus.differentiate(theta, 'Hybrid height')
    # bv_freq = np.sqrt(np.abs((g/theta[:,:-1,:,:].data)*d_theta.data))
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())

    # x_wind = x_wind.collapsed(tcoord, iris.analysis.MEAN)
    # y_wind = y_wind.collapsed(tcoord, iris.analysis.MEAN)

    x_diff = iris.analysis.calculus.differentiate(x_wind, tcoord)
    y_diff = iris.analysis.calculus.differentiate(y_wind, tcoord)
    x_grad = iris.analysis.calculus.differentiate(x_wind, 'longitude')
    y_grad = iris.analysis.calculus.differentiate(y_wind, 'latitude')
    print(np.min(y_diff.data), np.max(y_diff.data))

    mag = np.sqrt(x_wind.data**2 + y_wind.data**2)
    lats = x_wind.coord('latitude').points*(np.pi/180)
    # u_bar = x_wind.collapsed(['longitude', hcoord, tcoord], iris.analysis.MEAN)
    # f = 2*np.sin(lats)*(u_bar.data/(radius*np.cos(lats)))
    f10 = 2*omega*np.sin(lats[25])
    # f10 = 2*np.sin(lats[49])*(u_bar[49].data/(radius*np.cos(lats[49])))
    f = 2*omega*np.sin(lats)
    for i in range(15, 25):
        f[i] = f10
    # print(f)

    f = np.repeat(f[:, np.newaxis], (end-start), axis=0)
    f = np.reshape(f, ((end-start), latitudes))
    f = np.repeat(f[:, np.newaxis, :], vertical, axis=1)
    f = np.repeat(f[:, :, :, np.newaxis], longitudes, axis=3)
    # f = 2*np.sin(lats)*(u_bar.data/(radius*np.cos(lats)))
    # f = np.repeat(f[:,:,:,np.newaxis], 144, axis=3)
    print(f.shape)

    # numerator = np.sqrt(x_diff[:,:,:-1,:].data**2 + y_diff[:,:,:-1,:].data**2 + (x_wind[:-1,:,:-1,:].data*x_grad[:-1,:,:-1,:].data)**2 + (y_wind[:-1,:,:-1,:].data*y_grad[:-1,:,:,:].data)**2)
    numerator = np.sqrt((x_wind[:-1, :, :-1, :].data*x_grad[:-1, :, :-1, :].data)
                        ** 2 + (y_wind[:-1, :, :-1, :].data*y_grad[:-1, :, :, :].data)**2)
    print(np.min(numerator), np.max(numerator))
    # numerator = np.sqrt(x_diff.data**2 + y_diff.data**2)
    denominator = np.abs(f[:-1, :, :-1, :]*mag[:-1, :, :-1, :])
    print(np.min(denominator), np.max(denominator))
    rol = numerator/denominator
    # rol = np.log(numerator/denominator)

    # rol = rol/np.mean(rol)
    # print(rol.shape)
    rol = np.mean(rol, axis=0)
    rol = rol/np.min(rol)
    # rol = (rol - np.min(rol))/(np.max(rol) - np.min(rol))
    # print(np.where(rol==1.0))
    print(np.min(rol), np.max(rol))

    plt.figure(figsize=(10, 5))
    plt.contourf(np.arange(-longitudes/2, longitudes/2), np.array(heights),
                 np.roll(rol[:, lat_slice, :], 72, axis=1), np.linspace(0, 60, 30), cmap=reds)
    plt.title('Mean Lagrangian Rossby Number at Equator, t=%s to %s days' %
              (start/4, end/4))
    # plt.title('Mean Lagrangian Rossby Number at Equator')
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72, -60, -48, -36, -24, -12, 0, 12, 24, 36, 48, 60, 72), ('180W', '150W',
               '120W', '90W', '60W', '30W', '0', '30E', '60E', '90E', '120E', '150E', '180E'))
    plt.ylabel('Height [km]')
    cbar = plt.colorbar(pad=0.1)
    cbar.locator = ticker.AutoLocator()
    cbar.update_ticks()
    # cbar.set_ticklabels(['0','10','20','30','40','50','60'])
    # plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/laso/epsfigs/lrn_equator_normed.eps', format='eps')
    plt.show()

    plt.figure(figsize=(10, 5))
    plt.contourf(np.arange(-longitudes/2, longitudes/2),  y_axis[:-1], np.roll(
        rol[level, :, :], 72, axis=1), np.linspace(0, 60, 30), cmap=reds)
    plt.title('Mean Lagrangian Rossby Number at h=%s km, t=%s to %s days' %
              (heights[level], start/4, end/4))
    # plt.title('Mean Lagrangian Rossby Number at h=%s km' %(heights[level]))
    plt.xlabel('Longitude [degrees]')
    plt.xticks((-72, -60, -48, -36, -24, -12, 0, 12, 24, 36, 48, 60, 72), ('180W', '150W',
               '120W', '90W', '60W', '30W', '0', '30E', '60E', '90E', '120E', '150E', '180E'))
    plt.ylabel('Latitude [degrees]')
    plt.yticks((-40, -20, 0, 20, 40), ('40S', '20S', '0', '20N', '40N'))
    mbar = plt.colorbar(pad=0.1)
    mbar.locator = ticker.AutoLocator()
    mbar.update_ticks()
    # mbar.set_ticks([1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0])
    # mbar.set_ticklabels(['1','0.9', '0.8', '0.7', '0.6', '0.5', '0.4', '0.3', '0.2', '0.1', '0.0'])
    # plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/laso/epsfigs/lrn_%s_normed.eps' %(heights[level]), format='eps')
    plt.show()


def wavelets(cubes, radius=7160000, scales=512, wavelet='mexh', x=(106, 110), y=(43, 47), level=47, sampling=0.25):
    """ Performs wavelet transform for 1-D data
    Inputs: Iris CubeList, radius of planet, max scale for wavelet, mother wavelet shape,
    range of longitude columns to average over, range of latitude columns to average over,
    atmospheric level, number of samples per day

    Outputs: Plot of data, plot of scaleogram"""

    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube.copy()
    # Extract zonal wind data

    x_wind.coord('latitude').coord_system = GeogCS(radius)
    x_wind.coord('longitude').coord_system = GeogCS(radius)
    # Sets planet radius in m for area-weighted average. Default is radius of Proxima Centauri b

    run_length, height = x_wind.shape[0], x_wind.shape[1]
    time = np.arange(0, run_length)*sampling
    # Gives time in (Earth) days if the sampling rate is in samples per day
    heights = np.round(x_wind.coord('level_height').points*1e-03, 0)
    lats, lat_points = x_wind.coord(
        'latitude'), x_wind.coord('latitude').points
    longs, long_points = x_wind.coord(
        'longitude'), x_wind.coord('longitude').points
    # Extract info for labelling plots
    print(x_wind.shape)

    if lats.bounds == None:
        x_wind.coord('latitude').guess_bounds()
    if longs.bounds == None:
        x_wind.coord('longitude').guess_bounds()

    # Set grid box bounds if there are none

    patch = x_wind[:, level, y[0]:y[1], x[0]:x[1]].copy()
    patch_grid = iris.analysis.cartography.area_weights(patch)
    patch_mean = patch.collapsed(
        ['latitude', 'longitude'], iris.analysis.MEAN, weights=patch_grid)
    # Extract data from a small patch of the atmosphere and find the mean

    data = patch_mean.data
    coeffs, freqs = pywt.cwt(data, np.arange(
        1, scales), wavelet, 1./run_length)

    power = (abs(coeffs))**2
    # Find power from coefficients
    period = (1./freqs)*run_length*sampling
    # Gives period in days

    lat_index, long_index = int((y[0]+y[1])/2), int((x[0]+x[1])/2)
    # Indices for central lat and long coordinates, for the plot titles

    plt.plot(np.arange(0, run_length)*sampling,
             patch_mean.data, linestyle='-', color='b')
    plt.title('Zonal Wind at lat=%s, long=%s, h=%s km' %
              (lat_points[lat_index], long_points[long_index], heights[level]))
    plt.xlabel('Time [days]')
    plt.ylabel('Velocity [m s-1]')
    plt.show()

    plt.pcolormesh(time, period, power, cmap=magma)
    plt.title('Scaleogram for lat=%s, long=%s, h=%s km' %
              (lat_points[lat_index], long_points[long_index], heights[level]))
    plt.xlabel('Time [days]')
    plt.ylabel('Period [days]')
    plt.show()
    
def vapour_series(cubes, radius=7160000, level=47, x=(106,110), y=(43,47)):
    
    for cube in cubes:
        if cube.standard_name == 'specific_humidity':
            vapour = cube.copy()
            
    vapour.coord('latitude').coord_system = GeogCS(radius)
    vapour.coord('longitude').coord_system = GeogCS(radius)
    # Sets planet radius in m for area-weighted average. Default is radius of Proxima Centauri b
    
    run_length, height = vapour.shape[0], vapour.shape[1]
    # Gives time in (Earth) days if the sampling rate is in samples per day
    heights = np.round(vapour.coord('level_height').points*1e-03,0)
    lats, lat_points = vapour.coord('latitude'), vapour.coord('latitude').points
    longs, long_points = vapour.coord('longitude'), vapour.coord('longitude').points
    
    if lats.bounds == None:
        vapour.coord('latitude').guess_bounds()
    if longs.bounds == None:
        vapour.coord('longitude').guess_bounds()       
    # Set grid box bounds if there are none
        
    patch = vapour[:,level,y[0]:y[1],x[0]:x[1]].copy()    
    patch_grid = iris.analysis.cartography.area_weights(patch)
    patch_mean = patch.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights=patch_grid)
    data = patch_mean.data
    
    lat_index, long_index = int((y[0]+y[1])/2), int((x[0]+x[1])/2)
    
    avg = np.mean(data)
    smallest = np.min(data)
    biggest = np.max(data)
    print('The average specific humidity is %s kg/kg' %avg)
    print('The min and max are %s and %s' %(smallest, biggest))

    plt.figure(figsize=(10,5))   
    plt.plot(np.arange(0,run_length), data, linestyle='-', color='b')
    plt.title('Specific Humidity at Eastern Terminator, h=%s km' %(heights[level]))    
    # plt.title('Specific humidity at lat=%s, long=%s, h=%s km' %(lat_points[lat_index], long_points[long_index], heights[level]))
    plt.xlabel('Time [days]') 
    plt.ylabel('Water vapour [kg/kg]')
    # plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/laso/epsfigs/vapour_qbo.eps', format='eps')

    plt.show()
    
