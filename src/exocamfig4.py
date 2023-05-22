#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 11:12:46 2023

@author: Mo Cohen
"""

import iris
import warnings
import matplotlib.pyplot as plt
import iris.coords
import numpy as np
from numpy import unravel_index
import windspharm
from iris.analysis import calculus


warnings.filterwarnings('ignore')


def fig4a(datalist, start=400, end=600, nlat=90, nlon=144, level=8,
                   omega=1.19e-05, g=9.12, radius=5797818, lat=80, meaning=3,
                   savedir='/exports/csce/datastore/geos/users/s1144983/papers/cloudproject/epsfigs_v2/',
                   save=False):
    """ This function calculates the Rossby wave phase speed (including 
    zonal wind) over time. At certain latitudes, the phase velocity alternates 
    between positive and negative. This is the region where the cyclonic 
    structure appears to travel back and forth periodically. 
    
    Inputs:
        datalist - List of Iris CubeLists. Should include Control and Dry
        TRAP-1e simulations, in that order.
        start - Start of period to be plotted (int)
        end - End of period to be plotted (int)
        nlat - Number of latitudes (int)
        nlon - Number of longitudes (int)
        level - Level being plotted (int)
        omega - Rotation rate in rads/sec. Default TRAPPIST-1e (float)
        g - Gravitational constant. Default TRAPPIST-1e (float)
        lat - Latitude where gyre is being tracked. Default 71N (row 80) (int)
        meaning - Meaning period for rolling mean. Default 5 days. (int)
        savedir - Directory to save output plot (str)
        save - Save plot or no (Boolean)
        
    Outputs:
        Plot of Rossby wave phase velocities at 71N on Control and Dry 
        TRAPPIST-1e showing oscillation across 0 for the former only.
        
        Plot of northeast gyre location (in longitude east) and Rossby wave
        phase speed for Control TRAP-1e showing matching periods.
        
        """
    time_axis = np.arange(start, end-meaning+1)
    # Create numbered time axis
    
    cphase_list = []
    cphase_hlist = []
    names = ['Control TRAP-1e', 'Dry TRAP-1e']
    # Names of the two simulations we are comparing

    for cubes in datalist:
        # Take CubeList for each simulation separately
        for cube in cubes:
            if cube.standard_name == 'x_wind':
                x_wind = cube[start:end, :, :, :].copy()
                longterm_x_wind = cube.copy()
            if cube.standard_name == 'y_wind':
                y_wind = cube[start:end, :, :, :].copy()
            if cube.standard_name == 'air_potential_temperature':
                theta = cube[start:end, :, :, :].copy()
            if cube.standard_name == 'air_pressure':
                pressure = cube[start:end, :, :, :].copy()
        # Extract datacubes we need for the calculation
        y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
        # Regrid y-wind onto coordinates of x-wind cube
        km_heights = np.round(x_wind.coord('level_height').points*1e-03, 2)
        # Extract heights in km for labelling plots
        latitudes = x_wind.coord('latitude').points
        longitudes = x_wind.coord('longitude').points
        # Extract latitudes and longitudes

        lat_deg = int(latitudes[lat])
        # Convert input row number to latitude in degrees north
        print('Calculating at ' + str(lat_deg))


        shortterm_zmzw = x_wind[:, level, lat, :].collapsed('longitude',
                                                      iris.analysis.MEAN)
        st_zmzw = np.array(shortterm_zmzw.data)
        longterm_zmzw = longterm_x_wind[:, level, lat, :].collapsed('longitude',
                                                            iris.analysis.MEAN)
        lt_zmzw = np.array(np.mean(longterm_zmzw.data))

        lat_rad = lat_deg*(np.pi/180)
        # Convert input latitude to radians
        beta = 2*omega*np.cos(lat_rad)/radius
        # Beta factor
        circum = 2*np.pi*radius*np.cos(lat_rad)
        # Circumference in meters at input latitude
        x_num = 2*np.pi/circum
        # Zonal wavenumber in units of m^-1 at input latitde

        d_theta = iris.analysis.calculus.differentiate(theta, 'level_height')
        # Change in potential temperature with height
        bv_freq = np.mean(
            np.sqrt(np.abs((g/theta[:, :-1, :, :].data)*d_theta.data)), axis=-1)
        # Zonal mean Brunt-Vaisala frequency
        Ld = (bv_freq[:, level, lat]*6800)/(2*omega*np.sin(lat_rad))
        # Calculate Rossby radius of deformation using the BV frequency
        # for the height and latitude we are plotting
        # Scale height is fixed at 6800 for all simulations

        c_phase = (st_zmzw - lt_zmzw) - ((beta +
                                          ((st_zmzw - lt_zmzw)/(Ld**2)))/(x_num**2 + (1/Ld)**2))
        c_phase_meaned = np.convolve(np.array(c_phase),
                                     np.ones(meaning), 'valid')/meaning
        # Phase velocity as per Vallis 6.65
        cphase_list.append(c_phase_meaned)
        # Append to list
        c_phase_h = ((st_zmzw - lt_zmzw) - (beta/(x_num**2)))
        # Phase velocity not accounting for vertical effects
        c_phase_h_meaned = np.convolve(np.array(c_phase_h),
                                     np.ones(meaning), 'valid')/meaning
        cphase_hlist.append(c_phase_h_meaned)
        # Append to its own list

    # Figure 4a):

    fig, ax = plt.subplots(figsize=(6,4))
    ax.plot(time_axis, cphase_list[0], color='b', label=names[0])
    # Control TRAP-1e
    ax.plot(time_axis, cphase_list[1], color='r', label=names[1])
    # Dry TRAP-1e
    
    # plt.plot(cphase_hlist[0], color='r',
    #          linestyle='dashed', label=names[0] + ', no kd')
    # # Control TRAP-1e
    # plt.plot(cphase_hlist[1], color='b',
    #          linestyle='dashed', label=names[1] + ', no kd')
    # # Dry TRAP-1e

    ax.set_title('Rossby wave phase velocity at %s N' % (lat_deg), fontsize=14)
    ax.set_xlabel('Time [days]', fontsize=14)
    ax.set_ylabel('Velocity [m/s]', fontsize=14)
    ax.set_ylim(-20, 20)
    plt.legend()

    if save == True:
        plt.savefig(savedir + 'cphase_%s_to_%s.eps' %
                    (start, end), format='eps')
    else:
        pass
    plt.show()


def fig4b(cubes, start=500, end=600, nlat=90, nlon=144, level=8,
                   omega=1.19e-05, g=9.12, radius=5797818, lat=80, meaning=3,
                   savedir='/exports/csce/datastore/geos/users/s1144983/papers/cloudproject/epsfigs_v2/',
                   save=False):
    
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube[start:end, :, :, :].copy()
            longterm_x_wind = cube.copy()
        if cube.standard_name == 'y_wind':
            y_wind = cube[start:end, :, :, :].copy()
        if cube.standard_name == 'air_potential_temperature':
            theta = cube[start:end, :, :, :].copy()
        if cube.standard_name == 'air_pressure':
            pressure = cube[start:end, :, :, :].copy()
    # Extract datacubes we need for the calculation
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
    # Regrid y-wind onto coordinates of x-wind cube
    km_heights = np.round(x_wind.coord('level_height').points*1e-03, 2)
    # Extract heights in km for labelling plots
    latitudes = x_wind.coord('latitude').points
    longitudes = x_wind.coord('longitude').points
    # Extract latitudes and longitudes

    lat_deg = int(latitudes[lat])
    # Convert input row number to latitude in degrees north
    print('Calculating at ' + str(lat_deg))
    
    ## Code commented out below calculates the global area-weighted mean
    ## of the zonal wind. It can be used as U in the calculation instead of the
    ## latitude-specific zonal mean wind. This leads to a smaller offset
    ## between the two curves, but a poorer correlation in the amplitudes.

    # short_lats = x_wind.coord('latitude')
    # short_lons = x_wind.coord('longitude')

    # if short_lats.bounds == None:
    #     x_wind.coord('latitude').guess_bounds()
    # if short_lons.bounds == None:
    #     x_wind.coord('longitude').guess_bounds()

    # long_lats = longterm_x_wind.coord('latitude')
    # long_lons = longterm_x_wind.coord('longitude')

    # if long_lats.bounds == None:
    #     longterm_x_wind.coord('latitude').guess_bounds()
    # if long_lons.bounds == None:
    #     longterm_x_wind.coord('longitude').guess_bounds()

    # shortgrid = iris.analysis.cartography.area_weights(x_wind)
    # longgrid = iris.analysis.cartography.area_weights(longterm_x_wind)
    # shortterm_zmzw = x_wind[:, :, :, :].collapsed(['longitude', 'latitude'],
    #                                 iris.analysis.MEAN, weights=shortgrid)
    # st_zmzw = np.array(shortterm_zmzw[:, level].data)
    # longterm_zmzw = longterm_x_wind[:, :, :, :].collapsed(['longitude',
    #                         'latitude'], iris.analysis.MEAN, weights=longgrid)
    # lt_zmzw = np.array(np.mean(longterm_zmzw[:, level].data))
    
    shortterm_zmzw = x_wind[:,level,lat,:].collapsed('longitude', 
                                                     iris.analysis.MEAN)
    st_zmzw = np.array(shortterm_zmzw.data)
    longterm_zmzw = longterm_x_wind[:,level,lat,:].collapsed('longitude',
                                                    iris.analysis.MEAN)
    lt_zmzw = np.array(np.mean(longterm_zmzw.data))

    lat_rad = lat_deg*(np.pi/180)
    # Convert input latitude to radians
    beta = 2*omega*np.cos(lat_rad)/radius
    # Beta factor
    circum = 2*np.pi*radius*np.cos(lat_rad)
    # Circumference in meters at input latitude
    x_num = 2*np.pi/circum
    # Zonal wavenumber in units of m^-1 at input latitde

    d_theta = iris.analysis.calculus.differentiate(theta, 'level_height')
    # Change in potential temperature with height
    bv_freq = np.mean(
        np.sqrt(np.abs((g/theta[:, :-1, :, :].data)*d_theta.data)), axis=-1)
    # Zonal mean Brunt-Vaisala frequency
    Ld = (bv_freq[:, level, lat]*6800)/(2*omega*np.sin(lat_rad))
    # Calculate Rossby radius of deformation using the BV frequency
    # for the height and latitude we are plotting
    # Scale height is fixed at 6800 for all simulations

    c_phase = (st_zmzw - lt_zmzw) - ((beta +
                    ((st_zmzw - lt_zmzw)/(Ld**2)))/(x_num**2 + (1/Ld)**2))
    # Phase velocity as per Vallis 6.65
       

    time_axis = np.arange(start, end-meaning+1)
    # Create numbered time axis
    time_length = np.arange(0, time_axis.shape[0])
    # Just a list from 0 to the time axis length (used in loop below)
    c_phase_meaned = np.convolve(np.array(c_phase),
                                 np.ones(meaning), 'valid')/meaning

    v = y_wind[:, level, lat, 0:72].data

    lon_deg = []
    lon_ind = []
    for ytime in range(0, v.shape[0]):
        y_ind = np.where(np.diff(np.sign(v[ytime, :])) == 2.)[0]
        print(y_ind, longitudes[y_ind])

        if len(y_ind) == 1:
            lon_deg.append(longitudes[y_ind])
            lon_ind.append(y_ind)
        else:
            lon_deg.append(0)
            lon_ind.append(0)

    lon_deg_meaned = np.convolve(np.array(lon_deg).flatten(),
                                 np.ones(meaning), 'valid')/meaning

    zero_ind = np.where(np.diff(np.sign(c_phase_meaned)))[0]
    # Find indices where the sign of the phase velocity changes

    zeroes = []
    for zc in zero_ind:
        t1 = time_length[zc]
        t2 = time_length[zc+1]
        p1 = c_phase_meaned[zc]
        p2 = c_phase_meaned[zc+1]
        interpolated_zero = t1 + (0-p1) * ((t2-t1)/(p2-p1))
        zeroes.append(interpolated_zero)
        # This code block interpolates points where the
        # phase velocity is 0

    fig, ax1 = plt.subplots(figsize=(6,4))
    # Plot with gyre longitude location and phase velocity on the two y-axes
    ax1.set_xlabel('Time [days]', fontsize=14)
    ax1.set_ylabel('Velocity [m/s]', fontsize=14)
    ax1.plot(time_axis, c_phase_meaned, color='b', label='Phase vel')
    ax1.tick_params(axis='y', labelcolor='b')
    ax1.set_ylim(-20, 20)

    ax2 = ax1.twinx()
    ax2.set_ylabel('Longitude [deg E]', fontsize=14)
    ax2.plot(time_axis, lon_deg_meaned, color='k', label='Longitude')
    ax2.plot(time_axis, np.ones_like(lon_deg_meaned)
             * 85, color='k', linestyle='dashed')
    ax2.tick_params(axis='y', labelcolor='k')
    ax2.set_ylim(0, 170)

    plt.title('Phase velocity at %sN and gyre longitude' %lat_deg, fontsize=14)
    fig.tight_layout()

    if save == True:
        plt.savefig(savedir + 'gyre_lon_cphase_%s_to_%s.eps' % (start, end),
                    format='eps')
    else:
        pass
    plt.show()

    # Li'l extra

    plt.plot(time_axis, lon_deg_meaned)
    plt.title('Path travelled by northeast gyre')
    plt.xlabel('Time [days]')
    plt.ylabel('Longitude [deg E]')
    plt.show()
