#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 16:32:45 2023

@author: Mo Cohen
"""
import iris, windspharm
import warnings
from windspharm.iris import VectorWind
import iris.plot as iplt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
from matplotlib.colors import TwoSlopeNorm
from iris.analysis import calculus
from numpy import unravel_index

period = (506, 526)
warnings.filterwarnings('ignore')



def extract_core(winds, time_slice=-1, level=8):
   
    rel_vort = winds.vorticity()   

    rv_min = unravel_index(np.argmax(rel_vort[level,:20,0:72].data, axis=None), rel_vort[level,:20,0:72].shape)
    print(rv_min[0], rv_min[1])
    print(rel_vort[level,rv_min[0],rv_min[1]].data)      
            
    return rv_min


def model_rwave(cubes,startlon=30,start=500,end=600,nlat=90,nlon=144,level=8,
                omega=1.19e-05,g=9.12,radius=5797818,lat=80,meaning=5,save='no'):
    
    
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube[start:end,:,:,:].copy()
            longterm_x_wind = cube.copy()
        if cube.standard_name == 'y_wind':
            y_wind = cube[start:end,:,:,:].copy()
        if cube.standard_name =='air_potential_temperature':
            theta = cube[start:end,:,:,:].copy()
        if cube.standard_name =='air_pressure':
            pressure = cube[start:end,:,:,:].copy()           

    time_axis = np.arange(start, end-meaning+1)
    time_length = np.arange(0, time_axis.shape[0])
    
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
    km_heights = np.round(x_wind.coord('level_height').points*1e-03,2)
    latitudes = x_wind.coord('latitude').points
    longitudes = x_wind.coord('longitude').points
    
#     winds = windspharm.iris.VectorWind(x_wind[:,level,:,:],y_wind[:,level,:,:])
#     # Create a VectorWind data object from the x and y wind cubes
#     rel_vort = winds.vorticity() 
#     rel_vort = rel_vort.data
#     print(rel_vort.shape)
#     rel_vort = np.flip(rel_vort,axis=1)
    
#     core_lons = []
#     for time in range(0,rel_vort.shape[0]):
#         rossby_core = unravel_index(np.argmax(rel_vort[time,70:,12:72], axis=None), rel_vort[time,70:,12:72].shape)
#         print(rossby_core[1]+12, longitudes[rossby_core[1]+12])
#         core_lons.append(rossby_core[1]+12)
    
# #    lat_deg = [np.abs(int(latitudes[item])) for item in core_lats] # Convert input row number to latitude in degrees north
#     lon_deg = [longitudes[item] for item in core_lons]
#     print(lon_deg)
#     lon_deg_meaned = np.convolve(np.array(lon_deg).flatten(),np.ones(meaning),'valid')/meaning
#     lon_diff = np.diff(lon_deg_meaned)
    
    v = y_wind[:,level,lat,0:72].data
    
    lon_deg = []
    for ytime in range(0,v.shape[0]):   
        y_ind = np.where(np.diff(np.sign(v[ytime,:])) == 2.)[0]       
        print(y_ind, longitudes[y_ind])
        
        if len(y_ind) == 1:
            lon_deg.append(longitudes[y_ind])
        else:
            lon_deg.append(0)
    
#    lon_deg = [longitudes[item] for item in core_lons]
    lon_deg_meaned = np.convolve(np.array(lon_deg).flatten(),np.ones(meaning),'valid')/meaning
    lon_diff = np.diff(lon_deg_meaned)
    
    shortterm_zmzw = x_wind[:,level,lat,:].collapsed('longitude',iris.analysis.MEAN)
    st_zmzw = shortterm_zmzw.data
    longterm_zmzw = longterm_x_wind[:,level,lat,:].collapsed(['longitude','time'], iris.analysis.MEAN)
    lt_zmzw = longterm_zmzw.data

    lat_deg = latitudes[lat]
    lat_rad = np.array(lat_deg)*(np.pi/180) # Convert input latitude to radians
    beta = 2*omega*np.cos(lat_rad)/radius # Beta factor    
    circum = 2*np.pi*radius*np.cos(lat_rad) # Circumference in meters at input latitude 
    x_num = 2*np.pi/circum
    
    d_theta = iris.analysis.calculus.differentiate(theta, 'level_height')
    bv_freq = np.mean(np.sqrt(np.abs((g/theta[:,:-1,:,:].data)*d_theta.data)), axis=-1)
    Ld = (bv_freq[:,level,lat]*6800)/(2*omega*np.sin(lat_rad))    

        
    wave_num = beta + (np.array(st_zmzw) - np.array(lt_zmzw))*((1./np.array(Ld))**2)
    wave_denom = x_num**2 + (1./np.array(Ld))**2
    wave_component = wave_num/wave_denom
    
    c_phase = (np.array(st_zmzw) - np.array(lt_zmzw)) - wave_component 
    c_phase_meaned = np.convolve(c_phase.flatten(),np.ones(meaning),'valid')/meaning
    zero_ind = np.where(np.diff(np.sign(c_phase_meaned)))[0]
    
    zeroes = []
    for zc in zero_ind:
        t1 = time_length[zc]
        t2 = time_length[zc+1]
        p1 = c_phase_meaned[zc]
        p2 = c_phase_meaned[zc+1]
        interpolated_zero = t1 + (0-p1)* ((t2-t1)/(p2-p1))
        zeroes.append(interpolated_zero)
        
                            
    distance = np.cumsum(c_phase*60*60*24*(360/circum)) # Convert m/s to m/day = distance travelled in a day
    deriv = lon_diff*(1/(60*60*24))*(circum/360)

    
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('Time [days]')
    ax1.set_ylabel('Velocity [m/s]')
    ax1.plot(time_axis,c_phase_meaned, color='b', label='Phase vel')
#    ax1.plot(time_axis, np.zeros_like(c_phase_meaned), color='b', linestyle='dashed')
#    ax1.plot(time_axis[:-1], deriv, color='g', label='Deriv')
    ax1.plot([item+start for item in zeroes], np.zeros_like(zeroes), 'o', color='r')
    ax1.tick_params(axis='y', labelcolor='b')
    ax1.set_ylim(-30, 30)
    
    ax2 = ax1.twinx()
    ax2.set_ylabel('Longitude [deg E]')
    ax2.plot(time_axis, lon_deg_meaned-85, color='k', label='Longitude')
    ax2.plot(time_axis, np.ones_like(lon_deg_meaned)*0, color='k', linestyle='dashed')
    ax2.tick_params(axis='y', labelcolor='k')
    ax2.set_ylim(-85,85)
    
    plt.title('Rossby wave phase velocity lat %s and gyre longitude' %lat_deg) 
    fig.tight_layout()
    plt.show()
    
    
#     markers_on = [0, 50, 10, 15]
#     plt.plot(time_axis,c_phase, color='b', label='Phase vel')
#     plt.plot([item+500 for item in zeroes], np.zeros_like(zeroes), 'o', color='r')
# #    plt.plot(time_axis, c_phase, 'o', color='r', markevery=markers_on)
# #    plt.plot(time_axis, zmzw, color='r', label='ZMZW')
# #    plt.plot(time_axis, wave_component, color='r', linestyle='dashed', label='Wave comp')
# #    plt.plot(time_axis, np.array(lt_zmzw), color='g', label='Longterm ZMZW')
#     plt.title('Rossby wave phase velocity')
#     plt.xlabel('Time [days]')
#     plt.ylabel('Velocity [m/s]')
#     plt.legend(fontsize='small')
#     plt.show()
    
    plt.plot(time_axis, lon_deg_meaned)
#    plt.plot([item+500 for item in zeroes], [lon_deg[item] for item in zeroes], 'o', color='r')
#    plt.plot(time_axis, lon_deg,'o', color='r', markevery=markers_on)
    plt.title('Path travelled by northeast gyre')
    plt.xlabel('Time [days]')
    plt.ylabel('Longitude [deg E]')
    plt.show()
    
    
    return