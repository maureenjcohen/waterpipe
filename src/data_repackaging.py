#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 11:09:06 2023

@author: Mo Cohen
"""
import iris

dataset = iris.load('/exports/csce/datastore/geos/users/s1144983/um_data/cloudproject/runawaytrap.nc')

tstart = 400
tend = 600

cube_names = ['x_wind', 'y_wind', 'air_potential_temperature', 'air_pressure',
              'surface_temperature','upward_air_velocity',
              'mass_fraction_of_cloud_liquid_water_in_air' ,'mass_fraction_of_cloud_ice_in_air',
              'specific_humidity',
              'surface_net_downward_shortwave_flux'] 


subsetted_list = iris.cube.CubeList()
for cube in dataset:
    if (cube.standard_name in cube_names or cube.long_name in cube_names):
        print(cube.standard_name)
        if len(cube.shape) == 4:
            subsetted_cube = cube[tstart:tend,:,:,:]
            subsetted_list.append(subsetted_cube)
        elif len(cube.shape) == 3:
            subsetted_cube = cube[tstart:tend,:,:]
            subsetted_list.append(subsetted_cube)
        else:
            print('Something went wrong')

print(subsetted_list)

iris.save(subsetted_list,'/exports/csce/datastore/geos/users/s1144983/um_data/cloudproject/warmtrap.nc')
