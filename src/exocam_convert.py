#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 12:46:56 2023

@author: Mo Cohen
"""
import iris, os

exocamdir = '/exports/csce/datastore/geos/users/s1144983/exocam_data/'

dfiles = os.listdir(exocamdir)
dlist = []
for dfile in dfiles:
    data = iris.load(exocamdir + dfile)
    dlist.append(data)

fileorder = [0,2,1,3]
timesteps = iris.cube.CubeList()
x_wind = iris.cube.CubeList()
y_wind = iris.cube.CubeList()
cl_liq = iris.cube.CubeList()
cl_ice = iris.cube.CubeList()
temp = iris.cube.CubeList()
z_wind = iris.cube.CubeList()
for n in fileorder:
    for cube in dlist[n]:
        if cube.long_name == 'current timestep':
            timesteps.append(cube)
        if cube.long_name == 'Zonal wind':
            x_wind.append(cube)
        if cube.long_name == 'Meridional wind':
            y_wind.append(cube)
        if cube.long_name == 'Vertical velocity (pressure)':
            z_wind.append(cube)
        if cube.long_name == 'Temperature':
            temp.append(cube)
        if cube.long_name == 'Grid box averaged cloud liquid amount':
            cl_liq.append(cube)
        if cube.long_name == 'Grid box averaged cloud ice amount':
            cl_ice.append(cube)

x_wind = x_wind.concatenate()     
y_wind = y_wind.concatenate()
z_wind = z_wind.concatenate()
temp = temp.concatenate()
cl_liq = cl_liq.concatenate()
cl_ice = cl_ice.concatenate()

outlist = iris.cube.CubeList()
outlist.append(x_wind[0])
outlist.append(y_wind[0])
outlist.append(z_wind[0])
outlist.append(temp[0])
outlist.append(cl_liq[0])
outlist.append(cl_ice[0])

print(outlist) 

iris.save(outlist, exocamdir + 'exocam.nc')      
            