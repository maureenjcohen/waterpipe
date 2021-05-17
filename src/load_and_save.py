#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 14:12:21 2021

@author:Mo Cohen

"""

import os, iris

directory = '/exports/csce/datastore/geos/users/s1144983/um_data/ch1_control/raw/'

names = os.listdir(directory)
# Create list of names of files in source directory

allfiles = []
for entry in names:
    fullpath = os.path.join(directory, entry)
    allfiles.append(fullpath)
    
cubes = iris.load(allfiles)
print(cubes)

for cube in cubes:
    filename = cube.name()
    if cube.has_lazy_data() == True:
        cube.data
    else:
        pass

    iris.save(cube, directory+filename)
    
    print('Completed save of ' + filename)