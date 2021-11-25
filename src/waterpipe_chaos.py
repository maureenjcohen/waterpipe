#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 13:30:00 2021

@author: Mo Cohen
"""

import iris
import matplotlib.pyplot as plt
import numpy as np


def plot_chaos(cube1, cube2, cube3, long=0, lat=0, level=35):
    
    x_data = cube1[:, level, long, lat].data
    x_label = cube1.standard_name
    
    y_data = cube2[:, level, long, lat].data
    y_label = cube2.standard_name
    
    z_data = cube3[:, level, long, lat].data
    z_label = cube3.standard_name
    
    ax = plt.figure().add_subplot(projection='3d')
    
    ax.plot(x_data, y_data, z_data, lw=0.5)
    ax.set_xlabel('%s' %x_label)
    ax.set_ylabel('%s' %y_label)
    ax.set_zlabel('%s' %z_label)
    
    ax.view_init(-140,60)
    plt.show()