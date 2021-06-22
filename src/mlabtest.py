#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 18 13:24:18 2021

@author: s1144983
"""
import numpy as np
import iris
from mayavi import mlab 
import matplotlib.cm as mpl_cm


x = iris.Constraint('eastward_wind')
y = iris.Constraint('northward_wind')
z = iris.Constraint('upward_air_velocity')
u = iris.load('/exports/csce/datastore/geos/users/s1144983/um_data/ch1_control/qbo.nc', constraints=x)
v = iris.load('/exports/csce/datastore/geos/users/s1144983/um_data/ch1_control/qbo.nc', y)
w = iris.load('/exports/csce/datastore/geos/users/s1144983/um_data/ch1_control/qbo.nc', z)
u = u[0][-1,:,:,:].data
v = v[0][-1,:,1:,:].data
w = w[0][-1,1:,:,:].data


src = mlab.pipeline.vector_field(u,v,w)
# magnitude = mlab.pipeline.extract_vector_norm(src)
# flow = mlab.flow(u, v, w, seed_scale=1, seed_resolution=0.5, integration_direction='both', colormap='jet')
vcp = mlab.pipeline.vector_cut_plane(src, mask_points=50, scale_factor=10., colormap='jet')
# mlab.quiver3d(u,v,w,colormap='jet', mask_points=75)
mlab.show()
