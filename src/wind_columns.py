#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 16:16:28 2022

@author: Mo Cohen
"""
import pyvista as pv
import iris, windspharm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
from mpl_toolkits.mplot3d import axes3d
from pyvista import transform_vectors_sph_to_cart, grid_from_sph_coords
import plotly.graph_objects as go
import plotly.io as pio
pio.renderers.default='png'


def wind_columns(cubes,time_slice=-1,radius=7160000,xstride=10,ystride=10,wscale=1,vscale=0.05,vangle=270.0):
    
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube[time_slice,:,:,:].copy()
        if cube.standard_name == 'y_wind':
            y_wind = cube[time_slice,:,:,:].copy()
        if cube.standard_name == 'upward_air_velocity':
            z_wind = cube[time_slice,:,:,:].copy()
            
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
    z_wind = z_wind.regrid(x_wind, iris.analysis.Linear())
    
    verts = [('level_height', x_wind.coord('level_height').points)]
    z_wind = z_wind.interpolate(verts, iris.analysis.Linear())
    
    wind_level = radius*1.2
    
    heights = x_wind.coord('level_height').points + wind_level

    
    lon = x_wind.coord('longitude').points
    lat = 90.0 - x_wind.coord('latitude').points
    
    lons = lon[::xstride]
    lats = lat[::ystride]
    
    inv_axes = [*range(x_wind.ndim - 1, -1, -1)]
    
    u_data = x_wind.data[..., ::ystride, ::xstride].transpose(inv_axes)
    v_data = y_wind.data[..., ::ystride, ::xstride].transpose(inv_axes)
    w_data = z_wind.data[..., ::ystride, ::xstride].transpose(inv_axes)
    
    w_data *= wscale
    
    vectors = np.stack([i.transpose(inv_axes).swapaxes(x_wind.ndim - 2, x_wind.ndim - 1).ravel("C")
                        for i in transform_vectors_sph_to_cart(lons,lats,heights,u_data,-v_data,w_data)],axis=1)
    
    vectors *= radius*vscale
    grid = grid_from_sph_coords(lons,lats,heights)
    grid.point_data['vector3d'] = vectors
    
    arrow = pv.Arrow()
    glyphs = grid.glyph(orient='vector3d',scale='vector3d',geom=arrow,tolerance=0.005)
       
    
    pl = pv.Plotter()
    pl.add_mesh(pv.Sphere(radius=radius),color='blue',opacity=1.0)
    pl.add_mesh(glyphs,cmap='magma')
    pl.camera_position = 'yz'
    pl.camera.azimuth = vangle
    pl.show()
    
    
def vortex(cubes,time_slice=-1):
    
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube[time_slice,:,:,:].copy()
        if cube.standard_name == 'y_wind':
            y_wind = cube[time_slice,:,:,:].copy()
        if cube.standard_name == 'upward_air_velocity':
            z_wind = cube[time_slice,:,:,:].copy()
            
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
    z_wind = z_wind.regrid(x_wind, iris.analysis.Linear())
    
    verts = [('level_height', x_wind.coord('level_height').points)]
    z_wind = z_wind.interpolate(verts, iris.analysis.Linear())
    
    height = x_wind.coord('level_height').points    
    lon = x_wind.coord('longitude').points
    lat = x_wind.coord('latitude').points
    

    vert_u = []
    vert_v = []
    
    for level in range(0,x_wind.shape[0]):

        winds = windspharm.iris.VectorWind(x_wind[level,:,:],y_wind[level,:,:])
        # Create a VectorWind data object from the x and y wind cubes
        uchi,vchi,upsi,vpsi = winds.helmholtz(truncation=21)
        # Calculate the Helmholtz decomposition. Truncation is set to 21 because 
        # this is what Hammond and Lewis 2021 used.
        
        zonal_upsi = upsi.collapsed('longitude', iris.analysis.MEAN)
        zonal_vpsi = vpsi.collapsed('longitude', iris.analysis.MEAN)
        # Calculate zonal means of the x and y components of the rotational component
        eddy_upsi = upsi - zonal_upsi
        eddy_vpsi = vpsi - zonal_vpsi 
        # Subtract zonal means from the original cubes to get the x and y eddy rotational components
        vert_u.append(eddy_upsi.data)
        vert_v.append(eddy_vpsi.data)
        
    vert_u, vert_v = np.array(vert_u), np.array(vert_v)
    fig = go.Figure(data = go.Cone(
        x=lon,
        y=lat,
        z=height,
        u=vert_u,
        v=vert_v,
        w=z_wind.data,
        colorscale='Blues',
        sizemode='scaled',
        sizeref=0.5))
    
    # fig.update_layout(scene=dict(aspectratio=dict(x=1,y=1,z=1),
    #                               camera_eye=dict(x=1.2,y=1.2,z=0.6)))
    fig.show()


def mplvortex(cubes,time_slice=-1):
    
    for cube in cubes:
        if cube.standard_name == 'x_wind':
            x_wind = cube[time_slice,:,:,:].copy()
        if cube.standard_name == 'y_wind':
            y_wind = cube[time_slice,:,:,:].copy()
        if cube.standard_name == 'upward_air_velocity':
            z_wind = cube[time_slice,:,:,:].copy()
            
    y_wind = y_wind.regrid(x_wind, iris.analysis.Linear())
    z_wind = z_wind.regrid(x_wind, iris.analysis.Linear())
    
    verts = [('level_height', x_wind.coord('level_height').points)]
    z_wind = z_wind.interpolate(verts, iris.analysis.Linear())
    
    height = x_wind.coord('level_height').points    
    lon = x_wind.coord('longitude').points
    lat = x_wind.coord('latitude').points
    

    vert_u = []
    vert_v = []
    
    for level in range(0,x_wind.shape[0]):

        winds = windspharm.iris.VectorWind(x_wind[level,:,:],y_wind[level,:,:])
        # Create a VectorWind data object from the x and y wind cubes
        uchi,vchi,upsi,vpsi = winds.helmholtz(truncation=21)
        # Calculate the Helmholtz decomposition. Truncation is set to 21 because 
        # this is what Hammond and Lewis 2021 used.
        
        zonal_upsi = upsi.collapsed('longitude', iris.analysis.MEAN)
        zonal_vpsi = vpsi.collapsed('longitude', iris.analysis.MEAN)
        # Calculate zonal means of the x and y components of the rotational component
        eddy_upsi = upsi - zonal_upsi
        eddy_vpsi = vpsi - zonal_vpsi 
        # Subtract zonal means from the original cubes to get the x and y eddy rotational components
        vert_u.append(eddy_upsi.data)
        vert_v.append(eddy_vpsi.data)
        
    vert_u, vert_v = np.array(vert_u), np.array(vert_v)
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    
    x,y,z = np.meshgrid(lat,height,lon)
    
    q = ax.quiver(x,y,z,vert_u,vert_v,z_wind.data,length=0.1,cmap='magma',lw=2)
    plt.show()