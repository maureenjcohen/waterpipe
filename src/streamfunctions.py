#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 13:25:50 2023

@author: Mo Cohen
"""
import iris, aeolus
from aeolus.core import AtmoSim
from aeolus.const import init_const
from aeolus.calc import diag
from aeolus.model import um
from aeolus.coord import replace_z_coord
import matplotlib.pyplot as plt

#trap_const = init_const('trap1e')

drytrap = iris.load('/exports/csce/datastore/geos/users/s1144983/um_data/cloudproject/dry.nc')
trap = iris.load('/exports/csce/datastore/geos/users/s1144983/um_data/cloudproject/trap.nc')


drytrap = [replace_z_coord(cube) for cube in drytrap]

drytrap_obj = AtmoSim(cubes=drytrap, name='Dry TRAP-1e', description='Control Trappist-1e with dry atmosphere', planet='trap1e')

dry_zonal_sf = diag.zonal_mass_streamfunction(drytrap,const=drytrap_obj.const, model=um)

print(dry_zonal_sf)