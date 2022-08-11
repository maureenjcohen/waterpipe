#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 15:55:11 2022

@author: Mo Cohen
"""
import glob
import numpy as np
import matplotlib.pyplot as plt

path = '/exports/csce/datastore/geos/users/s1144983/psg_files/trapclose/spectra/'

files = sorted(glob.glob(str(path) + '*.txt'))   

time_axis = np.arange(0,len(files))
print(time_axis)

transit_depths = []
for file_number in range(0,len(files)):
  open_file = open(files[file_number],'r')
  txt_line = open_file.readlines()
  txt_line = txt_line[756]
  datapoints = txt_line.split(' ')
  transit_depth = datapoints[10].strip()
  transit_depths.append(transit_depth)
  open_file.close()
  
transit_floats = []
for string in transit_depths:
  data = float(string)
  transit_floats.append(data)
print(np.array(transit_floats))

fig,ax = plt.subplots(figsize=(10,5))
plt.plot(time_axis[:100],np.array(transit_floats[:100]))
plt.ylabel('Transit depth [ppm]')
plt.xlabel('Time [days]')
plt.title('Time series of 2.7 $\mu$m CO$_2$ feature')


plt.savefig('/exports/csce/datastore/geos/users/s1144983/papers/cloudproject/epsfigs/warmtrap_2p7CO2.eps', format='eps')

plt.show()