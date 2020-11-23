# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 14:54:07 2020

@author: Mo Cohen

Pipeline for post-processing UM output data

Experiment: Sensitivity of water cycle to variation of ___
Model: University of Exeter group tidally locked Proxima Centauri b, UM vn11.7

"""
# Import packages
import os, iris
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import iris.plot as iplot


# def main(): 

brewer_red = mpl_cm.get_cmap('brewer_YlOrRd_09')
brewer_bg = mpl_cm.get_cmap('brewer_PuBu_09')
brewer_reds = mpl_cm.get_cmap('brewer_Reds_09')
#Load colormaps to use in plots

# Block 1: Load UM output PP files into an Iris cube, convert data, and save
# to nc

filenames = os.listdir(directory)
print(filenames)
# Create list of names of files in source directory

allfiles = []
for entry in filenames:
    fullpath = os.path.join(directory, entry)
    allfiles.append(fullpath)
# Add full address to file name and append to a new list
    
cubes = iris.load(allfiles)
# Load files into an iris cube

for cube in cubes:
 if cube.has_lazy_data() == True:
     cube.data
 else:
     pass
 # Convert lazy data to real data
 
# Block 2: Generate and save plots
 
surf_temp = # Enter cube index
iplt.contourf(surf_temp, brewer_red.N, cmap=brewer_red)
ax = plt.gca()
ax.gridlines(draw_labels=True)
plt.title('Mean Surface Temperature [K]', y=1.20)
plt.colorbar(pad=0.1)
plt.show()
# Surface temperature

# Pressure-temperature profiles

spec_humid = # Enter cube index
iplt.contourf(spec_humid, brewer_bg.N, cmap=brewer_bg)
ax = plt.gca()
ax.gridlines(draw_labels=True)
plt.title('Specific Humidity at h= [kg kg-1]', y=1.20)
plt.colorbar(pad=0.1)
plt.show()

iplt.contourf(spec_humid[:,:,0], brewer_bg.N, cmap=brewer_bg)
plt.title('Specific Humidity at Eastern Terminator [kg kg-1]', y=1.05)
plt.ylabel('Height [m]')
plt.xlabel('Latitude')
plt.colorbar(pad=0.1)
plt.show()
# Specific humidity

x_wind = # Enter cube index
y_wind = # Enter cube inex
y = np.arange(-45,45)
x = np.arange(-72,72)
xlabels = np.arange(-180,181,40)
ylabels = np.arange(-90,91,30)

X,Y = np.meshgrid(x,y)
fig = plt.figure(figsize = (12, 7)) 
strm = plt.streamplot(X, Y, x_wind, y_wind, density = 0.5, color=speed, cmap=brewer_reds)
fig.colorbar(strm.lines)
plt.title('Wind speed and direction [m s-1], h=')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.xticks((-72,-52,-32,-12,0,12,32,52,72),('180W','140W','100W','60W','0','60E','100E','140E','180E'))
plt.yticks((-45,-30,-15,0,15,30,45),('90S','60S','30S','0','30N','60N','90N'))    
plt.show() 
# Wind stream plots



# if __name__ == '__main__':
#     print('Source folder:')
#     directory = input()

