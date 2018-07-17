import matplotlib.pyplot as plt
import numpy as np
from numpy import cos, pi
import h5py
from scipy.io import loadmat
from mpl_toolkits.basemap import Basemap
import netCDF4
from scipy import interpolate
import os
import pickle
from func_plot_parti import lambert_map_parti, pickle_load


plot_path_1993_2017='/home/ebent/plots/1993_2017/'
plot_path_2006_2011='/home/ebent/plots/2006_2011/'
plot_path_jup='/home/ebent/plots/2006_2011/jup2/'
load_path='/data/mmazloff/AVISO/'
load_path2='/data/SO12/runs/RUN_BLING_Dec2017/SO12_RUN/DIAGNOSTICS/'
load_path3='/data/soccom/GRID_12/'

mean_Salt_surf = pickle_load('mean_Salt_surf', '/data/ebent')
mean_Salt_30 = pickle_load('mean_Salt_30', '/data/ebent')
mean_Salt_100 = pickle_load('mean_Salt_100', '/data/ebent')
mean_Salt_200 = pickle_load('mean_Salt_200', '/data/ebent')
mean_Salt_500 = pickle_load('mean_Salt_500', '/data/ebent')

file1 = netCDF4.Dataset(load_path2+'so12_i0_year2006_5day_Salt.nc','r')

lon_min = 1950
lon_max = 2520
lat_min = 0
lat_max = 541

yc = file1.variables['lat'][lat_min:lat_max]
xc = file1.variables['lon'][lon_min:lon_max]
XC, YC = np.meshgrid(xc,yc)

os.system('python p_xy.py')
#run /home/ebent/Octopus/Octopus-master/scripts/p_xy.py
print lat, lon

lambert_map_parti(np.linspace(32.8,35.3,50), XC, YC, mean_Salt_30, lon, lat, 'Model : Salinity at 30 m', 'Salinity', plt.cm.jet)
