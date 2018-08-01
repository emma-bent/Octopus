import sys
import netCDF4
import numpy as np

load_path2='/data/SO12/runs/RUN_BLING_Dec2017/SO12_RUN/DIAGNOSTICS/'
file1 = netCDF4.Dataset(load_path2+'so12_i0_year2006_5day_Salt.nc','r')
lon_min = 1440 #1950
lon_max = 3241 #2520
lat_min = 0
lat_max = 1024 #541

# XC is the longitude vector of the domain of the model output, we need it to have lon0 to convert initial longitudes (in degrees) to indexes x
XC = file1.variables['lon'][lon_min:lon_max]
# YC is the latitude vector of the domain of the model output, we need it convert initial latitudes (in degrees) to indexes y
YC = file1.variables['lat'][lat_min:lat_max]

# Initial pos of the particles IN DEGREES in a box for npts = scale*scale = 100
min_lon, max_lon = 180., 230.
min_lat, max_lat = -69., -67.
scale = 100.
range_lon = (max_lon-min_lon)/scale
range_lat = (max_lat-min_lat)/scale
lon = np.arange(min_lon, max_lon, range_lon)
lat = np.arange(min_lat, max_lat, range_lat)
lon,lat = np.meshgrid(lon,lat)

#lon = np.ones(2)*170 use for npts=2
#lat = np.arange(-77, -75) use for npts=2
#lat = np.ones(1)*(-75) use for when npts=1

lon0  = XC[0]
res   = 12 # horizontal resolution of the model

# Conversion to indexes x and y
x     = (lon - lon0)*res

y     = np.zeros((lat.shape[0],lat.shape[1]))
# Prend la valeur absolue de lat - YC et cherche ou cette valeur est la plus petite, ainsi cela repere l'index de YC ou il y a la latitude que l'on cherche et la rempli dans y
# Takes the absolute value of lat - YC and search for the min value, this finds the index of YC where there is the lat that we are looking for and fills it in y
for i in range(lat.shape[0]):
	for j in range(lat.shape[1]):
			y[i,j] = np.array(np.where(np.abs(lat[i,j]-YC)==np.min(np.abs(lat[i,j]-YC))))

# Choose the level of depth for initial position IN INDEX
k     = 9

# This makes sure that the z values have the same dimension as x and y, no conversion from meters to indexes here
z     = np.ones_like(x)*k

# Number of particles 
n  = x.size

# Create a variable containing x, y and z values
xyz=np.r_[x.reshape(1,-1),y.reshape(1,-1),z.reshape(1,-1)]

# Save the file
path2run = '../run_work/'
xyz.astype('>f8').tofile(path2run + 'parti_init.bin')


'''
def test():
    center=r_[1500,107,25].reshape(-1,3)
    delta=r_[10,10,0.1].reshape(-1,3)
    xyz = random.random((npts,3)) - 0.5
    xyz = xyz*2*delta+center
    fn = 'particle_initial_xyz.bin'
    xyz.T.astype('>f8').tofile(fn)
    print "write particle initial positions to file "+fn
    return


if __name__=='__main__':
    test()
    '''
