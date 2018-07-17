""" Convert indexes to degrees/meters and plot some trajectories
Emma Bent
"""

import os
import numpy as np
import netCDF4

load_path2= '/data/SO12/runs/RUN_BLING_Dec2017/SO12_RUN/DIAGNOSTICS/'
file1     = netCDF4.Dataset(load_path2+'so12_i0_year2006_5day_Salt.nc','r')
lon_min   = 1440 #1950
lon_max   = 3241 #2520
lat_min   = 0 
lat_max   = 1024 #1260 # This is bigger than the YC for the domain but otherwise we can't convert from yround to lat, trick is below in "Really need to remember this !!"

# YC is the latitude vector of the domain of the model output, we need it convert indexes y to latitudes (in degrees)
YC        = file1.variables['lat'][lat_min:lat_max]
# XC is the longitude vector of the domain of the model output, we need it to have lon0 to convert indexes x to longitudes (in degrees)
XC        = file1.variables['lon'][lon_min:lon_max]

# RF and DRF are grid information from the model for the conversion of z to dep (IN METERS)
folder2   = '/home/ebent/Octopus/Octopus-master/grid/'
RF        = np.fromfile(os.path.join(folder2,'RF.data'),'>f4')
DRF       = np.fromfile(os.path.join(folder2,'DRF.data'),'>f4')

res       = 12 # horizontal resolution of the model
lon0      = XC[0]
npts      = 10000 # number of particles

# Load Octopus output
fn 	  = '12_10_4MLD_0006.XYZ.0000000001.0000001801.data'
folder    = '/data/ebent/Octopus/output/12_10_4mld/'
opt       = np.fromfile(os.path.join(folder,fn),'>f4').reshape(-1,3,npts)

print "data has %i records"%(opt.shape[0])
print 'glued data : ' + fn
print 'location of data : ' + folder

x,y,z     = opt[:,0,:],opt[:,1,:],opt[:,2,:] # selects x, y and z in the output

######################## Really need to remember this !! ##########################
# Trick to get all the pos in latitude that go out of Southern Ocean to be equal to northern boudary of SO in latitude YC[-1]=2.5 degrees, which is index 1259

x = np.ma.masked_where(x>1800., x)
x = np.ma.masked_where(x<0., x)

y = np.ma.masked_array(y)
y.mask = x.mask
y = np.ma.masked_where(y>1023., y)
y = np.ma.masked_where(y<0., y)

z = np.ma.masked_array(z)
z.mask = y.mask
z = np.ma.masked_where(z>103., z)
z = np.ma.masked_where(z<0., z)

x.mask = z.mask

'''
for i in range(y.shape[0]):
	for j in range(y.shape[1]):
		if y[i,j]>1259:
			y[i,j]=1259

for i in range(y.shape[0]):
	for j in range(y.shape[1]):
		if z[i,j]>103:
			z[i,j]=103

'''
##############################################################################


# Makes sure z has only positive values but not sure why Isa coded this for z ...
#z[z<0.]   = np.abs(z[z<0.])

dep       = np.zeros((z.shape[0],z.shape[1]), dtype='>f4')
zround    = np.zeros((z.shape[0],z.shape[1]), dtype='int')
yround    = np.zeros((y.shape[0],y.shape[1]), dtype='int')
xround    = np.zeros((x.shape[0],x.shape[1]), dtype='int')

# Rounds up all the z, y and x to integers (zround, yround and xround) and calculates dep IN METERS
for t in range(z.shape[0]):
	for n in range(z.shape[1]):
		xround[t,n]   = np.int_(np.floor(x[t,n])) 
		yround[t,n]   = np.int_(np.floor(y[t,n])) 
		zround[t,n]   = np.int_(np.floor(z[t,n])) 

		tmp           = zround[t,n]
		dep[t,n]      = abs(RF[tmp]-DRF[tmp]*(z[t,n]-tmp))
		
# Conversion to degrees instead of indices
lon       = x/float(res) + lon0

# This selects in YC all the elements of index yround and fills a matrix lat of same size as yround.shape
lat       = YC[yround]




'''
from numpy import *
from pylab import *

fn='DIMES_0004_0033.XYZ.0000000166.0000002161.data' #put a glued XYZ data filename 
npts=35000 #specify particle numbers
opt=fromfile(fn,'>f4').reshape(-1,3,npts)

print "data has %i records"%(opt.shape[0])

#plot some trajectories
x,y=opt[:,0,:10],opt[:,1,:10] #this is in model grid index coordinate, convert to lat-lon using x=x/6.0;y=y/6.0-77.875
plot(x,y,'-')
xlabel('x')
ylabel('y')
show()
'''
