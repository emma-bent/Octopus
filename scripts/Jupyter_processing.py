""" SECOND file containing all the functions need to process an Octopus output (first one is Octopus_functions.py).

Written by Emma Bent
Date : September 13, 2018
"""


#_______________________________Function to process the output____________________________________#

def count_in_zones(location, NPP, npts, xround, yround, zround):
    '''Function to select the index p of particles for each zone.

    Selects the index p of the particles for each zone (gyre, north, west, east, east_west, masked at all time) by looking at the size of np.ma.where with conditions on xround (exit the Ross Sea to the E or W) or xround.mask (particles masked at all time) or south_front (exit to the N).

    Args : 
        NPP (str)     : number from 01 to 12 for the ensemble of particles being processed
        xround (array): the longitude positions in indexes of the particles. shape[0] = time, shape[1] = the number of the particle
        yround (array): the latitude positions in indexes of the particles. shape[0] = time, shape[1] = the number of the particle
        zround (array): the depth positions in indexes of the particles. shape[0] = time, shape[1] = the number of the particle  
    
    Returns : 
        in_RG, to_north, to_west, to_east, to_west_east, are_masked : index of the particles that stay in the Ross Sea/Gyre, that exit to the N, W, E, W then E and that are masked for all of the time steps.
    '''
    import numpy as np
    from func_pickle import pickle_load, pickle_save
    
    southern_front = pickle_load('ACC_southern_front', '/data/ebent/', verbose=False)
    southern_front = np.ma.masked_array(southern_front, mask=False)
    
    path = '/data/ebent/Octopus/output/' + location + '/' 
    yround = pickle_load('NPP' + NPP + '_yround', path, verbose=False)
    xround = pickle_load('NPP' + NPP + '_xround', path, verbose=False)
    
    W_boundary = 12 # 121,04 degrees
    E_boundary = 1560 # 250,04 degrees

    to_east = []
    to_west = []
    to_north = []
    in_RG = []
    to_west_east = []
    are_masked = []
    
    for p in range(npts):
        south_front = southern_front[yround[:,p],xround[:,p]] # select southern_front of one particle traj
        south_front.mask = xround[:,p].mask
        
        north_outside_RG = np.squeeze(np.array(np.ma.where(south_front==0))) # don't forget to use np.ma when dealing with masks
        west_outside_RG = np.squeeze(np.array(np.ma.where(xround[:,p]<=W_boundary))) 
        east_outside_RG = np.squeeze(np.array(np.ma.where(xround[:,p]>=E_boundary)))
        parti_are_masked = np.squeeze(np.array(np.ma.where(xround.mask[:,p]==False)))
        
        if north_outside_RG.size==0. and west_outside_RG.size==0. and east_outside_RG.size==0. and parti_are_masked.size>0.:
            in_RG.append(p)
        if north_outside_RG.size>0. and parti_are_masked.size>0.:
            to_north.append(p)
        if north_outside_RG.size==0. and west_outside_RG.size>0. and east_outside_RG.size==0. and parti_are_masked.size>0.:
            to_west.append(p)
        if north_outside_RG.size==0. and west_outside_RG.size==0. and east_outside_RG.size>0. and parti_are_masked.size>0.:
            to_east.append(p)
        if north_outside_RG.size==0. and west_outside_RG.size>0. and east_outside_RG.size>0. and parti_are_masked.size>0.:
            to_west_east.append(p)
        if parti_are_masked.size==0.:
            are_masked.append(p)
    
    return in_RG, to_north, to_west, to_east, to_west_east, are_masked








def post_processing(casename, NPP, location):
    '''Function to post process an Octopus output.
    
    This function converts indexes to degrees/meters for lon, lat and dep using the original script p_xy.py. It masks the particles that go out of the domain.

    Args :
	casename (str): casename used to save Octopus opt
	NPP (str):      number from 01 to 12 for the ensemble of particles being processed
	location (str): path to the file with the saved opt such as path = '/data/ebent/Octopus/output/' + location + '/'
    
    Returns : 
	lon (array):    the longitude positions in degrees of the particles. shape[0] = time, shape[1] = the number of the particle
	lat (array):    the latitude positions in degrees of the particles. shape[0] = time, shape[1] = the number of the particle
        dep (array):    the depth positions in meters of the particles. shape[0] = time, shape[1] = the number of the particle
        xround (array): the longitude positions in indexes of the particles. shape[0] = time, shape[1] = the number of the particle
        yround (array): the latitude positions in indexes of the particles. shape[0] = time, shape[1] = the number of the particle
        zround (array): the depth positions in indexes of the particles. shape[0] = time, shape[1] = the number of the particle  

    '''

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
    fn 	      = casename + '_00' + NPP + '.XYZ.0000000001.0000001801.data'
    folder    = '/data/ebent/Octopus/output/' + location + '/'
    opt       = np.fromfile(os.path.join(folder,fn),'>f4').reshape(-1,3,npts)
    
    print ''
    print npp
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
    
    ################### Get number of parti that go wrong on z ######################
    
    # Counts particles with negative z values
    a = sorted(np.ma.where(z<0.)[1])
    b = []
    for i in range(len(a)):
    	if a[i]!=a[i-1]:
    		b.append(a[i])
    count_inf0 = len(b)
    print 'Number of parti with z < 0 :', count_inf0
    
    # Counts particles with z values greater than the bottom of model
    a = sorted(np.ma.where(z>103.)[1])
    bb = []
    for i in range(len(a)):
    	if a[i]!=a[i-1]:
    		bb.append(a[i])
    count_inf = len(bb) 
    print 'Number of parti with z > 103 (bottom of ocean) :', count_inf
    
    #################################################################################
    
    z = np.ma.masked_where(z>103., z)
    z = np.ma.masked_where(z<0., z)
    
    x.mask = z.mask
    
    
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
    	
    xround = np.ma.masked_array(xround, mask=False)
    yround = np.ma.masked_array(yround, mask=False)
    zround = np.ma.masked_array(zround, mask=False)
    
    xround.mask = x.mask
    yround.mask = y.mask
    zround.mask = z.mask
    
    # Conversion to degrees instead of indices
    lon       = x/float(res) + lon0
    
    # This selects in YC all the elements of index yround and fills a matrix lat of same size as yround.shape
    lat       = YC[yround]
    
    
    return(lon, lat, dep, xround, yround, zround)    


#_______________________________Function to mask the output____________________________________#

def func_mask(casename, NPP, location, lon_W, lon_E, lat_S, lat_N):
    '''Function to mask Octopus output.
    
    Masks the particles that go on land and particles with wrong initial position (because if particles were initialized on bathymetry, they are randomly shuffled somewere else. This function calls the function post_processing(casename, NPP, location) which is an automated version of p_xy.py

    Args :
	casename (str): casename used to save Octopus opt
	NPP (str):      number from 01 to 12 for the ensemble of particles being masked
	location (str): path to the file with the saved opt such as path = '/data/ebent/Octopus/output/' + location + '/'
        lon_W (float):  longitude west of box of initialisation for the particles as described in init_parti_xyz.py
	lon_E (float):  longitude east of box of initialisation for the particles as described in init_parti_xyz.py
	lat_S (float):  latitude south of box of initialisation for the particles as described in init_parti_xyz.py
	lat_N (float):  latitude north of box of initialisation for the particles as described in init_parti_xyz.py

    Saves made : 
    	LON, LAT, DEP, xround, yround, zround

    '''
    from func_pickle import pickle_load, pickle_save
    #from func_p_xy import post_processing
    from Octopus_functions import post_processing

    import numpy as np
    import h5py
    import netCDF4
    
    load_path2='/data/SO12/runs/RUN_BLING_Dec2017/SO12_RUN/DIAGNOSTICS/'
    load_path3='/data/soccom/GRID_12/'
    
    # Load files
    file1 = netCDF4.Dataset(load_path2+'so12_i0_year2006_5day_Salt.nc','r')
    
    # This is hFacC for the SOUTHERN HEMISPHERE
    file_h = h5py.File(load_path3 + 'grid.mat','r')
    hFacC = file_h.get('hFacC')
    hFacC = np.array(hFacC)
    Xf = file_h.get('XC')
    Xf = np.array(Xf)
    Yf = file_h.get('YC')
    Yf = np.array(Yf)
    
    # On selectionne la bonne taille de hFacC == taille de "bigger domain"
    lon_min   = 1440 
    lon_max   = 3241
    lat_min   = 0 
    lat_max   = 1024
    
    YC        = file1.variables['lat'][lat_min:lat_max] # de -77,99 a -40,05
    XC        = file1.variables['lon'][lon_min:lon_max] # de 120,04 a 270,04
    XC, YC    = np.meshgrid(XC, YC)
    
    hfacc = hFacC[:, lat_min:lat_max, lon_min:lon_max]
    
    
    ################### Extract the output of Octopus ###########################
    
    # Call function post_processing that I wrote    
    lon, lat, dep, xround, yround, zround = post_processing(casename=casename, NPP=NPP, location=location)
    
    LON = np.ma.masked_array(lon, mask=False)
    LAT = np.ma.masked_array(lat, mask=False)
    DEP = np.ma.masked_array(dep, mask=False)
    
    # Put a mask on hfacc==0
    for p in range(LON.shape[1]):
        for t in range(LON.shape[0]):
            if LON.mask[t,p]==True: 
                continue
            if hfacc[zround[t,p], yround[t,p], xround[t,p]]==0.:
                LON.mask[t:,p]=True
    
    # Put a mask on all time steps after the first LON.mask[t,p]==True for each parti
    for p in range(LON.shape[1]):
        for t in range(LON.shape[0]):
            if LON.mask[t,p]==True:
                LON.mask[t:,p]=True
                break
    
    # Mask particles that have an initial pos out of the box I chose (this happens because of the pb of bathymetry in initial pos)
    parti = []
    for p in range(LON.shape[1]):
        if LON[0,p]<lon_W or LON[0,p]>lon_E or LAT[0,p]<lat_S or LAT[0,p]>lat_N:
            LON.mask[:,p]=True
            parti.append(p)
    
    print ''
    print 'Number of particles that are out of the initial box on first time step :', len(parti)
    print ''

    # Make sure all pos have the same mask
    LAT.mask = LON.mask
    DEP.mask = LON.mask
    
    xround.mask = LON.mask
    yround.mask = LON.mask
    zround.mask = LON.mask
    
    # Save the variables
    path = '/data/ebent/Octopus/output/' + location + '/'
    pickle_save('NPP' + NPP + '_DEP', path, DEP, verbose=False)
    pickle_save('NPP' + NPP + '_LAT', path, LAT, verbose=False)
    pickle_save('NPP' + NPP + '_LON', path, LON, verbose=False)
    
    pickle_save('NPP' + NPP + '_zround', path, zround, verbose=False)
    pickle_save('NPP' + NPP + '_yround', path, yround, verbose=False)
    pickle_save('NPP' + NPP + '_xround', path, xround, verbose=False)



#_______________________________Function to count the parti in each zone____________________________________#

def count_zones(location, NPP, npts):
    '''Function to count how many particles of the Octopus output are in each 5 zones of the domain of study.
    
    Function first counts how many particles are "alive" for at least one time step then counts for each zone at each time step how many particles are in. Also counts the particles that "die" which is equivalent of one a particle gets masked.

    Args :
        location (str): path to the file with the saved opt such as path = '/data/ebent/Octopus/output/' + location + '/'
        NPP (str):      number from 01 to 12 for the ensemble of particles being studied
        npts (int):     total number of particles in this ensemble
    
    Loaded variables :
        LON, LAT, DEP, xround, yround, zround, southern_front

    Saves made : 
        Z1, Z2, Z3, Z4, Z5, dead
    '''

    import numpy as np
    from func_pickle import pickle_load, pickle_save
    
    southern_front = pickle_load('ACC_southern_front', '/data/ebent/', verbose=False)
    southern_front = np.ma.masked_array(southern_front, mask=False)
    
    path = '/data/ebent/Octopus/output/' + location + '/'
    DEP = pickle_load('NPP' + NPP + '_DEP', path, verbose=False)
    LAT = pickle_load('NPP' + NPP + '_LAT', path, verbose=False)
    LON = pickle_load('NPP' + NPP + '_LON', path, verbose=False)
    
    zround = pickle_load('NPP' + NPP + '_zround', path, verbose=False)
    yround = pickle_load('NPP' + NPP + '_yround', path, verbose=False)
    xround = pickle_load('NPP' + NPP + '_xround', path, verbose=False)
    
    W_boundary = 12 # 121,04 degrees
    E_boundary = 1560 # 250,04 degrees
    
    # Count the particles that are "alive" during part of the time serie 
    are_alive = []
    for p in range(npts):
        south_front = southern_front[yround[:,p],xround[:,p]] # select southern_front of one particle traj
        south_front.mask = xround[:,p].mask
        
        parti_are_dead = np.squeeze(np.array(np.ma.where(xround.mask[:,p]==False)))
        if parti_are_dead.size>0.:
            are_alive.append(p)
    
    print ''
    print 'Number of particles "alive" :', len(are_alive)
    
    
    # Particle analysis : count the particles in each zone according to time
    Z1 = np.zeros(xround.shape[0])
    Z2 = np.zeros(xround.shape[0])
    Z3 = np.zeros(xround.shape[0])
    Z4 = np.zeros(xround.shape[0])
    Z5 = np.zeros(xround.shape[0])
    dead = np.zeros(xround.shape[0])
    
    for p in are_alive: #range(npts)
        south_front = southern_front[yround[:,p],xround[:,p]] # select southern_front of one particle traj
        south_front.mask = xround[:,p].mask # mask elements of south_front that are irrelevant
        #print ''
        #print 'nb of parti :', p
        for t in range(xround.shape[0]):
            if south_front[t]!=0 and xround[t,p]>W_boundary and xround[t,p]<E_boundary:
                #print t, 'in RG'
                Z1[t]+=1
                #print Z1[t], Z2[t], Z3[t], Z4[t], Z5[t], dead[t]
            
            elif south_front[t]==0:
                #print t, 'north'
                Z2[t:]+=1
                #print Z1[t], Z2[t], Z3[t], Z4[t], Z5[t], dead[t]
                break # after the break, the particle is considered forever in the zone
            elif xround[t,p]<=W_boundary:
                #print t,'out to the west'
                Z3[t:]+=1
                #print Z1[t], Z2[t], Z3[t], Z4[t], Z5[t], dead[t]
                break
            elif xround[t,p]>=E_boundary:
                #print t, 'east'
                Z4[t:]+=1
                #print Z1[t], Z2[t], Z3[t], Z4[t], Z5[t], dead[t]
                break
            elif xround[t,p]<=W_boundary and xround[t,p]>=E_boundary:
                #print t, 'east and west'
                Z5[t:]+=1
                #print Z1[t], Z2[t], Z3[t], Z4[t], Z5[t], dead[t]
                break
                
            elif south_front[t-1]!=0 and xround[t-1,p]>W_boundary and xround[t-1,p]<E_boundary and xround.mask[t,p]==True:
                #print t, 'parti is dead !'
                dead[t:]+=1
                #print Z1[t], Z2[t], Z3[t], Z4[t], Z5[t], dead[t]
                #print Z1[t-1], Z2[t-1], Z3[t-1], Z4[t-1], Z5[t-1], dead[t-1]
                break
            ##print Z1[t], Z2[t], Z3[t], Z4[t], Z5[t], dead[t]
    
    print ''
    print 'Different zones at t = 0 :', Z1[0], Z2[0], Z3[0], Z4[0], Z5[0], dead[0]
    print 'Different zones at t =', xround.shape[0], ':', Z1[-1], Z2[-1], Z3[-1], Z4[-1], Z5[-1], dead[-1] 
    print 'check at t = 0 :', Z1[0]+Z2[0]+Z3[0]+Z4[0]+Z5[0]+dead[0]
    print 'check at t = 44 :', Z1[44]+Z2[44]+Z3[44]+Z4[44]+Z5[44]+dead[44]
    print 'check at t = 100 :', Z1[100]+Z2[100]+Z3[100]+Z4[100]+Z5[100]+dead[100]
    print 'check at t =', xround.shape[0], ':', Z1[-1]+Z2[-1]+Z3[-1]+Z4[-1]+Z5[-1]+dead[-1]
    print '' 
    
    pickle_save('NPP' + NPP + '_Z1', path, Z1, verbose=False)
    pickle_save('NPP' + NPP + '_Z2', path, Z2, verbose=False)
    pickle_save('NPP' + NPP + '_Z3', path, Z3, verbose=False)
    pickle_save('NPP' + NPP + '_Z4', path, Z4, verbose=False)
    pickle_save('NPP' + NPP + '_Z5', path, Z5, verbose=False)
    pickle_save('NPP' + NPP + '_dead', path, dead, verbose=False)


def add_zones(location, NPP, max_time_experiment):
        '''Function used to add zones of several ensembles together.

        Loaded variables :
            NPP01_Z1 or NPP11_dead for example, depends on provided arg NPP

        Saves made : 
            Sum_Zones (dict) : contains sums of Z1, Z2, Z3, Z4, Z5, dead

        Args : 
            location (str): path to the file with the saved opt such as path = '/data/ebent/Octopus/output/' + location + '/'
            NPP (str):      number from 01 to 12 for the ensemble of particles being studied
            max_time_experiment (int) : used to crop each zone from the end to have same sizes
        '''
        from func_pickle import pickle_load, pickle_save
        
        Ensembles = {}   
        Z = ['Z1','Z2','Z3','Z4','Z5','dead']
        npp = NPP
        path = '/data/ebent/Octopus/output/' + location + '/'
        max_time_experiment = max_time_experiment     

        # Start by loading each zone and saving it in a dictionary (Zones) of a dictionary (Ensembles)
        for NPP in npp:
            Zones = {}
            
            for zone in Z:
                Zones[zone] = pickle_load('NPP' + NPP + '_' + zone, path, verbose=False)[max_time_experiment:]
                
            Zones['dead']= Zones['dead'] - Zones['dead'][0]
            Ensembles[NPP] = Zones

        pickle_save('Ensembles', path, Ensembles)

        # Then sum up each different zone together with the other ones of each ensemble and save it in a dictionary (Sum_Zones)
        Sum_Zones = {}
        
        for NPP in npp:
            for zone in Z:
                if NPP == '01':
                    Sum_Zones[zone] = Ensembles[NPP][zone].copy()
                else:
                    Sum_Zones[zone] += Ensembles[NPP][zone]

        pickle_save('Sum_Zones', path, Sum_Zones)





