def func_mask(casename, NPP, location, lon_W, lon_E, lat_S, lat_N):
    '''Function mask Octopus output.
    
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
    from func_p_xy import post_processing
    
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
                #print t,p
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
    
    # Make sure all pos have the same mask
    LAT.mask = LON.mask
    DEP.mask = LON.mask
    
    xround.mask = LON.mask
    yround.mask = LON.mask
    zround.mask = LON.mask
    
    # Save the variables
    path = '/data/ebent/Octopus/output/' + location + '/'
    pickle_save('NPP' + NPP + '_DEP', path, DEP)
    pickle_save('NPP' + NPP + '_LAT', path, LAT)
    pickle_save('NPP' + NPP + '_LON', path, LON)
    
    pickle_save('NPP' + NPP + '_zround', path, zround)
    pickle_save('NPP' + NPP + '_yround', path, yround)
    pickle_save('NPP' + NPP + '_xround', path, xround)
    
    
    
    
