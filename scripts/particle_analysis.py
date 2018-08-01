def routine(casename, NPP, location, npts, lon_W, lon_E, lat_S, lat_N):

	'''This function is a routine that deals with the Octopus opt for 1 NPP.

	Calls three other functions : to process the output, to mask it, to count how many particles are in each zone.
	
	Args :

            casename (str): casename used to save Octopus opt
            NPP (str):      number from 01 to 12 for the ensemble of particles being processed
            location (str): path to the file with the saved opt such as path = '/data/ebent/Octopus/output/' + location + '/'
            npts (int):     total number of particles in the ensemble (NPP)
            lon_W (float):  longitude west of box of initialisation for the particles as described in init_parti_xyz.py
            lon_E (float):  longitude east of box of initialisation for the particles as described in init_parti_xyz.py
            lat_S (float):  latitude south of box of initialisation for the particles as described in init_parti_xyz.py
            lat_N (float):  latitude north of box of initialisation for the particles as described in init_parti_xyz.py

	'''
        import warnings
	warnings.filterwarnings('ignore')

        #from func_p_xy import post_processing
        #from func_mask_post_process import func_mask
        #from func_count_zones import count_zones
        from func_pickle import pickle_load, pickle_save
        from Octopus_functions import func_mask, count_zones

        #post_processing(casename, NPP, location)
        func_mask(casename, NPP, location, lon_W, lon_E, lat_S, lat_N) # this calls post_processing
        count_zones(location, NPP, npts)



def loop_routine(casename, NPP, location, npts, lon_W, lon_E, lat_S, lat_N):

        '''This function is a routine that deals with the Octopus opt, loops through several NPPs.

         Calls three other functions : to process the output, to mask it, to count how many particles are in each zone.
         
         Args :

             casename (str): casename used to save Octopus opt
             NPP (str):      number from 01 to 12 for the ensemble of particles being processed
             location (str): path to the file with the saved opt such as path = '/data/ebent/Octopus/output/' + location + '/'
             npts (int):     total number of particles in the ensemble (NPP)
             lon_W (float):  longitude west of box of initialisation for the particles as described in init_parti_xyz.py
             lon_E (float):  longitude east of box of initialisation for the particles as described in init_parti_xyz.py
             lat_S (float):  latitude south of box of initialisation for the particles as described in init_parti_xyz.py
             lat_N (float):  latitude north of box of initialisation for the particles as described in init_parti_xyz.py

         '''
        from particle_analysis import routine
        
        NPP = NPP
        
        for npp in NPP:
                print npp 
                routine(casename=casename, NPP=npp, location=location, npts=npts, lon_W=lon_W, lon_E=lon_E, lat_S=lat_S, lat_N=lat_N)
