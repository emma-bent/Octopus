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
                
            elif xround.mask[t,p]==True:
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
    
    pickle_save('NPP' + NPP + '_Z1', path, Z1)
    pickle_save('NPP' + NPP + '_Z2', path, Z2)
    pickle_save('NPP' + NPP + '_Z3', path, Z3)
    pickle_save('NPP' + NPP + '_Z4', path, Z4)
    pickle_save('NPP' + NPP + '_Z5', path, Z5)
    pickle_save('NPP' + NPP + '_dead', path, dead)
