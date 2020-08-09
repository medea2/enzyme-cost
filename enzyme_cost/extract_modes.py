""" extract directions (signs and counts) and list of reversible reactions from thermodynamic sampling 
    remove the two phosphate transfer reactions and remove duplicative modes         """

# goal: rewrite this as a function that takes as an input a directory and a filename and that gives as an output a list of mode objects


import os
import numpy as np
import h5py

import mode


def extract_modes(modes_dir, modes_file, curr_dir):
    
    # go to directory 
    os.chdir(modes_dir)
    
    #import matlab file using h5py
    # open the file
    f = h5py.File(modes_file,'r')
    # get to the 'subfile' directions
    directions = f.get('directions')
    
    # extract the signs and convert to numpy array
    signs = directions.get('signs')
    signs = np.array(signs)
    
    # extract the counts and convert to numpy array
    counts = directions.get('counts')
    counts = [float(el) for el in counts]
    counts = np.array(counts)
    
    # extract the reversible Reactions and convert to string array
    # revRcts = f.get('revRxns')
    # revRcts = np.array(revRcts)
    # print(revRcts)
    
    # how can I read out strings from the matlab file???   -> for now I've put in here manually the reaction names
    # original revRcts -> all
    revRcts = ['R_ATPS4rpp','R_ENO','R_FBA','R_FUM','R_GAPD','R_MDH','R_PGI',
               'R_PGK','R_PGM','R_PIt2r','R_RPE','R_TALA','R_TKT1','R_TKT2',
               'R_TPI','R_PIt2r_H2PO4']
    
    # reduced revRcts -> took away phosphate transporters
    # revRcts = ['R_ATPS4rpp', 'R_ENO', 'R_FBA', 'R_FUM', 'R_GAPD', 'R_MDH', 'R_PGI', 'R_PGK',
    #          'R_PGM','R_RPE', 'R_TALA','R_TKT1','R_TKT2','R_TPI']  

#original: 'ATPS4r','ENO','FBA','FUM','GAPD','MDH','PGI',
#'PGK','PGM','PIt2r','RPE','TALA','TKT1','TKT2',
#'TPI','PIt2r_H2PO4'

    # create a list of mode objects
    modes = list()
    
    #iterate over all modes and create a mode object for each mode with the corresponding values
    for i in range(0, len(counts)):
        newmode = mode.mode(revRcts, signs[i,:], counts[i], 0, 0, 0)    
        modes.append(newmode)
        
    
    #get indices of phosphate reactions to be deleted
    P1 = revRcts.index("R_PIt2r")
    P2 = revRcts.index("R_PIt2r_H2PO4")
    P2 -= 1 #because as soon as P1 is deleted P2 is shifted
    
    # I don't understand why, but this removes the two phosphate transfer reactions from all the revRcts arrays of all modes 
    # I thought it would do that only for the mode 0, but no that's not the case... ???
    modes[0].revRcts.pop(P1)
    modes[0].revRcts.pop(P2)
    
    # go through all modes
    for i in range(0, len(modes)):
        
        # delete the sign of the phosphate transfer reactions in each mode
        modes[i].signs = np.delete(modes[i].signs, P1)
        modes[i].signs = np.delete(modes[i].signs, P2)
         
    #check for now identical modes and merge them 
    for i in range(0, len(modes)):
        for j in range (i+1, len(modes)):
            # check if indices are not longer than the list of modes and check for equality of the modes (signs)
            if( i < len(modes) and  j < len(modes) and np.array_equal(modes[i].signs, modes[j].signs)):
                # add the counts of the two identical modes and remove one of them
                modes[i].count += modes[j].count
                modes.pop(j) 
        
    # calculate the probability of each mode
    sum_counts = 0
    for i in range(0,len(modes)):
        sum_counts += modes[i].count
    
    for i in range(0, len(modes)):
        modes[i].prob = 100*(modes[i].count /sum_counts)
        

    # go back to the original/ enzyme-cost directory
    os.chdir(curr_dir) 

    return modes
        



