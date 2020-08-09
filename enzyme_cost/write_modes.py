""" writing the modes and results into a file """

# currently -> test writing to files
import numpy as np
import csv
import os.path
from extract_modes import extract_modes

def write_modes(modes):
    
    # output results into a file         
    with open('overview-modes.csv', 'w', newline='') as csvfile:
        fieldnames =['mode','counts','prob','AA cost/gr','AA cost','growth rate']
        writer = csv.DictWriter(csvfile, fieldnames = fieldnames)
        writer.writeheader()
        
        for i in range(0,len(modes)):
            if modes[i].opt_GR != 0:
                cost_per_gr = modes[i].cost/modes[i].opt_GR
            else:
                cost_per_gr = 'undef'
                
            writer.writerow({'mode':str(i), 'counts': str(modes[i].count), 
                             'prob': str(np.round(modes[i].prob,2)), 'AA cost/gr': str(cost_per_gr),
                             'AA cost': str(modes[i].cost), 'growth rate': str(modes[i].opt_GR)})
            
if __name__ == '__main__':
    main()

