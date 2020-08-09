"""  run once the unconstrained RBA and save the results as a pickle file  """

# python 2/3 compatibility
from __future__ import division, print_function
from rba.xml.targets import TargetGroup, TargetReaction

# global imports
import sys
import os.path

# package imports
import rba
import pickle


def main():
    
    # set input and output directory for bacterial model 
    xml_dir = '/Users/medeafux/Desktop/ETH/Masterarbeit/Bacterial-RBA-models/Escherichia-coli-K12-WT'
    output_dir = xml_dir 
    # set current directory
    curr_dir = '/Users/medeafux/Desktop/ETH/Masterarbeit/enzyme_cost'
 
    # run the unconstrained model once to get the reference reaction.out file to check for directions
    print('Model building from XML files ...')
    model = rba.RbaModel.from_xml(xml_dir) 
    print('Starting iterative RBA resolution (default model)...')
    results = model.solve()
    print('Optimal growth rate is {}.'.format(results.mu_opt))
    results.write(output_dir)
        
    # save the results as a pickle file  
    os.chdir(curr_dir)     
    pickle.dump(results, open( "rba-unconstrained-results.p", "wb" ))
    
if __name__ == '__main__':
    main()
    
    
    

