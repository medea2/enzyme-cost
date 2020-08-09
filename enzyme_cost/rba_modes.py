"""  run RBA for a set of modes and compare their costs and growth rates   """

# python 2/3 compatibility
from __future__ import division, print_function
from rba.xml.targets import TargetGroup, TargetReaction

# global imports
import sys
import os.path
import csv

# package imports
import rba
from cost_from_AAfluxes import parse_enzyme_cost
from extract_modes import extract_modes
from write_modes import write_modes
#from unconstrained import run_unc
from mode import mode

def main():
    
    # set input and output directory for bacterial model 
    xml_dir = '/Users/medeafux/Desktop/ETH/Masterarbeit/Bacterial-RBA-models/Escherichia-coli-K12-WT'
    output_dir = xml_dir
    # set directory and file of modes origin (should in the end be used for the extract_modes function)
    modes_dir = '/Users/medeafux/Desktop/ETH/Masterarbeit/code-Mattia/hyperflux/data/results/core_example/glc/samples'
    modes_file = 'state_e_coli_gerosa2015_glc.mat'   
    # set current directory
    curr_dir = '/Users/medeafux/Desktop/ETH/Masterarbeit/enzyme_cost'

    # import modes, delete phosphate transfer and remove duplicate modes
    modes = extract_modes(modes_dir, modes_file, curr_dir) 
    os.chdir(curr_dir)
    
    plot()
    
    # run the unconstrained model
    print('Model building from XML files ...')
    model = rba.RbaModel.from_xml(xml_dir) 
    print('Starting iterative RBA resolution (default model)...')
    results = model.solve()
    print('Optimal growth rate is {}.'.format(results.mu_opt))
    results.write(output_dir)
          
    gr = results.mu_opt
    AA_cost = parse_enzyme_cost(output_dir,curr_dir)
    AA_cost_per_gr = AA_cost/gr
    
    mode_unc = mode(0, 0, 'unc', AA_cost, gr, 0)
    modes.append(mode_unc)
    
    ######################################
    # go through modes and run for each mode: 
    # constrain targets, run constraint RBA, save growth rate, calculate cost
    # loop over the modes
    for i in range(0, len(modes)-1):
            
        # create a new TargetGroup object = group of target fluxes
        target_group = TargetGroup('mode_targets')
       
        # loop over the reactions of each mode (are actually the same for all of them)
        for j in range(0, len(modes[i].revRcts)):
                
            # create a new TargetReaction = Association of a target value with a chemical reaction (use id of wanted reaction)
            target_reaction = TargetReaction(modes[i].revRcts[j])
            # constrain the upper or lower bound according to the sign of the rev reaction inside the mode
            if modes[i].signs[j] == 1 :
                target_reaction.lower_bound = 'zero'
            elif modes[i].signs[j] == -1 :
                target_reaction.upper_bound = 'zero' 
            # add the target flux value of the target reaction to the list of reaction fluxes in the target group object
            target_group.reaction_fluxes.append(target_reaction)

        # add the target reaction flux to the model (why is this step needed?)
        model.targets.target_groups.append(target_group)
        results = model.solve()        
        print('Optimal growth rate is {}.'.format(results.mu_opt))
        results.write(output_dir)
         
        #save optimal growth rate for each mode
        modes[i].opt_GR = results.mu_opt
        # calculate cost and save cost 
        modes[i].cost = parse_enzyme_cost(output_dir,curr_dir) 
        
        #remove target group from the model = reset model
        model.targets.target_groups.remove(target_group)
    
    # write results into a file
    write_modes(modes)

    
if __name__ == '__main__':
    main()
    
    
    

