"""  get all reversible reactions from the model and constrain their fluxes and run RBA 
     ->  to check with which reactions the direction change has the most influence """

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
# from extract_modes import extract_modes
from write_modes import write_modes
#from unconstrained import run_unc
from mode_plots import plot
from mode import mode

from rev_reaction import rev_reaction

import pickle

import cobra
from shutil import copyfile
from docplex.mp.utils import PublishResultAsDf
from test.support import LARGEST

def main():
    
    # set input and output directory for bacterial model 
    xml_dir = '/Users/medeafux/Desktop/ETH/Masterarbeit/Bacterial-RBA-models/Escherichia-coli-K12-WT'
    output_dir = xml_dir
    # set current directory
    curr_dir = '/Users/medeafux/Desktop/ETH/Masterarbeit/enzyme_cost'
    
    ## get a list of all reversible reactions of the model  
    # read the model
    model = cobra.io.read_sbml_model(xml_dir + "/data/iJO1366.xml")
    r_list = [r for r in model.reactions if 'c' in r.compartments and r.reversibility and not r.boundary]

    # Do FVA to check for real reversibility
    variability = cobra.flux_analysis.variability.flux_variability_analysis(model, reaction_list=r_list, fraction_of_optimum=0)
    bool_series = (variability.minimum < 0) & (variability.maximum > 0)
    
    # get out all reversible reactions and save their names as a list
    true_reversible = variability[bool_series]
    reactionNames = true_reversible.index.values
    listOfReactionNames = list(reactionNames) # to do: change to dictionary -> can also save original flux
 
    # run the unconstrained model once to get the reference reaction.out file to check for directions
    print('Model building from XML files ...')
    model = rba.RbaModel.from_xml(xml_dir) 
    print('Starting iterative RBA resolution (default model)...')
    results = model.solve()
    print('Optimal growth rate is {}.'.format(results.mu_opt))
    results.write(output_dir)
 
    # unpickle the unconstrained RBA results
    #with open('rba-unconstrained-results.p', 'rb') as pickle_file:
    #    results = pickle.load(pickle_file)
 
    # save the unconstrained fluxes we get into a dictionary, to be used later
    unconstrainedFluxes = results.reaction_fluxes()
    # to check if there's a difference to the not pickled one, seems not to be the case  ?????
    # pickledFluxes = unc_results.reaction_fluxes() 
    
    # get growth rate and calculate overall enzyme cost
    unc_gr = results.mu_opt
    unc_AA_cost = parse_enzyme_cost(output_dir, curr_dir)
    unc_AA_cost_per_gr = unc_AA_cost/unc_gr

    # create a dictionary to store all the reversible reactions
    rev_reactions = {} 
    # save the values for the unconstrained model as the 'unconstrained' reversible reaction 
    unconstrained_reaction = rev_reaction("unconstrained", 0, unc_AA_cost, 0.0, unc_gr, unc_AA_cost_per_gr, 1, 'undef')
    rev_reactions[unconstrained_reaction.reaction_name] = unconstrained_reaction
         
    # go through all reversible reactions
    for reaction in listOfReactionNames:
        
        # make a dictionary of reactions that belong to one reaction, e.g. incl. duplicates
        rct_all = {}
        # calculate original flux of main reaction and add it to the dictionary
        originalFlux = unconstrainedFluxes.get("R_" + reaction)
        rct_all['R_' + reaction] = originalFlux # add the reaction itself
        
        # go through all fluxes and find the ones with the name or duplicates and add them to the dictionary
        for key in unconstrainedFluxes: 
            rct_string =   'R' + '_' + reaction + '_'
            
            if (not isinstance(key, str)) or (not isinstance(rct_string, str)):
                print("something")
            
            if key.startswith(rct_string):
                # save found reaction incl. flux to the dictionary
                rct_all[key] = unconstrainedFluxes.get(key)
            
            # set a largest flux in case there will be a TypeError
            largest_flux = originalFlux
        
        try: 
            # sort the reactions according to their absolute flux, take largest 
            sorted_rct_all = sorted(rct_all, key=lambda dict_key: abs(rct_all[dict_key]), reverse = True)
            largest_flux = rct_all.get(sorted_rct_all[0])
        except TypeError:
            print("A flux was of Type 'None'.")
        
        if largest_flux == 0:
            # if no flux -> reversing reaction will not change anything
            new_rev_reaction = rev_reaction("R_" + reaction, largest_flux, unc_AA_cost, 0.0, unc_gr, unc_AA_cost_per_gr, len(rct_all), '0')
            # add new reversible reaction to the dictionary
            rev_reactions[new_rev_reaction.reaction_name] = new_rev_reaction
        elif largest_flux == None: # reaction couldn't be found 
            new_rev_reaction = rev_reaction("R_" + reaction, largest_flux, None , 0.0, None, None, len(rct_all), None)
            rev_reactions[new_rev_reaction.reaction_name] = new_rev_reaction
        else:
            # create a new reversible reaction, yet unknown parameters initiated to 0
            new_rev_reaction = rev_reaction("R_" + reaction, largest_flux, 0, 0, 0, 0, len(rct_all), None)
            # add new reversible reaction to the dictionary
            rev_reactions[new_rev_reaction.reaction_name] = new_rev_reaction
        
            # save the direction of the largest flux to adapt all the reactions 
            if largest_flux > 0: 
                positive = True
                # save the original direction of the flux
                new_rev_reaction.orig_direction = 'forward'
            else: 
                positive = False
                # save the original direction of the flux
                new_rev_reaction.orig_direction = 'backward'
            # change the reaction fluxes into the opposite direction
            # create a new TargetGroup object = group of target fluxes
            target_group = TargetGroup('reaction_targets')
           
            # go through all duplicate reactions
            for rctname in sorted_rct_all: 
                # create a new TargetReaction = Association of a target value with a chemical reaction (use id of wanted reaction)
                target_reaction = TargetReaction(rctname)
                # constrain the reaction flux accordingly
                if positive: 
                    target_reaction.upper_bound = 'zero' 
                else: 
                    target_reaction.lower_bound = 'zero'
                # add the target flux value of the target reaction to the list of reaction fluxes in the target group object
                target_group.reaction_fluxes.append(target_reaction)
        
            # add the target reaction group to the model 
            model.targets.target_groups.append(target_group)
            print('Model building from XML files with reversed reaction ...')
            results = model.solve() # solve the new model
            print('Optimal growth rate is {}.'.format(results.mu_opt))
            results.write(output_dir) # needed in order to calculate costs later
        
            #save optimal growth rate 
            new_rev_reaction.gr_reversed = results.mu_opt
            # calculate cost and save cost 
            new_rev_reaction.cost_rev_flux = parse_enzyme_cost(output_dir, curr_dir)
            # calculate cost delta
            new_rev_reaction.cost_delta =  new_rev_reaction.cost_rev_flux - unc_AA_cost
            # calculate cost per growth rate
            if new_rev_reaction.gr_reversed > 0:
                new_rev_reaction.cost_per_gr = new_rev_reaction.cost_rev_flux/new_rev_reaction.gr_reversed
            else: 
                new_rev_reaction.cost_per_gr = float("inf")  
        
            #remove target group from the model = reset model
            model.targets.target_groups.remove(target_group) 
    
        # save the dictionary of reversible reactions as a pickle file, to be unpickled later
        # im Loop, dass es auch etwas speichert wenn es irgendwann nicht mehr klappt
        pickle.dump(rev_reactions, open( "save_rev_rcts.p", "wb" ))

    
if __name__ == '__main__':
    main()
    
    
    

