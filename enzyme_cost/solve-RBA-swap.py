"""Solve standard RBA problem."""

# python 2/3 compatibility
from __future__ import division, print_function

# global imports
import sys
import os.path
import time

# package imports
import rba
from docplex.cp.modeler import true


def main():
    if len(sys.argv) < 2:
        print('Please provide path to directory containing xml files.')
    else:
        xml_dir = sys.argv[1]      
        if len(sys.argv) >= 3:
            # output_dir = sys.argv[2]      #this has to be used if not run on a specific model
            # I defined here the output directory manually, so that the program can be debugged easily
            output_dir = '/Users/medeafux/Desktop/ETH/Masterarbeit/Bacterial-RBA-models-master/Escherichia-coli-K12-WT'
        else:
            output_dir = xml_dir  #why should there ever be such an input that this is needed? why should the programme then continue? 

        # load model, build matrices and solve
        print('Model building from XML files ...')
        model = rba.RbaModel.from_xml(xml_dir)
        print('Starting iterative RBA resolution...')

        ##########
        #adapt the model according to the modes from hyperflux for the TKT2 test case 
        wantedRct = model.metabolism.reactions.get_by_id('R_GHMT2r')
        # set revesible boolean to True (if False) 
        wantedRct.reversible = False
        
        print(wantedRct.reversible) # just to check 
        
        #swap the reactants and the products
        temp = wantedRct.products
        wantedRct.products = wantedRct.reactants
        wantedRct.reactants = temp

        
        # start timer to track how long the rba calculation takes
        start = time.process_time()
        #solve the actual model 
        results = model.solve()
        # print the time it has taken to solve model 
        print(time.process_time()-start)

        print('Optimal growth rate is {}.'.format(results.mu_opt))
        results.write(output_dir)
        
        # used to write the results in such a way that they can be displayed as an Escher Map
        results.write_fluxes(os.path.join(output_dir, 'Escher'), file_type='json', merge_isozyme_reactions=True,only_nonzero=True,remove_prefix=True)
        # used to write the results in such a way that they can be displayed als a Proteomap
        results.write_proteins(os.path.join(output_dir, 'Proteomaps'), file_type='csv')
        
        
        results.export_matlab(output_dir)

if __name__ == '__main__':
    main()
