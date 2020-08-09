"""Extract the tRNA-AA fluxes to define the overall cost of enzyme production"""

import os

def parse_enzyme_cost(output_dir, curr_dir):

    # go into the right directory
    os.chdir(output_dir)
    
    #create an empty list for the AA reactions
    AAfluxes = list()
    # list of AA that we're looking for (manually done), TRS -> flux from TRNA-AA to the protein
    targetAA = ['ALATRS', 'ARGTRS', 'ASNTRS', 'ASPTRS', 'CYSTRS', 'GLUTRS', 'GLNTRS', 
                'GLYTRS','HISTRS','ILETRS', 'LEUTRS','LYSTRS','METTRS','PHETRS','PROTRS',
                'SERTRS', 'THRTRS','TRPTRS','TYRTRS','VALTRS']
    
    #open reactions.out file to get the fluxes
    with open('reactions.out', 'r') as reactions: 
       # go through reactions.out line by line
       for line in reactions:
           # check for each reaction if it contains the name of an amino acid
           if any(ele in line for ele in targetAA):
               #add found amino acid + flux to AAfluxes
               AAfluxes.append(line)          
    
    #go through AAfluxes and add up all the fluxes to get the overall cost
    cost = 0;
    #go through all found fluxes, split them up into name and flux and then add the cost to the total cost
    for AA in AAfluxes: 
        #print(AA)
        splitted = AA.split('\t')
        cost += float(splitted[1])
   
    # go back to original directors
    os.chdir(curr_dir)
    
    return cost 
    
 
  
  
    
  
   
  