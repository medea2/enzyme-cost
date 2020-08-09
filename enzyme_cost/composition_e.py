""" calculate the cost of an enzyme via it's composition, get composition from rba """


from lxml import etree
import os.path
from macrocomponent import macrocomponent
from enzyme import enzyme
import csv
import libsbml
import re
from statistics import median 

import sys 
# to get information out of brenda:
from brendapy import BrendaParser

# to calculate the geometric mean in an easy way
from scipy.stats.mstats import gmean

import pickle
from ipywidgets.widgets.interaction import empty

BRENDA_PARSER = BrendaParser()  # reuse parser     

def main():
     
    # go to directory containing the xml files 
    # TODO: make path relative
    project_dir = '/Users/medeafux/Desktop/ETH/Masterarbeit/enzyme_cost'
    model_dir = '/Users/medeafux/Desktop/ETH/Masterarbeit//Bacterial-RBA-models/Escherichia-coli-K12-WT'
    
    os.chdir(model_dir)
    
    # extract composition of all macrocomponents and save them as a list (?) of macromolecules (class)
    tree = etree.parse("proteins.xml")
    root = tree.getroot()
    # get to the list of Macromolecules
    listMM = root[1]
    
    #create a dictionary to store the macromolecules as instances of the class macrocomponent
    macrocomponents = {} 
    
    #iterate over all  macromolecules in list
    for MM in listMM: 
        id = MM.get('id') # get the id of the macromolecule
        comp_cost = 0   # initially set its cost to zero 
        # MM[0] = composition -> use this to calculate the cost of the macromolecule
        for compRef in MM[0]:
            comp_cost += float(compRef.get('stoichiometry'))
        # create a new macromolecule with the found id and calculated cost (= composition) and add to dictionary
        macrocomponents[id] = comp_cost
        
    # extract composition of enzymes, create list of enzymes 
    tree = etree.parse("enzymes.xml")
    root = tree.getroot()
    listEnzymes = root[0]
    
    os.chdir(project_dir)
    
    # make a dictionary to store all enzymes 
    enzymes = {}
    # go through all enzymes
    for enz in listEnzymes: 
        # get the name of the reaction that is catalyzed by a specific enzyme
        reaction_name = enz.get('reaction')
        
        try:
            # get the id of the macromolecule that constitutes the enzyme and it's corresponding stoichiometry
            enzy_macrocom = [enz[0][0][0].get('species'), float(enz[0][0][0].get('stoichiometry'))]
            #create new enzyme and append it to the list
            new_enzyme = enzyme(reaction_name, enzy_macrocom, 0, 0, 0)
            # add new enzyme to the dictionary of enzymes
            enzymes[reaction_name] = new_enzyme
           
        # catch the error that appears when enzymes don't have a machinery composition and list them into one file
        except IndexError:
            # print('Enzyme of ' + reaction_name + ' does not have machinery composition')
            missing = open("missing_enzymes.txt", 'a')
            missing.write(reaction_name + ', ')
    
    # calculate the cost of each enzyme according to macromolecule composition
    # go over all enzymes
    for enz in enzymes:
        # search for id of enzyme component in the dictionary of macrocomponents and get according cost
        if macrocomponents.__contains__(enzymes[enz].macrocomponents[0]):
            # calculate the cost of the enzyme according to its composition and stoichiometry
            enzymes[enz].cost = enzymes[enz].macrocomponents[1] * macrocomponents.get(enzymes[enz].macrocomponents[0])
        
    # get the EC number of each enzyme and save it into the enzymes
    # got to the model file 
    os.chdir(model_dir + '/data')
    reader = libsbml.SBMLReader()
    document = reader.readSBML('iJO1366.xml')
    sbml_model = document.getModel()
    
    # create empty dictionary to store reaction name and EC number
    Dict = {} 
    # go through all reactions of the model 
    for r in sbml_model.getListOfReactions():
        # extract the ec number if the reaction has one
        if "ec-code" in r.annotation_string:
            text = r.annotation_string
            ec = re.findall("/ec-code/\S.\S.\S.\S", text)
           
            ec = str(ec)
            ec = ec[11:]
            ec = ec[:-2]
            
            # take the last EC number if we have more then one EC number in one enzyme
            if len(ec) > 7:
                ec = ec[-7:]
                
            # if an EC number contains 'l' replace it with '1'
            ec.replace('l','1')
            
            # turn ec number into string and remove everything except the number itself
            # unsolved problem: some reactions have several ec numbers!
            Dict[r.id] = ec
            
    # go over all enzymes and store the EC number for the enzymes we have it
    # and search in brendapy and get out the turnover number 
    for enz in enzymes:
        # get EC number
        if Dict.__contains__(enzymes[enz].reaction_name):
            enzymes[enz].EC_number = Dict.get(enzymes[enz].reaction_name)
            
            # search on BRENDA for the EC number
            try:
                proteins = BRENDA_PARSER.get_proteins(enzymes[enz].EC_number)
            except KeyError:
                print('EC number invalid')  # what to do with EC numbers that contain '-' ? 
                
                
            # calculate the median of all available TN numbers for this EC number 
            # TODO: change somehow so that this only needs to be calculated when there is no TN found for E. Coli
            allorg_tn_values = [] # list to store all the TN numbers calculated for all organism
            for p in proteins.values():
                oneorg_tn_values = [] #list for all tn values
                
                if p.TN is not None:
                    # try for each condition to get a value (not all of them have)
                    for condition in p.TN:
                        # add the value of each condition, if there is one
                        try:
                            oneorg_tn_values.append(condition['value'])
                        except KeyError:
                            print("No value for this TN condition was found.")
                            
                    # calculate geometric mean of all found values to get final turnover number
                    final_tn = gmean(oneorg_tn_values)
                    
                    allorg_tn_values.append(final_tn)
                        
            # calculate the median of all tn values found for all organisms 
            if len(allorg_tn_values) > 0:
                backup_tn = median(allorg_tn_values)     
            
        
            # go through all the proteins found for that EC number
            for p in proteins.values(): 
                if p.organism == "Escherichia coli":
                    # list to store all turnover numbers of all conditions of the enzyme
                    tn_values = []
            
                    # try to get to the turnover numbers in general -> several conditions
                    if p.TN is not None:
                        # try for each condition to get a value (not all of them have)
                        for condition in p.TN:
                            # add the value of each condition, if there is one
                            try:
                                print(condition['value'])
                                tn_values.append(condition['value'])
                            except KeyError:
                                print("No value for this TN condition was found.")
                         
                            # calculate geometric mean of all found values to get final turnover number
                            final_tn = gmean(tn_values)
            
                    # there is no TN 
                    else:
                        print("TN not found")
                        # TODO: add estimation of this value -> = median over all organisms
                        final_tn = backup_tn  

                    #add tn to the enzyme
                    print(final_tn)
                    enzymes[enz].turnover = final_tn
        
        # no EC number found           
        else:
            missingEC = open("missing_EC.txt", 'a')
            missing.write(reaction_name)
            enzymes[enz].EC_number = 'not found'
    
    # save the list of enzymes as a pickle file, to be unpickled later
    # im Loop, dass es auch etwas speichert wenn es irgendwann nicht mehr klappt
    os.chdir(project_dir + '/enzyme_cost')
    pickle.dump(enzymes, open( "save_enzymes.p", "wb" ))
    
 
main()
    
    
    
    