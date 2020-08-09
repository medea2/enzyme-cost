""" writing the modes and results into a file """

# currently -> test writing to files
import pickle
import csv
import os.path

def main():
    
    os.chdir('/Users/medeafux/Desktop/ETH/Masterarbeit/enzyme_cost')
    
    with open('save_rev_rcts.p', 'rb') as pickle_file:
        rev_reactions = pickle.load(pickle_file)
    
    
    # output results into a file         
    with open('reversible-reactions.csv', 'w', newline='') as csvfile:
        fieldnames =['reaction name','original flux','AA costs','AA cost delta','growth rate','AA cost per growth rate', '#duplicate reactions', 'original direction']
        writer = csv.DictWriter(csvfile, fieldnames = fieldnames)
        writer.writeheader()
       
        # go through all reversible reactions and write their content into a csv file
        for reaction in rev_reactions:  
            writer.writerow({'reaction name': reaction, 'original flux': rev_reactions[reaction].original_flux,
                             'AA costs': rev_reactions[reaction].cost_rev_flux, 'AA cost delta': rev_reactions[reaction].cost_delta, 
                             'growth rate': rev_reactions[reaction].gr_reversed, 'AA cost per growth rate': rev_reactions[reaction].cost_per_gr,
                             '#duplicate reactions': rev_reactions[reaction].duplicates,
                             'original direction': rev_reactions[reaction].orig_direction})
            
if __name__ == '__main__':
    main()

