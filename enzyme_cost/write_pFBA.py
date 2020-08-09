""" writing the modes and results into a file """

# currently -> test writing to files
import pickle
import csv
import os.path

def main():
    
    os.chdir('/Users/medeafux/Desktop/ETH/Masterarbeit/enzyme_cost')
    
    with open('pfba_solutions.p', 'rb') as pickle_file:
        pfba_sol = pickle.load(pickle_file)
    with open('weighted_solutions.p', 'rb') as pickle_file2:
        weighted_sol = pickle.load(pickle_file2)
        
    # output results into a file: pFBA
    with open('pfba-sol.csv', 'a', newline='') as csvfile:
        fieldnames =['Reaction Name','cost forward', 'cost backward']
        writer = csv.DictWriter(csvfile, fieldnames = fieldnames)
        writer.writeheader()    
       
        # go through all reactions and write down
        for reaction in pfba_sol:  
           writer.writerow({'Reaction Name': reaction,
                            'cost forward':  pfba_sol[reaction]['lower_bound_cost'],
                            'cost backward': pfba_sol[reaction]['upper_bound_cost']})
    
    # output results into a file: weighted pFBA
    with open('weighted-pfba-sol.csv', 'a', newline='') as csvfile:
        fieldnames =['Reaction Name','cost forward', 'cost backward']
        writer = csv.DictWriter(csvfile, fieldnames = fieldnames)
        writer.writeheader()    
       
        # go through all reactions and write down
        for reaction in weighted_sol:  
            writer.writerow({'Reaction Name': reaction,
                             'cost forward': weighted_sol[reaction]['lower_bound_cost'],
                            'cost backward': weighted_sol[reaction]['upper_bound_cost']})
            
if __name__ == '__main__':
    main()

