
import csv
import matplotlib.pyplot as plt
import os.path
from scipy import stats

def main():
    
    project_dir = '/Users/medeafux/Desktop/ETH/Masterarbeit/enzyme_cost'
    os.chdir(project_dir)
    
    # defining lists to save the data
    modes = []
    counts = []
    prob = []
    AA_cost_pgr = []
    AA_cost = []
    gr = []
    
    # parse the data from the result file
    # TODO: make this different and use pickle to save the class as class and then plot from there
    with open('overview-modes.csv') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        next(reader)
    
        for row in reader:
            modes.append(float(row[0]))
            # counts.append(float(row[1]))   -> unconstrained has value 'unc' in it = string!
            prob.append(float(row[2]))
            AA_cost_pgr.append(float(row[3]))
            AA_cost.append(float(row[4]))
            gr.append(float(row[5]))
            
    # Set the font dictionaries (for plot title and axis titles)
    title_font = {'fontname':'Palatino', 'size':'24', 'color':'black', 'weight':'normal',
              'verticalalignment':'bottom'} # Bottom vertical alignment for more space
    axis_font = {'fontname':'Palatino', 'size':'20'}

    plt.scatter(gr, prob)
    plt.xlabel('Maximal growth rate of mode [1/s]', **axis_font)
    plt.ylabel('Probability of mode [%]', **axis_font)
    #plt.title('RBA growth rate vs. probability of mode', **title_font)
    plt.ylim(0,100)
    plt.show()
    
    #calculate the correlation
    prob_gr_corr = stats.pearsonr(gr, prob)
    print(prob_gr_corr)
    
    plt.scatter(AA_cost_pgr, prob)
    plt.xlabel('AA cost per growth rate')
    plt.ylabel('probability of mode')
    plt.title('Cost per mode vs. probability of mode')
    plt.ylim(0,100)
    plt.show()
    
    # plot the difference
    a = AA_cost_pgr[len(AA_cost_pgr)-1] # get the unconstrained value
    AAcost_diff_val = [x - a for x in AA_cost_pgr] # substract it from each AA cost per growth rate value 
    
    plt.scatter(AAcost_diff_val, prob)
    plt.xlabel('AA cost per growth rate relative to unconstrained model')
    plt.ylabel('probability of mode')
    plt.title('relative cost per mode vs. probability of mode')
    plt.ylim(0,100)
    plt.show()
        
main()