""" extracting the cost calculated by all three methods (RBA, pFBA, weighted pFBA) and correlate the results + plots """

# currently -> test writing to files
import pickle
import csv
import os.path
import rba
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sn
from heatmap import heatmap, corrplot
import scipy

def main():
        
    os.chdir('/Users/medeafux/Desktop/ETH/Masterarbeit/enzyme_cost')
    
    with open('pfba_solutions.p', 'rb') as pickle_file:
        pfba_sol = pickle.load(pickle_file)
    with open('weighted_solutions.p', 'rb') as pickle_file2:
        weighted_sol = pickle.load(pickle_file2)
    with open('save_rev_rcts.p', 'rb') as pickle_file3:
        rba_sol = pickle.load(pickle_file3)
        
    
    #save always the constraint version 
    reaction_names =[]
    pfba = []
    w_pfba = []
    rba = []
    
    # extract values into lists 
    for reaction in pfba_sol: 
        # get the higher value = constrained version
        reaction_names.append(reaction)
        pfba.append(max(pfba_sol[reaction]['lower_bound_cost'], pfba_sol[reaction]['upper_bound_cost']))
        w_pfba.append(max(weighted_sol[reaction]['lower_bound_cost'], weighted_sol[reaction]['upper_bound_cost']))
       
        
        # rba.append(rba_sol['R_' + reaction].cost_rev_flux)
        # changed to growth rate -> we compare the growth rate of RBA
        rba.append(rba_sol['R_' + reaction].gr_reversed)
          
    # store values to file       
    with open('costs-comparison.csv', 'w', newline='') as csvfile:
        fieldnames =['reaction','rba','pfba', 'weighted pfba']
        writer = csv.DictWriter(csvfile, fieldnames = fieldnames)
        writer.writeheader()  
        for i in range(len(rba)):
            writer.writerow({'reaction': reaction_names[i], 'rba': rba[i] , 'pfba': pfba[i], 'weighted pfba': w_pfba[i]})
    
    # Set the font dictionaries (for plot title and axis titles)
    title_font = {'fontname':'Palatino', 'size':'24', 'color':'black', 'weight':'normal',
              'verticalalignment':'bottom'} # Bottom vertical alignment for more space
    axis_font = {'fontname':'Palatino', 'size':'20'}
    
    # creating dataframe
    data = {'rba':rba, 'pfba': pfba, 'w_pfba': w_pfba}
    df = pd.DataFrame(data,columns=['rba','pfba','w_pfba']) 
    #remove 'None' values (otherwise later problems with the calculation of correlation)
    df = df[df.rba.notnull()]
    # sort the data according to rba growth rate in ascending order
    sorted_df = df.sort_values(by ='rba' )
    # remove all values with growth rate lower than 0.58
    zoom_df = sorted_df[sorted_df.rba > 0.58]

    # plots
    #plt.figure(1)   
    # correlation matrix            
    corrMatrix = df.corr(method='pearson', min_periods=1)
    print(corrMatrix)
    #sn.heatmap(corrMatrix, annot=True)
    #plt.title('Pearson Correlation Coefficients')

    #plt.figure(2)   
    # correlation matrix             
    zoom_corrMatrix = zoom_df.corr(method='pearson', min_periods=1)
    print(zoom_corrMatrix)
    #sn.heatmap(corrMatrix, annot=True)
    #plt.title('Spearman Correlation Coefficients')
    
    rba_pfba_corr = scipy.stats.pearsonr(df['rba'], df['pfba'])
    rba_wpfba_corr = scipy.stats.pearsonr(df['rba'], df['w_pfba'])
    zoom_rba_pfba_corr = scipy.stats.pearsonr(zoom_df['rba'], zoom_df['pfba'])
    zoom_rba_wpfba_corr = scipy.stats.pearsonr(zoom_df['rba'], zoom_df['pfba'])
    
    print('rba vs pfba') 
    print(rba_pfba_corr)
    print('rba vs wpfba')
    print(rba_wpfba_corr)
    print('ZOOM rba vs pfba')
    print(zoom_rba_pfba_corr)
    print('ZOOM rba vs wpfba')
    print(zoom_rba_wpfba_corr)
        
    plt.figure(1)
    # plot pfba vs rba
    plt.scatter(df['rba'], df['pfba'])
    plt.xlabel('RBA growth rate [1/h]', **axis_font)
    plt.ylabel('pFBA calculated costs [a.u.]', **axis_font)
    #plt.title('RBA vs. pFBA', **title_font)
    
    plt.figure(2)
    # plot weighted pfba vs rba
    plt.scatter(df['rba'], df['w_pfba'])
    plt.xlabel('RBA growth rate [1/h]', **axis_font)
    plt.ylabel('weighted pFBA calculated costs [a.u.]', **axis_font)
    #plt.title('RBA vs. weighted pFBA', **title_font)
    
    plt.figure(3)
    # plot pfba vs rba with growth rate cutoff at 0.58
    plt.scatter(zoom_df['rba'], zoom_df['pfba'])
    plt.xlabel('RBA growth rate [1/h]', **axis_font)
    plt.ylabel('pFBA calculated costs [a.u.]', **axis_font)
    #plt.title('RBA vs. pFBA - growth rate cutoff at > 0.58', **title_font)
    
    plt.figure(4)
    # plot weighted pfba vs rba with growth rate cutoff at 0.58
    plt.scatter(zoom_df['rba'], zoom_df['w_pfba'])
    plt.xlabel('RBA growth rate [1/h]', **axis_font)
    plt.ylabel('weighted pFBA calculated costs [a.u.]', **axis_font)
    #plt.title('RBA vs. weighted pFBA - growth rate cutoff at > 0.58', **title_font)
    
    
    
    plt.show()
    

    
    # trying fancy heatmap plot
#    plt.figure(figsize=(8, 8))
#    corrplot(df.corr(), size_scale=300);
            
if __name__ == '__main__':
    main()

