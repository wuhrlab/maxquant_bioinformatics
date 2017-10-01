from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
import seaborn as sns
sns.set_style('ticks')

def calculate_missed_cleavage_rate(enzyme = 'trypsin'):
    # read in the data
    data = pd.read_csv('modificationSpecificPeptides.txt', '\t')
    
    # count the total number of peptides in the data
    total_num_of_peptides = float(len(data))
    
    # for a trypsin digest:
    if enzyme.lower() == 'trypsin':
        # count number of internal K and R
        data['has_missed_cleavage'] = data['Sequence'].apply(lambda x: ('K' in x[:-1]) or ('R' in x[:-1]) )
        
        pep_w_missed_cleavage = float(len(data[data['has_missed_cleavage']]))
        
        print('Total number of peptides in modificationSpecificPeptides.txt: ' + str(total_num_of_peptides))
        print('Number of peptides with a missed cleavage (Trypsin) is: ' + str(pep_w_missed_cleavage))
        print('Estimated missed cleavage rate: '+  str(np.around(100. *pep_w_missed_cleavage / total_num_of_peptides ,1)) + '%')
        
    if enzyme.lower() == 'lysc':
        data['has_missed_cleavage'] = data['Sequence'].apply(lambda x: 'K' in x[:-1] )
        pep_w_missed_cleavage = float(len(data[data['has_missed_cleavage']]))
        
        print('Total number of peptides in modificationSpecificPeptides.txt: ' + str(total_num_of_peptides))
        print('Number of peptides with a missed cleavage (LysC) is: ' + str(pep_w_missed_cleavage))
        print('Estimated missed cleavage rate: '+  str(np.around(100. *pep_w_missed_cleavage / total_num_of_peptides ,1)) + '%')

calculate_missed_cleavage_rate(enzyme = 'lysc')