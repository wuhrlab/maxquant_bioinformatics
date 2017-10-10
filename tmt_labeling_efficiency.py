from __future__ import print_function
import pandas as pd
import numpy as np


def calculate_tmt_labeling_efficiency():
    # read in the data
    data = pd.read_csv('modificationSpecificPeptides.txt', '\t')
    
    # remove peptides with an n-terminal cysteine 
    data['has_n-term_cysteine'] = data['Sequence'].apply(lambda x: x[0] == 'C')
    data = data[~data['has_n-term_cysteine']]
    
    # determine whether TMT was on the n-terminus or not for each peptide
    data['has_TMT_on_n-term'] = data['Modifications'].apply(lambda x: 'dynamicTMTnTerm' in x)
    
    # count the total number of peptides in the data
    total_num_of_peptides = float(len(data))
    
    # count the number of peptides with an n-terminal TMT
    num_of_peptides_w_n_term_TMT = float(len(data[data['has_TMT_on_n-term']]))
    
    # print the total number of peptides identified
    print('Total number of peptides in modificationSpecificPeptides.txt: ' + str(total_num_of_peptides))
    
    # print the number of peptides with a TMT on the n-terminus
    print('Number of peptides in modificationSpecificPeptides.txt with a TMT on the n-terminus: ' + str(num_of_peptides_w_n_term_TMT))
    
    # print the fraction of all peptides with a TMT on the n-terminus
    print('Estimated labeling efficiency is: ' + str(np.around(100. * num_of_peptides_w_n_term_TMT / total_num_of_peptides, 1)) + str('%'))
    
calculate_tmt_labeling_efficiency()