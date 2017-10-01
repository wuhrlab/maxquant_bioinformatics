from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
import seaborn as sns
sns.set_style('ticks')

def calculate_percent_cysteines():
    data = pd.read_csv('modificationSpecificPeptides.txt', '\t')
    total_num_of_peptides = float(len(data))
    data['pep_contains_C'] = data['Sequence'].apply(lambda x: 'C' in x)
    num_peptides_with_cysteine = float(len(data[data['pep_contains_C']]))
    print('Total number of peptides in modificationSpecificPeptides.txt: ' + str(total_num_of_peptides))
    print('Number of peptides with a Cysteine (C) is: ' + str(num_peptides_with_cysteine))
    print('Fraction of peptides with a Cysteine (C) is: '+  str(np.around(100. *num_peptides_with_cysteine / total_num_of_peptides ,1)) + '%')

calculate_percent_cysteines()