from __future__ import print_function
import pandas as pd
import numpy as np

# read in protein data from MaxQuant
protein_data = pd.read_csv('proteinGroups.txt', '\t')

# read in spreadsheet with mappings from X. laevis contigs to human gene symbols
# using the tool available at http://kirschner.med.harvard.edu/tools/genesym_assignment.html 
phrog_annotation = pd.read_csv('phrog_annotation.csv')

print ('Patience! This script takes approximately 1 sec / 50 proteins in your proteinGroups.txt file.')

def map_gene_symbol_for_each_row(protein_id_row):
    gene_symbol_group_list = np.array([])
    xen_laevis_ids = protein_id_row.split(';')
    
    all_phrog_contig_ids = phrog_annotation['ID'].values
    
    for i in xen_laevis_ids:
        if i in all_phrog_contig_ids:
            gs = phrog_annotation[phrog_annotation['ID'] == i]['Gene Symbol'].values[0]
            if isinstance(gs, str):
                gene_symbol_group_list = np.append(gene_symbol_group_list, np.array([gs]))
        else:
            gene_symbol_group_list = []
    gene_symbol_group_list = np.unique(gene_symbol_group_list)
    count = 0
    gene_symb_string = ''
    for j in gene_symbol_group_list:
        if count == 0:
            gene_symb_string = gene_symb_string + j
        else:
            gene_symb_string = gene_symb_string + ';' + j
        count = count + 1
    return gene_symb_string

# for each row map the frog reference id to the human gene symbol if such a mapping exists
# in the phrog_annotation file
protein_data['mapped_gene_symbols'] = protein_data['Protein IDs'].apply(lambda x: map_gene_symbol_for_each_row(x))
print ('Number of proteins in your file: ' + str(len(protein_data)))
print ('Number of proteins matched to at least one human gene symbol: ' + str(len(protein_data[~(protein_data['mapped_gene_symbols'] == '')]))) 

# create a new file containing the human gene symbol mapping in the last column (will be blank if no map)
protein_data.to_csv('proteins_mapped_to_gene_symbols.csv', index = False)
