from __future__ import print_function
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
import seaborn as sns
sns.set_style('ticks')

def map_gene_symbol_for_each_row(protein_id_row):
    gene_symbol_group_list = np.array([])
    xen_laevis_ids = protein_id_row.split(';')
    print (xen_laevis_ids)
    for i in xen_laevis_ids:
        gs = phrog_annotation[phrog_annotation['ID'] == i]['Gene Symbol'].values[0]
        if isinstance(gs, str):
            gene_symbol_group_list = np.append(gene_symbol_group_list, np.array([gs]))
#             if len(gene_symbol_group_list) >= 1:
# #                 gene_symbol_group_list = gene_symbol_group_list.append(';' + gs)
#                 gene_symbol_group_list = np.append(gene_symbol_group_list, np.array([';' + gs]))
#             else:
# #                 gene_symbol_group_list = gene_symbol_group_list.append(gs)
#                 gene_symbol_group_list = np.append(gene_symbol_group_list, np.array([gs]))
    gene_symbol_group_list = np.unique(gene_symbol_group_list)
    return gene_symbol_group_list