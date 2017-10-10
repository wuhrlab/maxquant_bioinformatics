from __future__ import print_function
import pandas as pd
import numpy as np
import argparse
from sys import exit, argv

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

parser = argparse.ArgumentParser(description='Missed cleavage rate command line')
parser.add_argument('--num_clusters', type = int)
parser.add_argument('--line_opacity', type = float)
parser.add_argument('--folder_name')

args = parser.parse_args()
num_clusters = args.num_clusters
opacity = args.line_opacity
folder_name = args.folder_name

import os
cwd = os.getcwd()

import os
import errno

def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

make_sure_path_exists(cwd + '\\' + folder_name)


# read in protein data w gene symbosl matched
protein_data = pd.read_csv('proteins_mapped_to_gene_symbols.csv')

# calculate the sum of all the reporter intensities from all 10 channels
protein_data['sum_all_reporter_intensities'] = protein_data[[u'Reporter intensity corrected 0', u'Reporter intensity corrected 1',
       u'Reporter intensity corrected 2', u'Reporter intensity corrected 3',
       u'Reporter intensity corrected 4', u'Reporter intensity corrected 5',
       u'Reporter intensity corrected 6', u'Reporter intensity corrected 7',
       u'Reporter intensity corrected 8', u'Reporter intensity corrected 9']].sum(axis = 1)

# drop any rows that don't have any signal in any channel
protein_data = protein_data[protein_data['sum_all_reporter_intensities'] > 0]

# normalize the reporter ion counts so that the total sum for each row (protein) is 1.0
protein_data[[u'Reporter intensity corrected 0', u'Reporter intensity corrected 1',
       u'Reporter intensity corrected 2', u'Reporter intensity corrected 3',
       u'Reporter intensity corrected 4', u'Reporter intensity corrected 5',
       u'Reporter intensity corrected 6', u'Reporter intensity corrected 7',
       u'Reporter intensity corrected 8', u'Reporter intensity corrected 9']] = protein_data[[u'Reporter intensity corrected 0', u'Reporter intensity corrected 1',
       u'Reporter intensity corrected 2', u'Reporter intensity corrected 3',
       u'Reporter intensity corrected 4', u'Reporter intensity corrected 5',
       u'Reporter intensity corrected 6', u'Reporter intensity corrected 7',
       u'Reporter intensity corrected 8', u'Reporter intensity corrected 9']].apply(lambda row: row / np.sum(row), axis = 1)

# perform kmeans clustering on the data
kmeans = KMeans(n_clusters=num_clusters)
kmeans.fit(protein_data[[u'Reporter intensity corrected 0', u'Reporter intensity corrected 1',
       u'Reporter intensity corrected 2', u'Reporter intensity corrected 3',
       u'Reporter intensity corrected 4', u'Reporter intensity corrected 5',
       u'Reporter intensity corrected 6', u'Reporter intensity corrected 7',
       u'Reporter intensity corrected 8', u'Reporter intensity corrected 9']].values)

protein_data['k_means_label'] = kmeans.labels_
protein_data.to_csv(cwd + '\\' + folder_name + '\\' + 'proteins_mapped_to_gene_symbols_w_kmeans_clusters.csv', index = False)

# make a plot for each K Means cluster
for cluster in np.arange(num_clusters):
	# plot the line of each protein
	for j in protein_data[protein_data['k_means_label'] == cluster].T:
		plt.plot(np.arange(10), protein_data.ix[j][[u'Reporter intensity corrected 0', u'Reporter intensity corrected 1',
       u'Reporter intensity corrected 2', u'Reporter intensity corrected 3',
       u'Reporter intensity corrected 4', u'Reporter intensity corrected 5',
       u'Reporter intensity corrected 6', u'Reporter intensity corrected 7',
       u'Reporter intensity corrected 8', u'Reporter intensity corrected 9']], color = 'grey', alpha = opacity)
	plt.ylim((-.05, 1))
	plt.xticks(size = 12)
	plt.yticks(size = 12)
	plt.xlabel('Reporter channel #', size = 12, weight = 'semibold')
	plt.ylabel('Relative abundance', size = 12, weight = 'semibold')
	plt.title('Cluster number: ' + str(num_clusters) + ' ,N = ' + str(len(protein_data[protein_data['k_means_label'] == cluster])) + ' proteins')
	plt.savefig(cwd + '\\' + folder_name + '\\' + 'cluster_' + str(cluster) + '_protein_relative_abundances.png', dpi = 300)
	plt.savefig(cwd + '\\' + folder_name + '\\' + 'cluster_' + str(cluster) + '_protein_relative_abundances.svg')
	plt.close()

	cluster_gs_list = []
	for k in protein_data[protein_data['k_means_label'] == cluster]['mapped_gene_symbols']:
		if k is not np.nan:
			for l in k.split(';'):
				cluster_gs_list.append(l)
	cluster_gs_list = set(cluster_gs_list)
	pd.DataFrame({'human_gene_symbols' : list(cluster_gs_list)}).to_csv(cwd + '\\' + folder_name + '\\' +'cluster_num_' + str(cluster) + '.csv', index = False)

# determine the background (all human gene symbols identified)
background = []
for i in protein_data['mapped_gene_symbols']:
	if i is not np.nan:
		for j in i.split(';'):
			background.append(j)
# get rid of duplicates
background = set(background)
pd.DataFrame({'human_gene_symbols' : list(background)}).to_csv(cwd + '\\' + folder_name + '\\' +'protein_background.csv', index = False)
