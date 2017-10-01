import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
import seaborn as sns]
from sklearn.cluster import KMeans
sns.set_style('ticks')

phrog = pd.read_csv('lysC_meera_24_fracs_X_laevis_developmental_time_series_PHROG_peptides_after_protein_assembler.csv')

human_gs_mapping_jgi = pd.read_csv('JGI_9_1_mapped_human_gs.csv')
human_gs_mapping_phrog = pd.read_csv('klab_08r2_rm_w_annotation.csv')

phrog[[ u'126 Sn', u'127N Sn', u'127C Sn',
       u'128N Sn', u'128C Sn', u'129N Sn', u'129C Sn', u'130N Sn', u'130C Sn',
       u'131 Sn']] = phrog[[ u'126 Sn', u'127N Sn', u'127C Sn',
       u'128N Sn', u'128C Sn', u'129N Sn', u'129C Sn', u'130N Sn', u'130C Sn',
       u'131 Sn']].apply(lambda row: row / np.sum(row), axis = 1)

phrog[[ u'126 Sn', u'127N Sn', u'127C Sn',
       u'128N Sn', u'128C Sn', u'129N Sn', u'129C Sn', u'130N Sn', u'130C Sn',
       u'131 Sn']].head()



kmeans = KMeans(n_clusters=4, random_state=0).fit(phrog[[ u'126 Sn', u'127N Sn', u'127C Sn',
       u'128N Sn', u'128C Sn', u'129N Sn', u'129C Sn', u'130N Sn', u'130C Sn',
       u'131 Sn']].values)

phrog['kmeans_cluster_label'] = kmeans.labels_



plt.figure(figsize = (10, 7))
for i in phrog.T:
    plt.plot(np.arange(10), phrog[phrog['protein_id'] == i][[u'126 Sn', u'127N Sn', u'127C Sn',
       u'128N Sn', u'128C Sn', u'129N Sn', u'129C Sn', u'130N Sn', u'130C Sn',
       u'131 Sn']].values.reshape(10, 1), color = 'black', alpha = 0.005)
plt.ylim((-.02, 1.02))
plt.ylim((-.02, 1.02))
plt.legend(loc = 'best', fontsize = 'xx-large')
plt.xticks(np.arange(10), [1, 2, 8, 10, 12, 17, 20, 25, 30, 35], size = 26)
plt.yticks(size = 26)
plt.xlabel('Developmental stage', size = 26, weight = 'semibold')
plt.ylabel('Relative abundance', size = 26, weight = 'semibold')
plt.title('Relative protein abundance\n of all X. laevis proteins, N = 9,345 proteins', size = 20 ,weight = 'semibold')
plt.savefig('all_proteins_relative_abundance_correct_stages_black.svg', transparent = True)


plt.figure(figsize = (10, 7))
plt.plot(np.arange(10), phrog[phrog['kmeans_cluster_label'] == 0][[u'126 Sn', u'127N Sn', u'127C Sn',
       u'128N Sn', u'128C Sn', u'129N Sn', u'129C Sn', u'130N Sn', u'130C Sn',
       u'131 Sn']].median(), markersize = 20, marker = 'o', color = '#1b9e77', lw = 4, label = 'Cluster 0 median\nN = 5,451 X. laevis proteins')
plt.plot(np.arange(10), phrog[phrog['kmeans_cluster_label'] == 2][[u'126 Sn', u'127N Sn', u'127C Sn',
       u'128N Sn', u'128C Sn', u'129N Sn', u'129C Sn', u'130N Sn', u'130C Sn',
       u'131 Sn']].median(), markersize = 20, marker = 'o', color = 'grey', lw = 4, label = 'Cluster 2 median\nN = 2,697 X. laevis proteins')
plt.plot(np.arange(10), phrog[phrog['kmeans_cluster_label'] == 1][[u'126 Sn', u'127N Sn', u'127C Sn',
       u'128N Sn', u'128C Sn', u'129N Sn', u'129C Sn', u'130N Sn', u'130C Sn',
       u'131 Sn']].median(), markersize = 20, marker = 'o', color = '#7570b3', lw = 4, label = 'Cluster 1 median\nN = 860 X. laevis proteins')
plt.plot(np.arange(10), phrog[phrog['kmeans_cluster_label'] == 3][[u'126 Sn', u'127N Sn', u'127C Sn',
       u'128N Sn', u'128C Sn', u'129N Sn', u'129C Sn', u'130N Sn', u'130C Sn',
       u'131 Sn']].median(), markersize = 20, marker = 'o', color = '#e7298a', lw = 4, label = 'Cluster 3 median\nN = 338 X. laevis proteins')
plt.ylim((-.02, 1.02))
plt.legend(loc = 'best', fontsize = 'xx-large')
plt.xticks(np.arange(10), [1, 2, 8, 10, 12, 17, 20, 25, 30, 35], size = 26)
plt.yticks(size = 26)
plt.xlabel('Developmental stage', size = 26, weight = 'semibold')
plt.ylabel('Relative abundance', size = 26, weight = 'semibold')
plt.title('Median relative protein abundance\n of each K Means Cluster', size = 20 ,weight = 'semibold')
plt.savefig('median_Kmeans_clustering_4clusters_triceratops_processing_phrog.svg', transparent = True)

plt.figure(figsize = (10, 7))
for i in phrog[phrog['kmeans_cluster_label'] == 0].T:
    plt.plot(np.arange(10), phrog[phrog['protein_id'] == i][[u'126 Sn', u'127N Sn', u'127C Sn',
       u'128N Sn', u'128C Sn', u'129N Sn', u'129C Sn', u'130N Sn', u'130C Sn',
       u'131 Sn']].values.reshape(10, 1), color = '#1b9e77', alpha = 0.05)
plt.ylim((-.02, 1.02))
plt.ylim((-.02, 1.02))
plt.legend(loc = 'best', fontsize = 'xx-large')
plt.xticks(np.arange(10), [1, 2, 8, 10, 12, 17, 20, 25, 30, 35], size = 26)
plt.yticks(size = 26)
plt.xlabel('Developmental stage', size = 26, weight = 'semibold')
plt.ylabel('Relative abundance', size = 26, weight = 'semibold')
plt.title('Relative protein abundance\n of each protein in cluster 0, N = 5,451 proteins', size = 20 ,weight = 'semibold')
plt.savefig('colored_Kmeans_cluster_0_protein_data_triceratops_processing_phrog.svg', transparent = True)

plt.figure(figsize = (10, 7))
for i in phrog[phrog['kmeans_cluster_label'] == 1].T:
    plt.plot(np.arange(10), phrog[phrog['protein_id'] == i][[u'126 Sn', u'127N Sn', u'127C Sn',
       u'128N Sn', u'128C Sn', u'129N Sn', u'129C Sn', u'130N Sn', u'130C Sn',
       u'131 Sn']].values.reshape(10, 1), color = '#7570b3', alpha = 0.05)
plt.ylim((-.02, 1.02))
plt.ylim((-.02, 1.02))
plt.legend(loc = 'best', fontsize = 'xx-large')
plt.xticks(np.arange(10), [1, 2, 8, 10, 12, 17, 20, 25, 30, 35], size = 26)
plt.yticks(size = 26)
plt.xlabel('Developmental stage', size = 26, weight = 'semibold')
plt.ylabel('Relative abundance', size = 26, weight = 'semibold')
plt.title('Relative protein abundance\n of each protein in cluster 1, N = 860 proteins', size = 20 ,weight = 'semibold')
plt.savefig('colored_Kmeans_cluster_1_protein_data_triceratops_processing_phrog.svg', transparent = True)

plt.figure(figsize = (10, 7))
for i in phrog[phrog['kmeans_cluster_label'] == 2].T:
    plt.plot(np.arange(10), phrog[phrog['protein_id'] == i][[u'126 Sn', u'127N Sn', u'127C Sn',
       u'128N Sn', u'128C Sn', u'129N Sn', u'129C Sn', u'130N Sn', u'130C Sn',
       u'131 Sn']].values.reshape(10, 1), color = 'grey', alpha = 0.05)
plt.ylim((-.02, 1.02))
plt.ylim((-.02, 1.02))
plt.legend(loc = 'best', fontsize = 'xx-large')
plt.xticks(np.arange(10), [1, 2, 8, 10, 12, 17, 20, 25, 30, 35], size = 26)
plt.yticks(size = 26)
plt.xlabel('Developmental stage', size = 26, weight = 'semibold')
plt.ylabel('Relative abundance', size = 26, weight = 'semibold')
plt.title('Relative protein abundance\n of each protein in cluster 2, N = 2,697 proteins', size = 20 ,weight = 'semibold')
plt.savefig('colored_Kmeans_cluster_2_protein_data_triceratops_processing_phrog.svg', transparent = True)

plt.figure(figsize = (10, 7))
for i in phrog[phrog['kmeans_cluster_label'] == 3].T:
    plt.plot(np.arange(10), phrog[phrog['protein_id'] == i][[u'126 Sn', u'127N Sn', u'127C Sn',
       u'128N Sn', u'128C Sn', u'129N Sn', u'129C Sn', u'130N Sn', u'130C Sn',
       u'131 Sn']].values.reshape(10, 1), color = '#e7298a', alpha = 0.05)
plt.ylim((-.02, 1.02))
plt.ylim((-.02, 1.02))
plt.legend(loc = 'best', fontsize = 'xx-large')
plt.xticks(np.arange(10), [1, 2, 8, 10, 12, 17, 20, 25, 30, 35], size = 26)
plt.yticks(size = 26)
plt.xlabel('Developmental stage', size = 26, weight = 'semibold')
plt.ylabel('Relative abundance', size = 26, weight = 'semibold')
plt.title('Relative protein abundance\n of each protein in cluster 3, N = 338 proteins', size = 20 ,weight = 'semibold')
plt.savefig('colored_Kmeans_cluster_3_protein_data_triceratops_processing_phrog.svg', transparent = True)

phrog[phrog['kmeans_cluster_label'] == 0]['gene_symbol'].to_csv('cluster_0_human_gene_symbols_for_GO.csv', index = False)
phrog[phrog['kmeans_cluster_label'] == 2]['gene_symbol'].to_csv('cluster_2_human_gene_symbols_for_GO.csv', index = False)
phrog[phrog['kmeans_cluster_label'] == 1]['gene_symbol'].to_csv('cluster_1_human_gene_symbols_for_GO.csv', index = False)
phrog[phrog['kmeans_cluster_label'] == 3]['gene_symbol'].to_csv('cluster_3_human_gene_symbols_for_GO.csv', index = False)
phrog['gene_symbol'].to_csv('all_human_gene_symbols_background.csv', index = False)