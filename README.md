# maxquant_bioinformatics
Analysis of TMT-MS3 proteomics data with MaxQuant

missed_cleavage_rate.py , tmt_labeling_efficiency.py and percent_cysteines.py are for use with the modificationSpecificPeptides.txt output file from MaxQuant from the sample quality control step

map_gene_symbol.py and k_means_clustering.py are for use with the proteinGroups.txt output file from MaxQuant 

USAGE:
Use of these scripts assumed you have installed and run MaxQuant on your desktop. If you have installed it elsewhere please modify the directory accordingly.

tmt_labeling_efficiency.py

Download the script and place it in Desktop\MaxQuant_1.6.0.16\MaxQuant\combined\txt
Open the terminal by simultaneously pressing windows button + R
type cmd and press enter
type cd Desktop\MaxQuant_1.6.0.16\MaxQuant\combined\txt and press enter
type python tmt_labeling_efficiency.py and press enter

missed_cleavage_rate.py
Download the script and place it in Desktop\MaxQuant_1.6.0.16\MaxQuant\combined\txt
Open the terminal by simultaneously pressing windows button + R
type cmd and press enter
type cd Desktop\MaxQuant_1.6.0.16\MaxQuant\combined\txt and press enter
type python missed_cleavage_rate.py --protease <protease_name> where <protease_name> is lysc or trypsin
eg: python missed_cleavage_rate.py --protease lysc
eg: python missed_cleavage_rate.py --protease trypsin

percent_cysteines.py
Download the script and place it in Desktop\MaxQuant_1.6.0.16\MaxQuant\combined\txt
Open the terminal by simultaneously pressing windows button + R
type cmd and press enter
type cd Desktop\MaxQuant_1.6.0.16\MaxQuant\combined\txt and press enter
type python percent_cysteines.py and press enter

map_gene_symbol.py
Download PHROG_annotation.csv from the PHROG_annotation.zip file on https://scholar.princeton.edu/wuehr/sample_prep
and place into the txt MaxQuant folder (eg. Desktop\MaxQuant_1.6.0.16\MaxQuant\24fracs\combined\txt)
Download map_gene_symbol.py and place it in Desktop\MaxQuant_1.6.0.16\MaxQuant\combined\txt
Open the terminal by simultaneously pressing windows button + R
type cmd and press enter
type cd Desktop\MaxQuant_1.6.0.16\MaxQuant\combined\txt and press enter
type python map_gene_symbol.py and press enter

k_means_clustering.py
Download the script and place it in the appropriate txt MaxQuant folder as above
This must be in the same folder as where map_gene_symbol.py was placed and executed
Open the terminal by simultaneously pressing windows button + R
type cmd and press enter
type (eg.) cd Desktop\MaxQuant_1.6.0.16\MaxQuant\combined\txt and press enter
type python k_means_clustering.py --num_clusters <number> --line_opacity <opacity> --folder_name <name>
  where you replace <number> with an integer such as 4 which will be the number of clusters used by the k means clustering algorithm
  replace <opacity> with a float such as .01 to tune the opacity of the individual lines of the protein plots generated
  replace <name> with the desired folder name you wish to create
  
  eg. python k_means_clustering.py --num_clusters 4 --line_opacity .01 --folder_name clustering_results
  
  This will output all results into the newly created folder clustering_results
