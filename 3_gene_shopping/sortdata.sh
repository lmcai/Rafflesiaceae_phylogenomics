#path: /n/holyscratch01/davis_lab/lmcai/sapria_phylogeny/sortdata

######################
#Get the root-to-tip variance
python ~/programs/SortaDate/src/get_var_length.py gene_trees_sortdata/ --flend .inclade1.ortho1.tre --outg Crossopetalum,Oxalis,Elaeocarpus --outf gene_tree_root2tip_var.tsv

#There wu=ill be error message: 
#Error: no matching tip labels. Returning original tree.
#Error: this really only works with nexus or newick. Exiting.
#/bin/sh: line 1: 344793 Segmentation fault      (core dumped) pxrmt -t gene_trees_sortdata//2151.inclade1.ortho1.tre -n Crossopetalum,Oxalis,Elaeocarpus
#     344794                       (core dumped) | pxlstr

#This may due to missing data in the gene trees

######################
#Get the bipartition support

python ~/programs/SortaDate/src/get_bp_genetrees.py gene_trees_sortdata/ G2999.astral.rooted.tre --flend .inclade1.ortho1.tre --outf gene_tree_bp.tsv

######################
#Combine the results from these two runs

python ~/programs/SortaDate/src/combine_results.py gene_tree_root2tip_var.tsv gene_tree_bp.tsv --outf gene_tree_combined_statistics.tsv