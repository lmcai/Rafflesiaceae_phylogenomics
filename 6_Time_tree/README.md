# Divergence time estimation

The divergence time of Malpighiales, Rafflesiaceae and Apodanthaceae in particular, was estimated based on the penalized likelihood method implemented in treePL and the under the Bayesian framework implemented in PhyloBayes v4.1.

## TreePL analysis

1. The input phylogeny has the topology estimated from the coalescent method using 2135 genes with ASTRAL4 (analysis 15 in Table S3). We then use this fixed topology to infer branch length in mutational units using concatenated DNA sequence filtered by HmmCleaner with a threshold of 10, 3rd codon removed, and sites with >70% gap removed (analysis 5 in Table S3). The IQTREE command used for branch length estimation is as follows:

```
iqtree2 -s G2135.cds.Hmm10_no3rd.fas -p G2135.cds.Hmm10_no3rd.partition.best_scheme.nex -te analysis15.astral4.MFP.tre -T 4 --prefix analysis15_topo.analysis5_aln
```

2. This step resulted in `analysis15_topo.analysis5_aln.tre`, which was used to run treePL with the control file `treePL.in.txt`. treePL was run three times to prime, cross-validate, and search thoroughly for an optimum tree. The best smooth value was determined to be 0.01 based on the chi-square test.

The final ultrametric tree is `malpig_treePL.tre`.


## Phylobayes analysis

Phylobayes has the advantage of both Bayesian analysis and complex mixture models to better model substitutional processes. It is also computationally intense, and we thus subsampled the input data to select the most clock-like genes

1. Rank the genes based on clock-likeness using Sortadata. The input includes rooted gene trees and one rooted species tree. The scripts are provided in `sortadate.sh`.

```
######################
#Get the root-to-tip variance
python ~/programs/SortaDate/src/get_var_length.py rooted_gene_trees/ --flend .tre --outg Crossopetalum,Oxalis,Elaeocarpus --outf gene_tree_root2tip_var.tsv

#There wu=ill be error message: 
#Error: no matching tip labels. Returning original tree.
#Error: this really only works with nexus or newick. Exiting.
#/bin/sh: line 1: 344793 Segmentation fault      (core dumped) pxrmt -t gene_trees_sortdata//2151.inclade1.ortho1.tre -n Crossopetalum,Oxalis,Elaeocarpus
#     344794                       (core dumped) | pxlstr

#This may be due to missing data in the gene trees

######################
#Get the bipartition support

python ~/programs/SortaDate/src/get_bp_genetrees.py rooted_gene_trees/ ../23_add_apo_na_aln_geneTr/sp.tre --flend .tre --outf gene_tree_bp.tsv

######################
#Combine the results from these two runs

python ~/programs/SortaDate/src/combine_results.py gene_tree_root2tip_var.tsv gene_tree_bp.tsv --outf gene_tree_combined_statistics.tsv

#######################
#Sort and get the list of the good genes

#remove rows without four columns
awk 'NF==4{print}{}' gene_tree_combined_statistics.tsv >gene_tree_combined_statistics_filtered.tsv

python ~/programs/SortaDate/src/get_good_genes.py gene_tree_combined_statistics.tsv --max 500 --order 1,2,3 --outf best_gene.list
```

2. Based on the result of SortaDate, select ## genes meeting the following criteria (1) All genes with >31 species, containing both Rafflesiaceae and Apodanthaceae, and (2) ranked in the top 50% based on clock-likeness. 

This resulted in 829 genes listed in `G829.list`.

3. Divergence time estimation in PhyloBayes v4.1.

With a fixed topology (same as the treePL topology above, required by PhyloBayes), we used the 829-gene DNA alignment, uniformly distributed fossil calibration points, and CAT-GTR mixture model to infer divergence time in PhyloBayes. The following command was used:

```
~/programs/phylobayes/bin/pb -d G829.sortadate.no3rd.hmm10.gt03.fas.phy -T analysis15_topo.analysis5_aln.treefile -
r malp.outgroup -cal malp.clib -ln -gtr -cat -f malp.pb1
```

To summarize the result with confidence intervals with first 1000 (25%) trees as burn-in and sampling every 10 trees, use the following command:
```
~/programs/phylobayes/bin/readdiv -x 1000 5 malp.pb1
```
The resulting time tree with 95% HPD interval is `malp3_sample.chronogram`
