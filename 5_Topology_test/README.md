# gCF and sCF characterization

We used the species tree, gene trees, and DNA alignments from the coalescent-based Analysis 15 (table S3) to assess gCF and sCF. 

For gCF:
```
iqtree2 -te G2141.wastral.hybrid.tre --gcf G2141.MFP.trees --prefix G2135.gCF
```
This result in `G2135.gCF.cf.tree.nex`

For sCF with the complete DNA alignment:
```
iqtree2 -te ../coalescent/G2141.wastral.hybrid.tre --scfl 200 -s G2135.cds.fas --prefix G2135.cds.sCF
```
This result in `G2135.cds.sCF.cf.tree.nex`

For sCF with the complete protein alignment:
```
iqtree2 -te ../coalescent/G2141.wastral.hybrid.tre --scfl 200 -s G2135.aa.fas --prefix G2135.aa.sCF
```
This result in `G2135.aa.sCF.cf.tree.nex`


# Gene- and site-wise support for alternative topologies

We assessed the gene and site-specific likelihood (LL) for eight alternative placements of Rafflesiaceae+Apodanthaceae with the following steps.

1. Optimize branch length with constrained topology in IQTREE using the concatenated alignment. The eight alternative placements are indicated in H*.ref.tre
```
#For DNA
iqtree2 -s G2135.cds.fas -redo -p G2135.cds.partition -m MFP -g H${SLURM_ARRAY_TASK_ID}.ref.tre --prefix cds.H${SLURM_ARRAY_TASK_ID} -T AUTO
#For protein
iqtree2 -s G2135.aa.mask_trim.fas -redo -p G2135.aa.mask_trim.partition -m MFP -g H${SLURM_ARRAY_TASK_ID}.ref.tre --prefix aa.H${SLURM_ARRAY_TASK_ID} -T AUTO
```
2. Calculate per partition (gene) LL using the `-wpl` flag:
```
# For DNA
cat cds.H*.treefile >cds.tree_candidates.trees 
iqtree2 -s G2135.cds.fas -p RateperS_treesearch.best_model.nex -wpl --prefix cds.LLperG -z cds.tree_candidates.trees -n 0 -redo

#For protein
cat aa.H*.treefile >aa.tree_candidates.trees 
iqtree2 -s G2135.aa.mask_trim.fas -p RateperS_treesearch.best_model.nex -wpl --prefix aa.LLperG -z aa.tree_candidates.trees -n 0 -redo
```
In the log files of these runs, there will be a line starting with "Partition-specific rates:  1.067 1.019 1.569 1.175 0.778...". These data were used as a measurement of rate for each gene.

3. Calculate per site rate and LL:

This is only done with DNA. We first estimate site-specific rate using the `-wsr` flag
```
iqtree2 -s G2135.cds.fas -p G2135.cds.partition -wsr --prefix RateperS_treesearch
```
Then we calculate per site LL for each of the eight alternative topologies:
```
iqtree2 -s G2135.cds.fas -p RateperS_treesearch.best_model.nex -wsl --prefix cds.LLperS -z csd.tree_candidates.trees -n 0 -redo
```

# AU test

We used AU tests to statisticlly evaluate whether a phylogenetic hypothesis can be accepted or rejected. The process of which largely follows the gene- and site-wise likelihood calculation. Briefly, we optimize the branch length for alternative topologies (e.g., H1, H2, H3...) using corresponding alignments and then concatenate all trees into a single file and use the following command in IQTREE to perform AU test:
```
cat codon1.H1.tre codon1.H2.tre codon1.H3.tre > codon1.all.trees

iqtree -s codon1.fas -p codon1.best_model.nex -z codon1.all.trees -n 0 -zb 10000 -au
```
