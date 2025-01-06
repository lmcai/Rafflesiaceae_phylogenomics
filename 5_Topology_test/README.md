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
