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

1. Rank the genes based on clock-likeness using Sortadata.