# Multispecies coalescent model simulation

1. Gene tree topology simulation under the MSC model

We subsampled the ASTRAL species tree from Analysis 15 to include only five species plus an outgroup. The branch legnth of this species tree was optimized in MPEST first with theh following commands:
```
mpest -i G2135.MFP.rooted.trees -u analysis15.tre -B
```
We then used the function ‘sim.coal.mpest’ in the R package Phybase v2.0 (Liu and Yu 2010) to simulate 10,000 gene trees using the reference species tree. 
```
#R
library(Phybase)
number_geneTr=10000
SpTr_text=readLines("analysis15.mpest.tre)
sim_geneTrees=sim.coal.mpest(SpTr_text,number_geneTr)

```

2. Modeling gene tree branch lengths in mutation units

We first extract the empirical branch length from individual gene trees estimated from IQTREE using the python script `modify_MSC_sim_tree_branch_length.py`. This resulted in `*.brlen.txt` files containing empirical terminal branch length for each clade and a `int.brlen.txt` file for all internal branches.

Then we fit a gaussian distribution to each branch length using the `fitdistr` function from the R package MASS. 

Finally, we randomly sample from these gaussian distributions to assign branch length to each gene tree simulated from the MSC model above, using the script `modify_MSC_sim_tree_branch_length.py`.

3. Simulating DNA alignments

With the GTR model parameters inferred from the empirical G2135 dataset in IQTREE, we used the following command in seqgen to simulate DNA alignments of different length (500, 1000, 1500, and 2000 bp) based on individual gene trees:
```
#generate DNA alignments under the GTR model
seq-gen -m GTR -l 500 -a 0.5 -r 1.56 3.91 1.31 1.38 5.01 1.00 -f 0.28 0.20 0.25 0.27 -of [input tree file] > [output]
```

4. After obtaining the simulated DNA alignment, we reconstruct individual gene trees and species tree using both coalescent and concatenation method to evaluate the effects of long branch attraction and incomplete lineage sorting.
