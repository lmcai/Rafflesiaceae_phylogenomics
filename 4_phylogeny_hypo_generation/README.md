# Coalescent analyses

Coalescent analyses were conducted on the complete dataset of 2135 genes, inferred from the five combinations of partition and subsitution models described in `2_Alignment_geneTree`.

## I. ASTRAL

ASTRAL IV offers several modes to weigh branch length and branch support, we tested all of these modes in an initial run:
```
ASTER-Linux/bin/astral4 -i G2141.aster.genetree.trees -o G2141.aster4.tre
ASTER-Linux/bin/wastral -i G2141.aster.genetree.trees -o G2141.wastral.support.tre -r 16 -s 16 -x 100 -n 0 --mode 2
ASTER-Linux/bin/wastral -i G2141.aster.genetree.trees -o G2141.wastral.brlen.tre -r 16 -s 16 -x 100 -n 0 --mode 3
ASTER-Linux/bin/wastral -i G2141.aster.genetree.trees -o G2141.wastral.hybrid.tre -r 16 -s 16 -x 100 -n 0
#full node annotation and polytomy test
java -jar ../../Astral/astral.5.7.8.jar -q G2141.aster4.tre -i G2141.aster.genetree.trees -t 2 -o G2141.astral_anno.tre
java -jar ../../Astral/astral.5.7.8.jar -q G2141.aster4.tre -i G2141.aster.genetree.trees -t 10 -o G2141.astral_polytomy.tre
```

The resulting species trees were similar to each other. We decided to use the ASTRAL4 algorithm without any branch collapsing.

The commands we used for final analyses were as follows:
```
ASTER-Linux/bin/astral4 -i G2135.GTR_codon.trees -o G2135.GTR_codon.aster4.tre
ASTER-Linux/bin/astral4 -i G2135.MFP_codon.trees -o G2135.MFP_codon.aster4.tre
ASTER-Linux/bin/astral4 -i G2135.msetGTR_codon.trees -o G2135.msetGTR_codon.aster4.tre
ASTER-Linux/bin/astral4 -i G2135.msetGTR.trees -o G2135.msetGTR.aster4.tre
ASTER-Linux/bin/astral4 -i G2135.MFP.trees -o G2135.MFP.aster4.tre
ASTER-Linux/bin/astral4 -i G434.mono.trees -o G434.mono.aster4.tre
```
This will output a species tree and branch support was evaluated with local posterial probability.

## MP-ESP

MP-EST v3.0 was conducted similarly, using five sets of gene trees. A sample MP-EST command is as follows: 
```
mpest -i G2135.GTR_codon.rooted.trees -n 1 -s 432567
```
The branch support was evaluated using the non-paramatric bootstrap gene trees. Where 1 BP gene tree was sampled for 2135 genes, in 100 replications. Then 100 BP species tree was inferred and summarized onto the best species tree to generate BP support. 

# Concatenation

We applied different partition, profile mixture, and site-spesific frequency models to the DNA and protein alignments.
