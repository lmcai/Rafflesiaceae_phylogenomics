# I. Concatenation

We applied different partition, profile mixture, and site-spesific frequency models to the DNA and protein alignments. All analyses were conducted in IQTREE, branch support was evaluated using 1000 ultrafast bootstrap.

## 1. DNA concatenation analyses

### Analysis 1
Dataset = G2135; Site masking= HmmCleaner threshold 50 + TrimAL remove sites >70% gap; Partition = gene-based, ModelFinder best fitting
```
iqtree2 -s G2135.cds.fas -p G2135.cds.partition -m MFP+MERGE -T AUTO -B 1000
```

### Analysis 2
Dataset = G2135; Site masking= HmmCleaner threshold 50 + TrimAL remove sites >70% gap; Partition = **codon-based**, ModelFinder best fitting
```
iqtree2 -s G2135.cds.fas -p G2135.cds.partition_codon -m MFP+MERGE -T AUTO -B 1000
```
### Analysis 3
Dataset = G2135; Site masking= **No 3rd codon** + HmmCleaner threshold 50 + TrimAL remove sites >70% gap; Partition = codon-based, ModelFinder best fitting
```
iqtree2 -s G2135.cds.fas -p G2135.cds.partition_codonno3rd -m MFP+MERGE -T AUTO -B 1000
```
### Analysis 4
Dataset = G2135; Site masking= No 3rd codon + **HmmCleaner threshold 10** + TrimAL remove sites >70% gap; Partition = codon-based partition, ModelFinder best fitting
```
iqtree2 -s G2135.cds.Hmm10_no3rd.fas -p G2135.cds.Hmm10_no3rd.partition -m MFP+MERGE -T AUTO -B 1000
```
### Analysis 5
Dataset = G2135; Site masking= HmmCleaner threshold 50 + TrimAL remove sites >70% gap; **GHOST heterotachy model with 4 classes**
```
iqtree2 -s G2135.cds.fas -m GTR+FO*H4 -wspm -T AUTO --prefix G2135.cds_codon.GHOST -B 1000
```
### Analysis 6
Dataset = **G434 (monophyletic Rafflesiaceae+Apodanthaceae)**; Site masking= HmmCleaner threshold 50 + TrimAL remove sites >70% gap; Partition = gene-based partition, ModelFinder best fitting
```
iqtree2 -s G434.cds.fas -p G434.cds.partition -m MFP+MERGE -T AUTO -B 1000
```

## 2. Protein concatenation analyses
### Analysis 7
Dataset = G2135; Site masking= HmmCleaner threshold 50 + TrimAL remove sites >70% gap; Partition = gene-based, ModelFinder best fitting
```
iqtree2 -s G2135.aa.mask_trim.fas -p G2135.aa.mask_trim.partition -m MFP+MERGE -T AUTO -B 1000
```
### Analysis 8
Dataset = G2135; Site masking= **HmmCleaner threshold 10** + TrimAL remove sites >70% gap; Partition = gene-based, ModelFinder best fitting
```
iqtree2 -s G2135.aa.Hmm10_trim03.fas -B 1000 -p G2135.aa.Hmm10_trim03.partition -m MFP+MERGE -T AUTO
```
### Analysis 9
Dataset = **G829 (clock-like+>32sp)**; Site masking= HmmCleaner threshold 50 + TrimAL remove sites >70% gap; **site-specific profile mixture model**
```
iqtree2 -s G829.aa.fas -m LG+C30+F+G -ft G2135.aa.mask_trim.partition.treefile -b 100 -T AUTO --prefix G829.aa.PMSF
```
### Analysis 10
Dataset = **G434 (monophyletic Rafflesiaceae+Apodanthaceae)**; Site masking= HmmCleaner threshold 50 + TrimAL remove sites >70% gap; Partition = gene-based partition, ModelFinder best fitting
```
iqtree2 -s G434.aa.fas -p G434.aa.partition -m MFP+MERGE -T AUTO -B 1000
```

# II. Coalescent analyses

Coalescent analyses were conducted on the complete dataset of 2135 genes, inferred from the five combinations of partition and subsitution models described in `2_Alignment_geneTree`.

## 1. ASTRAL

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

## 2. MP-ESP

MP-EST v3.0 was conducted similarly, using five sets of gene trees. A sample MP-EST command is as follows: 
```
mpest -i G2135.GTR_codon.rooted.trees -n 1 -s 432567
```
The branch support was evaluated using the non-paramatric bootstrap gene trees. Where 1 BP gene tree was sampled for 2135 genes, in 100 replications. Then 100 BP species tree was inferred and summarized onto the best species tree to generate BP support. 


