# Merge phylogenomic data from two studies

We combined the phylogenomic data from Cai et al. 2019 (Widespread ancient whole‐genome duplications in Malpighiales coincide with Eocene global climatic upheaval) and Cai et al. 2021 (The perfect storm: gene tree estimation error, incomplete lineage sorting, and ancient gene flow explain the most recalcitrant ancient angiosperm clade, Malpighiales) with the following steps.

1. Based on the gene IDs of three species commonly sampled in the two studies (Manihot esculenta, Populus trichocarpa, and Jatropha curcas), a preliminary one-2-one ortholog mapping was generated. This one-2-one mapping is provided in `malpighiales_sapria_orthogroup.list`.

2. To extract amino acid sequences and output to fasta file, use `aa_aln_prep.py`.

To extract DNA sequences and output to fasta file, use `na_aln_prep.py`.
