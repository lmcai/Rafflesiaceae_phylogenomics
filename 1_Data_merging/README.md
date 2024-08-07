# Merge phylogenomic data from two studies

We combined the phylogenomic data from Cai et al. 2019 (Widespread ancient whole‐genome duplications in Malpighiales coincide with Eocene global climatic upheaval) and Cai et al. 2021 (The perfect storm: gene tree estimation error, incomplete lineage sorting, and ancient gene flow explain the most recalcitrant ancient angiosperm clade, Malpighiales) with the following steps.

1. Based on the gene IDs of three species commonly sampled in the two studies (Manihot esculenta, Populus trichocarpa, and Jatropha curcas), a preliminary one-2-one ortholog ID mapping was generated. This one-2-one ID mapping is provided in `malpighiales_sapria_orthogroup.list`.

   A total of 5113 genes in WGD paper, 4498 have either Populus or Manihot. 3086/4498 malpighiales orthogroup overlap with 3038/11008 sapria VGT orthogroup. After requiring these orthogroups to have at least one species from Rafflesiaceae and at least 10 total species, a total of 2221 orthogroup meet the criteria (2141 after the application of yang et al, see 2_Alignment_geneTree).

3. To extract amino acid sequences and output to fasta file, use `aa_aln_prep.py`.

To extract DNA sequences and output to fasta file, use `na_aln_prep.py`.
