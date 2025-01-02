# Merge phylogenomic data from three studies

We combined the phylogenomic data from Cai et al. 2019 (Widespread ancient whole‐genome duplications in Malpighiales coincide with Eocene global climatic upheaval), Cai et al. 2021 (The perfect storm: gene tree estimation error, incomplete lineage sorting, and ancient gene flow explain the most recalcitrant ancient angiosperm clade, Malpighiales), and a newly generated Apodanthaceae RNA sequencing data with the following steps.

1. Based on the gene IDs of three species commonly sampled in Cai et al. 2019 and Cai et al. 2021 (Manihot esculenta, Populus trichocarpa, and Jatropha curcas), a preliminary one-2-one ortholog ID mapping was generated. This one-2-one ID mapping is provided in `malpighiales_sapria_orthogroup.list`.

2. A among the 5113 genes in the WGD paper, 4498 have either Populus or Manihot sequences. 3086/4498 orthogroups in the WGD paper overlap with 3038/11008 sapria VGT orthogroup. After requiring these orthogroups to have at least one species from Rafflesiaceae and at least 10 total species, a total of 2221 orthogroup were retained (2135 after the application of yang et al, see 2_Alignment_geneTree).

3. To extract amino acid sequences and output to fasta file, use `aa_aln_prep.py`. To extract DNA sequences and output to fasta file, use `na_aln_prep.py`.

4. To add orthologous Apodanthaceae sequences (Apodanthes and Pilostyles) to this dataset, we used HMMER to identify potential orthologs.

   An HMM profile as generated for each alignment using
   ```
   hmmbuild 1.hmm 1.aln.fas
   ```
   To search ortholog in Apodanthaceae, use the following commanad:
   ```
   hmmsearch --domtblout 1.tabfile -o 1.hmmsearch 1.hmm Apodanthes_CDS.fa
   ```
