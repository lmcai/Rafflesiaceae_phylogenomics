# Sequence alignment, gene tree inference, data cleaning

Sequence data cleaning, alignment, and gene tree inference were conducted with the following steps:

1. First-round alignment + tree inference
  
   Use the mafft-linsi algorithm to align DNA sequences, trim with trimAL (-gt 0.10), remove outgroup sequences (non-Malpighiales) with `get_outgroup.py`, then infer a preliminary maximum likelihood gene trees with IQTREE (1000 ultrafast bootstrap replication and optimal models determined by ModelFinder).

The bash script to run these analyses altogether for each gene is provided in `mafft_trimal_fasttree_iqtree.sh`.

2. Use the pipeline from Yang and Smith (2016) to remove paralogs and create one-to-one ortholog for the cleaned dataset with at least 10 species.

```
python yangya-phylogenomic_dataset_construction-489685700c2a/prune_paralogs_RT.py test/ .treefile na_tree_yang_pruned/ 10 malp.taxa
```
This resulted in 2141 orthogroups for final phylogenetic investigation. Newly added orthologous Apodanthaceae for each orthogroup is included in `orthologous_apodanthaceae.tsv`.

3. Second round alignment

   Use the mafft-linsi algorithm to align **RAW** AA sequences (untrimmed), back translate to DNA codon alignment
```
mafft --localpair --maxiterate 1000 aa.fas > aa.aln.fas
pal2nal.pl aa.aln.fas na.fas -output fasta > na.aln.fas

```

4. Alignment masking with HmmCleaner

Use the `HmmCleaner.pl` from [MACSE_ALFIX_V01](https://github.com/ranwez/MACSE_V2_PIPELINES/tree/master) to remove non-homologous regions in the alignment. I have difficulty installing HmmCleaner from cpan or source, but it was contained in the MACSE pipeline. A threshold of 50 was used after trying this parameter from 10 to 50. HmmCleaner directly deletes the regions, which will break the codon position. So I write a custom Python script `hmmcleaner_codon_aware_masking.py` to take in the log file from HmmCleaner and mask the corresponding position in the alignment.

```
perl MACSE_ALFIX_V01/HMMcleanerV1_8_VR2/HMMcleanNuc_VR.pl 2675.na.aln.fas 50
python hmmcleaner_codon_aware_masking.py 2675.na.aln.fas 2675.na.aln_Hmm50.log 

```
or for protein alignments
```
perl MACSE_ALFIX_V01/HMMcleanerV1_8_VR2/HMMcleanNuc_VR.pl 2675.aa.aln.fas 50
python hmmcleaner_codon_aware_masking.py 2675.aa.aln.fas 2675.aa.aln_Hmm50.log 
```

This will generate *.masked.fas for each fasta file, with all sites flagged by HmmCleaner masked ('N' for DNA and '?' for protein).

5. Infer a final maximum likelihood gene trees with IQTREE (1000 ultrafast bootstrap replication, -bnni to reduce the risk of overestimating branch supports with UFBoot, and optimal models determined by ModelFinder).
```
iqtree2 -s $ID.na.mask.fas -o $OUTGROUP -st CODON -T AUTO -B 1000 -bnni -redo
```
