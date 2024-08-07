# Sequence alignment, gene tree inference, data cleaning

Sequence alignment, gene tree inference, and data cleaning was conducted with the following steps:

1. First round alignment + tree inference
  
   Use the mafft-linsi algorithm to align DNA sequences, trim with trimAL (-gt 0.10), remove outgroup sequences (non-Malpighiales) with `get_outgroup.py`, then infer a preliminary maximum likelihood gene trees with IQTREE (1000 ultrafast bootstrap replication and optimal models determined by ModelFinder).

The bash script to run these analyses altogether for each gene is provided in `mafft_trimal_fasttree_iqtree.sh`.

2. Use the pipeline from Yang and Smith (2016) to remove paralogs and create one-to-one ortholog for the cleaned dataset with at least 10 species.

```
python yangya-phylogenomic_dataset_construction-489685700c2a/prune_paralogs_RT.py test/ .treefile na_tree_yang_pruned/ 10 malp.taxa
```
This resulted in 2141 orthogroups for final phylogenetic investigation.

3. Second round alignment + tree inference

   Use the mafft-linsi algorithm to align **RAW** AA sequences (untrimmed), back translate to DNA codon alignment, then infer a final maximum likelihood gene trees with IQTREE (3000 ultrafast bootstrap replication and optimal models determined by ModelFinder).
```
mafft --localpair --maxiterate 1000 aa.fas > aa.aln.fas
pal2nal.pl aa.aln.fas na.fas -output fasta > na.aln.fas

```
