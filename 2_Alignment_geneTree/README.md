# Sequence alignment, gene tree inference, data cleaning

Sequence alignment, gene tree inference, and data cleaning was conducted with the following steps:

1. Use the mafft-linsi algorithm to align AA sequence, then convert to DNA codon alignment with pal2nal, trim with trimAL (-gt 0.10), remove outgroup sequences (non-Malpighiales) with `get_outgroup.py`, then infer maximum likelihood gene trees with IQTREE (3000 ultrafast bootstrap replication and optimal models determined by ModelFinder).

The bash script to run these analyses altogether for each gene is provided in `mafft_trimal_fasttree_iqtree.sh`.
