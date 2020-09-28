#!/bin/bash
#
#SBATCH -n 1                 # Number of cores
#SBATCH -N 1                 # Number of nodes for the cores
#SBATCH -t 3-20:05           # Runtime in D-HH:MM format
#SBATCH -p shared    # Partition to submit to
#SBATCH --mem=20000            # Memory pool for all CPUs
#SBATCH -o lmcai.out      # File to which standard out will be written
#SBATCH -e lmcai.err      # File to which standard err will be written


ID=$1
#ID=${SLURM_ARRAY_TASK_ID}
echo $ID
cd /n/holyscratch01/davis_lab/lmcai/sapria_phylogeny/na_aln

module load mafft/7.407-fasrc01
#if there are less than 800 gene copies, use the accurate linsi; otherwise use mafft
gene_num=$((grep '>' $ID.na.fas| wc -l) 2>&1)
if [ "$gene_num" -lt 500 ]
then
mafft-linsi --adjustdirection $ID.na.fas | sed 's/_R_//g' > $ID.na.aln.fas
else
mafft --adjustdirection $ID.na.fas | sed 's/_R_//g' > $ID.na.aln.fas
fi
echo 'na alignment done...'


#mildly trim the alignment
trimal -in $ID.na.aln.fas -out $ID.na.aln.trimmed.fas -gt 0.10

#fastree
FastTree -nt -gtr < $ID.na.aln.trimmed.fas >$ID.fast.tre

echo 'fasttree done...'

module load gcc/7.1.0-fasrc01 openmpi/3.1.3-fasrc01 iqtree/1.6.10-fasrc02

OUT=$((python ../get_outgroup.py $ID.na.aln.trimmed.fas) 2>&1)
#model selection currentlt cannot run on multiple mpi
mpirun -np 1 iqtree-mpi -s $ID.na.aln.trimmed.fas -o $OUT -pre $ID -nt AUTO -bb 3000 -bnni -alrt 2000 -nm 3000 -redo

rm $ID.bionj
rm $ID.ckp.gz
rm $ID.uniqueseq.phy
rm $ID.mldist

rm $ID.contree
rm $ID.log
rm $ID.model.gz
rm $ID.splits.nex
rm $ID.ufboot
