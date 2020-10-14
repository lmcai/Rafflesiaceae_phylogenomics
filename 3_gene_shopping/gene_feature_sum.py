from ete3 import Tree
from Bio import SeqIO
import os
from numpy import median

gene_tree_dir='/n/holyscratch01/davis_lab/lmcai/sapria_phylogeny/8_GGI/G2141_1to1_trees/'
genes=os.listdir(gene_tree_dir)

gene_feature={}
for f in genes:
	gene=f.split('.')[0]
	rec={}
	t=Tree(gene_tree_dir+f)
	tips=[leaf.name for leaf in t]
	g1=[i for i in tips if i.startswith(('Ricinus','Hevea','Manihot','Endospermum','Jatropha'))]
	g2=[i for i in tips if i.startswith('Clutia')]
	g3=[i for i in tips if i.startswith('Drypetes')]
	g4=[i for i in tips if i.startswith('Galearia')]
	g5=[i for i in tips if i.startswith(('Erythroxylum','Rhizophora'))]
	g6=[i for i in tips if i.startswith('Ixonanthes')]
	#g7=[i for i in tips if i.startswith('Passiflora')]
	rec['sp_num']=str(len(tips))
	rec['Pandaceae']='N'
	if len(g4)>0:
		rec['Pandaceae']='Y'
	rec['focal_clades']=str(len([i for i in [g1,g2,g3,g4,g5] if len(i)>0]))
	rec['Ixonanthaceae']='N'
	if len(g6)>0:
		rec['Ixonanthaceae']='Y'
	aa_seqs=SeqIO.parse('../2_aa_aln_geneTr/'+gene+'.aa.fas','fasta')
	aa_seq_lens=[]
	for i in aa_seqs:
		aa_seq_lens.append(len(i.seq.ungap("-")))
	rec['aa_length']=str(median(aa_seq_lens))
	gene_feature[gene]=rec

#tsv header
out=open('G2141_features.tsv','a')
out.write('\t'.join(['Gene']+[i for i in gene_feature['0'].keys()])+'\n')
for k in gene_feature.keys():
	out.write(k+'\t'+'\t'.join([gene_feature[k][i] for i in gene_feature[k].keys()])+'\n')

out.close()