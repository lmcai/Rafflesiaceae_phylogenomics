#get chi-square test of compositional bias in the alignments
#use iqtree to output alignment information
#iqtree2 -s 4_na_aln_1to1/0.aln -m GTR+G -n 0 -alninfo

import os
from ete3 import Tree
from numpy import mean
import ete3

files=os.listdir('.')
files=[i for i in files if i.endswith('.log')]
out=open('na_aln_chisquare_compositional_bias.tsv','w')
out.write('GeneID\tSp_num\taln_length\taverage_node_support\tRafflesiaceae_best_compositional_homogeneity_pvalue\tApodanthaceae_best_compositional_homogeneity_pvalue\n')
for i in files:
	#print(i)
	x=open(i).readlines()
	#aln_len=x[15].split()[5]
	num_sp=x[13].split()[2]
	aln_len=x[13].split()[5]
	node_bp=[]
	#t=Tree('../3_na_tree_1to1_yang_subtree/'+i.split('.')[0]+'.inclade1.ortho1.tre',format=1)
	try:t=Tree(i.split('.')[0]+'.na.mask.fas.treefile',format=1)
	except ete3.parser.newick.NewickError:
		print(i)
		continue
	for node in t.traverse("postorder"):
		if not node.is_leaf():
			try:node_bp.append(int(node.name))
			except ValueError:pass
			#try:node_bp.append(int(node.name.split('_')[1]))
			#except IndexError:pass
	#j=18
	j=17
	raff_cur_best=0
	apo_cur_best=0
	while True:
		ll=x[j].split()
		if ll[0]=='WARNING:' or ll[0]=='****':break
		if ll[1].startswith(('Sap','Rca','Rtu','Rhi')):
			if float(ll[-1][:-1])>raff_cur_best:
				raff_cur_best=float(ll[-1][:-1])
		elif ll[1].startswith(('Apodan','Pilo')):
			if float(ll[-1][:-1])>raff_cur_best:
				apo_cur_best=float(ll[-1][:-1])
		j=j+1
	sp=[node.name for node in t]
	if len([k for k in sp if k.startswith(('Apodan','Pilo'))])==0:apo_cur_best='NA'
	d=out.write(i.split('.')[0]+'\t'+num_sp+'\t'+aln_len+'\t'+str(mean(node_bp))+'\t'+str(raff_cur_best)+'\t'+str(apo_cur_best)+'\n')

out.close()