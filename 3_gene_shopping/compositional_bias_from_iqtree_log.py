#get chi-square test of compositional bias in the alignments
#use iqtree to output alignment information
#iqtree2 -s 4_na_aln_1to1/0.aln -m GTR+G -n 0 -alninfo

import os

files=os.listdir('.')
files=[i for i in files if i.endswith('.log')]
out=open('na_aln_chisquare_compositional_bias.tsv')
out.write('GeneID\tRafflesiaceae_best_compositional_homogeneity_pvalue\n')
for i in files:
	x=open(i).readlines()
	j=18
	while True:
		cur_best=0
		ll=x[j].split()
		if ll[0]=='WARNING:':break
		if ll[1] in ['Sap','Rca','Rtu','Rhi']:
			if float(ll.split()[-1][:-1])>cur_best:
				cur_best=float(ll.split()[-1][:-1])
	out.write(i.split('.')[0]+'\t'+str(cur_best)+'\n')

out.close()