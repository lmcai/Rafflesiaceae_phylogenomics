#get chi-square test of compositional bias in the alignments
#use iqtree to output alignment information
#iqtree2 -s 4_na_aln_1to1/0.aln -m GTR+G -n 0 -alninfo

import os

files=os.listdir('.')
files=[i for i in files if i.endswith('.log')]
out=open('na_aln_chisquare_compositional_bias.tsv','a')
out.write('GeneID\tRafflesiaceae_best_compositional_homogeneity_pvalue\n')
for i in files:
	x=open(i).readlines()
	j=18
	cur_best=0
	while True:
		ll=x[j].split()
		if ll[0]=='WARNING:' or ll[0]=='****':break
		if ll[1] in ['Sap','Rca','Rtu','Rhi']:
			if float(ll[-1][:-1])>cur_best:
				print j
				cur_best=float(ll[-1][:-1])
				print cur_best
		j=j+1
	out.write(i.split('.')[0]+'\t'+str(cur_best)+'\n')

out.close()