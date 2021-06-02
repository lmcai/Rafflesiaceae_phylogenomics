#get branch length distribution from empirical gene trees
from ete3 import Tree
import os
import ete3
genes=[i for i in os.listdir('.') if i.endswith('.treefile')]
internal=[]
cro=[]
gal=[]
clu=[]
sap=[]
rhi=[]
dry=[]

for g in genes:
	t=Tree(g)
	sp2keep=[i.name for i in t if i.name in ['Sap','Crossopetalum','Clutia','Galearia','Rhizophora','Drypetes']]
	if len(sp2keep)>2:
		t.prune(sp2keep,preserve_branch_length=True)
		for n in t.traverse("postorder"):
			if not n.is_leaf() and not n.is_root():internal.append(n.dist)
			elif n.is_leaf():
				if n.name=='Crossopetalum':cro.append(n.dist)
				elif n.name=='Galearia':gal.append(n.dist)
				elif n.name=='Clutia':clu.append(n.dist)
				elif n.name=='Sap':sap.append(n.dist)
				elif n.name=='Rhizophora':rhi.append(n.dist)
				elif n.name=='Drypetes':dry.append(n.dist)

out=open('int.brlen.txt','a')
out.write('\n'.join([str(i) for i in internal]))
out.close()

out=open('cro.brlen.txt','a')
out.write('\n'.join([str(i) for i in cro]))
out.close()

out=open('gal.brlen.txt','a')
out.write('\n'.join([str(i) for i in gal]))
out.close()

out=open('clu.brlen.txt','a')
out.write('\n'.join([str(i) for i in clu]))
out.close()

out=open('sap.brlen.txt','a')
out.write('\n'.join([str(i) for i in sap]))
out.close()

out=open('rhi.brlen.txt','a')
out.write('\n'.join([str(i) for i in rhi]))
out.close()

out=open('dry.brlen.txt','a')
out.write('\n'.join([str(i) for i in dry]))
out.close()

from numpy import random
out=open('G2141.mpest.5clades.10000sim.mubrlen.trees','a')
x=open('../phybase_MSC_sim/G2141.mpest.5clades.10000sim.trees').readlines()
for l in x:
	t=Tree(l)
	for n in t.traverse("postorder"):
		if not n.is_leaf() and not n.is_root():
			n.dist=random.uniform(0.002,0.006,1)[0]
		elif n.is_leaf():
			if n.name=='Crossopetalum':n.dist=random.normal(0.261,0.214,1)[0]
			elif n.name=='Galearia':n.dist=random.normal(0.0949,0.0443,1)[0]
			elif n.name=='Clutia':n.dist=random.normal(0.186,0.0628,1)[0]
			elif n.name=='Sap':n.dist=random.normal(0.475,0.129,1)[0]
			elif n.name=='Rhizophora':n.dist=random.normal(0.162,0.0489,1)[0]
			elif n.name=='Drypetes':n.dist=random.normal(0.117,0.0429,1)[0]
	for n in t.get_children():
		if not n.is_leaf():n.dist=random.uniform(0.04,0.06,1)[0]
	out.write(t.write()+'\n')
		
#generate DNA alignments under the GTR model
#seq-gen -m GTR -l 500 -a 0.5 -r 1.56 3.91 1.31 1.38 5.01 1.00 -f 0.28 0.20 0.25 0.27 -of [input tree file] > [output]

#calculate sCF
#(((((Clutia,Drypetes),Galearia),Sap),Rhizophora),Crossopetalum)
#iqtree2 -t seqgen.ref.tre -s concatenated_aln.nex --scf 100 --prefix seqgen.sCF