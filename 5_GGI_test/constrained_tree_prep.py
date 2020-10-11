from ete3 import Tree
import sys

h1=Tree('../H1.ref.tre')
h2=Tree('../H2.ref.tre')
h3=Tree('../H3.ref.tre')
h4=Tree('../H4.ref.tre')
h5=Tree('../H5.ref.tre')
h6=Tree('../H6.ref.tre')
h7=Tree('../H7.ref.tre')

total_sp_name=[leaf.name for leaf in h1]

sp_num_out=open('../sp_num.list','a')

def generate_constrained_tree(ID):
	t=Tree('../G2141_1to1_trees/'+ID+'.inclade1.ortho1.tre')
	sp_names=[leaf.name for leaf in t]
	d=sp_num_out.write(ID+'\t'+str(len(sp_names))+'\n')
	#define groups
	g1=[i for i in sp_names if i.startswith(('Ricinus','Hevea','Manihot','Endospermum','Jatropha'))]
	g2=[i for i in sp_names if i.startswith('Clutia')]
	g3=[i for i in sp_names if i.startswith('Drypetes')]
	g4=[i for i in sp_names if i.startswith('Galearia')]
	g5=[i for i in sp_names if i.startswith(('Erythroxylum','Rhizophora'))]
	g6=[i for i in sp_names if i.startswith('Ixonanthes')]
	g7=[i for i in sp_names if i.startswith('Passiflora')]
	#sp2prune=[i for i in total_sp_name if not i in sp_names]
	if len(g1)>0:
		h1.prune(sp_names)
		h1.write(outfile=ID+'.H1.tre',format=9)
		if len(g2)>0:
			h2.prune(sp_names)
			h2.write(outfile=ID+'.H2.tre',format=9)
			if len(g3)>0:
				h3.prune(sp_names)
				h3.write(outfile=ID+'.H3.tre',format=9)
		if len(g6)>0:
			h6.prune(sp_names)
			h6.write(outfile=ID+'.H6.tre',format=9)
		if len(g7)>0:
			h7.prune(sp_names)
			h7.write(outfile=ID+'.H7.tre',format=9)


generate_constrained_tree(sys.argv[1])
