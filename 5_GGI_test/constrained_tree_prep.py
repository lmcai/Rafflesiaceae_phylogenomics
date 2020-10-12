from ete3 import Tree
import sys

h1=Tree('../H1.ref.tre')
h2=Tree('../H2.ref.tre')
h4=Tree('../H4.ref.tre')
h5=Tree('../H5.ref.tre')
h7=Tree('../H7.ref.tre')
h9=Tree('../H9.ref.tre')
h12=Tree('../H12.ref.tre')
h15=Tree('../H15.ref.tre')
h16=Tree('../H16.ref.tre')
h17=Tree('../H17.ref.tre')

h3=Tree('../H3.ref.tre')
h6=Tree('../H6.ref.tre')
h8=Tree('../H8.ref.tre')
h10=Tree('../H10.ref.tre')
h11=Tree('../H11.ref.tre')
h13=Tree('../H13.ref.tre')
h14=Tree('../H14.ref.tre')

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
	if len(g1)>0 and len(g2)>0 and len(g3)>0 and len(g5)>0:
		h1.prune(sp_names)
		h1.write(outfile=ID+'.H1.tre',format=9)
		h2.prune(sp_names)
		h2.write(outfile=ID+'.H2.tre',format=9)
		h4.prune(sp_names)
		h4.write(outfile=ID+'.H4.tre',format=9)
		h5.prune(sp_names)
		h5.write(outfile=ID+'.H5.tre',format=9)
		h7.prune(sp_names)
		h7.write(outfile=ID+'.H7.tre',format=9)
		h9.prune(sp_names)
		h9.write(outfile=ID+'.H9.tre',format=9)
		h12.prune(sp_names)
		h12.write(outfile=ID+'.H12.tre',format=9)
		h15.prune(sp_names)
		h15.write(outfile=ID+'.H15.tre',format=9)
		if len(g4)>0:
			h3.prune(sp_names)
			h3.write(outfile=ID+'.H3.tre',format=9)
			h6.prune(sp_names)
			h6.write(outfile=ID+'.H6.tre',format=9)
			h8.prune(sp_names)
			h8.write(outfile=ID+'.H8.tre',format=9)
			h10.prune(sp_names)
			h10.write(outfile=ID+'.H10.tre',format=9)
			h11.prune(sp_names)
			h11.write(outfile=ID+'.H11.tre',format=9)
			h13.prune(sp_names)
			h13.write(outfile=ID+'.H13.tre',format=9)
			h14.prune(sp_names)
			h14.write(outfile=ID+'.H14.tre',format=9)
	if len(g1)>0 and len(g6)>0:
		h16.prune(sp_names)
		h16.write(outfile=ID+'.H16.tre',format=9)
	if len(g1)>0 and len(g7)>0:
		h17.prune(sp_names)
		h17.write(outfile=ID+'.H17.tre',format=9)


generate_constrained_tree(sys.argv[1])
