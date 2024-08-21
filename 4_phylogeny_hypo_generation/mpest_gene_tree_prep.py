from ete3 import Tree
import ete3
x=open('G2141.aster.genetree.trees').readlines()

out=open('G2141.mpest.gene.trees','a')
trees=[]
for l in x:
	t=Tree(l.strip())
	sp=[node.name for node in t]
	try:
		raff=set(sp) & {'Sap','Rhi','Rca','Rtu'}
		if len(raff)>1:
			temroot=t.get_common_ancestor(raff)
			t.set_outgroup(temroot)
		else:
			t.set_outgroup(list(raff)[0])
		if 'Oxalis' in sp and 'Elaeocarpus' in sp:
			root=t.get_common_ancestor('Oxalis','Elaeocarpus')
			t.set_outgroup(root)
		elif 'Oxalis' in sp:
			t.set_outgroup(t&"Oxalis")
		elif 'Elaeocarpus' in sp:
			t.set_outgroup(t&"Elaeocarpus")
		elif 'Crossopetalum' in sp:
			t.set_outgroup(t&"Crossopetalum")
		trees.append(t.write(format=9))
	except:
		print(t.write())


out.write('\n'.join(trees))
out.close()