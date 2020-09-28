from ete3 import Tree
x=open('G2299.gene.trees').readlines()

out=open('G2299.mpest.gene.trees','a')
trees=[]
for l in x:
	t=Tree(l.strip())
	trees.append(t.write(format=9))

out.write('\n'.join(trees))
out.close()