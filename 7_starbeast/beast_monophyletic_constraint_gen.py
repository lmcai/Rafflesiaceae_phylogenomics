#to set outgroups in beast, the mrca for all gene trees, species tree, and lod need to be set
#a initial tree in newick format also needs to be provided


import os
out=open('beast_monophyletic_constraint.txt','a')

#add mrca constraint for two trees in beauti to see where the following text should be inserted
genes=os.listdir('.')
genes=[g for g in genes if g.endswith('fasta')]
j=1
for g in genes:
	x=open(g).readlines()
	sp=[i.strip()[1:] for i in x if i.startswith('>')]
	ingroup=[i for i in sp if not i.startswith(('Oxa','Cro'))]
	out.write("            <distribution id=\"ingroup"+`j`+".prior\" spec=\"beast.math.distributions.MRCAPrior\" monophyletic=\"true\" tree=\"@Tree.t:"+g.split('.')[0]+"\">\n")
	out.write("                <taxonset id=\"ingroup"+`j`+"\" spec=\"TaxonSet\">\n")
	for i in ingroup:
		out.write("                    <taxon idref=\""+i+"\"/>\n")
	out.write("                </taxonset>\n")
	out.write("            </distribution>\n")
	j=j+1

j=1

out.write('\n\n#add the following after \"<log idref=\"ExtinctionFraction.t:Species\"/>\"..\n')
for g in genes:
	out.write("        <log idref=\"ingroup"+`j`+".prior\"/>\n")
	j=j+1

out.close()

#prepare the starting tree
#
import os
import ete3
x=os.listdir('.')
x=[i.split('.')[0] for i in x if i.startswith(('1','2','3','4','5','8'))]
out=open('beast_monophyletic_constraint.txt','a')
from ete3 import Tree
for i in x:
	t=Tree('../3_na_tree_1to1_yang_subtree/'+i+'.inclade1.ortho1.tre',format=1)
	sp=open(i+'.trees').readlines()
	sp=[j.strip() for j in sp[5:20] if j.startswith('\t\t\t')]
	sp2out=[leaf.name for leaf in t if leaf.name.startswith(tuple(sp))]
	t.prune(sp2out)
	ingroup=[leaf.name for leaf in t if leaf.name.startswith(('Sap','Rhi','Rca','Rtu','Ery', 'Galearia', 'Drypetes', 'Clutia', 'Ricinus', 'Ixonanthes'))]
	#monophyletic ingroup tree
	if t.check_monophyly(values=ingroup, target_attr="name")[0]:
		out.write("    <init spec='beast.util.TreeParser' id='NewickTree.t:"+i+"' initial=\"@Tree.t:"+i+"\" taxa='@"+i+"' IsLabelledNewick=\"true\" newick=\""+t.write(format=5)+"\"/>\n")
	else:
		t=Tree('misc/'+i+'.tre',format=1)
		out.write("    <init spec='beast.util.TreeParser' id='NewickTree.t:"+i+"' initial=\"@Tree.t:"+i+"\" taxa='@"+i+"' IsLabelledNewick=\"true\" newick=\""+t.write(format=5)+"\"/>\n")

#non-monophyletic gene trees:

#replace the following bock to newick trees
<init estimate="false" id="RandomTree.t:XYZ26" initial="@Tree.t:XYZ26" spec="beast.evolution.tree.RandomTree" taxa="@XYZ26">
    <populationModel id="ConstantPopulation0.t:XYZ26" spec="ConstantPopulation">
        <parameter id="randomPopSize.t:XYZ26" name="popSize">1.0</parameter>
    </populationModel>
</init>



