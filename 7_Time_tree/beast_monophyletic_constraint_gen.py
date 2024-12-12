#to set outgroups in beast, the mrca for all gene trees, species tree, and lod need to be set
#a initial tree in newick format also needs to be provided


import os
out=open('beast_monophyletic_constraint.txt','a')

#add mrca constraint for two trees in beauti to see where the following text should be inserted
genes=os.listdir('.')
genes=[g for g in genes if g.endswith('fasta')]
j=1
#j=56
#j=82
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
#j=56
#j=82
out.write('\n\n#add the following after \"<log idref=\"ExtinctionFraction.t:Species\"/>\"..\n')
for g in genes:
	out.write("        <log idref=\"ingroup"+`j`+".prior\"/>\n")
	j=j+1

out.close()

#there are two 1670 and two 1653, add '1' for one of the genes

#prepare the starting tree
#
import os
x=os.listdir('.')
x=[i.split('.')[0] for i in x if i.startswith(('1','2','3','4','5','8'))]
from ete3 import Tree
for i in x:
	#t=Tree('../3_na_tree_1to1_yang_subtree/'+i+'.inclade1.ortho1.tre',format=1)
	#sp=open(i+'.trees').readlines()
	#sp=[j.strip() for j in sp[5:20] if j.startswith('\t\t\t')]
	#sp2out=[leaf.name for leaf in t if leaf.name.startswith(tuple(sp))]
	#i=i.strip()
	t=Tree('../3_na_tree_1to1_yang_subtree/'+i.split('.')[0]+'.inclade1.ortho1.tre',format=1)
	#sp=open('/n/home08/lmcai/EX26_short300bp/'+i).readlines()
	sp=open('G28_Pandaceae_manual_curate/'+i).readlines()
	t
	sp2out=[leaf.name for leaf in t if leaf.name.startswith(tuple(sp))]
	t.prune(sp2out)
	ingroup=[leaf.name for leaf in t if leaf.name.startswith(('Sap','Rhi','Rca','Rtu','Ery', 'Galearia', 'Drypetes', 'Clutia', 'Ricinus', 'Ixonanthes'))]
	#monophyletic ingroup tree
	if t.check_monophyly(values=ingroup, target_attr="name")[0]:
		t.write(outfile=i.split('.')[0]+'.tre',format=1)
	else:
		print i
		#t=Tree('misc/'+i+'.tre',format=1)
		
#These trees then need to be reformated to ultrametric tree with a relatively large root age so that all gene trees coalesce at the root of the species tree
#in R

library(ape)
files=list.files(path = ".", pattern = '*.tre')
for (i in files){
	tr=read.tree(i)
	chr=chronos(tr, lambda = 1, model = "correlated", calibration = makeChronosCalib(tr,node = "root", age.min = 100,age.max=120))
	write.tree(chr,paste(i,'.chr',sep=''))
}

#in python
import os 
timetrees=[i for i in os.listdir('.') if i.endswith('.chr')]
out=open('beast_monophyletic_constraint.txt','a')

for i in timetrees:
	out.write("    <init spec='beast.util.TreeParser' id='NewickTree.t:"+i.split('.')[0]+"' initial=\"@Tree.t:"+i.split('.')[0]+"\" taxa='@"+i.split('.')[0]+"' IsLabelledNewick=\"true\" newick=\""+open(i).readline().strip()+"\"/>\n")

out.close()
#non-monophyletic gene trees:

#replace the following bock to newick trees
<init estimate="false" id="RandomTree.t:XYZ26" initial="@Tree.t:XYZ26" spec="beast.evolution.tree.RandomTree" taxa="@XYZ26">
    <populationModel id="ConstantPopulation0.t:XYZ26" spec="ConstantPopulation">
        <parameter id="randomPopSize.t:XYZ26" name="popSize">1.0</parameter>
    </populationModel>
</init>

#move the SBI section to right above 'mcmc'
#replace the SBI section with a defined species tree and make sure to add scale="0.001" in the species tree.
#This way species tree will have a very recent divergence time and all gene tree will coalesce at the base of the tree.
#Previous errors of can;t find a proper starting place is cause by the violation of coalescent model: some gene trees have more recent coalescent time than the species tree.
#It works!!! YAY!!!

