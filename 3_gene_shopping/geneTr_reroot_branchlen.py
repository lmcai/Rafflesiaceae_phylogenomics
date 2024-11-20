#/usr/bin/python
from ete3 import Tree
import os
import ete3
x=os.listdir('.')
x=[i for i in x if i.endswith('treefile')]
for j in x:
	t=Tree(j,format=1)
	outgroup=[]
	out_sp=[]
	out_sp=[i.name for i in t if i.name.startswith(('Elaeocarpus','Oxalis'))]
	if len(out_sp)>1:
		outgroup = t.get_common_ancestor(out_sp)
	elif len(out_sp)==1:
		outgroup = t&out_sp[0]
	else:
		out_sp=[i.name for i in t if i.name.startswith('Crossopetalum')]
		if len(out_sp)>0:
			outgroup = t&out_sp[0]
	if outgroup:
		try:
			t.set_outgroup(outgroup)
		except ete3.coretype.tree.TreeError:
			t.set_outgroup(out_sp[0])
		t.write(format=1, outfile=j.split('.')[0]+'.rooted.tre')
	else:
		#choose a non-raff to be outgroup
		out_sp=[i.name for i in t if i.name.startswith(('Rinorea','Viola','Malesherbia','Passiflora','Casearia','Flacourtia','Populus','SalixSu','Salix'))]
		out_sp=out_sp+[i.name for i in t if i.name.startswith(('Tristellateia','Galphimia','Elatine','Bergia','Sauropus','Bischofia'))]
		out_sp=out_sp+[i.name for i in t if i.name.startswith(('Calophyllum','Mammea','Hypericum','Podostemum','Clusia','Garcinia','Ochna','Linum','Chrysobalanus'))]
		if len(out_sp)>0:
			t.set_outgroup(out_sp[0])
			t.write(format=1, outfile=j.split('.')[0]+'.rooted.tre')
		else:
			print j

x=open('../G2141.branch_len.tsv').readlines()
x=[i.split()[0] for i in x[1:]]
out=open('G2141.branch_len.tsv','a')
out.write('ID\tstem_raff_branch_len\tancestor_of_stem_raff_bran_len\tBP\tdist2root\n')
for j in x:
	try:
		t=Tree(j+'.rooted.tre',format=1)
		raff=[i.name for i in t if i.name.startswith(('Sap','Rca','Rtu','Rhi')) and not i.name.startswith('Rhizo')]
		#if raff are monophyletic
		if t.check_monophyly(values=raff, target_attr="name")[0]:
			if len(raff)>1:
				raff_node=t.get_common_ancestor(raff)
				out.write(j+'\t'+str(raff_node.dist)+'\t'+str(raff_node.get_ancestors()[0].dist)+'\t'+raff_node.get_ancestors()[0].name+'\t'+str(len(raff_node.get_ancestors()))+'\n')
			elif len(raff)==1:
				raff_node=t&raff[0]
				out.write(j+'\tNA\tNA\t'+raff_node.get_ancestors()[0].name+'\t'+str(len(raff_node.get_ancestors()))+'\n')
		#choose the node with most species
		else:
			for leaf in t:
				if leaf.name in raff:leaf.add_features(raff='Y')
			curr_sp_num=0
			for node in t.get_monophyletic(values='Y', target_attr="raff"):
				if len(node.get_leaf_names())>curr_sp_num:
					curr_sp_num=len(node.get_leaf_names())
					raff_node=node
			if curr_sp_num>1:
				out.write(j+'\t'+str(raff_node.dist)+'\t'+str(raff_node.get_ancestors()[0].dist)+'\t'+raff_node.get_ancestors()[0].name+'\t'+str(len(raff_node.get_ancestors()))+'\n')
			elif curr_sp_num==1:
				out.write(j+'\tNA\tNA\t'+raff_node.get_ancestors()[0].name+'\t'+str(len(raff_node.get_ancestors()))+'\n')
	except IOError:
		print j
	except IndexError:
		print j
	except ete3.parser.newick.NewickError:
		print j

################
##IMPORTANT!!! ete3 generate incorrect node labels when rerooting, see this paper https://academic.oup.com/mbe/article/34/6/1535/3077051
#To deal with this problem, I reroot the trees in Notung, using the following command. Then I get the BP support values using ete3
java -jar ~/programs/Notung-2.8.1.3-beta.jar -s sp.tre -b gene.list --root --treeoutput newick --speciestag prefix --silent --progressbar --outputdir rooted_geneTr --nolosses

#/usr/bin/python
from ete3 import Tree
import os
x=os.listdir('.')

out=open('G2141.BP.tsv','a')
out.write('ID\tBP\n')
for j in x:
	t=Tree(j,format=1)
	BP='NA'
	raff=[i.name for i in t if i.name.startswith(('Sap','Rca','Rtu','Rhi')) and not i.name.startswith('Rhizo')]
	if len(raff)>1:
		#if rafflesiaceae species are monophyletic
		if t.check_monophyly(values=raff, target_attr="name")[0]:
			raff_node=t.get_common_ancestor(raff)
			BP=raff_node.get_ancestors()[0].name
		else:
			#print('nonmonophyletic raff '+j)
			for leaf in t:
					if leaf.name in raff:leaf.add_features(raff='Y')
			curr_sp_num=0
			for node in t.get_monophyletic(values='Y', target_attr="raff"):
				if len(node.get_leaf_names())>curr_sp_num:
					BP=node.get_ancestors()[0].name
					curr_sp_num=len(node.get_leaf_names())
	elif len(raff)==1:
		raff_node=t&raff[0]
		BP=raff_node.get_ancestors()[0].name
	else:
		print j
	out.write(j.split('.')[0]+'\t'+BP+'\n')



###############

#The following gene trees do not have outgroups from the 100 bootstrap run in IQTREE
3024.aln.treefile
2874.aln.treefile
2871.aln.treefile
764.aln.treefile
2708.aln.treefile
2321.aln.treefile
2732.aln.treefile
284.aln.treefile
2448.aln.treefile
1288.aln.treefile
270.aln.treefile
965.aln.treefile
2712.aln.treefile
1715.aln.treefile
3055.aln.treefile
2547.aln.treefile
1871.aln.treefile
2049.aln.treefile
1519.aln.treefile
139.aln.treefile
937.aln.treefile
2068.aln.treefile
2818.aln.treefile
1737.aln.treefile
2946.aln.treefile
1327.aln.treefile
2159.aln.treefile
2445.aln.treefile
1274.aln.treefile
2765.aln.treefile
390.aln.treefile
810.aln.treefile
2735.aln.treefile
1585.aln.treefile
2937.aln.treefile
952.aln.treefile
1807.aln.treefile
594.aln.treefile
2647.aln.treefile
2305.aln.treefile
3072.aln.treefile
2654.aln.treefile
2836.aln.treefile
2901.aln.treefile
2867.aln.treefile
1785.aln.treefile
2840.aln.treefile
2509.aln.treefile
2412.aln.treefile
2552.aln.treefile
2366.aln.treefile
547.aln.treefile
2968.aln.treefile
911.aln.treefile
2800.aln.treefile
1562.aln.treefile
1141.aln.treefile
2895.aln.treefile
2845.aln.treefile
366.aln.treefile
1891.aln.treefile
1173.aln.treefile
2688.aln.treefile
3000.aln.treefile
2062.aln.treefile
2576.aln.treefile
1860.aln.treefile
2969.aln.treefile
3006.aln.treefile
2717.aln.treefile
1864.aln.treefile
2887.aln.treefile
1340.aln.treefile
1550.aln.treefile
1467.aln.treefile
2521.aln.treefile
2808.aln.treefile
2545.aln.treefile
1769.aln.treefile
2706.aln.treefile
280.aln.treefile
2217.aln.treefile
2294.aln.treefile
2926.aln.treefile
2832.aln.treefile
2425.aln.treefile
3009.aln.treefile
1998.aln.treefile
57.aln.treefile
2376.aln.treefile
1187.aln.treefile
3022.aln.treefile
802.aln.treefile
2512.aln.treefile
806.aln.treefile