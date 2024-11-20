###################
#Get nodal support

################
##IMPORTANT!!! ete3 generate incorrect node labels when rerooting, see this paper https://academic.oup.com/mbe/article/34/6/1535/3077051
#To deal with this problem, I reroot the trees in Notung, using the following command. Then I get the BP support values using ete3
java -jar ~/programs/Notung-2.8.1.3-beta.jar -s sp.tre -b gene.list --root --treeoutput newick --speciestag prefix --silent --progressbar --outputdir rooted_geneTr --nolosses

#These rerooted gene trees are good for getting the node labels, but internal branch length are all modified by Notung
#/usr/bin/python
from ete3 import Tree
import os
files=os.listdir('.')
files=[i for i in files if i.endswith('rooting.0.nwk')]

out=open('G2141.node_support.tsv','w')
out.write('ID\tRaff_UFBP\tApo_UFBP\n')

for j in files:
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
		print(j+':No raff')
	apo_BP='NA'
	apo=[i.name for i in t if i.name.startswith(('Apodan','Pilo'))]
	if len(apo)>1:
		if t.check_monophyly(values=apo, target_attr="name")[0]:
			apo_node=t.get_common_ancestor(apo)
			apo_BP=apo_node.get_ancestors()[0].name
		else:
			apo_node=t&apo[0]
			apo_BP=apo_node.get_ancestors()[0].name
	elif len(apo)==1:
		apo_node=t&apo[0]
		apo_BP=apo_node.get_ancestors()[0].name
	d=out.write(j.split('.')[0]+'\t'+BP+'\t'+apo_BP+'\n')

####################################
#To get branch length distributions
#/usr/bin/python
from ete3 import Tree
import os
import ete3
from numpy import median, min

x=open('G2141.node_support.tsv').readlines()
out=open('G2141.node_support_branlen.tsv','w')
out.write('ID\tRaff_UFBP\tApo_UFBP\tmedian_root2tip\tmin_raff_root2tip\traff_brlen\traff_dist2root\tmin_apo_root2tip\tapo_brlen\traff_apo_nomophyly\n')

for j in x[1:]:
	jj=j
	j=j.split()[0]
	#print(j)
	t=Tree(j+'.na.mask.fas.treefile',format=1)
	rooted_notung=Tree('rooted_geneTr/'+j+'.na.mask.fas.treefile.rooting.0.nwk',format=1)
	root_candidate1=[node.name for node in rooted_notung.children[1]]
	root_candidate2=[node.name for node in rooted_notung.children[0]]
	if len(root_candidate1)<len(root_candidate2):
		outgroup=root_candidate1
		backup=root_candidate2
	else:
		outgroup=root_candidate2
		backup=root_candidate1
	if len(outgroup)>1:outgroup_node=t.get_common_ancestor(outgroup)
	else:outgroup_node=t&outgroup[0]
	try:t.set_outgroup(outgroup_node)
	except ete3.coretype.tree.TreeError:
		t.set_outgroup(backup[0])
		if len(outgroup)>1:outgroup_node=t.get_common_ancestor(outgroup)
		else:outgroup_node=t&outgroup[0]
		t.set_outgroup(outgroup_node)
	t.write(format=1, outfile=j+'.rooted.tre')
	median_root2tip=median([outgroup_node.get_distance(node.name) for node in t if not node.name.startswith(('Sap','Rca','Rtu','Rhi_','Apod','Pilo'))])
	min_raff_root2tip='NA'
	raff_nodedist2root='NA'
	raff_brlen='NA'
	min_apo_root2tip='NA'
	apo_brlen='NA'
	raff=[i.name for i in t if i.name.startswith(('Sap','Rca','Rtu','Rhi')) and not i.name.startswith('Rhizo')]
	apo=[i.name for i in t if i.name.startswith(('Apod','Pilo'))]
	#rafflesiaceae related features
	if len(raff)>0:
		min_raff_root2tip=min([outgroup_node.get_distance(node) for node in raff])
		#if raff are monophyletic
		if t.check_monophyly(values=raff, target_attr="name")[0]:
			if len(raff)>1:
				raff_node=t.get_common_ancestor(raff)
			elif len(raff)==1:
				raff_node=t&raff[0]
		#if not monophyletic, choose the node with most species
		else:
			for leaf in t:
				if leaf.name in raff:leaf.add_features(raff='Y')
			curr_sp_num=0
			for node in t.get_monophyletic(values='Y', target_attr="raff"):
				if len(node.get_leaf_names())>curr_sp_num:
					curr_sp_num=len(node.get_leaf_names())
					raff_node=node
		raff_nodedist2root=outgroup_node.get_distance(raff_node, topology_only=True)
		raff_brlen=raff_node.get_ancestors()[0].dist
	#apodanthaceae related
	if len(apo)>0:
		min_apo_root2tip=min([outgroup_node.get_distance(node) for node in apo])
		if t.check_monophyly(values=apo, target_attr="name")[0]:
			if len(apo)>1:
				apo_node=t.get_common_ancestor(apo)
			elif len(apo)==1:
				apo_node=t&apo[0]
		else:apo_node=t&apo[0]
		apo_brlen=apo_node.get_ancestors()[0].dist
	mono='NA'
	if len(raff)>0 and len(apo)>0:
		if t.check_monophyly(values=raff+apo, target_attr="name")[0]:mono='Y'
		else:mono='N'
	print(mono)
	d=out.write(jj.strip()+'\t'+str(median_root2tip)+'\t'+str(min_raff_root2tip)+'\t'+str(raff_brlen)+'\t'+str(raff_nodedist2root)+'\t'+str(min_apo_root2tip)+'\t'+str(apo_brlen)+'\t'+mono+'\n')			


out.close()