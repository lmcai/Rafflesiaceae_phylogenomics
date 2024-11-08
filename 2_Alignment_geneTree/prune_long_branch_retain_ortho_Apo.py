from ete3 import Tree
import os
from scipy import stats

#remove long internal branches (potential paralogs)
#reroot the tree
#identify orthologs
treefiles=os.listdir('.')
treefiles=[i for i in treefiles if i.endswith('.treefile')]

def prune_long_br(tr,length_limit):
	#calculate trimmed mean branch length
	br_len=[node.dist for node in tr.traverse("postorder")]
	mean_br = stats.trim_mean(br_len, 0.05)
	#first round set to >15 times, second time set to >10 times
	tips2remove=[]
	for node in tr.traverse("postorder"):
		if node.dist>length_limit*mean_br:tips2remove=tips2remove+[j.name for j in node]
	tips2remove=list(set(tips2remove))
	#do not remove anything in our data yet
	#tips2remove=[j for j in tips2remove if not j in jdata]
	if len(tips2remove)>0:
		taxa=[node.name for node in tr if not node.name in tips2remove]
		tr.prune(taxa,preserve_branch_length =True)
	return(tr)

raff_headers = (">Sap", ">Rhi", ">Rca",">Rtu") 
malp_headers=('Ptr','Mes','Jat','Sap','Rhi','Rtu','Rca','Apoda','Pilo')
for i in treefiles:
	tree=Tree(i)
	#num_taxa_before=len(tree)
	#reroot the tree
	sp=[node.name for node in tree]
	root_candidate=[i for i in sp if i.startswith('Amb')]
	root_candidate=root_candidate+[i for i in sp if i.startswith('Cin')]
	root_candidate=root_candidate+[i for i in sp if i.startswith('Osa')]
	root_candidate=root_candidate+[i for i in sp if i.startswith('Sor')]
	root_candidate=root_candidate+[i for i in sp if i.startswith('Aqu')]
	root_candidate=root_candidate+[i for i in sp if i.startswith('Nel')]
	try:
		tree.set_outgroup(root_candidate[0])
		#two rounds of long branch pruning; first remove branches longer than 15 times then average; second round 10 times
		#tree = prune_long_br(tree,10)
		tree = prune_long_br(tree,15)
		tree.write(outfile=i.split('.')[0]+'.nolongbr_rooted.tre')
		#label Malpighiales
		for node in tree:
			if node.name.startswith(malp_headers):node.add_features(label="Malp")
		#get the ortholog Apodanthaceae for each Rafflesiaceae clade
		raff=open('../add_apo_rna_geneTr_rnd1/'+i.split('.')[0]+'.combined.fas').readlines()
		raff=[i[1:].strip() for i in raff if i.startswith(raff_headers)]
		for node in tree.get_monophyletic(values=["Malp"], target_attr="label"):
			malp_sp=[leaf.name for leaf in node]
			if len(set(malp_sp) & set(raff))>0:
				apo_output=[j for j in malp_sp if (j.startswith('Apodan') or j.startswith('Pilo'))]
		if len(apo_output):print(i+': '+', '.join(apo_output))
		else:print(i+': No orthologous Apodanthaceae')
	except:
		print(i+': Failed rooting')


#################
#extract seq according to 
from Bio import SeqIO

treefiles=os.listdir('.')
#treefiles=[i for i in treefiles if i.endswith('.nolongbr.tre')]
treefiles=[i for i in treefiles if i.endswith('.nolongbrnodup.tre')]

for i in treefiles:
	tr=Tree(i)
	taxa=[j.name for j in tr]
	#recs=SeqIO.parse(i.split('.')[0]+'.aln.fas','fasta')
	recs=SeqIO.parse(i.split('.')[0]+'.nolongbr.Rtrim.aln.fas','fasta')
	#out=open(i.split('.')[0]+'.nolongbr.aln.fas','a')
	out=open(i.split('.')[0]+'.nolongbr.Rtrim.nodup.fas','a')
	for rec in recs:
		if rec.id in taxa:
			d=SeqIO.write(rec,out,'fasta')
	out.close()

