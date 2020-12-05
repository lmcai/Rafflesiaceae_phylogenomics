#filter gene for starbeast
#taxon sampling Sap, Rhi, Rca, Rtu, Ery, Rhizo, Galearia, Drypetes, Clutia, Ricinus, Ixonanthes, Crossopetalum, Oxalis
from ete3 import Tree

x=open('../6_sortdata/gene_tree_combined_statistics.tsv').readlines()

out=open('gene_tree_SortDataStatistics_TaxonSampling.tsv','a')

for l in x:
	t=Tree('../3_na_tree_1to1_yang_subtree/'+l.split()[0],format=1)
	sap=[leaf.name for leaf in t if leaf.name.startswith(('Sap','Rhi','Rca','Rtu')) and not leaf.name.startswith('Rhizo')]
	focal_sp=[leaf.name for leaf in t if leaf.name.startswith(('Sap','Rhi','Rca','Rtu','Ery', 'Galearia', 'Drypetes', 'Clutia', 'Ricinus', 'Ixonanthes', 'Crossopetalum', 'Oxalis'))]
	Ixon=[leaf.name for leaf in t if leaf.name.startswith('Ixo')]
	Pan=[leaf.name for leaf in t if leaf.name.startswith('Galea')]
	dry=[leaf.name for leaf in t if leaf.name.startswith('Drype')]
	eupho=[leaf.name for leaf in t if leaf.name.startswith(('Clutia','Ricinus'))]
	ER=[leaf.name for leaf in t if leaf.name.startswith(('Ery','Rhizo'))]
	outgr=[leaf.name for leaf in t if leaf.name.startswith(('Crossopetalum','Oxalis'))]
	out.write(l.strip()+'\t'+'\t'.join([str(i) for i in [len(sap),len(Ixon),len(Pan),len(dry),len(eupho),len(ER),len(outgr),len(focal_sp)]])+'\n')

out.close()

from Bio import SeqIO
x=open('G81_no_missing_data.list').readlines()
for l in x:
	recs=SeqIO.parse('../4_na_aln_1to1/'+l.strip()+'.aln','fasta')
	out=open(l.strip()+'.aln','a')
	for rec in recs:
		if rec.id.startswith(('Sap','Rhi','Rca','Rtu','Ery', 'Galearia', 'Drypetes', 'Clutia', 'Ricinus', 'Ixonanthes', 'Crossopetalum', 'Oxalis')):
			d=SeqIO.write(rec,out,'fasta')
		