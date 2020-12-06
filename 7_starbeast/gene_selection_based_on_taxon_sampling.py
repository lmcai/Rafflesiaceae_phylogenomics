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

#for the 377 orthogroups that have only one missing clade, get their sapria and ricinus ID
x=open('G377_one_missing.list').readlines()
out=open('G377_one_missing_sap_ricinus.list','a')
for l in x:
	t=Tree('3_na_tree_1to1_yang_subtree/'+l.strip(),format=1)
	sap=[leaf.name for leaf in t if leaf.name.startswith('Sap')]
	ricinus=[leaf.name for leaf in t if leaf.name.startswith('Rici')]
	out.write(l.strip()+'\t'+'\t'.join(sap+ricinus)+'\n')

out.close()
#add intron and gene length info
x=open('Sapria_longintron.rnd2_old_snap_evidence_based.GAAS.filtered.csv').readlines()
a={}
for l in x:
	a[l.split(',')[0]]='\t'.join(l.split(',')[1:])

z=open('G377_one_missing_sap_ricinus.list').readlines()
out=open('G377_one_missing_sap_ricinus_geneInfo.list','a')
for l in z:
	try:
		out.write(l.strip()+'\t'+a[l.split()[1]])
	except:
		out.write(l)

#get position of the longest CDS in sapria gene
from gff3 import Gff3
gff = Gff3('Sapria_longintron.rnd2_old_snap_evidence_based.GAAS.gff')
CDS_length={}
CDS_pos={}
for l in gff.lines:
	if l['type']=='mRNA':
	#if l['type']=='mRNA' and float(l['attributes']['_AED']) <0.5:
		rna_nam=l['attributes']['Name']
		CDS_end_pos=0
		CDS_max=0
		for rec in l['children']:
			if rec['type']=='CDS':
				CDS_len=rec['end']-rec['start']+1
				CDS_end_pos=CDS_end_pos+CDS_len
				if CDS_len>CDS_max:
					CDS_max=CDS_len
					CDS_max_pos=[CDS_end_pos-CDS_len,CDS_end_pos]
		CDS_length[rna_nam]=CDS_max
		CDS_pos[rna_nam]=CDS_max_pos
		
x=open('/n/holyscratch01/davis_lab/lmcai/sapria_phylogeny/12_starbeast2/G377_one_missing_sap_ricinus_geneInfo.list').readlines()
out=open('G377_one_missing_sap_ricinus_geneInfo.list','a')
out.write(x[0])
for l in x[1:]:
	try:
		out.write(l.strip()+'\t'+str(CDS_length[l.split()[1]])+'\t'+str(CDS_pos[l.split()[1]])+'\n')
	except:
		out.write(l)


###############
#generate alignment
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
#x=open('G81_no_missing_data.list').readlines()
x=open('G377_one_missing_sap_ricinus_geneInfo.list').readlines()

for l in x:
	try:
		exon_pos=l.split()[-1].split(',')
		recs=SeqIO.parse('../4_na_aln_1to1/'+l.split('.')[0]+'.aln','fasta')
		recs=MultipleSeqAlignment(recs)
		#get exon position
		sap_rec=[i for i in recs if i.id=='Sap']
		sap_rec=sap_rec[0]
		#get start and end position
		cur_sap_na=0
		aln_exon=[]
		for i in range(0,len(sap_rec.seq)):
			if not sap_rec.seq[i]=='-':
				cur_sap_na=cur_sap_na+1
				if cur_sap_na>int(exon_pos[0]) and cur_sap_na<int(exon_pos[1]):aln_exon.append(i)
		recs=recs[:,aln_exon[0]:aln_exon[-1]]
	#recs=SeqIO.parse('../4_na_aln_1to1/'+l.strip()+'.aln','fasta')
		out=open(l.split('.')[0]+'.aln','a')
		for rec in recs:
			if rec.id.startswith(('Sap','Rhi','Rca','Rtu','Ery', 'Galearia', 'Drypetes', 'Clutia', 'Ricinus', 'Ixonanthes', 'Crossopetalum', 'Oxalis')):
				d=SeqIO.write(rec,out,'fasta')
	except:print l.split('.')[0]
		