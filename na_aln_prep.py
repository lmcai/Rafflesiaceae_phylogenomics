from ete3 import Tree
from Bio import SeqIO

orthos=open('../malpighiales_sapria_orthogroup.list').readlines()

for i in range(0,len(orthos)):
	malp_ortho=orthos[i].split()[0]
	malp_tr=Tree('/n/davis_lab/Users/lmcai/Malpighiales/WGD/2_round/homo2_22sp/'+malp_ortho)
	malp_sp=[leaf.name for leaf in malp_tr]
	try:
		malp_seqs=SeqIO.index('/n/davis_lab/Users/lmcai/Malpighiales/WGD/2_round/NA_aln/'+malp_ortho.split('.')[0]+'.na.fas','fasta')
	except ValueError:
		malp_seqs=SeqIO.index('./tem/'+malp_ortho.split('.')[0]+'.na.fas','fasta')
		
	outfile=open(`i`+'.na.fas','a')
	for sp in malp_sp:
		d=SeqIO.write(malp_seqs[sp],outfile,'fasta')
		#sap_ortho=orthos[i].split()[1]
		#sap_tr=Tree('../yang_subtree/'+sap_ortho)
		#sap_sp=[leaf.name for leaf in sap_tr if leaf.name.startswith(('Rca','Rtu','Rhi','Sap','pSHI'))]
		#sap_seqs=SeqIO.index('/scratch/lmcai/15_orthofinder_seq/na_aln/'+sap_ortho.split('.')[0]+'.na.fas','fasta')
		#for sp in sap_sp:
		#	d=SeqIO.write(malp_seqs[sp],outfile,'fasta')
	outfile.close()
	
	