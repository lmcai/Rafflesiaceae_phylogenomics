from ete3 import Tree
from Bio import SeqIO

orthos=open('../malpighiales_sapria_orthogroup.list').readlines()

#malp sequences
for i in range(0,len(orthos)):
	malp_ortho=orthos[i].split()[0]
	malp_tr=Tree('/n/davis_lab/Users/lmcai/Malpighiales/WGD/2_round/homo2_22sp/'+malp_ortho)
	malp_sp=[leaf.name for leaf in malp_tr if not leaf.name.startswith('Podostemum')]
	malp_seqs=SeqIO.parse('./tem/'+malp_ortho.split('.')[0]+'.aa.fas','fasta')
	outfile=open(`i`+'.aa.fas','a')
	for seq in seqs:
		if seq.id.startswith(('Vitis','Elaeocarpus','Crossopetalum','Oxalis')) or seq.id in malp_sp:
			d=SeqIO.write(seq,outfile,'fasta')
	outfile.close()

#sap sequences
for i in range(0,len(orthos)):
	sap_ortho=orthos[i].split()[1]
	sap_tr=Tree('../yang_subtree/'+sap_ortho)
	sap_sp=[leaf.name for leaf in sap_tr if leaf.name.startswith(('Rca','Rtu','Rhi','Sap'))]
	if len(sap_sp)>0:
		outfile=open(`i`+'.na.fas','a')
		sap_seqs=SeqIO.index('/scratch/lmcai/15_orthofinder_seq/aa_aln/'+sap_ortho.split('.')[0]+'.aa.fas','fasta')
		for sp in sap_sp:
			d=SeqIO.write(sap_seqs[sp],outfile,'fasta')
		outfile.close()
	else:
		print('rm '+`i`+'.aa.fas')

