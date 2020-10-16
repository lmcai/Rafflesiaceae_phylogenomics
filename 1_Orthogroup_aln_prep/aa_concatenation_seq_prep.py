from Bio import SeqIO
from ete3 import Tree

orthos=open('../malpighiales_sapria_orthogroup.list').readlines()

#prepare 1 to 1 orthogroup sequences
for i in range(0,len(orthos)):
	try:
		malp_ortho=orthos[i].split()[0]
		malp_tr=Tree('../na_tree_yang_pruned/'+`i`+'.inclade1.ortho1.tre',format=1)
		malp_sp=[leaf.name for leaf in malp_tr]
		malp_seqs=SeqIO.parse(`i`+'.aa.aln.trimmed.fas','fasta')
		outfile=open(`i`+'.aa.aln.trimmed.1to1.fas','a')
		for seq in malp_seqs:
			if seq.id in malp_sp:
				d=SeqIO.write(seq,outfile,'fasta')
		outfile.close()
	except:
		print i

#prepare concatenation sequences
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
import fnmatch, os


filename=[]
for fn in os.listdir('.'):
    if fnmatch.fnmatch(fn,'*fas'):filename.append(fn)

for fn in filename:
    x=AlignIO.read(fn,'fasta',alphabet=Gapped(IUPAC.protein))
    new_filename='.'.join(fn.split('.')[:-1])+'.nex'
    g = open(new_filename, "w")
    d=g.write(x.format("nexus"))
    g.close()


#################
#concatenate nexus to super matrix in a order
#
from Bio.Nexus import Nexus

file_list = open('loci.order').readlines()
file_list =[l.strip() for l in file_list]
nexi =  [(fname, Nexus.Nexus(fname)) for fname in file_list]


combined = Nexus.combine(nexi)
output=open('malp_500G.nex', 'w')
combined.write_nexus_data(output)
output.close()