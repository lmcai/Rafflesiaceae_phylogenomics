from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
import fnmatch, os, argparse


parser = argparse.ArgumentParser(description='Concatenate fasta alignments and generate gene partition file (nexus).')
parser.add_argument('--input_dir', help='path to the alignments', required=True)
parser.add_argument('--suffix',  help='suffix of alignment files', required=True)
parser.add_argument('--output',  help='output directory', required=True)
parser.add_argument('--loci_order',  help='(optional) a file containing the order of the loci in the concatenation')

args = parser.parse_args()


filename=[]
new_filenames=[]
for fn in os.listdir(args.input_dir):
    if fnmatch.fnmatch(fn,'*'+args.suffix):filename.append(fn)

for fn in filename:
    x=AlignIO.read(args.input_dir+'/'+fn,'fasta',alphabet=Gapped(IUPAC.protein))
    new_filename=args.output+'/'+'.'.join(fn.split('.')[:-1])+'.nex'
    new_filenames.append(new_filename)
    g = open(new_filename, "w")
    d=g.write(x.format("nexus"))
    g.close()


#################
#concatenate nexus to super matrix in a order
#
from Bio.Nexus import Nexus

try:
	file_list = open(args.loci_order).readlines()
	file_list =[l.strip() for l in file_list]
except:
	file_list =new_filenames

nexi =  [(fname, Nexus.Nexus(fname)) for fname in file_list]

combined = Nexus.combine(nexi)
output=open('concatenated_aln.nex', 'w')
combined.write_nexus_data(output)
output.close()