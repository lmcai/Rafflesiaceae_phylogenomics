from Bio import SeqIO
import sys, os

def parse_log_file(log_file):
    mask_positions = {}
    with open(log_file, 'r') as file:
        current_id = None
        for line in file:
            if '\t' not in line:
                current_id = line.split()[0]
                if current_id not in mask_positions:
                    mask_positions[current_id] = []
            else:
                if current_id:
                    positions = line.strip().split('-')
                    start = int(positions[0]) - 1  # convert to 0-based index
                    end = int(positions[1])        # end is exclusive
                    mask_positions[current_id].append((start, end))
    return mask_positions

def mask_fasta(fasta_file, mask_positions, output_file):
    sequences = SeqIO.index(fasta_file, "fasta")
    output_handle=open(output_file, "w")
    for seq_id, positions in mask_positions.items():
        seq = list(str(sequences[seq_id].seq))
        for start, end in positions:
            seq[start:end] = 'N' * (end - start)
        #sequences[seq_id].seq = "".join(seq)
        #print(seq_id,sequences[seq_id].seq)
        output_handle.write('>'+seq_id+'\n'+"".join(seq)+'\n')

# execute
log_file = sys.argv[2]
fasta_file = sys.argv[1]
output_file = sys.argv[1].split('.')[0]

mask_positions = parse_log_file(log_file)
mask_fasta(fasta_file, mask_positions, output_file+'.tem.fas')

##############################
#remove completely ambigous codons
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

# Function to check if a site is ambiguous
def is_ambiguous(site):
    ambiguous_chars = {"N", "-","?"}
    return all(base in ambiguous_chars for base in site)

# Function to remove codons with ambiguous sites
def remove_ambiguous_codons(alignment):
    num_seqs = len(alignment)
    codon_size = 3
    #new_alignment = []
    new_sequences = [""] * num_seqs 
    for i in range(0, alignment.get_alignment_length(), codon_size):
        codon_sites = alignment[:, i:i+codon_size]
        ambiguous_count = sum(is_ambiguous(codon) for codon in codon_sites)
        if (any(is_ambiguous(codon_sites[:, j]) for j in range(codon_size))) or (ambiguous_count>0.8*num_seqs):
            continue  # Skip this codon if any site is ambiguous
        #new_alignment.append(codon_sites)
        for j in range(num_seqs):
            new_sequences[j] = new_sequences[j] + codon_sites[j]
    # Combine the selected codons back into sequences
    #new_alignment_strs = ["".join([new_alignment[j][i] for j in range(len(new_alignment))]) for i in range(num_seqs)]
    return new_sequences

# Parse the alignment and remove ambiguous codons
alignment = AlignIO.read(output_file+'.tem.fas', "fasta")
cleaned_alignment = remove_ambiguous_codons(alignment)

# Convert the cleaned alignment back to a string and save it to a file
SeqIO.write(cleaned_alignment, output_file+'.na.mask.fas', "fasta")

os.system("rm "+output_file+'.tem.fas')  
