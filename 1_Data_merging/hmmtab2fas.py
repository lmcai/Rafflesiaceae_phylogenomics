import os
import re
from Bio import SeqIO

# Define directory paths and file names
tabfile_directory = './'  # Update with your directory path
fasta_file = 'Apo_Pil.cdhit1.trans.fas'
output_directory = './output_sequences'

# Ensure the output directory exists
os.makedirs(output_directory, exist_ok=True)

# Function to parse .tabfile to find best hits
def parse_tabfile(filepath):
    best_hits = {"Apodanthes": None, "Pilostyles": None}
    with open(filepath, 'r') as file:
        for line in file:
            # Skip comment lines
            if line.startswith("#"):
                continue           
            # Parse columns (target name, e-value)
            columns = line.split()
            target_name = columns[0]
            evalue = float(columns[5])           
            # Check for Apodanthes and Pilostyles records and update best hits
            for genus in best_hits:
                if target_name.startswith(genus):
                    if best_hits[genus] is None or evalue < best_hits[genus][1]:
                        best_hits[genus] = (target_name, evalue)   
    return best_hits

# Function to extract sequences from FASTA file
def extract_sequences(best_hits, output_file):
    with open(output_file, 'w') as fasta_out:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            if record.id in (best_hits["Apodanthes"][0], best_hits["Pilostyles"][0]):
                SeqIO.write(record, fasta_out, 'fasta')

# Main processing loop
for tabfile in os.listdir(tabfile_directory):
    if tabfile.endswith('.tabfile'):
        filepath = os.path.join(tabfile_directory, tabfile)       
        # Parse .tabfile to get best hits
        best_hits = parse_tabfile(filepath)       
        # Skip files without valid hits for both genera
        if not best_hits["Apodanthes"] or not best_hits["Pilostyles"]:
            continue     
        # Define output filename using the .tabfile prefix
        file_prefix = os.path.splitext(tabfile)[0]
        output_file = os.path.join(output_directory, f"{file_prefix}.fasta")      
        # Extract and save sequences to output FASTA file
        extract_sequences(best_hits, output_file)

