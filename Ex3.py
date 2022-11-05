import sys
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqIO
from Bio.SubsMat import MatrixInfo 

from constants import FASTA_TYPES
from file_helper import valid_fasta_file, valid_output_file

def validate_params(input_file, species_list, output_file):
  if not valid_fasta_file(input_file):
    print(f"Input file path must exist and must be FASTA (.{FASTA_TYPES[0]} o .{FASTA_TYPES[1]})")
    exit(1)
  
  for specie_file in species_list:
    if not valid_fasta_file(specie_file):
        print(specie_file)
        print(f"Species files path must exist and must be FASTA (.{FASTA_TYPES[0]} o .{FASTA_TYPES[1]})")
        exit(1)
  
  if not valid_output_file(output_file):
    print(f"Output file must be .txt | .output | .out")
    exit(1)

def perform_sequence_alignment(origin_file, species_list, output_file):
    # Use the Blosum62 matrix
    matrix = MatrixInfo.blosum62

    original_sequence = list(SeqIO.parse(open(origin_file,'r'), 'fasta'))[0].seq
    sequences = []

    for file in species_list:
        sequences.append(list(SeqIO.parse(open(file,'r'), 'fasta'))[0].seq)
    
    results = []
    for seq in sequences:
        results.append(pairwise2.align.globaldx(original_sequence, seq, matrix))
    
    output = open(f'{output_file}', 'w')
    for i in range(len(results)):
        result = results[i]
        file_name = species_list[i]
        output.write("MSA made with > {f}\n".format(f = file_name))
        output.write(format_alignment(*result[-1]))
        output.write('\n')

def run_exercise_3(origin_file, species_list, output_file):
    validate_params(origin_file, species_list, output_file)
    perform_sequence_alignment(origin_file, species_list, output_file)
    print("MSA successfully completed. ")