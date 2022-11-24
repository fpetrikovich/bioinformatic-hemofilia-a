from Bio.Align.Applications import ClustalwCommandline

from error_helper import exit_with_error
from constants import FASTA_TYPES, MSA_OUTPUT_DIR
from file_helper import valid_fasta_file, valid_output_file, create_bash_file, run_bash_file


def validate_params(input_file, output_file):
  if not valid_fasta_file(input_file):
    exit_with_error("Input file path must exist and must have a FASTA extension.")
  
  if not valid_output_file(output_file):
    exit_with_error("Output file must be .txt | .output | .out")

def perform_sequence_alignment(input_file, output_file):
    bash_file = MSA_OUTPUT_DIR + "command.sh"
    clustalw_cline = ClustalwCommandline("clustalw", infile=input_file, outfile=output_file, pwgapext=0.1, pwgapopen=10, gapext=0.2, gapopen=10)
    print(clustalw_cline)
    create_bash_file(bash_file, str(clustalw_cline))
    run_bash_file(bash_file)
    
def run_exercise_3(origin_file, output_file):
    validate_params(origin_file, output_file)
    perform_sequence_alignment(origin_file, output_file)
    print("MSA successfully completed. ")