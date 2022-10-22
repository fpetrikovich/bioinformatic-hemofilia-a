import os.path
from constants import FASTA_FILES, mRNA_GENBANK_FILES
from subprocess import call
from constants import FASTA_TYPES, OUTPUT_TYPES, EX1_OUTPUT_DIR

def file_param_to_file_name(param):
    file_names = {
        "1A": mRNA_GENBANK_FILES.ISOFORM_A_PREPROTEIN.value,
        "1B": mRNA_GENBANK_FILES.ISOFORM_B.value,
        "2A": FASTA_FILES.ORFS_ISOFORM_A_PREPROTEIN.value,
        "2B": FASTA_FILES.ORFS_ISOFORM_B.value,
        "2AProtein": FASTA_FILES.ISOFORM_A_PREPROTEIN.value,
        "2BProtein": FASTA_FILES.ISOFORM_B.value
    }

    folder_names = {
        "1A": "inputs/",
        "1B": "inputs/",
        "2A": "/",
        "2B": EX1_OUTPUT_DIR,
        "2AProtein": EX1_OUTPUT_DIR,
        "2BProtein": EX1_OUTPUT_DIR,
    }
 
    return folder_names.get(param, "") + file_names.get(param, "")

def generate_output_path(filename):
    return EX1_OUTPUT_DIR + filename

def generate_report_path(filename, is_local = True):
    return "reports/" + ("local_" if is_local else "") + filename

def save_file(file_name, content): 
	f = open(file_name, "w")
	f.write(content)
	f.close()

def create_bash_file(file_name, command):
    f = open(file_name, "w")
    f.write("#!/bin/sh\n")
    f.write(command)
    f.close()

def run_bash_file(file_name):
    rc = call(file_name, shell=True)

def file_exists(file):
    return os.path.exists(file)

def valid_fasta_file(file):
  return file_exists(file) and file.split('.')[-1] in FASTA_TYPES

def valid_output_file(file):
  return file.split('.')[-1] in OUTPUT_TYPES
