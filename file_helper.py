from constants import mRNA_GENBANK_FILES
from subprocess import call

def file_param_to_file_name(param):
    file_names = {
        "1A": mRNA_GENBANK_FILES.ISOFORM_A_PREPROTEIN.value,
        "1B": mRNA_GENBANK_FILES.ISOFORM_B.value,
        "2A": mRNA_GENBANK_FILES.FASTA_FILE_A.value,
        "2B": mRNA_GENBANK_FILES.FASTA_FILE_B.value
    }

    folder_names = {
        "1A": "inputs/",
        "1B": "inputs/",
        "2A": "outputs/",
        "2B": "outputs/",
    }
 
    return folder_names.get(param, "") + file_names.get(param, "")

def generate_output_path(filename):
    return "outputs/" + filename


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
