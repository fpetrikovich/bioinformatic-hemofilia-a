import os.path
from subprocess import call
from constants import FASTA_TYPES, OUTPUT_TYPES, EX1_OUTPUT_DIR, FILE_IDENTIFIERS, FOLDER_IDENTIFIERS
from error_helper import exit_with_error

def file_param_to_file_name(param):
    return FOLDER_IDENTIFIERS.get(param, "") + FILE_IDENTIFIERS.get(param, "")

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
    try:
        call(file_name, shell=True)
    except Exception as e:
        exit_with_error("Error while running bash script.", e)
        
def run_bash_file_with_arguments(file_name, args):
    try: 
        call(['bash', file_name, *args])
    except Exception as e:
        exit_with_error("Error while running bash script.", e)

def file_exists(file):
    return os.path.exists(file)

def valid_fasta_file(file):
  return file_exists(file) and file.split('.')[-1] in FASTA_TYPES

def valid_output_file(file):
  return file.split('.')[-1] in OUTPUT_TYPES

def delete_file(path):
    try:
        os.remove(path)
    except Exception as e:
        exit_with_error("Error while deleting file.", e)
