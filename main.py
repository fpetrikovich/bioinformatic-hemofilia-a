import argparse
import time

from Ex1 import run_exercise_1
from Ex2 import run_exercise_2
from Ex3 import run_exercise_3
from Ex4 import run_exercise_4
from Ex5 import run_exercise_5
from argument_helper import print_curr_arguments, handle_arguments
from file_helper import check_file_is_not_dir, file_param_to_file_name, generate_output_path, generate_report_path
from constants import CONFIG_FILE, MSA_OUTPUT_DIR, ORFS_FILE_SUFFIX, FASTA_EXTENSION, CORRECT_ORF_FILE_SUFFIX, BLAST_REPORTS_DIR
from constants import EX1_GB, EX2_DB, EX2_QUERY, EX2_LOCAL, EX2_REPORT, EX3_SEQS, EX3_OUT, EX4_BLAST, EX4_PATTERN, EX5_OUT, EX5_SEQ, EX5_SIZE
from error_helper import exit_with_error


def main():

    # Parse arguments
    parser = argparse.ArgumentParser(description="Bioinformatics Sequencing")

    # Add arguments
    parser.add_argument('-e', '--exercise', required=True, help = 'exercise to run, can be 1 to 5.')          # Ejercicio para correr
    parser.add_argument('-c', '--use_config', action='store_true', help='add this flag to use the config file.')   # Si se debe usar el archivo de config o los parametros
    parser.add_argument('-cf', '--config_file', help='name of the config file you want to use, configuration.ini by default.',
                        type=str, default = CONFIG_FILE, required=False)
    ######################### Exercise 1 params #########################
    parser.add_argument('-gb', '--' + EX1_GB, help='EX1: identifier of genbank input file.',
                        type=str, required=False)
    ######################### Exercise 2 params #########################
    parser.add_argument('-db', '--' + EX2_DB, help='EX2: database to use for remote consults. Can be swissprot (default) or nr.',
                        type=str, default='swissprot', required=False)
    parser.add_argument('-q', '--' + EX2_QUERY, help='EX2: Identifier of fasta file to query.',
                        type=str, required=False)
    parser.add_argument('-l', '--' + EX2_LOCAL, action='store_true', help = 'EX2: Add if BLAST consult should be local. Will only use swissprot if so.')
    parser.add_argument('-r', '--' + EX2_REPORT, help='EX2: Report output file, \'myblast\' by default.',
                        type=str, default='myblast', required=False)
    ######################### Exercise 3 params #########################
    parser.add_argument('-ss', '--' + EX3_SEQS, help='EX3: File with fasta sequences to do MSA.',
                        type=str, required=False)
    parser.add_argument('-out', '--' + EX3_OUT, help='EX3: Output file name, \'msa_output.out\' by default.',
                        default = 'msa_output.out', type=str, required=False)
    ######################### Exercise 4 params #########################
    parser.add_argument('-b', '--' + EX4_BLAST, help='EX4: Blast report file name to use as input.',
                        type=str, required=False)
    parser.add_argument('-p', '--' + EX4_PATTERN, help='EX4: Pattern to find in description of blast report.',
                        type=str, required=False)    
    ######################### Exercise 5 params #########################
    parser.add_argument('-seq', '--' + EX5_SEQ, help='EX5: File with one or more nucleotide sequences.',
                        type=str, required=False)
    parser.add_argument('-outseq', '--' + EX5_OUT, help='EX5: File where possible AA ORFs will be printed, \'protein_orf\' by default.',
                        default = 'protein_orf', type=str, required=False)
    parser.add_argument('-msize', '--' + EX5_SIZE, help='EX5: Min size ORF can have, 2700 by default.',
                        default=2700, required=False)

    try:
        args = parser.parse_args()

        item = int(args.exercise)

        input_file = ""
        output_file = ""
        prosite_file = ""
        database = "swissprot"
        args_to_use = handle_arguments(item, args)

    except Exception as e:
        exit_with_error("Invalid arguments. Check if the input values are correct.", e)

    # Param parsing and setup
    try:
        ###### EXERCISE 1 ARGUMENT HANDLING ######
        if item == 1 and args_to_use[EX1_GB] != None: 
            # Check only the valid identifiers are used
            if not args_to_use[EX1_GB] in ["1A", "1B"]:
                exit_with_error("Invalid genbank identifier. Must be 1A or 1B.")

            # Generate the input file name and the output file name
            input_file = file_param_to_file_name(str(args_to_use[EX1_GB]))
            output_file = generate_output_path(str(args_to_use[EX1_GB]))

        ###### EXERCISE 2 ARGUMENT HANDLING ######
        elif item == 2 and args_to_use[EX2_QUERY] != None and args_to_use[EX2_REPORT] != None and args_to_use[EX2_DB] != None:
            # Check only the valid identifiers are used
            if not args_to_use[EX2_QUERY] in ["2A", "2B", "2AProtein", "2BProtein"]:
                exit_with_error("Invalid query identifier. Must be 2A, 2AProtein, 2B, or 2BProtein.")

            if not check_file_is_not_dir(str(args_to_use[EX2_REPORT])):
                exit_with_error("Report should be a filename, not a directory.")

            # Generate the input file name and the output file name
            input_file = file_param_to_file_name(str(args_to_use[EX2_QUERY]))
            output_file = generate_report_path(str(args_to_use[EX2_REPORT]).split(".")[0], bool(args_to_use[EX2_LOCAL]))
            # handle database
            database = str(args_to_use[EX2_DB])
            if not database in ["swissprot", "nr"]:
                exit_with_error("Invalid database option. Must be swissprot or nr.")
        
        ###### EXERCISE 3 ARGUMENT HANDLING ######
        elif item == 3 and args_to_use[EX3_SEQS] != None and args_to_use[EX3_OUT] != None:
            if not check_file_is_not_dir(str(args_to_use[EX3_OUT])):
                exit_with_error("Output file should be a filename, not a directory.")

            sequences = args_to_use[EX3_SEQS]
            output_file = MSA_OUTPUT_DIR + args_to_use[EX3_OUT]

        ###### EXERCISE 4 ARGUMENT HANDLING ######
        elif item == 4 and args_to_use[EX4_BLAST] != None and args_to_use[EX4_PATTERN] != None:
            input_file = args_to_use[EX4_BLAST]
        
        ###### EXERCISE 5 ARGUMENT HANDLING ######
        elif item == 5 and args_to_use[EX5_SEQ] != None and args_to_use[EX5_OUT] != None and args_to_use[EX5_SIZE] != None:
            if not check_file_is_not_dir(str(args_to_use[EX5_OUT])):
                exit_with_error("Output file should be a filename, not a directory.")

            output_file = args_to_use[EX5_OUT].split(".")[0]
            input_file = args_to_use[EX5_SEQ]

        else:
            exit_with_error("Missing or invalid combination of params. Check the manual.")

    except Exception as e:
        exit_with_error("Error while parsing arguments.", e)


    # Run the exercise with the parsed params
    print("[INFO] Running exercise", item, "...")
    print_curr_arguments(item, args_to_use)
    start_time = time.time()

    try:
        if item == 1:
            nucleotide_file = output_file + '_nucleotides' + FASTA_EXTENSION
            proteins_file = output_file + ORFS_FILE_SUFFIX + FASTA_EXTENSION
            final_protein_file = output_file + CORRECT_ORF_FILE_SUFFIX + FASTA_EXTENSION
            run_exercise_1(input_file, nucleotide_file, proteins_file, final_protein_file)
        
        elif item == 2:
            run_exercise_2(input_file, output_file, bool(args_to_use[EX2_LOCAL]), database)
        
        elif item == 3:
            run_exercise_3(sequences, output_file)

        elif item == 4:
            run_exercise_4(input_file, args_to_use[EX4_PATTERN])
        
        elif item == 5:
            orf_file = output_file + ".orf"
            prosite_file = output_file + ".patmatmotifs"
            run_exercise_5(input_file, orf_file, prosite_file, int(args_to_use[EX5_SIZE]))

    except Exception as e:
        exit_with_error("Error when running the exercise.", e)
    
    execution_time = int(time.time() - start_time)
    print("[DONE] Execution took %i mins %i seconds" % (int(execution_time / 60), execution_time % 60))


if __name__ == '__main__':
    main()
