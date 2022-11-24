import argparse
import time
import configparser

from Ex1 import run_exercise_1
from Ex2 import run_exercise_2
from Ex3 import run_exercise_3
from Ex4 import run_exercise_4
from Ex5 import run_exercise_5
from file_helper import file_exists, file_param_to_file_name, generate_output_path, generate_report_path
from constants import ARGUMENTS, CONFIG_FILE, CONFIG_DICT, MSA_OUTPUT_DIR, ORFS_FILE_SUFFIX, FASTA_EXTENSION, CORRECT_ORF_FILE_SUFFIX, BLAST_REPORTS_DIR, EMBOSS_DIR
from constants import EX1_GB, EX2_DB, EX2_QUERY, EX2_LOCAL, EX2_REPORT, EX3_SEQS, EX3_OUT, EX4_BLAST, EX4_PATTERN, EX5_OUT, EX5_SEQ

def create_config(config_file, exercise):
    if not file_exists(config_file):
        print("[ERROR] File %s does not exists in the root directory." % (config_file))
        exit(1)
    config = configparser.ConfigParser()
    config.read(config_file)
    exercise_config = config[CONFIG_DICT[exercise]]
    print(exercise_config)
    return exercise_config

def handle_arguments(exercise, args):
    use_config = bool(args.use_config)
    arguments = ARGUMENTS
    if (use_config):
        config_file = str(args.config_file)
        # Get the arguments from the config file
        config = create_config(config_file, exercise)
        for key in config: 
            value = config.get(key)
            # Convert the value too boolean if necessary
            if value == "yes" or value == "no": value = config.getboolean(key) 
            arguments[key] = value
    else:
        # Get the arguments from the command line
        arguments[EX1_GB] = args.genbank
        arguments[EX2_DB] = args.database
        arguments[EX2_LOCAL] = args.local
        arguments[EX2_QUERY] = args.query
        arguments[EX2_REPORT] = args.report
        arguments[EX3_SEQS] = args.sequences
        arguments[EX3_OUT] = args.output
        arguments[EX4_BLAST] = args.blast
        arguments[EX4_PATTERN] = args.pattern
        arguments[EX5_SEQ] = args.sequence
        arguments[EX5_OUT] = args.outputseq
    print(arguments)
    return arguments

def main():
    # Get the current time for the output files
    str_time = time.strftime("%Y%m%d-%H%M%S");

    # Parse arguments
    parser = argparse.ArgumentParser(description="Bioinformatics Sequencing")

    # Add arguments
    parser.add_argument('-e', '--exercise', required=True)          # Ejercicio para correr
    parser.add_argument('-c', '--use_config', action='store_true')   # Si se debe usar el archivo de config o los parametros
    parser.add_argument('-cf', '--config_file', help='name of the config file you want to use',
                        type=str, default = CONFIG_FILE, required=False)
    ######################### Exercise 1 params #########################
    parser.add_argument('-gb', '--' + EX1_GB, help='identifier of genbank input file',
                        type=str, required=False)
    ######################### Exercise 2 params #########################
    parser.add_argument('-db', '--' + EX2_DB, help='database to use for remote consults (swissprot or nr)',
                        type=str, default='swissprot', required=False)
    parser.add_argument('-q', '--' + EX2_QUERY, help='Identifier of fasta file to query',
                        type=str, required=False)
    parser.add_argument('-l', '--' + EX2_LOCAL, action='store_true')
    parser.add_argument('-r', '--' + EX2_REPORT, help='Report output file',
                        type=str, default='myblast', required=False)
    ######################### Exercise 3 params #########################
    parser.add_argument('-ss', '--' + EX3_SEQS, help='File with fasta sequences to do MSA',
                        type=str, required=False)
    parser.add_argument('-out', '--' + EX3_OUT, help='Output file name',
                        type=str, required=False)
    ######################### Exercise 4 params #########################
    parser.add_argument('-b', '--' + EX4_BLAST, help='Blast report file name to use as input',
                        type=str, required=False)
    parser.add_argument('-p', '--' + EX4_PATTERN, help='Pattern to find in description of blast report',
                        type=str, required=False)    
    ######################### Exercise 5 params #########################
    parser.add_argument('-seq', '--' + EX5_SEQ, help='File with one or more nucleotide sequences',
                        type=str, required=False)
    parser.add_argument('-outseq', '--' + EX5_OUT, help='File where possible AA ORFs will be printed',
                        type=str, required=False)

    try:
        args = parser.parse_args()

        item = int(args.exercise)

        input_file = ""
        output_file = ""
        prosite_file = ""
        database = "swissprot"
        args_to_use = handle_arguments(item, args)

    except Exception as e:
        print("[ERROR] Invalid arguments. Check if the input values are correct.\n[MESSAGE] " + str(e))
        exit(1)

    # Param parsing and setup
    try:
        ###### EXERCISE 1 ARGUMENT HANDLING ######
        if item == 1 and args_to_use[EX1_GB] != None: 
            # Check only the valid identifiers are used
            if not args_to_use[EX1_GB] in ["1A", "1B"]:
                print("[ERROR] Invalid genbank identifier. Must be 1A or 1B.")
                exit(1)
            # Generate the input file name and the output file name
            input_file = file_param_to_file_name(str(args_to_use[EX1_GB]))
            output_file = generate_output_path(str(args_to_use[EX1_GB]))

        ###### EXERCISE 2 ARGUMENT HANDLING ######
        elif item == 2 and args_to_use[EX2_QUERY] != None:
            # Check only the valid identifiers are used
            if not args_to_use[EX2_QUERY] in ["2A", "2B", "2AProtein", "2BProtein"]:
                print("[ERROR] Invalid query identifier. Must be 2A, 2AProtein, 2B, or 2BProtein.")
                exit(1)
            # Generate the input file name and the output file name
            input_file = file_param_to_file_name(str(args_to_use[EX2_QUERY]))
            output_file = generate_report_path(str(args_to_use[EX2_REPORT]), bool(args_to_use[EX2_LOCAL]))
            # handle database if parameter was added
            if args_to_use[EX2_DB] != None:
                database = str(args_to_use[EX2_DB])
                if not database in ["swissprot", "nr"]:
                    print("[ERROR] Invalid database option. Must be swissprot or nr.")
                    exit(1)
        
        ###### EXERCISE 3 ARGUMENT HANDLING ######
        elif item == 3 and args_to_use[EX3_SEQS] != None and args_to_use[EX3_OUT] != None:
            sequences = args_to_use[EX3_SEQS]
            output_file = MSA_OUTPUT_DIR + args_to_use[EX3_OUT]

        ###### EXERCISE 4 ARGUMENT HANDLING ######
        elif item == 4 and args_to_use[EX4_BLAST] != None and args_to_use[EX4_PATTERN] != None:
            input_file = BLAST_REPORTS_DIR + args_to_use[EX4_BLAST]
        
        ###### EXERCISE 5 ARGUMENT HANDLING ######
        elif item == 5 and args_to_use[EX5_SEQ] != None:
            output_file = "protein_orfs" if args_to_use[EX5_OUT] == None else args_to_use[EX5_OUT].split(".")[0]
            input_file = args_to_use[EX5_SEQ]

        else:
            print("[ERROR] Missing or invalid combination of params. Check the manual.")
            exit(1)

    except Exception as e:
        print("[ERROR] " + str(e))
        exit(1)

    # Run the exercise with the parsed params
    print("[INFO] Running exercise", item, "...")
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
            run_exercise_5(input_file, orf_file, prosite_file)

    except Exception as e:
        print("[ERROR] Unknown error when running the exercise.")
        print("[ERROR][MESSAGE] " + str(e))
        exit(1)
    
    execution_time = int(time.time() - start_time)
    print("[DONE] Execution took %i mins %i seconds" % (int(execution_time / 60), execution_time % 60))


if __name__ == '__main__':
    main()
