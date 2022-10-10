import argparse
import time

from Ex1 import run_exercise_1
from Ex2 import run_exercise_2
from Ex3 import run_exercise_3
from file_helper import file_param_to_file_name, generate_output_path
from constants import ORFS_FILE_SUFFIX, FASTA_EXTENSION, CORRECT_ORF_FILE_SUFFIX


def main():
    # Get the current time for the output files
    str_time = time.strftime("%Y%m%d-%H%M%S");

    # Parse arguments
    parser = argparse.ArgumentParser(description="Bioinformatics Sequencing")

    # Add arguments
    parser.add_argument('-e', '--exercise', required=True)   # Ejercicio para correr
    parser.add_argument('-gb', '--genbank', help='identifier of genbank input file',
                        type=str, required=False)
    parser.add_argument('-q', '--query', help='Identifier of fasta file to query',
                        type=str, required=False)
    parser.add_argument('-r', '--report', help='Report output file',
                        type=str, default='myblast', required=False)
    parser.add_argument('-origin', '--origin', help='With origin file to compare',
                        type=str, required=False)
    parser.add_argument('-compare', '--compare', help='Files to compare',
                        type=str, required=False, nargs='+')
    parser.add_argument('-out', '--output', help='Output file name',
                        type=str, required=False)
    
    args = parser.parse_args()

    input_file = ""
    output_file = ""
    output_file_local = ""
    item = 0

    # Param parsing and setup
    try:
        item = int(args.exercise)
        if args.genbank != None: 
            input_file = file_param_to_file_name(str(args.genbank))
            output_file = generate_output_path(str(args.genbank))

        elif args.query != None and args.report != None:
            input_file = file_param_to_file_name(str(args.query))
            output_file = 'reports/' + args.report
            output_file_local = 'reports/local_' + args.report

        elif args.origin != None and args.compare != None and args.output != None:
            origin_sequence = args.origin
            species_list = args.compare
            output_file = args.output


    except Exception as e:
        print("[ERROR] " + str(e))
        exit(0)

    # Run the exercise with the parsed params
    print("[INFO] Running exercise", item, "...")
    if item == 1:
        nucleotide_file = output_file + '_nucleotides.faa'
        proteins_file = output_file + ORFS_FILE_SUFFIX + FASTA_EXTENSION
        final_protein_file = output_file + CORRECT_ORF_FILE_SUFFIX + FASTA_EXTENSION
        run_exercise_1(input_file, nucleotide_file, proteins_file, final_protein_file)
    
    elif item == 2:
        run_exercise_2(input_file, output_file, output_file_local)

    elif item == 3:
        run_exercise_3(origin_sequence, species_list, output_file)

if __name__ == '__main__':
    main()
