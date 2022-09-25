import argparse
import time

from Ex1 import run_exercise_1
from Ex2 import run_exercise_2
from file_helper import file_param_to_file_name, generate_output_path


def main():
    # Get the current time for the output files
    str_time = time.strftime("%Y%m%d-%H%M%S");

    # Parse arguments
    parser = argparse.ArgumentParser(description="Bioinformatics Sequencing")

    # Add arguments
    parser.add_argument('-e', dest='exercise', required=True)   # Ejercicio para correr
    parser.add_argument('-gb', '--genbank', help='identifier of genbank input file',
                        type=str, required=False)
    parser.add_argument('-q', '--query', help='FASTA file for BLAST query',
                        type=str, required=False)
    parser.add_argument('-r', '--report', help='Report output file',
                        type=str, default='myblast.report', required=False)
    
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
            output_file = generate_output_path(str(args.genbank) + "_" + str_time)

        elif args.query != None and args.report != None:
            input_file = file_param_to_file_name(str(args.query))
            output_file = 'reports/' + args.report
            output_file_local = 'reports/local_' + args.report
    
    except:
        print("[ERROR] Invalid option input")
        exit(0)

    # Run the exercise with the parsed params
    print("[INFO] Running exercise", item, "...")
    if item == 1:
        nucleotide_file = output_file + '_nucleotides.fasta'
        proteins_file = output_file + '_proteins.fasta'
        run_exercise_1(input_file, nucleotide_file, proteins_file)
    
    elif item == 2:
        run_exercise_2(input_file, output_file, output_file_local)

if __name__ == '__main__':
    main()
