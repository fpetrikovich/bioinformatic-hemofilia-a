import argparse
import time

from Ex1 import run_exercise_1
from file_helper import file_param_to_file_name, generate_output_path

def main():
    # Get the current time for the output files
    str_time = time.strftime("%Y%m%d-%H%M%S");

    # Parse arguments
    parser = argparse.ArgumentParser(description="Bioinformatics Sequencing")

    # Add arguments
    parser.add_argument('-e', dest='exercise', required=True)   # Ejercicio para correr
    parser.add_argument('-i', '--input', help='identifier of genbank input file',
                        type=str, default='1A', required=False)
    
    args = parser.parse_args()

    input_file = ""
    output_file = ""
    item = 0

    # Param parsing and setup
    try:
        item = int(args.exercise)
        if args.input != None:
            input_file = file_param_to_file_name(str(args.input))
            output_file = generate_output_path(str(args.input) + "_" + str_time)

    except:
        print("[ERROR] Invalid option input")
        exit(0)

    # Run the exercise with the parsed params
    print("[INFO] Running exercise", item, "...")
    if item == 1:
        nucleotide_file = output_file + '_nucleotides.fasta'
        proteins_file = output_file + '_proteins.fasta'
        run_exercise_1(input_file, nucleotide_file, proteins_file)

if __name__ == '__main__':
    main()
