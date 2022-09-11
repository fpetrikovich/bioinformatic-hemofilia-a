import argparse

from exercise1 import run_exercise_1
from file_helper import file_param_to_file_name

def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Bioinformatics Sequencing")

    # Add arguments
    parser.add_argument('-f', dest='file', required=True)   # Archivo para usar
    parser.add_argument('-e', dest='exercise', required=True)   # Archivo para usar

    # The following two are if using Python 3.9 and up
    #parser.add_argument('-v', dest='verbose', action=argparse.BooleanOptionalAction, default=False)  # Verbose, print or not
    #parser.add_argument('-vv', dest='veryVerbose', action=argparse.BooleanOptionalAction, default=False)  # Verbose, print or not

    # The following are for Python 3.8 and under
    parser.add_argument('-v', dest='verbose', action='store_true')  # Verbose, print or not
    parser.add_argument('-vv', dest='veryVerbose', action='store_true')  # Verbose, print or not
    args = parser.parse_args()

    file_name = ""

    # Param parsing and setup
    try:
        item = int(args.exercise)
        if args.file != None:
            file_name = file_param_to_file_name(str(args.file))
    except:
        print("[ERROR] Invalid option input")
        exit(0)

    # Run the exercise with the parsed params
    print("[INFO] Running exercise", item, "...")
    if item == 1:
        run_exercise_1(file_name)

if __name__ == '__main__':
    main()
